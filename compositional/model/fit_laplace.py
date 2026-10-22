"""Fit NBSR at its posterior mode and draw from the Laplace approximation.

Usage: fit_laplace.py --input DIR --out DIR --seed N [options]

Called by run_fit.R for the `laplace` algorithm. The input directory holds
the Y.csv / X.csv that R wrote from its shared input preparation (samples
with no reads already dropped, `trt` a numeric 0/1 column with the control
group at 0), so NBSR is fit to exactly the data the Stan methods see, under
the same design [intercept, trt] and the same reference feature (--pivot:
the last feature's coefficients are fixed at zero, as in nbsr.stan).

NBSR is imported from the repository's submodule (external/nbsr) and not
from whatever copy happens to be installed: run_fit.R puts the submodule
first on PYTHONPATH, and this script records which file was imported so the
completion marker proves it.

Outputs, in --out:
  beta_draws.bin     S x P x (J-1) float64, C order, little-endian
  laplace_fit.json   dims, layout, seed, settings, loss diagnostics, versions
plus NBSR's own files (nbsr_beta.csv, checkpoint.pth, hessian.npy, ...).
"""
import json
import os
import shutil
import sys
import time
from pathlib import Path

import click
import numpy as np
import torch

REPO_ROOT = Path(__file__).resolve().parent.parent
NBSR_SUBMODULE = REPO_ROOT / "external" / "nbsr"

DRAWS_FILE = "beta_draws.bin"
FIT_INFO_FILE = "laplace_fit.json"


def import_nbsr():
    if not (NBSR_SUBMODULE / "nbsr" / "__init__.py").is_file():
        raise click.ClickException(
            f"NBSR submodule is empty at {NBSR_SUBMODULE}. Run: git submodule update --init")
    try:
        import nbsr
        from nbsr.fit import fit_once, hessian_filename
        from nbsr.nbsr_config import NBSRConfig
        from nbsr.prior import apply_prior, elicit_prior
    except ImportError as error:
        raise click.ClickException(f"Cannot import NBSR: {error}") from error
    imported_from = Path(nbsr.__file__).resolve()
    if NBSR_SUBMODULE.resolve() not in imported_from.parents:
        raise click.ClickException(
            f"nbsr was imported from {imported_from}, not from the submodule {NBSR_SUBMODULE}. "
            "run_fit.R puts the submodule first on PYTHONPATH; do not run this script directly "
            "with another copy of NBSR ahead of it.")
    return nbsr, fit_once, hessian_filename, NBSRConfig, elicit_prior, apply_prior


def package_version(name):
    try:
        from importlib.metadata import version
        return version(name)
    except Exception:  # noqa: BLE001 - provenance only
        return None


@click.command()
@click.option("--input", "input_dir", type=click.Path(exists=True, file_okay=False), required=True,
              help="Directory with Y.csv and X.csv as written by run_fit.R.")
@click.option("--out", "out_dir", type=click.Path(file_okay=False), required=True)
@click.option("--seed", type=int, required=True,
              help="Seeds the Adam initialisation and the Laplace draws (the unit's fit seed).")
@click.option("--iterations", type=int, default=20000, show_default=True, help="Adam iterations.")
@click.option("--lr", type=float, default=0.05, show_default=True, help="Adam learning rate.")
@click.option("--draws", type=int, default=16000, show_default=True,
              help="Draws of beta from N(beta_hat, -H^-1).")
@click.option("--runs", type=int, default=1, show_default=True,
              help="Independent initialisations (seed, seed+1, ...); the best log posterior is kept.")
@click.option("--beta-prior-sd", type=float, required=True,
              help="Prior sd of the treatment row of beta: sqrt(sigma_beta2) of the Stan model.")
@click.option("--intercept-prior-sd", type=float, default=None,
              help="Prior sd of the intercept row of beta; default: --beta-prior-sd.")
@click.option("--beta-prior", type=click.Choice(["fixed", "empirical"]), default="fixed", show_default=True,
              help="fixed: --beta-prior-sd on every covariate. empirical: NBSR's two-stage elicitation, a "
                   "diffuse fit (--stage1-prior-sd) then per-covariate sds matched to the --prior-quantile of "
                   "|beta| (nbsr.prior.elicit_prior); the main fit is warm-started from the diffuse one.")
@click.option("--prior-quantile", type=float, default=0.95, show_default=True,
              help="empirical only: |beta| quantile matched to the prior's tail.")
@click.option("--stage1-prior-sd", type=float, default=10.0, show_default=True,
              help="empirical only: prior sd of every covariate in the diffuse stage-1 fit.")
@click.option("--sigma-b", type=float, required=True,
              help="Prior sd of the coefficients on the dispersion covariates: sqrt(sigma_b2).")
@click.option("--b-pi-loc", type=float, required=True)
@click.option("--b-pi-scale", type=float, required=True)
@click.option("--sigma-bj-prior-sd", type=float, default=0.5, show_default=True,
              help="Half-normal scale of sigma_bj; 0.5 is the literal in nbsr.stan.")
@click.option("--z-column", "z_columns", multiple=True,
              help="Column of X.csv used as an external dispersion covariate (the method's dispersion_covariates).")
def main(input_dir, out_dir, seed, iterations, lr, draws, runs, beta_prior_sd, intercept_prior_sd, beta_prior,
         prior_quantile, stage1_prior_sd, sigma_b, b_pi_loc, b_pi_scale, sigma_bj_prior_sd, z_columns):
    nbsr, fit_once, hessian_filename, NBSRConfig, elicit_prior, apply_prior = import_nbsr()
    if intercept_prior_sd is None:
        intercept_prior_sd = beta_prior_sd
    if intercept_prior_sd <= 0:
        raise click.ClickException("--intercept-prior-sd must be positive.")
    input_dir, out_dir = Path(input_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if runs < 1 or draws < 1 or iterations < 1:
        raise click.ClickException("--runs, --draws and --iterations must be positive.")

    # One BLAS thread per task, as launch_jobs.sh sets for R. torch reads
    # OMP_NUM_THREADS itself; this makes the intent explicit.
    threads = int(os.environ.get("OMP_NUM_THREADS", "0") or 0)
    if threads > 0:
        torch.set_num_threads(threads)

    def config_for(path):
        return NBSRConfig(
            counts_path=input_dir / "Y.csv", coldata_path=input_dir / "X.csv", output_path=path,
            column_names=["trt"], z_columns=list(z_columns), z_log=False, z_total_counts=False,
            b_pi_prior_mean=b_pi_loc, b_pi_prior_sd=b_pi_scale, sigma_bj_prior_sd=sigma_bj_prior_sd,
            sigma_b=sigma_b, dispersion_link="logit", feature_offsets=True,
            lr=lr, iterations=iterations, trended_dispersion=True, pivot=True,
            beta_prior="fixed", beta_prior_sd=[intercept_prior_sd, beta_prior_sd],   # intercept, trt
            prior_quantile=prior_quantile, stage1_prior_sd=stage1_prior_sd,
            use_cuda_if_available=False,
        )

    # Best of `runs` initialisations, each seeded so a rerun reproduces it.
    started = time.perf_counter()
    best = None
    for run_no in range(runs):
        run_seed = seed + run_no
        run_dir = out_dir if runs == 1 else out_dir / f"run{run_no}"
        torch.manual_seed(run_seed)
        np.random.seed(run_seed % (2 ** 32))
        config = config_for(run_dir)
        prior_sd_used = [intercept_prior_sd, beta_prior_sd]
        if beta_prior == "empirical":
            # Stage 1: diffuse fit, then per-covariate sds by quantile matching
            # (writes prior.json beside the fit); stage 2 is warm-started from it.
            prior = elicit_prior(config)
            config = apply_prior(config, prior)
            prior_sd_used = list(prior["beta_prior_sd"])
            torch.manual_seed(run_seed)
        loss_history, model = fit_once(config)
        best_loss = float(np.min(loss_history))
        print(f"run {run_no} (seed {run_seed}): best loss {best_loss:.4f}, "
              f"final loss {loss_history[-1]:.4f}")
        if best is None or best_loss < best["best_loss"]:
            best = dict(run=run_no, seed=run_seed, dir=run_dir, model=model,
                        best_loss=best_loss, loss_history=loss_history, prior_sd_used=prior_sd_used)
    fit_seconds = time.perf_counter() - started
    if runs > 1:
        for f in best["dir"].iterdir():
            if f.is_file():
                shutil.copy2(f, out_dir / f.name)

    model = best["model"]
    loss_history = np.asarray(best["loss_history"], dtype=np.float64)
    if not np.isfinite(loss_history[-1]):
        raise click.ClickException("The fit ended with a non-finite loss.")

    # Draws from N(beta_hat, I^-1), I the negative Hessian of the log posterior
    # at the mode (over beta, with the dispersion model held at its mode). The
    # flat beta is covariate-major: index p * (J-1) + j.
    beta_hat = model.beta.detach().double().flatten()
    P, dim = model.covariate_count, model.dim
    I = torch.from_numpy(np.load(out_dir / hessian_filename)).double()
    if I.shape != (P * dim, P * dim):
        raise click.ClickException(f"Hessian has shape {tuple(I.shape)}, expected {(P * dim, P * dim)}.")
    try:
        L = torch.linalg.cholesky(I)
    except RuntimeError as error:
        raise click.ClickException(
            f"The negative Hessian is not positive definite ({error}); the fit did not reach a mode.") from error
    gen = torch.Generator().manual_seed(seed)
    Z = torch.randn(P * dim, draws, generator=gen, dtype=torch.float64)
    B = torch.linalg.solve_triangular(L.T, Z, upper=True)          # L^-T z ~ N(0, I^-1)
    beta_draws = (beta_hat.unsqueeze(1) + B).T.reshape(draws, P, dim)
    beta_draws.numpy().astype("<f8").tofile(out_dir / DRAWS_FILE)

    # The whole loss history is large and rarely read; keep its shape, plus
    # one convergence number: how much the best loss still improved during
    # the last tenth of the iterations. Adam at this learning rate bounces
    # around the mode rather than settling, so the final loss is not the
    # criterion; the best state is what was kept. A drop far above the
    # bounce (a few units) says the run was still descending.
    checkpoints = np.unique(np.linspace(0, len(loss_history) - 1, 21).astype(int))
    last_tenth = loss_history[int(0.9 * len(loss_history)):]
    best_before_last_tenth = float(np.min(loss_history[:int(0.9 * len(loss_history))]))
    best_loss_drop_last_tenth = best_before_last_tenth - float(np.min(last_tenth))
    info = {
        "draws_file": DRAWS_FILE,
        "dims": {"S": draws, "P": P, "J_minus_1": dim},
        "layout": "C order (S, P, J-1), float64 little-endian; beta flat index p*(J-1)+j",
        "seed": seed,
        "iterations": iterations, "lr": lr, "runs": runs, "best_run": best["run"],
        "best_run_seed": best["seed"],
        "best_loss": best["best_loss"], "final_loss": float(loss_history[-1]),
        "best_loss_drop_last_tenth": max(best_loss_drop_last_tenth, 0.0),
        "final_loss_minus_best": float(loss_history[-1]) - best["best_loss"],
        "loss_at": {int(i): float(loss_history[i]) for i in checkpoints},
        "beta_prior": beta_prior, "beta_prior_sd": beta_prior_sd, "intercept_prior_sd": intercept_prior_sd,
        "beta_prior_sd_used": best["prior_sd_used"],      # per covariate: intercept, trt
        "prior_quantile": prior_quantile, "stage1_prior_sd": stage1_prior_sd,
        "sigma_b": sigma_b,
        "b_pi_loc": b_pi_loc, "b_pi_scale": b_pi_scale, "sigma_bj_prior_sd": sigma_bj_prior_sd,
        "z_columns": list(z_columns),
        "fit_seconds": fit_seconds,
        "nbsr_file": str(Path(nbsr.__file__).resolve()),
        "python": sys.version.split()[0],
        "package_versions": {"torch": torch.__version__, "numpy": np.__version__,
                             "pandas": package_version("pandas"), "scipy": package_version("scipy")},
    }
    with open(out_dir / FIT_INFO_FILE, "w") as fh:
        json.dump(info, fh, indent=2)
    print(f"Wrote {draws} Laplace draws of beta ({P} x {dim}) to {out_dir / DRAWS_FILE}")


if __name__ == "__main__":
    main()
