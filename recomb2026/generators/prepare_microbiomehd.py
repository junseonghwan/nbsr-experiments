#!/usr/bin/env python3
"""Build Y.csv / X.csv for NBSR from a MicrobiomeHD dataset.

Y.csv: features x samples counts (first column = feature name, header = sample ids),
X.csv: samples x metadata with a 'trt' column (DiseaseState), same sample order as Y's columns,
features.csv: feature name -> full lineage, number of OTUs collapsed, prevalence.

Usage: prepare_microbiomehd.py <dataset> [--level genus|otu] [--min-samples N] [--otu-prevalence F]
  dataset in: crc_baxter cdi_schubert ibd_gevers_2014
"""
import argparse, os, re, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "microbiome")
RAW = os.path.join(DATA, "microbiomehd")
OUT = os.path.join(DATA, "processed")

CONFIG = {
    "crc_baxter": dict(id_col="Sample_Name_s", norm=lambda s: s,
        covariates=["DiseaseState", "diagnosis_s", "Age_s", "Gender_s", "BMI_s", "fit_result_s", "Smoke_s", "Diabetic_s"]),
    "cdi_schubert": dict(id_col="sample_id", norm=lambda s: s,
        covariates=["DiseaseState", "disease_stat", "age", "gender", "antibiotics >3mo", "protonpump", "race"]),
    "ibd_gevers_2014": dict(id_col="sample", norm=lambda s: s.replace(".", ""),
        covariates=["DiseaseState", "Diagnosis", "AB_exposure", "PCDAI", "race", "sample_location", "subject"]),
}
RANKS = ["k", "p", "c", "o", "f", "g", "s", "d"]

def parse_lineage(s):
    tok = dict()
    for part in s.split(";"):
        if "__" in part:
            r, name = part.split("__", 1)
            tok[r] = name
    return tok

def genus_key(tok):
    """Lineage up to genus; unclassified genera fall back to the lowest assigned rank."""
    if tok.get("g", ""):
        return ";".join(f"{r}__{tok.get(r,'')}" for r in RANKS[:6]), tok["g"]
    for r in ["f", "o", "c", "p", "k"]:
        if tok.get(r, ""):
            i = RANKS.index(r)
            return ";".join(f"{q}__{tok.get(q,'')}" for q in RANKS[:i+1]) + ";g__unclassified", f"{tok[r]}_unclassified"
    return "unclassified", "unclassified"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset", choices=CONFIG)
    ap.add_argument("--level", choices=["genus", "otu"], default="genus")
    ap.add_argument("--min-samples", type=int, default=2, help="keep features with counts in >= N samples")
    ap.add_argument("--otu-prevalence", type=float, default=0.25, help="otu level: keep OTUs present in >= this fraction of samples")
    a = ap.parse_args()
    cfg = CONFIG[a.dataset]
    d = os.path.join(RAW, f"{a.dataset}_results")
    rdp = pd.read_csv(os.path.join(d, "RDP", f"{a.dataset}.otu_table.100.denovo.rdp_assigned"), sep="\t", index_col=0)
    rdp.columns = [str(c) for c in rdp.columns]
    meta = pd.read_csv(os.path.join(d, f"{a.dataset}.metadata.txt"), sep="\t", encoding="latin1", dtype=str)
    meta.index = [cfg["norm"](str(x)) for x in meta[cfg["id_col"]]]
    meta = meta[~meta.index.duplicated()]
    common = [s for s in rdp.columns if s in meta.index]
    print(f"[{a.dataset}] OTU table {rdp.shape[0]} OTUs x {rdp.shape[1]} samples; metadata {meta.shape[0]}; matched {len(common)}")
    Y = rdp[common]
    meta = meta.loc[common, [c for c in cfg["covariates"] if c in meta.columns]].copy()

    if a.level == "genus":
        keys = [genus_key(parse_lineage(ix)) for ix in Y.index]
        lineage = pd.Series([k for k, _ in keys], index=Y.index)
        display = dict(keys)                      # lineage key -> genus name
        Yg = Y.groupby(lineage.values).sum()      # rows indexed by lineage key
        n_otu = lineage.value_counts()
        # unique readable names: append family, then order, ... until unique
        names = pd.Series({k: display[k] for k in Yg.index})
        for rank in ["f", "o", "c", "p"]:
            dup = names[names.duplicated(keep=False)].index
            if len(dup) == 0:
                break
            for k in dup:
                names[k] = f"{names[k]}__{parse_lineage(k).get(rank, '')}"
        assert names.is_unique, "genus names still not unique"
        feat = pd.DataFrame({"lineage": Yg.index, "n_otus": [int(n_otu[k]) for k in Yg.index]}, index=names.loc[Yg.index].values)
        Yf = Yg.copy(); Yf.index = feat.index
    else:
        prev = (Y > 0).mean(axis=1)
        for thr in [0.05, 0.10, 0.25, 0.50]:
            print(f"   OTUs present in >= {int(thr*100)}% of samples: {(prev >= thr).sum()}")
        Yf = Y[prev >= a.otu_prevalence].copy()
        Yf.index = [re.search(r"d__(\w+)", ix).group(1) for ix in Yf.index]
        feat = pd.DataFrame({"lineage": [ix for ix in Y.index[prev >= a.otu_prevalence]], "n_otus": 1}, index=Yf.index)

    # drop empty samples, then rare features
    lib = Yf.sum(axis=0)
    keep_s = lib > 0
    Yf = Yf.loc[:, keep_s]; meta = meta.loc[Yf.columns]
    nz = (Yf > 0).sum(axis=1)
    dropped = int((nz < a.min_samples).sum())
    Yf = Yf[nz >= a.min_samples]
    feat = feat.loc[Yf.index]
    feat["prevalence"] = (Yf > 0).mean(axis=1).round(3)
    feat["mean_proportion"] = (Yf / Yf.sum(axis=0)).mean(axis=1)
    feat = feat.sort_values("mean_proportion", ascending=False)

    X = meta.copy(); X.insert(0, "trt", X["DiseaseState"])
    out = os.path.join(OUT, f"{a.dataset}_{a.level}"); os.makedirs(out, exist_ok=True)
    Yf.astype(int).to_csv(os.path.join(out, "Y.csv"), index_label="")
    X.to_csv(os.path.join(out, "X.csv"), index_label="")
    feat.to_csv(os.path.join(out, "features.csv"), index_label="feature")

    lib = Yf.sum(axis=0); P = Yf / lib
    print(f"   -> {out}")
    print(f"   Y: {Yf.shape[0]} features x {Yf.shape[1]} samples ({dropped} features dropped for < {a.min_samples} nonzero samples)")
    print(f"   groups: {X['trt'].value_counts().to_dict()}")
    q = lib.quantile([0, .05, .5, .95, 1]).astype(int).to_dict()
    print(f"   library size quantiles (0,5,50,95,100%): {list(q.values())}")
    print(f"   zero fraction: {float((Yf == 0).mean().mean()):.3f}; features with prevalence >= 50%: {int((feat.prevalence >= .5).sum())}")
    top = feat.head(8)
    print("   top features (mean proportion):", ", ".join(f"{i} {v:.3f}" for i, v in top.mean_proportion.items()))
    print(f"   share of top 10 features: {feat.mean_proportion.head(10).sum():.3f}")

if __name__ == "__main__":
    main()
