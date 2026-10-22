#!/usr/bin/env python3
"""HMP V35 stool: collapse OTUs to genus, keep the first visit per subject, write Y.csv / X.csv / features.csv."""
import os, pandas as pd, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "..", "data", "microbiome")
Y = pd.read_csv(os.path.join(D, "hmp16s", "hmp16s_v35_stool_counts.csv"), index_col=0); Y.columns = Y.columns.astype(str)
tax = pd.read_csv(os.path.join(D, "hmp16s", "hmp16s_v35_taxonomy.csv"), index_col=0)
meta = pd.read_csv(os.path.join(D, "hmp16s", "hmp16s_v35_stool_metadata.csv"), index_col=0); meta.index = meta.index.astype(str)
tax = tax.loc[Y.index]
def key(r):
    for rank in ["GENUS", "FAMILY", "ORDER", "CLASS", "PHYLUM"]:
        v = r[rank]
        if isinstance(v, str) and v and v.lower() not in ("unclassified", "nan"):
            return (v if rank == "GENUS" else f"{v}_unclassified"), ";".join(f"{q[0].lower()}__{r[q] if isinstance(r[q], str) else ''}" for q in ["PHYLUM","CLASS","ORDER","FAMILY","GENUS"])
    return "unclassified", "unclassified"
kk = tax.apply(key, axis=1)
name = pd.Series([k[0] for k in kk], index=Y.index); lineage = pd.Series([k[1] for k in kk], index=Y.index)
Yg = Y.groupby(name.values).sum()
first = meta.sort_values("VISITNO").drop_duplicates("RSID")          # one sample per subject
keep = [s for s in Yg.columns if s in first.index]
Yg = Yg[keep]; X = first.loc[keep, ["RSID", "VISITNO", "SEX", "RUN_CENTER"]].copy(); X.insert(0, "trt", "H")
nz = (Yg > 0).sum(axis=1); Yg = Yg[nz >= 2]
feat = pd.DataFrame({"lineage": lineage.groupby(name.values).first().loc[Yg.index], "n_otus": name.value_counts().loc[Yg.index],
                     "prevalence": (Yg > 0).mean(axis=1).round(3), "mean_proportion": (Yg / Yg.sum()).mean(axis=1)}).sort_values("mean_proportion", ascending=False)
out = os.path.join(D, "processed", "hmp_stool_genus"); os.makedirs(out, exist_ok=True)
Yg.astype(int).to_csv(os.path.join(out, "Y.csv"), index_label=""); X.to_csv(os.path.join(out, "X.csv"), index_label=""); feat.to_csv(os.path.join(out, "features.csv"), index_label="feature")
lib = Yg.sum()
print(f"[hmp_stool] {Y.shape[0]} OTUs x {Y.shape[1]} samples -> {Yg.shape[0]} genera x {Yg.shape[1]} subjects (first visit); library median {int(lib.median())}, range {int(lib.min())}-{int(lib.max())}; zero fraction {float((Yg==0).mean().mean()):.3f}")
print("   top:", ", ".join(f"{i} {v:.3f}" for i, v in feat.mean_proportion.head(6).items()))
