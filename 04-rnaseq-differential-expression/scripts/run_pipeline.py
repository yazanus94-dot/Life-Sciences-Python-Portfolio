"""Run differential expression.

Usage:
  python scripts/run_pipeline.py --counts data/airway_counts.csv --coldata data/airway_coldata.csv \
      --factors cell dex --contrast dex --ref untrt --out results
counts: genes x samples CSV, first column = gene symbol (or ID).
coldata: samples x covariates CSV, first column = sample name.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rnaseq_de import model, normalize, plots, validate  # noqa: E402

np.random.seed(42)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--coldata", required=True)
    ap.add_argument("--factors", nargs="+", default=["cell", "dex"])
    ap.add_argument("--contrast", default="dex")
    ap.add_argument("--ref", default="untrt")
    ap.add_argument("--out", default="results")
    ap.add_argument("--skip-validation", action="store_true", help="skip the airway known-gene check (for non-airway data)")
    a = ap.parse_args()

    out = Path(a.out); (out / "figures").mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(a.counts, index_col=0)
    coldata = pd.read_csv(a.coldata, index_col=0).loc[counts.columns]
    counts = counts.groupby(level=0).sum()  # collapse duplicated gene labels
    print(f"{counts.shape[0]} genes x {counts.shape[1]} samples")

    counts = normalize.filter_low_counts(counts)
    sf = normalize.median_of_ratios(counts)
    print(f"{counts.shape[0]} genes after filtering; size factors:\n{sf.round(3).to_string()}")
    X, target = model.design_matrix(coldata, a.factors, a.contrast, a.ref)
    res = model.fit_genes(counts, sf, X, target)
    res.to_csv(out / "de_results.csv")
    sig = res[(res.padj < 0.05)]
    print(f"{len(sig)} genes with FDR < 0.05 ({(sig.log2FC > 0).sum()} up, {(sig.log2FC < 0).sum()} down)")

    if not a.skip_validation:
        val = validate.check_known_genes(res)
        val.to_csv(out / "known_gene_validation.csv", index=False)
        print(val.round(4).to_string(index=False))

    plots.volcano(res, out / "figures" / "volcano.png", label_genes=validate.KNOWN_UP)
    plots.ma_plot(res, out / "figures" / "ma_plot.png")
    plots.pca_plot(counts, sf, coldata, a.factors[0], a.contrast, out / "figures" / "pca.png")
    plots.top_gene_heatmap(counts, sf, res, out / "figures" / "top_genes_heatmap.png")


if __name__ == "__main__":
    main()
