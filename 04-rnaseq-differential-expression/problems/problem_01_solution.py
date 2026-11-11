"""Reference solution for problems/problem_01_normalisation.md.

A few very highly expressed genes are up-regulated in treated samples and NOTHING
else changes. Compare false discoveries when normalising by total counts (CPM-style)
versus median-of-ratios.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rnaseq_de import model, normalize  # noqa: E402
from rnaseq_de.simulate import simulate  # noqa: E402

counts, cd, _ = simulate(n_de=0, seed=7, dominant=(5, 8.0))  # 5 top genes x8, no other DE
counts = normalize.filter_low_counts(counts)
X, tgt = model.design_matrix(cd, ["subject", "treat"], "treat", "ctl")
dominant = set(counts.mean(axis=1).sort_values().index[-5:])

for name, sf in [("total-count (CPM-style)", normalize.library_size_factors(counts)),
                 ("median-of-ratios", normalize.median_of_ratios(counts))]:
    res = model.fit_genes(counts, sf, X, tgt)
    sig = res[res.padj < 0.05]
    other = sig[~sig.index.isin(dominant)]
    print(f"{name:26s} genes called: {len(sig):4d}; non-dominant genes falsely called: {len(other):4d} "
          f"of {len(res) - len(dominant)} ({len(other) / (len(res) - len(dominant)):.1%})")
