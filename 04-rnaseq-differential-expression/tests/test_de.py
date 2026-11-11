import numpy as np
import pandas as pd

from rnaseq_de import model, normalize, validate
from rnaseq_de.simulate import simulate


def run(counts, coldata, sf=None):
    c = normalize.filter_low_counts(counts)
    sf = normalize.median_of_ratios(c) if sf is None else sf
    X, tgt = model.design_matrix(coldata, ["subject", "treat"], "treat", "ctl")
    return model.fit_genes(c, sf, X, tgt), c


def test_size_factors_recover_depth_ratios():
    rng = np.random.default_rng(0)
    base = np.exp(rng.normal(6, 1, 2000))
    depth = np.array([0.5, 1.0, 2.0])
    counts = pd.DataFrame(rng.poisson(base[:, None] * depth), columns=list("abc"))
    sf = normalize.median_of_ratios(counts)
    np.testing.assert_allclose(sf / sf["b"], depth / depth[1], rtol=0.05)


def test_recall_and_fdr_on_simulated_truth():
    counts, cd, truth = simulate(seed=1)
    res, _ = run(counts, cd)
    called = set(res.index[res.padj < 0.05])
    recall = len(called & truth) / len(truth)
    fdp = len(called - truth) / max(len(called), 1)
    assert recall > 0.85
    assert fdp < 0.2


def test_null_data_gives_few_discoveries():
    counts, cd, _ = simulate(n_de=0, seed=2)
    res, _ = run(counts, cd)
    assert (res.padj < 0.05).mean() < 0.02


def test_direction_of_effect_is_correct():
    counts, cd, truth = simulate(seed=3, lfc=2.0)
    res, _ = run(counts, cd)
    sig = res[(res.padj < 0.05) & res.index.isin(truth)]
    assert len(sig) > 20
    assert (np.abs(sig.log2FC) > 1.0).mean() > 0.9


def test_known_gene_check_logic():
    res = pd.DataFrame({"log2FC": [2.0, -1.0], "padj": [0.001, 0.001]}, index=["DUSP1", "KLF15"])
    out = validate.check_known_genes(res, genes=["DUSP1", "KLF15", "PER1"]).set_index("gene")
    assert out.loc["DUSP1", "recovered"] and not out.loc["KLF15", "recovered"] and not out.loc["PER1", "recovered"]
