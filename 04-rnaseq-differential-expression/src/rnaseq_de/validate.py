"""Biological validation against published glucocorticoid-response genes.

Himes et al. 2014 (PLoS ONE) reported DUSP1, KLF15, PER1 and TSC22D3 among the
dexamethasone-induced genes in airway smooth muscle cells; FKBP5 and ZBTB16 are
classic glucocorticoid-receptor targets. The test is directional: a correct
pipeline should call these up-regulated and significant.
"""
from __future__ import annotations

import pandas as pd

KNOWN_UP = ["DUSP1", "KLF15", "PER1", "TSC22D3", "FKBP5", "ZBTB16"]


def check_known_genes(res: pd.DataFrame, genes=KNOWN_UP, alpha: float = 0.05) -> pd.DataFrame:
    rows = []
    for g in genes:
        if g in res.index:
            r = res.loc[g]
            rows.append(dict(gene=g, log2FC=r.log2FC, padj=r.padj,
                             recovered=bool(r.log2FC > 0 and r.padj < alpha)))
        else:
            rows.append(dict(gene=g, log2FC=float("nan"), padj=float("nan"), recovered=False))
    return pd.DataFrame(rows)
