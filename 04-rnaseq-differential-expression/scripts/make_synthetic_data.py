"""Write a synthetic paired count matrix (demo only; NOT real biology)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rnaseq_de.simulate import simulate  # noqa: E402

counts, cd, truth = simulate(n_genes=3000, n_de=300, seed=5)
Path("data").mkdir(exist_ok=True)
counts.to_csv("data/synthetic_counts.csv"); cd.to_csv("data/synthetic_coldata.csv")
print(f"wrote data/synthetic_counts.csv {counts.shape} and data/synthetic_coldata.csv; {len(truth)} true DE genes")
