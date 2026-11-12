"""The multi-capacity instances (#89) regenerate byte for byte and meet their spec."""

import csv
import hashlib
import statistics
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
GENERATOR = REPO_ROOT / "benchmarks" / "gen_multicap.py"
MULTICAP = REPO_ROOT / "benchmarks" / "instances" / "multicap"
CORR_BOUND = 0.5


def _manifest() -> dict[str, str]:
    lines = (MULTICAP / "manifest.sha256").read_text().splitlines()
    return {rel: digest for digest, rel in (line.split("  ", 1) for line in lines)}


def _instances() -> list[str]:
    return sorted(rel for rel in _manifest() if "/m" in rel)


def _rows(path: Path) -> tuple[list[int], list[list[int]]]:
    """Return (capacities, demand rows per node) of a multi-capacity file."""
    lines = path.read_text().splitlines()
    caps = next(
        [int(v) for v in s.split(":", 1)[1].split()] for s in lines if s.startswith("CAPACITY")
    )
    start = lines.index("DEMAND_SECTION") + 1
    dims = int(next(s.split(":", 1)[1] for s in lines if s.startswith("DIMENSION")))
    demands = [[int(v) for v in lines[start + j].split()[1:]] for j in range(dims)]
    return caps, demands


def test_regenerates_byte_for_byte(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(GENERATOR), "--out", str(tmp_path)],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    manifest = _manifest()
    generated = {str(p.relative_to(tmp_path)) for p in tmp_path.rglob("*") if p.is_file()} - {
        "manifest.sha256"
    }
    assert generated == set(manifest)
    committed = {str(p.relative_to(MULTICAP)) for p in MULTICAP.rglob("*") if p.is_file()}
    assert committed - {"manifest.sha256"} == set(manifest)  # no stale files
    for rel, digest in manifest.items():
        assert hashlib.sha256((tmp_path / rel).read_bytes()).hexdigest() == digest, rel
        assert (MULTICAP / rel).read_bytes() == (tmp_path / rel).read_bytes(), rel
    assert (MULTICAP / "manifest.sha256").read_bytes() == (
        tmp_path / "manifest.sha256"
    ).read_bytes()


def test_covers_both_sets_and_every_m() -> None:
    rels = _instances()
    for set_name, ext in (("spprclib", ".sppcc"), ("roberti", ".vrp")):
        bases = sorted(
            p.stem for p in (REPO_ROOT / "benchmarks/instances" / set_name).glob(f"*{ext}")
        )
        for m in (2, 5, 10):
            assert [f"{set_name}/m{m}/{b}-m{m}{ext}" for b in bases] == [
                r for r in rels if r.startswith(f"{set_name}/m{m}/")
            ]


@pytest.mark.parametrize("rel", _instances())
def test_rows_positive_and_every_customer_fits(rel: str) -> None:
    caps, demands = _rows(MULTICAP / rel)
    m = int(rel.split("/")[1][1:])
    assert len(caps) == m
    assert all(q > 0 for q in caps)
    assert all(len(d) == m for d in demands)
    # Node 0 is the depot in both base sets. Row 1 is the base's own demand
    # vector, which has zero-demand customers in G-n262-k25-316; every
    # generated row (k >= 2) is positive on every customer.
    assert demands[0] == [0] * m
    customers = demands[1:]
    assert all(d[0] >= 0 for d in customers)
    assert all(v > 0 for d in customers for v in d[1:])
    assert all(d[k] <= caps[k] for d in customers for k in range(m))
    cols = [[d[k] for d in customers] for k in range(m)]
    for j in range(m):
        for k in range(j + 1, m):
            assert abs(statistics.correlation(cols[j], cols[k])) <= CORR_BOUND + 1e-9, (j, k)


def test_nested_rows_and_original_row_one() -> None:
    for rel in (r for r in _instances() if "/m10/" in r):
        set_name, _, fname = rel.split("/")
        base = fname.rsplit("-m10", 1)[0]
        ext = Path(fname).suffix
        caps10, dem10 = _rows(MULTICAP / rel)
        for m in (2, 5):
            caps, dem = _rows(MULTICAP / f"{set_name}/m{m}/{base}-m{m}{ext}")
            assert caps == caps10[:m]
            assert dem == [d[:m] for d in dem10]
        orig = (REPO_ROOT / "benchmarks/instances" / set_name / f"{base}{ext}").read_text()
        q = next(int(s.split(":", 1)[1]) for s in orig.splitlines() if s.startswith("CAPACITY"))
        assert caps10[0] == q


def test_correlations_within_bound() -> None:
    with (MULTICAP / "correlations.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == len(_instances())
    assert all(float(r["max_abs_corr"]) <= CORR_BOUND for r in rows)
