#!/usr/bin/env python3
"""Generate multi-capacity versions (m = 2, 5, 10) of the benchmark instances.

Issue #89 (part of #86). Each base instance in ``benchmarks/instances/spprclib``
(``*.sppcc``) and ``benchmarks/instances/roberti`` (``*.vrp``) gets versions
with m capacity constraints ``sum_i d^k_i y_i <= Q_k``, k = 1..m, which a
solution must all satisfy. Graph, costs and profits are the base's.

File format (TSPLIB extension; every other line is copied byte for byte):

* ``CAPACITY : Q_1 Q_2 ... Q_m`` -- m values on the one CAPACITY line.
* ``DEMAND_SECTION`` lines are ``id d^1 d^2 ... d^m``.
* A ``MULTICAP : m=<m> seed=<seed> version=<v>`` header line after NAME.

Row 1 is the original (d, Q), unchanged, so it keeps the base's own
zero-demand customers (G-n262-k25-316 has twelve); rows k >= 2 are positive
on every customer. Today's single-capacity loader reads such a file as its
m = 1 restriction (first CAPACITY value, first demand column).

Parameters:

* ``MASTER_SEED`` below. Row k of instance b is drawn from SplitMix64 seeded
  with the first 8 bytes of SHA-256("<set>/<b>/<k>/<attempt>/<MASTER_SEED>"),
  so a row depends only on its base and on rows < k: rows are nested (m = 5 is
  m = 2 plus rows 3-5, m = 10 adds rows 6-10), and adding instances or rows
  never changes existing ones.
* Generator cycle for k = 2, 3, ...: ``GENERATORS`` (uniform, correlated,
  anti-correlated, spatial, heavy-tailed, repeat). maxd is the largest base
  customer demand; every value is an integer >= 1, depot demand stays 0.

  - uniform: U{1..maxd}.
  - correlated: round(rho d + (1 - rho) u), u ~ U{1..maxd}, rho = 3/10.
  - anti-correlated: as correlated with d replaced by maxd + 1 - d.
  - spatial: 3 random centre customers c with weights a_c ~ U{1..10};
    w_i = sum_c a_c s^2 / (s^2 + dist(c, i)^2), s = the mean depot distance;
    d_i = 1 + round((maxd - 1) w_i / max w). dist is the explicit matrix
    (SPPRCLIB) or Euclidean coordinates (Roberti).
  - heavy-tailed: with probability 17/20 U{1..max(1, maxd // 4)}, else
    U{maxd..3 maxd}.

* Tightness: Q_k = max(max_i d^k_i, round(tau_k sum_i d^k_i)) with
  tau_k = tau_1 f_k, tau_1 = Q / sum_i d_i (the original row's tightness),
  f_k ~ U{80..120} / 100: every row is about as tight as the original, drawn
  from [0.8, 1.2] tau_1. The max term keeps every customer feasible alone.
* Diversity: |Pearson corr(row k, row j)| <= ``CORR_BOUND`` for every j < k
  (the original included); otherwise the row is redrawn with attempt + 1, up
  to ``MAX_ATTEMPTS``.

Determinism: no library samplers and no libm exp/pow; integers are drawn by
rejection from SplitMix64, rounding is exact (integer or ``Fraction``), and
the only floating-point operations are + - * / and sqrt (correctly rounded
under IEEE 754), so output is identical on every platform.

Outputs, under ``benchmarks/instances/multicap/``:
``<set>/m<m>/<base>-m<m>.<ext>``, ``manifest.sha256`` (``sha256sum`` format,
paths relative to that directory) and ``correlations.csv`` (per instance and
m: the largest |corr| over all row pairs, and each row's Q_k and tau_k).

Usage: ``python benchmarks/gen_multicap.py [--out DIR]``; ``--check`` exits 1
when the files under DIR differ from what the script generates or its
``<set>/m<m>/`` directories hold files it does not generate; writing removes
such stale files there.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable

MASTER_SEED = 89
VERSION = 1
MS = (2, 5, 10)
CORR_BOUND = 0.5
MAX_ATTEMPTS = 100
RHO = Fraction(3, 10)
TAU_FACTOR_RANGE = (80, 120)  # f_k in hundredths
HEAVY_TAIL_PROB = (17, 20)  # probability of the small part
SPATIAL_CENTRES = 3

REPO = Path(__file__).resolve().parent.parent
BASE_SETS = {
    "spprclib": (REPO / "benchmarks" / "instances" / "spprclib", ".sppcc"),
    "roberti": (REPO / "benchmarks" / "instances" / "roberti", ".vrp"),
}
DEFAULT_OUT = REPO / "benchmarks" / "instances" / "multicap"

MASK64 = (1 << 64) - 1


class SplitMix64:
    """SplitMix64 (Steele, Lea, Flood 2014), written out for portability."""

    def __init__(self, seed: int) -> None:
        self.state = seed & MASK64

    def next(self) -> int:
        self.state = (self.state + 0x9E3779B97F4A7C15) & MASK64
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
        return z ^ (z >> 31)

    def randint(self, lo: int, hi: int) -> int:
        """Uniform integer in [lo, hi] by rejection (no modulo bias)."""
        span = hi - lo + 1
        if span <= 0:
            raise ValueError(f"empty range [{lo}, {hi}]")
        limit = (1 << 64) - ((1 << 64) % span)
        while True:
            v = self.next()
            if v < limit:
                return lo + v % span


def row_rng(set_name: str, base: str, k: int, attempt: int) -> SplitMix64:
    key = f"{set_name}/{base}/{k}/{attempt}/{MASTER_SEED}".encode()
    return SplitMix64(int.from_bytes(hashlib.sha256(key).digest()[:8], "big"))


def round_half_up(x: Fraction) -> int:
    return math.floor(x + Fraction(1, 2))


@dataclass
class Instance:
    set_name: str
    base: str
    ext: str
    lines: list[str]  # base file, line endings stripped
    name_line: int
    capacity_line: int
    demand_lines: list[int]  # line index of each node's demand, node order
    capacity: int
    demands: list[int]  # per node, depot 0
    depot: int
    dist: Callable[[int, int], float]

    @property
    def customers(self) -> list[int]:
        return [i for i in range(len(self.demands)) if i != self.depot]


def _is_keyword(line: str) -> bool:
    s = line.strip()
    return s == "EOF" or s.endswith("_SECTION") or ":" in s


def _matrix_dist(lines: list[str], start: int, n: int) -> Callable[[int, int], float]:
    flat: list[float] = []
    i = start
    while len(flat) < n * n:
        flat.extend(float(t) for t in lines[i].split())
        i += 1
    return lambda a, b: flat[a * n + b]


def _coord_dist(lines: list[str], start: int, n: int) -> Callable[[int, int], float]:
    xy = [tuple(float(t) for t in lines[start + j].split()[1:3]) for j in range(n)]

    def dist(a: int, b: int) -> float:
        dx = xy[a][0] - xy[b][0]
        dy = xy[a][1] - xy[b][1]
        return math.sqrt(dx * dx + dy * dy)

    return dist


def parse(set_name: str, path: Path) -> Instance:
    lines = path.read_text().splitlines()
    idx = {s.strip(): j for j, s in enumerate(lines) if s.strip().endswith("_SECTION")}
    hdr: dict[str, int] = {}
    for j, s in enumerate(lines):
        if ":" in s and not s.strip().endswith("_SECTION"):
            hdr.setdefault(s.split(":", 1)[0].strip(), j)
    n = int(lines[hdr["DIMENSION"]].split(":", 1)[1])
    demand_lines = list(range(idx["DEMAND_SECTION"] + 1, idx["DEMAND_SECTION"] + 1 + n))
    if any(_is_keyword(lines[j]) for j in demand_lines):
        raise ValueError(f"{path}: DEMAND_SECTION shorter than DIMENSION")
    demands = [int(lines[j].split()[1]) for j in demand_lines]
    depot = int(lines[idx["DEPOT_SECTION"] + 1]) - 1 if "DEPOT_SECTION" in idx else 0
    if demands[depot] != 0:
        raise ValueError(f"{path}: depot demand is not 0")
    if "EDGE_WEIGHT_SECTION" in idx:
        dist = _matrix_dist(lines, idx["EDGE_WEIGHT_SECTION"] + 1, n)
    else:
        dist = _coord_dist(lines, idx["NODE_COORD_SECTION"] + 1, n)
    return Instance(
        set_name=set_name,
        base=path.stem,
        ext=path.suffix,
        lines=lines,
        name_line=hdr["NAME"],
        capacity_line=hdr["CAPACITY"],
        demand_lines=demand_lines,
        capacity=int(lines[hdr["CAPACITY"]].split(":", 1)[1]),
        demands=demands,
        depot=depot,
        dist=dist,
    )


# --- demand generators: (instance, rng) -> demand per customer, in order ---


def gen_uniform(inst: Instance, rng: SplitMix64) -> list[int]:
    maxd = max(inst.demands)
    return [rng.randint(1, maxd) for _ in inst.customers]


def _mixed(rng: SplitMix64, base: int, maxd: int) -> int:
    v = round_half_up(RHO * base + (1 - RHO) * rng.randint(1, maxd))
    return max(1, v)


def gen_correlated(inst: Instance, rng: SplitMix64) -> list[int]:
    maxd = max(inst.demands)
    return [_mixed(rng, inst.demands[i], maxd) for i in inst.customers]


def gen_anticorrelated(inst: Instance, rng: SplitMix64) -> list[int]:
    maxd = max(inst.demands)
    return [_mixed(rng, maxd + 1 - inst.demands[i], maxd) for i in inst.customers]


def gen_spatial(inst: Instance, rng: SplitMix64) -> list[int]:
    maxd = max(inst.demands)
    cust = inst.customers
    centres = [cust[rng.randint(0, len(cust) - 1)] for _ in range(SPATIAL_CENTRES)]
    amps = [rng.randint(1, 10) for _ in centres]
    # math.fsum is correctly rounded on every Python version; the built-in
    # sum() of floats changed to compensated summation in 3.12.
    s = math.fsum(inst.dist(inst.depot, i) for i in cust) / len(cust)
    s2 = s * s
    w = [
        math.fsum(
            a * s2 / (s2 + inst.dist(c, i) * inst.dist(c, i))
            for c, a in zip(centres, amps, strict=True)
        )
        for i in cust
    ]
    wmax = max(w)
    return [1 + round_half_up(Fraction(maxd - 1) * Fraction(wi) / Fraction(wmax)) for wi in w]


def gen_heavy(inst: Instance, rng: SplitMix64) -> list[int]:
    maxd = max(inst.demands)
    small_num, small_den = HEAVY_TAIL_PROB
    out = []
    for _ in inst.customers:
        if rng.randint(1, small_den) <= small_num:
            out.append(rng.randint(1, max(1, maxd // 4)))
        else:
            out.append(rng.randint(maxd, 3 * maxd))
    return out


GENERATORS: list[tuple[str, Callable[[Instance, SplitMix64], list[int]]]] = [
    ("uniform", gen_uniform),
    ("correlated", gen_correlated),
    ("anticorrelated", gen_anticorrelated),
    ("spatial", gen_spatial),
    ("heavy", gen_heavy),
]


def pearson(a: list[int], b: list[int]) -> float:
    n = len(a)
    ma = Fraction(sum(a), n)
    mb = Fraction(sum(b), n)
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b, strict=True))
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((y - mb) ** 2 for y in b)
    if va == 0 or vb == 0:
        return 1.0  # a constant row: treat as fully correlated, forces a redraw
    return float(cov) / math.sqrt(float(va) * float(vb))


@dataclass
class Row:
    demands: list[int]  # per customer, in inst.customers order
    capacity: int
    tau: Fraction
    generator: str
    attempt: int


def make_row(inst: Instance, k: int, rows: list[Row]) -> Row:
    gname, gen = GENERATORS[(k - 2) % len(GENERATORS)]
    tau1 = Fraction(inst.capacity, sum(inst.demands))
    for attempt in range(MAX_ATTEMPTS):
        rng = row_rng(inst.set_name, inst.base, k, attempt)
        d = gen(inst, rng)
        if any(v < 1 for v in d):
            raise AssertionError(f"{inst.base} row {k}: non-positive demand")
        if all(abs(pearson(d, r.demands)) <= CORR_BOUND for r in rows):
            tau = tau1 * Fraction(rng.randint(*TAU_FACTOR_RANGE), 100)
            q = max(max(d), round_half_up(tau * sum(d)))
            return Row(d, q, tau, gname, attempt)
    raise RuntimeError(f"{inst.base} row {k}: no draw within |corr| <= {CORR_BOUND}")


def build_rows(inst: Instance, m: int) -> list[Row]:
    cust = inst.customers
    tau1 = Fraction(inst.capacity, sum(inst.demands))
    rows = [Row([inst.demands[i] for i in cust], inst.capacity, tau1, "original", 0)]
    for k in range(2, m + 1):
        rows.append(make_row(inst, k, rows))
    return rows


def render(inst: Instance, rows: list[Row]) -> str:
    m = len(rows)
    out = list(inst.lines)
    key = out[inst.capacity_line].split(":", 1)[0]
    out[inst.capacity_line] = f"{key}: " + " ".join(str(r.capacity) for r in rows)
    pos = {node: j for j, node in enumerate(inst.customers)}
    for node, li in enumerate(inst.demand_lines):
        extra = (
            ["0"] * (m - 1) if node == inst.depot else [str(r.demands[pos[node]]) for r in rows[1:]]
        )
        out[li] = out[li].rstrip() + " " + " ".join(extra)
    out.insert(inst.name_line + 1, f"MULTICAP : m={m} seed={MASTER_SEED} version={VERSION}")
    return "\n".join(out) + "\n"


def generate() -> dict[str, bytes]:
    """Return {relative path: file bytes} for every generated file."""
    files: dict[str, bytes] = {}
    csv = ["set,instance,m,max_abs_corr,capacities,taus,generators,attempts"]
    for set_name, (src, ext) in BASE_SETS.items():
        for path in sorted(src.glob(f"*{ext}")):
            inst = parse(set_name, path)
            rows = build_rows(inst, max(MS))
            for m in MS:
                sub = rows[:m]
                rel = f"{set_name}/m{m}/{inst.base}-m{m}{ext}"
                files[rel] = render(inst, sub).encode()
                corr = max(
                    abs(pearson(a.demands, b.demands))
                    for j, a in enumerate(sub)
                    for b in sub[j + 1 :]
                )
                csv.append(
                    f"{set_name},{inst.base},{m},{corr:.4f},"
                    + " ".join(str(r.capacity) for r in sub)
                    + ","
                    + " ".join(f"{float(r.tau):.4f}" for r in sub)
                    + ","
                    + " ".join(r.generator for r in sub)
                    + ","
                    + " ".join(str(r.attempt) for r in sub)
                )
    files["correlations.csv"] = ("\n".join(csv) + "\n").encode()
    manifest = "".join(
        f"{hashlib.sha256(b).hexdigest()}  {rel}\n" for rel, b in sorted(files.items())
    )
    files["manifest.sha256"] = manifest.encode()
    return files


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--check", action="store_true", help="compare instead of writing")
    args = ap.parse_args(argv)
    files = generate()
    # Stale files are looked for only in the <set>/m<m>/ directories this
    # script owns, so a wrong --out never loses unrelated files.
    owned = [args.out / s / f"m{m}" for s in BASE_SETS for m in MS]
    stale = sorted(
        str(p.relative_to(args.out))
        for d in owned
        for p in d.glob("*")
        if p.is_file() and str(p.relative_to(args.out)) not in files
    )
    if args.check:
        bad = [
            r
            for r, b in files.items()
            if not (args.out / r).is_file() or (args.out / r).read_bytes() != b
        ]
        for r in bad:
            print(f"differs: {r}", file=sys.stderr)
        for r in stale:
            print(f"not generated: {r}", file=sys.stderr)
        return 1 if bad or stale else 0
    for rel in stale:
        (args.out / rel).unlink()
    for rel, data in files.items():
        p = args.out / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)
    print(f"wrote {len(files)} files to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
