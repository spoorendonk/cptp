#!/usr/bin/env python3
"""Run PathWyse on an instance and verify against a known optimal value.

Usage:
    # Auto-lookup optimal from optimal.csv:
    python benchmarks/compare_pathwyse.py benchmarks/instances/spprclib/B-n45-k6-54.sppcc

    # Supply expected optimal explicitly:
    python benchmarks/compare_pathwyse.py \
        benchmarks/instances/spprclib/B-n45-k6-54.sppcc --expected -74278

    # Run all instances in a directory:
    python benchmarks/compare_pathwyse.py benchmarks/instances/spprclib/

    # Custom PathWyse binary:
    python benchmarks/compare_pathwyse.py --pathwyse-bin /path/to/pathwyse <instance>

    # PWBucket, with PathWyse's suggested settings for cyclic networks:
    python benchmarks/compare_pathwyse.py --algo PWBucket <instance>

    # Multiple capacities (CAPACITY : Q_1 .. Q_m) become m CAP resources:
    python benchmarks/compare_pathwyse.py benchmarks/instances/multicap/roberti/m5/

Environment:
    PATHWYSE_BIN  — override PathWyse binary location
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).parent
DEFAULT_PATHWYSE_BIN = SCRIPT_DIR / ".pathwyse" / "bin" / "pathwyse"


def _euc2d(x1: float, y1: float, x2: float, y2: float) -> int:
    """EUC_2D distance (TSPLIB convention: nint)."""
    return int(math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2) + 0.5)


def parse_tsplib(filepath: Path) -> dict[str, Any]:  # noqa: C901 - one linear pass per TSPLIB section
    """Parse a TSPLIB-format .sppcc or .vrp file.

    Handles two variants:
    - SPPRCLIB .sppcc: EDGE_WEIGHT_SECTION (full matrix) + NODE_WEIGHT_SECTION
    - Roberti .vrp: NODE_COORD_SECTION (EUC_2D) + PROFIT_SECTION
    """
    with open(filepath) as f:
        lines = f.readlines()

    # Parse into sections
    header: dict[str, str] = {}
    sections: dict[str, list[str]] = {}
    current_section = None

    section_names = {
        "NODE_COORD_SECTION",
        "DEMAND_SECTION",
        "DEPOT_SECTION",
        "EDGE_WEIGHT_SECTION",
        "NODE_WEIGHT_SECTION",
        "PROFIT_SECTION",
    }

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if line == "EOF":
            break

        # Check for section headers
        if line in section_names:
            current_section = line
            sections[current_section] = []
            continue

        # A "KEY : value" header line terminates any open section
        # Keys contain only uppercase letters, underscores, digits
        pos = line.find(":")
        if pos != -1:
            key_candidate = line[:pos].strip()
            if re.fullmatch(r"[A-Z][A-Z_0-9]*", key_candidate):
                header[key_candidate] = line[pos + 1 :].strip()
                current_section = None
                continue

        if current_section is not None:
            sections[current_section].append(line)

    n = int(header.get("DIMENSION", 0))
    # One bound per capacity resource: "CAPACITY : Q" or "CAPACITY : Q_1 .. Q_m"
    capacities = [float(x) for x in header.get("CAPACITY", "0").split()]
    m = len(capacities)
    name = header.get("NAME", filepath.stem)
    comment = header.get("COMMENT", "")

    # Build distance matrix
    dist: list[list[float]] = [[0.0] * n for _ in range(n)]

    if "EDGE_WEIGHT_SECTION" in sections:
        # Full matrix (SPPRCLIB .sppcc)
        numbers: list[float] = []
        for line in sections["EDGE_WEIGHT_SECTION"]:
            numbers.extend(float(x) for x in line.split())
        assert len(numbers) == n * n, f"Expected {n * n} edge weights, got {len(numbers)}"
        for i in range(n):
            for j in range(n):
                dist[i][j] = numbers[i * n + j]
    elif "NODE_COORD_SECTION" in sections:
        # EUC_2D coordinates
        coords: list[tuple[float, float]] = []
        for line in sections["NODE_COORD_SECTION"]:
            parts = line.split()
            coords.append((float(parts[1]), float(parts[2])))
        for i in range(n):
            for j in range(n):
                if i != j:
                    dist[i][j] = _euc2d(coords[i][0], coords[i][1], coords[j][0], coords[j][1])

    # Node weights / profits (PathWyse node_cost = -profit)
    # SPPRCLIB: NODE_WEIGHT_SECTION has costs (positive = cost, negative = profit)
    #           -> node_cost for PathWyse = node_weight (already cost)
    # Roberti:  PROFIT_SECTION has positive profits
    #           -> node_cost for PathWyse = -profit
    node_costs: list[float] = [0.0] * n

    if "NODE_WEIGHT_SECTION" in sections:
        values: list[float] = []
        for line in sections["NODE_WEIGHT_SECTION"]:
            values.extend(float(x) for x in line.split())
        for i in range(n):
            node_costs[i] = values[i]
    elif "PROFIT_SECTION" in sections:
        for line in sections["PROFIT_SECTION"]:
            parts = line.split()
            idx = int(parts[0]) - 1  # 1-indexed
            profit = float(parts[1])
            if 0 <= idx < n:
                node_costs[idx] = -profit  # PathWyse: cost = -profit

    # Demands (1-indexed in TSPLIB): "node d_1 .. d_m", one column per capacity
    demands: list[list[float]] = [[0.0] * m for _ in range(n)]
    if "DEMAND_SECTION" in sections:
        for line in sections["DEMAND_SECTION"]:
            parts = line.split()
            if len(parts) != m + 1:
                raise ValueError(
                    f"{filepath}: demand line {line!r} has {len(parts) - 1} values, "
                    f"expected {m} (one per CAPACITY)"
                )
            idx = int(parts[0]) - 1
            if 0 <= idx < n:
                demands[idx] = [float(x) for x in parts[1:]]

    # Depot (0-indexed)
    depot = 0
    if "DEPOT_SECTION" in sections:
        for line in sections["DEPOT_SECTION"]:
            d = int(line.strip())
            if d >= 1:
                depot = d - 1
                break

    # Depot demand = 0
    demands[depot] = [0.0] * m

    return {
        "name": name,
        "comment": comment,
        "dimension": n,
        "capacities": capacities,
        "dist": dist,
        "node_costs": node_costs,
        "demands": demands,
        "depot": depot,
    }


def _needs_scaling(data: dict[str, Any]) -> int:
    """Determine scale factor for PathWyse (which uses integer arithmetic).

    Returns 1000 if fractional values are present, 1 otherwise.
    """
    for v in data["node_costs"]:
        if v != int(v):
            return 1000
    for row in data["dist"]:
        for v in row:
            if v != int(v):
                return 1000
    return 1


def convert_to_pathwyse(data: dict[str, Any], output_path: Path) -> int:
    """Write PathWyse native format from parsed TSPLIB data.

    Returns the scale factor applied (1 if none). PathWyse uses int costs,
    so fractional values are scaled up. Divide PathWyse objective by this
    factor to get the real value.
    """
    n = data["dimension"]
    capacities: list[float] = data["capacities"]
    demands: list[list[float]] = data["demands"]
    m = len(capacities)
    scale = _needs_scaling(data)

    # Capacities and demands are resource values, which PathWyse reads as
    # integers and which are not scaled: refuse a fractional one.
    for v in [*capacities, *(d for row in demands for d in row)]:
        if v != int(v):
            raise ValueError(f"{data['name']}: fractional capacity or demand {v}")

    with open(output_path, "w") as f:
        depot = data.get("depot", 0)
        f.write(f"NAME : {data['name']}\n")
        f.write(f"COMMENT : {data['comment']}\n")
        f.write(f"SIZE : {n}\n")
        f.write("DIRECTED : 1\n")
        f.write("CYCLIC : 1\n")
        f.write(f"ORIGIN : {depot}\n")
        f.write(f"DESTINATION : {depot}\n")
        # One CAP resource per capacity row, each bounded by its own capacity
        f.write(f"RESOURCES : {m}\n")
        f.write(f"RES_NAMES : {' '.join(str(r) for r in range(m))}\n")
        f.write("RES_TYPE\n")
        for r in range(m):
            f.write(f"{r} CAP\n")
        f.write("END\n")
        f.write("RES_BOUND\n")
        for r in range(m):
            f.write(f"{r} 0 {int(capacities[r])}\n")
        f.write("END\n")

        # Edge costs (scaled to int)
        f.write("EDGE_COST\n")
        dist = data["dist"]
        for i in range(n):
            for j in range(n):
                f.write(f"{i} {j} {int(round(dist[i][j] * scale))}\n")
        f.write("END\n")

        # Node costs (scaled to int)
        f.write("NODE_COST\n")
        for i in range(n):
            f.write(f"{i} {int(round(data['node_costs'][i] * scale))}\n")
        f.write("END\n")

        # Node consumption (demands — resource values, not scaled)
        f.write("NODE_CONSUMPTION\n")
        for r in range(m):
            for i in range(n):
                f.write(f"{r} {i} {int(demands[i][r])}\n")
        f.write("END\n")

    return scale


def find_pathwyse_bin(override: str | None) -> Path:
    """Resolve PathWyse binary location."""
    if override:
        return Path(override)
    env = os.environ.get("PATHWYSE_BIN")
    if env:
        return Path(env)
    return DEFAULT_PATHWYSE_BIN


def load_optimal_csv(instance_path: Path) -> float | None:
    """Look up expected optimal value from optimal.csv in the instance's directory."""
    csv_path = instance_path.parent / "optimal.csv"
    if not csv_path.exists():
        return None
    stem = instance_path.stem
    with open(csv_path) as f:
        for row in csv.DictReader(filter(lambda line: not line.startswith("#"), f)):
            if row["instance"] == stem:
                return float(row["optimal"])
    return None


def pathwyse_revision(pw_bin: Path) -> str:
    """Commit of the PathWyse checkout the binary was built in, or "" if unknown.

    "" also when the binary's directory is not the top of a git checkout, so a
    binary inside another repository does not get that repository's commit.
    """
    root = pw_bin.resolve().parent.parent
    try:
        proc = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--show-toplevel", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return ""
    toplevel, commit = proc.stdout.split()
    return commit if Path(toplevel).resolve() == root else ""


# Overrides of the checkout's own pathwyse.set, per algorithm. PWDefault runs
# with pathwyse.set as shipped, as the paper's rows did; PWBucket with the
# settings PathWyse 1.1.0 suggests for cyclic networks
# (docs/documentation/SUGGESTED_SETTINGS.md).
ALGO_SETTINGS: dict[str, dict[str, str]] = {
    "PWDefault": {},
    "PWBucket": {
        "algo_names": "PWBucket",
        "algo_active": "1",
        "algo/default/autoconfig": "1",
        "algo/default/search": "bidirectional",
        "algo/default/parallel": "1",
        "algo/default/compare_unreachables": "1",
        "algo/default/ng": "standard",
        "algo/default/ng/set_size": "8",
        "algo/default/dssr": "restricted",
        "algo/default/join/type": "bucket_node",
        "algo/bucket/size_mode": "min_res_consumption",
    },
}


def needs_unit_buckets(algo: str, data: dict[str, Any] | None) -> bool:
    """True if PWBucket's min_res_consumption bucket size would be 0.

    PWBucket buckets on resource 0 (the first capacity row) with the smallest
    demand of a customer as bucket size, and divides by it: a customer with
    zero row-0 demand (G-n262-k25-316 has 12) makes it fail. Such instances
    run with unit_size buckets instead, and their rows say so.
    """
    if algo != "PWBucket" or data is None:
        return False
    depot = data.get("depot", 0)
    return any(d[0] == 0 for i, d in enumerate(data["demands"]) if i != depot)


def write_settings(pw_bin: Path, algo: str, output_path: Path, unit_buckets: bool = False) -> None:
    """Write a PathWyse parameter file: the checkout's pathwyse.set with algo's overrides.

    Without a pathwyse.set beside the binary's checkout (a binary installed
    elsewhere), only the overrides are written and PathWyse's built-in
    defaults stand for the rest, as when it finds no parameter file.
    """
    base = pw_bin.parent.parent / "pathwyse.set"
    lines = base.read_text().splitlines() if base.is_file() else []
    overrides = dict(ALGO_SETTINGS[algo])
    if unit_buckets:
        overrides["algo/bucket/size_mode"] = "unit_size"
    out: list[str] = []
    for line in lines:
        key = line.split("=", 1)[0].strip() if "=" in line and not line.startswith("#") else ""
        if key in overrides:
            out.append(f"{key} = {overrides.pop(key)}")
        else:
            out.append(line)
    out.extend(f"{key} = {value}" for key, value in overrides.items())
    output_path.write_text("\n".join(out) + "\n")


_ANSI = re.compile(r"\x1b\[[0-9;]*m")
_KEY_VALUE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_ ]*?)\s*:\s*(.*?)\s*$")


def run_pathwyse(
    pw_bin: Path,
    instance_path: Path,
    time_limit: int,
    algo: str = "PWDefault",
    settings_path: Path | None = None,
) -> dict[str, Any]:
    """Run PathWyse on an instance and parse the result.

    Reads both output styles: "Obj: v" and "PWDefault global time: t" (before
    1.0.0), and the aligned "Obj   : v" and "PWDefault Time   : t" (1.1.0).
    """
    result: dict[str, Any] = {"status": "ERROR", "obj": None, "tour": [], "time": None}
    cmd = [str(pw_bin), str(instance_path)]
    if settings_path is not None:
        cmd += ["-param", str(settings_path)]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=time_limit + 10,
            cwd=str(pw_bin.parent.parent),
        )
        output = _ANSI.sub("", proc.stdout + proc.stderr)
        for line in output.split("\n"):
            kv = _KEY_VALUE.match(line)
            if kv is None:
                continue
            key, value = kv.group(1), kv.group(2)
            if key == "Obj":
                result["obj"] = float(value)
            elif key == "Tour":
                result["tour"] = list(map(int, value.split()))
            elif key == "Solution Status":
                result["status"] = value
            elif key in (f"{algo} Time", f"{algo} global time"):
                result["time"] = float(value)
    except subprocess.TimeoutExpired:
        result["status"] = "TIMEOUT"
    except Exception as e:
        result["status"] = f"ERROR: {e}"
    return result


def instance_set(instance_path: Path) -> str:
    """Instance set from the path below instances/: e.g. "spprclib" or "multicap/roberti/m5"."""
    parts = instance_path.resolve().parent.parts
    if "instances" in parts:
        below = parts[len(parts) - parts[::-1].index("instances") :]
        if below:
            return "/".join(below)
    text = str(instance_path)
    for known in ("roberti", "spprclib"):
        if known in text:
            return known
    return "unknown"


def convert_instance(instance_path: Path, tmpdir: Path) -> tuple[Path, int, dict[str, Any] | None]:
    """Convert TSPLIB instance to PathWyse format if needed.

    Returns (path_to_use, scale_factor, parsed data or None for a native file).
    """
    suffix = instance_path.suffix.lower()
    if suffix in (".sppcc", ".vrp"):
        data = parse_tsplib(instance_path)
        pw_path = tmpdir / (instance_path.stem + ".pw")
        scale = convert_to_pathwyse(data, pw_path)
        return pw_path, scale, data
    return instance_path, 1, None


def print_result(result: dict[str, Any], expected: float | None) -> None:
    """Print one instance's result, as single-instance mode reports it."""
    print(f"\nInstance: {result['instance']}")
    obj_str = f"{result['obj']:.6g}" if result["obj"] is not None else "N/A"
    time_str = f"{result['time']:.2f}s" if result["time"] else "N/A"
    print(f"  PathWyse obj:    {obj_str}")
    print(f"  PathWyse time:   {time_str}")
    print(f"  PathWyse status: {result['status']}")
    if result["note"]:
        print(f"  Note:            {result['note']}")
    if result["tour"]:
        tour_str = " -> ".join(str(n) for n in result["tour"][:10])
        if len(result["tour"]) > 10:
            tour_str += " -> ..."
        print(f"  Tour ({len(result['tour'])} nodes): {tour_str}")
    if expected is not None:
        print(f"  Expected:        {expected:.6g}")
        if result["obj"] is not None:
            diff = abs(result["obj"] - expected)
            if diff < 0.01:
                print(f"  MATCH (diff={diff:.6f})")
            else:
                print(f"  MISMATCH (diff={diff:.4f})")


def run_single(
    pw_bin: Path,
    instance_path: Path,
    expected: float | None,
    time_limit: int,
    verbose: bool,
    algo: str = "PWDefault",
) -> dict[str, Any]:
    """Run PathWyse on one instance, return result dict.

    An instance that cannot be converted gives an ERROR row, not an exception,
    so one bad file does not lose the rest of a batch.
    """
    tmpdir = Path(tempfile.mkdtemp(prefix="cptp_pw_"))
    note = ""
    scale = 1
    try:
        pw_instance, scale, data = convert_instance(instance_path, tmpdir)
        unit_buckets = needs_unit_buckets(algo, data)
        if unit_buckets:
            note = "unit_size buckets"
        settings_path = tmpdir / "pathwyse.set"
        write_settings(pw_bin, algo, settings_path, unit_buckets=unit_buckets)
        result = run_pathwyse(
            pw_bin, pw_instance, time_limit=time_limit, algo=algo, settings_path=settings_path
        )
    except (OSError, ValueError, AssertionError, IndexError) as e:
        result = {"status": f"ERROR: {e}", "obj": None, "tour": [], "time": None}
    finally:
        import shutil

        shutil.rmtree(tmpdir, ignore_errors=True)

    # Unscale the objective
    if result["obj"] is not None and scale != 1:
        result["obj"] = result["obj"] / scale

    name = instance_path.stem
    result["instance"] = name
    result["set"] = instance_set(instance_path)
    result["note"] = note

    if expected is None:
        expected = load_optimal_csv(instance_path)

    if verbose:
        print_result(result, expected)

    result["expected"] = expected
    return result


CSV_HEADER = [
    "instance",
    "set",
    "algo",
    "obj",
    "time_s",
    "status",
    "timestamp",
    "machine",
    "pathwyse",
    "note",
]


def csv_needs_header(csv_path: Path) -> bool:
    """True if csv_path has no data lines yet; exits if its header is not CSV_HEADER.

    benchmarks/pathwyse.csv (the paper's rows) has the older header without
    the algo and pathwyse columns, so new rows go to a file of their own.
    """
    if not csv_path.exists():
        return True
    with open(csv_path) as check:
        for line in check:
            if line.strip() and not line.startswith("#"):
                header = next(csv.reader([line]))
                if header != CSV_HEADER:
                    print(
                        f"Error: {csv_path} has columns {header}, expected {CSV_HEADER}; "
                        "write to another file",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                return False
    return True


def main() -> None:  # noqa: C901 - argument handling and the report loop
    parser = argparse.ArgumentParser(
        description="Run PathWyse on instance(s) and verify against known optimal values."
    )
    parser.add_argument("instance", help="Instance file or directory of instances")
    parser.add_argument(
        "--expected", type=float, help="Expected optimal value (overrides optimal.csv lookup)"
    )
    parser.add_argument("--pathwyse-bin", help="Path to PathWyse binary")
    parser.add_argument(
        "--time-limit",
        type=int,
        default=120,
        help="Time limit per instance in seconds (default: 120)",
    )
    parser.add_argument(
        "--algo",
        choices=sorted(ALGO_SETTINGS),
        default="PWDefault",
        help="PathWyse algorithm (default: PWDefault, with the checkout's pathwyse.set)",
    )
    parser.add_argument("--csv", help="Append results to CSV file")
    parser.add_argument("--machine", help="Machine description (e.g. 'AMD Ryzen 9 3950X 32T')")
    args = parser.parse_args()

    pw_bin = find_pathwyse_bin(args.pathwyse_bin)
    if not pw_bin.exists():
        print(f"Error: PathWyse binary not found at {pw_bin}", file=sys.stderr)
        print("Run ./benchmarks/setup_pathwyse.sh first", file=sys.stderr)
        sys.exit(1)

    write_header = csv_needs_header(Path(args.csv)) if args.csv else False

    target = Path(args.instance)
    if target.is_dir():
        instances = sorted(target.glob("*.sppcc")) + sorted(target.glob("*.vrp"))
        if not instances:
            print(f"No .sppcc or .vrp files found in {target}", file=sys.stderr)
            sys.exit(1)
    elif target.is_file():
        instances = [target]
    else:
        print(f"Error: {target} not found", file=sys.stderr)
        sys.exit(1)

    # Table header for batch mode
    batch = len(instances) > 1
    if batch:
        print(
            f"Running PathWyse {args.algo} on {len(instances)} instances "
            f"(time limit: {args.time_limit}s)"
        )
        print(f"{'Instance':<30} {'PW Obj':>12} {'PW Time':>10} {'Expected':>12} {'Match':>7}")
        print("-" * 73)

    results = []
    matches = mismatches = errors = unchecked = 0

    for inst_path in instances:
        expected = args.expected  # only meaningful for single instance
        if batch:
            expected = None  # always use optimal.csv in batch mode

        r = run_single(
            pw_bin, inst_path, expected, args.time_limit, verbose=not batch, algo=args.algo
        )
        results.append(r)

        if batch:
            name = inst_path.stem
            obj_str = f"{r['obj']:.6g}" if r["obj"] is not None else "N/A"
            time_str = f"{r['time']:.2f}s" if r["time"] else "N/A"
            exp_str = f"{r['expected']:.6g}" if r["expected"] is not None else "N/A"

            match_str = "?"
            if r["obj"] is not None and r["expected"] is not None:
                if abs(r["obj"] - r["expected"]) < 0.01:
                    match_str = "OK"
                    matches += 1
                else:
                    match_str = "FAIL"
                    mismatches += 1
            elif r["obj"] is not None:
                match_str = "-"  # solved, no known optimum (e.g. multicap)
                unchecked += 1
            elif r["status"] == "TIMEOUT":
                match_str = "TLIM"
                errors += 1
            else:
                match_str = "ERR"
                errors += 1

            print(f"{name:<30} {obj_str:>12} {time_str:>10} {exp_str:>12} {match_str:>7}")

    if batch:
        print("-" * 73)
        print(
            f"{matches} match, {mismatches} mismatch, {unchecked} solved without a known "
            f"optimum, {errors} timeout/error out of {len(instances)}"
        )

    # Write CSV if requested
    if args.csv:
        timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
        machine = args.machine or ""
        revision = pathwyse_revision(pw_bin)
        with open(args.csv, "a", newline="") as f:
            writer = csv.writer(f)
            if write_header:
                writer.writerow(CSV_HEADER)
            for r in results:
                writer.writerow(
                    [
                        r["instance"],
                        r["set"],
                        args.algo,
                        f"{r['obj']:.12g}" if r["obj"] is not None else "",
                        f"{r['time']:.3f}" if r["time"] is not None else "",
                        r["status"],
                        timestamp,
                        machine,
                        revision,
                        r["note"],
                    ]
                )

    if mismatches > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
