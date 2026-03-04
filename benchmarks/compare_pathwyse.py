#!/usr/bin/env python3
"""Run PathWyse on an instance and verify against a known optimal value.

Usage:
    # Auto-lookup optimal from optimal.csv:
    python benchmarks/compare_pathwyse.py benchmarks/instances/spprclib/B-n45-k6-54.sppcc

    # Supply expected optimal explicitly:
    python benchmarks/compare_pathwyse.py benchmarks/instances/spprclib/B-n45-k6-54.sppcc --expected -74278

    # Run all instances in a directory:
    python benchmarks/compare_pathwyse.py benchmarks/instances/spprclib/

    # Custom PathWyse binary:
    python benchmarks/compare_pathwyse.py --pathwyse-bin /path/to/pathwyse <instance>

Environment:
    PATHWYSE_BIN  — override PathWyse binary location
"""

import argparse
import csv
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
DEFAULT_PATHWYSE_BIN = SCRIPT_DIR / ".pathwyse" / "bin" / "pathwyse"


def _euc2d(x1: float, y1: float, x2: float, y2: float) -> int:
    """EUC_2D distance (TSPLIB convention: nint)."""
    return int(math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2) + 0.5)


def parse_tsplib(filepath: Path) -> dict:
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
        "NODE_COORD_SECTION", "DEMAND_SECTION", "DEPOT_SECTION",
        "EDGE_WEIGHT_SECTION", "NODE_WEIGHT_SECTION", "PROFIT_SECTION",
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
                header[key_candidate] = line[pos + 1:].strip()
                current_section = None
                continue

        if current_section is not None:
            sections[current_section].append(line)

    n = int(header.get("DIMENSION", 0))
    capacity = float(header.get("CAPACITY", 0))
    name = header.get("NAME", filepath.stem)
    comment = header.get("COMMENT", "")

    # Build distance matrix
    dist: list[list[float]] = [[0.0] * n for _ in range(n)]

    if "EDGE_WEIGHT_SECTION" in sections:
        # Full matrix (SPPRCLIB .sppcc)
        numbers = []
        for line in sections["EDGE_WEIGHT_SECTION"]:
            numbers.extend(float(x) for x in line.split())
        assert len(numbers) == n * n, f"Expected {n*n} edge weights, got {len(numbers)}"
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
                    dist[i][j] = _euc2d(coords[i][0], coords[i][1],
                                         coords[j][0], coords[j][1])

    # Node weights / profits (PathWyse node_cost = -profit)
    # SPPRCLIB: NODE_WEIGHT_SECTION has costs (positive = cost, negative = profit)
    #           -> node_cost for PathWyse = node_weight (already cost)
    # Roberti:  PROFIT_SECTION has positive profits
    #           -> node_cost for PathWyse = -profit
    node_costs: list[float] = [0.0] * n

    if "NODE_WEIGHT_SECTION" in sections:
        values = []
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

    # Demands (1-indexed in TSPLIB)
    demands: list[float] = [0.0] * n
    if "DEMAND_SECTION" in sections:
        for line in sections["DEMAND_SECTION"]:
            parts = line.split()
            if len(parts) >= 2:
                idx = int(parts[0]) - 1
                if 0 <= idx < n:
                    demands[idx] = float(parts[1])

    # Depot (0-indexed)
    depot = 0
    if "DEPOT_SECTION" in sections:
        for line in sections["DEPOT_SECTION"]:
            d = int(line.strip())
            if d >= 1:
                depot = d - 1
                break

    # Depot demand = 0
    demands[depot] = 0.0

    return {
        "name": name,
        "comment": comment,
        "dimension": n,
        "capacity": capacity,
        "dist": dist,
        "node_costs": node_costs,
        "demands": demands,
        "depot": depot,
    }


def _needs_scaling(data: dict) -> int:
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


def convert_to_pathwyse(data: dict, output_path: Path) -> int:
    """Write PathWyse native format from parsed TSPLIB data.

    Returns the scale factor applied (1 if none). PathWyse uses int costs,
    so fractional values are scaled up. Divide PathWyse objective by this
    factor to get the real value.
    """
    n = data["dimension"]
    scale = _needs_scaling(data)

    with open(output_path, "w") as f:
        depot = data.get("depot", 0)
        f.write(f"NAME : {data['name']}\n")
        f.write(f"COMMENT : {data['comment']}\n")
        f.write(f"SIZE : {n}\n")
        f.write("DIRECTED : 1\n")
        f.write("CYCLIC : 1\n")
        f.write(f"ORIGIN : {depot}\n")
        f.write(f"DESTINATION : {depot}\n")
        f.write("RESOURCES : 1\n")
        f.write("RES_NAMES : 0\n")
        f.write("RES_TYPE\n")
        f.write("0 CAP\n")
        f.write("END\n")
        f.write("RES_BOUND\n")
        # Capacity is a resource bound, not scaled
        f.write(f"0 0 {int(data['capacity'])}\n")
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
        for i in range(n):
            f.write(f"0 {i} {int(data['demands'][i])}\n")
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


def run_pathwyse(pw_bin: Path, instance_path: Path, time_limit: int) -> dict:
    """Run PathWyse on an instance and parse the result."""
    result = {"status": "ERROR", "obj": None, "tour": [], "time": None}
    try:
        proc = subprocess.run(
            [str(pw_bin), str(instance_path)],
            capture_output=True, text=True, timeout=time_limit + 10,
            cwd=str(pw_bin.parent.parent),
        )
        output = proc.stdout + proc.stderr
        for line in output.split("\n"):
            line = line.strip()
            if line.startswith("Obj:"):
                result["obj"] = float(line.split(":")[1].strip())
            if line.startswith("Tour:"):
                result["tour"] = list(map(int, line.split(":")[1].strip().split()))
            if line.startswith("Solution Status:"):
                result["status"] = line.split(":")[1].strip()
            m = re.match(r"PWDefault global time:\s*([\d.]+)", line)
            if m:
                result["time"] = float(m.group(1))
    except subprocess.TimeoutExpired:
        result["status"] = "TIMEOUT"
    except Exception as e:
        result["status"] = f"ERROR: {e}"
    return result


def convert_instance(instance_path: Path, tmpdir: Path) -> tuple[Path, int]:
    """Convert TSPLIB instance to PathWyse format if needed.

    Returns (path_to_use, scale_factor).
    """
    suffix = instance_path.suffix.lower()
    if suffix in (".sppcc", ".vrp"):
        data = parse_tsplib(instance_path)
        pw_path = tmpdir / (instance_path.stem + ".pw")
        scale = convert_to_pathwyse(data, pw_path)
        return pw_path, scale
    return instance_path, 1


def run_single(pw_bin: Path, instance_path: Path, expected: float | None,
               time_limit: int, verbose: bool) -> dict:
    """Run PathWyse on one instance, return result dict."""
    tmpdir = Path(tempfile.mkdtemp(prefix="cptp_pw_"))
    try:
        pw_instance, scale = convert_instance(instance_path, tmpdir)
        result = run_pathwyse(pw_bin, pw_instance, time_limit=time_limit)
    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)

    # Unscale the objective
    if result["obj"] is not None and scale != 1:
        result["obj"] = result["obj"] / scale

    name = instance_path.stem
    result["instance"] = name

    if expected is None:
        expected = load_optimal_csv(instance_path)

    if verbose:
        print(f"\nInstance: {name}")
        obj_str = f"{result['obj']:.6g}" if result["obj"] is not None else "N/A"
        time_str = f"{result['time']:.2f}s" if result["time"] else "N/A"
        print(f"  PathWyse obj:    {obj_str}")
        print(f"  PathWyse time:   {time_str}")
        print(f"  PathWyse status: {result['status']}")
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

    result["expected"] = expected
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run PathWyse on instance(s) and verify against known optimal values."
    )
    parser.add_argument("instance", help="Instance file or directory of instances")
    parser.add_argument("--expected", type=float, help="Expected optimal value (overrides optimal.csv lookup)")
    parser.add_argument("--pathwyse-bin", help="Path to PathWyse binary")
    parser.add_argument("--time-limit", type=int, default=120, help="Time limit per instance in seconds (default: 120)")
    parser.add_argument("--csv", help="Append results to CSV file")
    args = parser.parse_args()

    pw_bin = find_pathwyse_bin(args.pathwyse_bin)
    if not pw_bin.exists():
        print(f"Error: PathWyse binary not found at {pw_bin}", file=sys.stderr)
        print("Run ./benchmarks/setup_pathwyse.sh first", file=sys.stderr)
        sys.exit(1)

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
        print(f"Running PathWyse on {len(instances)} instances (time limit: {args.time_limit}s)")
        print(f"{'Instance':<30} {'PW Obj':>12} {'PW Time':>10} {'Expected':>12} {'Match':>7}")
        print("-" * 73)

    results = []
    matches = mismatches = errors = 0

    for inst_path in instances:
        expected = args.expected  # only meaningful for single instance
        if batch:
            expected = None  # always use optimal.csv in batch mode

        r = run_single(pw_bin, inst_path, expected, args.time_limit, verbose=not batch)
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
            elif r["status"] == "TIMEOUT":
                match_str = "TLIM"
                errors += 1
            else:
                match_str = "ERR"
                errors += 1

            print(f"{name:<30} {obj_str:>12} {time_str:>10} {exp_str:>12} {match_str:>7}")

    if batch:
        print("-" * 73)
        print(f"{matches} match, {mismatches} mismatch, {errors} timeout/error out of {len(instances)}")

    # Write CSV if requested
    if args.csv:
        csv_path = Path(args.csv)
        write_header = not csv_path.exists()
        with open(csv_path, "a", newline="") as f:
            writer = csv.writer(f)
            if write_header:
                writer.writerow(["instance", "set", "obj", "time_s", "status"])
            for r in results:
                # Infer set from path
                inst_path = Path(args.instance)
                if "roberti" in str(inst_path):
                    inst_set = "roberti"
                elif "spprclib" in str(inst_path):
                    inst_set = "spprclib"
                else:
                    inst_set = "unknown"
                writer.writerow([
                    r["instance"],
                    inst_set,
                    f"{r['obj']:.6g}" if r["obj"] is not None else "",
                    f"{r['time']:.3f}" if r["time"] is not None else "",
                    r["status"],
                ])

    if mismatches > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
