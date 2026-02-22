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
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
DEFAULT_PATHWYSE_BIN = SCRIPT_DIR / ".pathwyse" / "build" / "pathwyse"


def parse_tsplib_sppcc(filepath: Path) -> dict:
    """Parse a TSPLIB-format .sppcc file."""
    data = {
        "name": "",
        "comment": "",
        "dimension": 0,
        "edge_weights": [],  # NxN matrix
        "node_weights": [],  # N values
        "capacity": 0,
        "demands": [],       # N values (0-indexed)
    }

    with open(filepath) as f:
        content = f.read()

    # Parse header fields
    m = re.search(r"NAME\s*:\s*(.+)", content)
    if m: data["name"] = m.group(1).strip()
    m = re.search(r"COMMENT\s*:\s*(.+)", content)
    if m: data["comment"] = m.group(1).strip()
    m = re.search(r"DIMENSION\s*:\s*(\d+)", content)
    if m: data["dimension"] = int(m.group(1))
    m = re.search(r"CAPACITY\s*:\s*(\d+)", content)
    if m: data["capacity"] = int(m.group(1))

    n = data["dimension"]

    # Parse EDGE_WEIGHT_SECTION (full matrix)
    idx = content.index("EDGE_WEIGHT_SECTION")
    rest = content[idx + len("EDGE_WEIGHT_SECTION"):]
    end_idx = rest.index("NODE_WEIGHT_SECTION")
    numbers = list(map(int, rest[:end_idx].split()))
    assert len(numbers) == n * n, f"Expected {n*n} edge weights, got {len(numbers)}"
    data["edge_weights"] = [numbers[i*n:(i+1)*n] for i in range(n)]

    # Parse NODE_WEIGHT_SECTION
    idx = content.index("NODE_WEIGHT_SECTION")
    rest = content[idx + len("NODE_WEIGHT_SECTION"):]
    end_idx = rest.index("CAPACITY")
    node_weights = list(map(int, rest[:end_idx].split()))
    assert len(node_weights) == n, f"Expected {n} node weights, got {len(node_weights)}"
    data["node_weights"] = node_weights

    # Parse DEMAND_SECTION (1-indexed)
    idx = content.index("DEMAND_SECTION")
    rest = content[idx + len("DEMAND_SECTION"):]
    end_marker = "EOF" if "EOF" in rest else "DEPOT_SECTION"
    end_idx = rest.index(end_marker) if end_marker in rest else len(rest)
    demand_text = rest[:end_idx].strip()
    demands = [0] * n
    for line in demand_text.split("\n"):
        parts = line.strip().split()
        if len(parts) >= 2:
            node_id = int(parts[0]) - 1  # 1-indexed -> 0-indexed
            demand = int(parts[1])
            if 0 <= node_id < n:
                demands[node_id] = demand
    data["demands"] = demands

    return data


def convert_to_pathwyse(data: dict, output_path: Path):
    """Convert parsed SPPCC data to PathWyse native format."""
    n = data["dimension"]
    with open(output_path, "w") as f:
        f.write(f"NAME : {data['name']}\n")
        f.write(f"COMMENT : {data['comment']}\n")
        f.write(f"SIZE : {n}\n")
        f.write("DIRECTED : 1\n")
        f.write("CYCLIC : 1\n")
        f.write("RESOURCES : 1\n")
        f.write("RES_NAMES : 0\n")
        f.write("RES_TYPE\n")
        f.write("0 CAP\n")
        f.write("END\n")
        f.write("RES_BOUND\n")
        f.write(f"0 0 {data['capacity']}\n")
        f.write("END\n")

        # Edge costs
        f.write("EDGE_COST\n")
        for i in range(n):
            for j in range(n):
                cost = data["edge_weights"][i][j]
                f.write(f"{i} {j} {cost}\n")
        f.write("END\n")

        # Node costs
        f.write("NODE_COST\n")
        for i in range(n):
            f.write(f"{i} {data['node_weights'][i]}\n")
        f.write("END\n")

        # Node consumption (demands)
        f.write("NODE_CONSUMPTION\n")
        for i in range(n):
            f.write(f"0 {i} {data['demands'][i]}\n")
        f.write("END\n")


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


def convert_instance(instance_path: Path, tmpdir: Path) -> Path:
    """Convert TSPLIB instance to PathWyse format if needed. Returns path to use."""
    suffix = instance_path.suffix.lower()
    if suffix in (".sppcc", ".vrp"):
        data = parse_tsplib_sppcc(instance_path)
        pw_path = tmpdir / (instance_path.stem + ".pw")
        convert_to_pathwyse(data, pw_path)
        return pw_path
    return instance_path


def run_single(pw_bin: Path, instance_path: Path, expected: float | None,
               time_limit: int, verbose: bool) -> dict:
    """Run PathWyse on one instance, return result dict."""
    tmpdir = Path(tempfile.mkdtemp(prefix="rcspp_pw_"))
    try:
        pw_instance = convert_instance(instance_path, tmpdir)
        result = run_pathwyse(pw_bin, pw_instance, time_limit=time_limit)
    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)

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
                if diff < 1.0:
                    print(f"  MATCH (diff={diff:.4f})")
                else:
                    print(f"  MISMATCH (diff={diff:.2f})")

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
                if abs(r["obj"] - r["expected"]) < 1.0:
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
