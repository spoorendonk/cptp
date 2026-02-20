#!/usr/bin/env python3
"""Convert SPPRCLIB instances from TSPLIB to PathWyse format, run both solvers,
and compare optimal objectives."""

import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

SPPRCLIB_DIR = Path(__file__).parent / "instances" / "spprclib"
PATHWYSE_BIN = Path.home() / "code" / "pathwyse" / "bin" / "pathwyse"
CPTP_BIN = Path(__file__).parent.parent / "build" / "cptp-solve"
PATHWYSE_DIR = Path.home() / "code" / "pathwyse"


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
    # Read all numbers until NODE_WEIGHT_SECTION
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


def run_pathwyse(instance_path: Path, time_limit: int = 120) -> dict:
    """Run PathWyse on an instance and parse the result."""
    result = {"status": "ERROR", "obj": None, "tour": [], "time": None}
    try:
        proc = subprocess.run(
            [str(PATHWYSE_BIN), str(instance_path)],
            capture_output=True, text=True, timeout=time_limit + 10,
            cwd=str(PATHWYSE_DIR),
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


def run_cptp(instance_path: Path, time_limit: int = 120) -> dict:
    """Run our CPTP solver on an instance and parse the result."""
    result = {"status": "ERROR", "obj": None, "tour": [], "time": None, "cuts": 0}
    try:
        proc = subprocess.run(
            [str(CPTP_BIN), str(instance_path),
             "--time_limit", str(time_limit), "--output_flag", "false"],
            capture_output=True, text=True, timeout=time_limit + 30,
        )
        output = proc.stdout + proc.stderr
        for line in output.split("\n"):
            line = line.strip()
            m = re.match(r"Objective:\s*([-\d.eE+]+)\s+Bound:\s*([-\d.eE+]+)\s+Gap:\s*([-\d.eE+]+)%\s+Time:\s*([\d.]+)s\s+Nodes:\s*(\d+)", line)
            if m:
                result["obj"] = float(m.group(1))
                result["bound"] = float(m.group(2))
                result["gap"] = float(m.group(3))
                result["time"] = float(m.group(4))
                result["nodes"] = int(m.group(5))
                result["status"] = "Optimal" if result["gap"] < 0.01 else "Feasible"
            if line.startswith("Solution:"):
                parts = line.replace("Solution:", "").strip().split(" -> ")
                result["tour"] = [int(x) for x in parts if x.strip()]
            m2 = re.match(r"User cuts:\s*(\d+)", line)
            if m2:
                result["cuts"] = int(m2.group(1))
    except subprocess.TimeoutExpired:
        result["status"] = "TIMEOUT"
    except Exception as e:
        result["status"] = f"ERROR: {e}"
    return result


def main():
    instances = sorted(SPPRCLIB_DIR.glob("*.sppcc"))
    if not instances:
        print(f"No .sppcc instances found in {SPPRCLIB_DIR}")
        sys.exit(1)

    # Parse arguments
    time_limit = 120
    max_nodes = 9999
    pattern = None
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == "--time-limit" and i + 1 < len(sys.argv):
            time_limit = int(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == "--max-nodes" and i + 1 < len(sys.argv):
            max_nodes = int(sys.argv[i + 1])
            i += 2
        else:
            pattern = sys.argv[i]
            i += 1
    if pattern:
        instances = [p for p in instances if pattern in p.name]
    if max_nodes < 9999:
        filtered = []
        for p in instances:
            with open(p) as fh:
                for line in fh:
                    m = re.match(r"DIMENSION\s*:\s*(\d+)", line.strip())
                    if m:
                        if int(m.group(1)) <= max_nodes:
                            filtered.append(p)
                        break
        instances = filtered

    print(f"Comparing {len(instances)} SPPRCLIB instances: our solver vs PathWyse")
    print(f"{'Instance':<30} {'PW Obj':>10} {'PW Time':>8} {'Our Obj':>12} {'Our Time':>9} {'Gap%':>7} {'Cuts':>5} {'Match':>6}")
    print("-" * 96)

    tmpdir = Path(tempfile.mkdtemp(prefix="spprc_"))
    matches = 0
    mismatches = 0
    errors = 0

    for inst_path in instances:
        name = inst_path.stem

        # Parse and convert to PathWyse format
        try:
            data = parse_tsplib_sppcc(inst_path)
        except Exception as e:
            print(f"{name:<30} PARSE ERROR: {e}")
            errors += 1
            continue

        pw_path = tmpdir / inst_path.name
        convert_to_pathwyse(data, pw_path)

        # Run PathWyse on converted file
        pw_result = run_pathwyse(pw_path, time_limit=time_limit)

        # Run our solver on original TSPLIB file
        cptp_result = run_cptp(inst_path, time_limit=time_limit)

        # Format output
        pw_obj_str = f"{pw_result['obj']:.0f}" if pw_result["obj"] is not None else "N/A"
        pw_time_str = f"{pw_result['time']:.2f}s" if pw_result["time"] else "N/A"
        our_obj_str = f"{cptp_result['obj']:.0f}" if cptp_result["obj"] is not None else "N/A"
        our_time_str = f"{cptp_result['time']:.2f}s" if cptp_result["time"] else "N/A"
        gap_str = f"{cptp_result.get('gap', -1):.2f}" if cptp_result.get("gap") is not None else "N/A"
        cuts_str = str(cptp_result.get("cuts", 0))

        # Compare objectives:
        # PathWyse obj = sum(arc_costs) + sum(node_weights) for visited nodes (including depot)
        # Our obj      = sum(edge_costs) - sum(profits) including depot profit
        # Both now include depot weight, so objectives should match directly.
        match_str = "?"
        our_gap = cptp_result.get("gap", 100.0)
        if pw_result["obj"] is not None and cptp_result["obj"] is not None:
            expected_our = pw_result["obj"]
            if pw_result["status"] == "Optimal" and our_gap is not None and our_gap < 0.01:
                # Both claim optimal — check match
                if abs(cptp_result["obj"] - expected_our) < 1.0:
                    match_str = "OK"
                    matches += 1
                else:
                    match_str = "FAIL"
                    mismatches += 1
                    print(f"  DEBUG: depot_w={depot_weight}, expected_our={expected_our:.0f}, actual_our={cptp_result['obj']:.0f}")
            elif pw_result["status"] == "Optimal":
                # Our solver not proven optimal — check if we found a good solution
                if cptp_result["obj"] <= expected_our + 1.0:
                    match_str = "~OK"
                    matches += 1
                else:
                    match_str = "GAP"
                    errors += 1
            else:
                match_str = "TLIM"
                errors += 1
        else:
            match_str = "ERR"
            errors += 1

        print(f"{name:<30} {pw_obj_str:>10} {pw_time_str:>8} {our_obj_str:>12} {our_time_str:>9} {gap_str:>7} {cuts_str:>5} {match_str:>6}")

    print("-" * 96)
    print(f"Results: {matches} match, {mismatches} mismatch, {errors} timeout/error out of {len(instances)} instances")

    # Clean up
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)

    sys.exit(0 if mismatches == 0 else 1)


if __name__ == "__main__":
    main()
