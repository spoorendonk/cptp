#!/usr/bin/env python3
"""Benchmark: compare normal solve vs prove-only mode (known UB cutoff + no heuristics).

Reads known optimal values from progress.csv, runs each instance in both modes,
and reports the time comparison.

Usage:
    python benchmarks/bench_ubs_cutoff.py [--time-limit 60] [--max-nodes 100] [pattern]
"""

import csv
import re
import subprocess
import sys
from pathlib import Path

CPTP_BIN = Path(__file__).parent.parent / "build" / "cptp-solve"
SPPRCLIB_DIR = Path(__file__).parent / "instances" / "spprclib"
PROGRESS_CSV = Path(__file__).parent / "progress.csv"


def load_known_optimals() -> dict[str, float]:
    """Load known optimal objectives from progress.csv (gap_pct == 0)."""
    optimals = {}
    with open(PROGRESS_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            gap = float(row["gap_pct"])
            if gap == 0:
                optimals[row["instance"]] = float(row["obj"])
    return optimals


def run_solver(instance_path: Path, time_limit: int,
               cutoff: float | None = None,
               skip_warm_start: bool = False) -> dict:
    """Run cptp-solve and parse output."""
    args = [str(CPTP_BIN), str(instance_path),
            "--time_limit", str(time_limit),
            "--output_flag", "false"]
    if cutoff is not None:
        args += ["--cutoff", str(cutoff)]
    if skip_warm_start:
        args += ["--heu_ws", "false"]

    result = {"status": "ERROR", "obj": None, "bound": None, "gap": None,
              "time": None, "nodes": 0, "cuts": 0}
    try:
        proc = subprocess.run(args, capture_output=True, text=True,
                              timeout=time_limit + 30)
        output = proc.stdout + proc.stderr
        for line in output.split("\n"):
            line = line.strip()
            m = re.match(
                r"Objective:\s*([-\d.eE+]+)\s+Bound:\s*([-\d.eE+]+)"
                r"\s+Gap:\s*([-\d.eE+]+)%\s+Time:\s*([\d.]+)s"
                r"\s+Nodes:\s*(\d+)", line)
            if m:
                result["obj"] = float(m.group(1))
                result["bound"] = float(m.group(2))
                result["gap"] = float(m.group(3))
                result["time"] = float(m.group(4))
                result["nodes"] = int(m.group(5))
                result["status"] = "Optimal" if result["gap"] < 0.01 else "Feasible"
            m2 = re.match(r"User cuts:\s*(\d+)", line)
            if m2:
                result["cuts"] = int(m2.group(1))
    except subprocess.TimeoutExpired:
        result["status"] = "TIMEOUT"
    except Exception as e:
        result["status"] = f"ERROR: {e}"
    return result


def main():
    time_limit = 60
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

    if not CPTP_BIN.exists():
        print(f"Error: solver binary not found at {CPTP_BIN}")
        print("Build first: cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)")
        sys.exit(1)

    optimals = load_known_optimals()
    instances = sorted(SPPRCLIB_DIR.glob("*.sppcc"))

    if pattern:
        instances = [p for p in instances if pattern in p.name]

    # Filter by node count
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

    # Only run instances with known optimals
    instances = [p for p in instances if p.stem in optimals]

    print(f"Benchmarking {len(instances)} instances: normal vs prove-only (cutoff + no heuristics)")
    print(f"Time limit: {time_limit}s per run")
    print()

    hdr = (f"{'Instance':<28} {'Optimal':>10} "
           f"{'Normal':>8} {'Nodes':>7} {'Cuts':>5} "
           f"{'Prove':>8} {'Nodes':>7} {'Cuts':>5} "
           f"{'Speedup':>8}")
    print(hdr)
    print("-" * len(hdr))

    total_normal = 0.0
    total_prove = 0.0
    count = 0

    for inst_path in instances:
        name = inst_path.stem
        opt = optimals[name]

        # Normal solve
        normal = run_solver(inst_path, time_limit)

        # Prove-only: feed known optimal as cutoff, disable heuristics
        prove = run_solver(inst_path, time_limit, cutoff=opt, skip_warm_start=True)

        # Format
        n_time = f"{normal['time']:.2f}s" if normal["time"] else "N/A"
        n_nodes = str(normal["nodes"])
        n_cuts = str(normal["cuts"])
        p_time = f"{prove['time']:.2f}s" if prove["time"] else "N/A"
        p_nodes = str(prove["nodes"])
        p_cuts = str(prove["cuts"])

        speedup = "N/A"
        if normal["time"] and prove["time"] and prove["time"] > 0:
            sp = normal["time"] / prove["time"]
            speedup = f"{sp:.2f}x"
            total_normal += normal["time"]
            total_prove += prove["time"]
            count += 1

        status = ""
        if prove["status"] != "Optimal":
            status = f"  [{prove['status']}]"

        print(f"{name:<28} {opt:>10.0f} "
              f"{n_time:>8} {n_nodes:>7} {n_cuts:>5} "
              f"{p_time:>8} {p_nodes:>7} {p_cuts:>5} "
              f"{speedup:>8}{status}")

    print("-" * len(hdr))
    if count > 0:
        print(f"{'Total':<28} {'':>10} "
              f"{total_normal:>7.1f}s {'':>7} {'':>5} "
              f"{total_prove:>7.1f}s {'':>7} {'':>5} "
              f"{total_normal/total_prove:>7.2f}x")
        print(f"\nAggregate speedup over {count} instances: "
              f"{(total_normal/total_prove):.2f}x")


if __name__ == "__main__":
    main()
