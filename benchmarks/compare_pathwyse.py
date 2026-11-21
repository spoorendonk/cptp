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
import hashlib
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
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
        return Path(override).resolve()
    env = os.environ.get("PATHWYSE_BIN")
    if env:
        return Path(env).resolve()
    return DEFAULT_PATHWYSE_BIN.resolve()


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


UNKNOWN_REVISION = "unknown"


def _read_stamp(stamp: Path) -> tuple[str, dict[str, str]] | None:
    """(commit, {path: sha256}) from a provenance stamp, or None if unreadable."""
    try:
        lines = stamp.read_text().splitlines()
    except (OSError, UnicodeDecodeError):
        return None
    if not lines or not re.fullmatch(r"commit [0-9a-f]{40}", lines[0]):
        return None
    digests: dict[str, str] = {}
    for line in lines[1:]:
        # sha256sum / shasum -a 256 lines: "<hex>  <path>" ("<hex> *<path>" in binary mode)
        match = re.fullmatch(r"([0-9a-f]{64}) [ *](.+)", line)
        if match is None:
            return None
        digests[match.group(2)] = match.group(1)
    return lines[0].removeprefix("commit "), digests


def pathwyse_revision(pw_bin: Path) -> str:
    """Commit the binary was built at, as its provenance stamp proves it, or "unknown".

    setup_pathwyse.sh writes the stamp, <binary>.provenance, after building at
    the pinned commit: the commit, and the sha256 of each file a run uses
    (paths relative to the binary's directory: the binary, the shared library
    it loads, ../pathwyse.set). The commit stands only while the binary is
    listed and every listed file still hashes as recorded, so a binary rebuilt
    or a file changed since then is "unknown", whatever the checkout's HEAD is.
    """
    pw_bin = pw_bin.resolve()
    stamp = _read_stamp(pw_bin.with_name(pw_bin.name + ".provenance"))
    if stamp is None:
        return UNKNOWN_REVISION
    commit, digests = stamp
    if pw_bin.name not in digests:
        return UNKNOWN_REVISION
    for rel, digest in digests.items():
        try:
            actual = hashlib.sha256((pw_bin.parent / rel).read_bytes()).hexdigest()
        except OSError:
            return UNKNOWN_REVISION
        if actual != digest:
            return UNKNOWN_REVISION
    return commit


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


def write_settings(pw_bin: Path, algo: str, output_path: Path, unit_buckets: bool = False) -> bool:
    """Write a PathWyse parameter file: the checkout's pathwyse.set with algo's overrides.

    Without a pathwyse.set beside the binary's checkout (a binary installed
    elsewhere), only the overrides are written and PathWyse's built-in
    defaults stand for the rest, as when it finds no parameter file.
    """
    base = pw_bin.parent.parent / "pathwyse.set"
    found = base.is_file()
    lines = base.read_text().splitlines() if found else []
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
    return found


_ANSI = re.compile(r"\x1b\[[0-9;]*m")
_KEY_VALUE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_ ]*?)\s*:\s*(.*?)\s*$")

# The only PathWyse solution status that proves optimality (src/core/data/path.cpp);
# "Feasible", "Super optimal", "NG-route" and the others are unproven.
PROVEN_STATUS = "Optimal"


def _no_result(status: str) -> dict[str, Any]:
    return {"status": status, "obj": None, "tour": [], "time": None}


def parse_output(output: str, algo: str) -> dict[str, Any]:
    """Result fields from PathWyse's output; status "ERROR" if it names no status.

    Reads both output styles: "Obj: v" and "PWDefault global time: t" (before
    1.0.0), and the aligned "Obj   : v" and "PWDefault Time   : t" (1.1.0).
    """
    result = _no_result("ERROR")
    for line in _ANSI.sub("", output).split("\n"):
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
    return result


def run_pathwyse(
    pw_bin: Path,
    instance_path: Path,
    time_limit: int,
    algo: str = "PWDefault",
    settings_path: Path | None = None,
) -> dict[str, Any]:
    """Run PathWyse on an instance and parse the result.

    A run that exits nonzero, or is killed by a signal, is an ERROR whatever it
    printed: its objective and time are not kept, so they cannot be counted.
    """
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
    except subprocess.TimeoutExpired:
        return _no_result("TIMEOUT")
    except (OSError, ValueError) as e:
        return _no_result(f"ERROR: {e}")
    if proc.returncode < 0:
        return _no_result(f"ERROR: signal {-proc.returncode}")
    if proc.returncode != 0:
        return _no_result(f"ERROR: exit {proc.returncode}")
    try:
        return parse_output(proc.stdout + proc.stderr, algo)
    except ValueError as e:  # an unreadable Obj, Tour or time value
        return _no_result(f"ERROR: {e}")


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


def classify(result: dict[str, Any]) -> str:
    """One run's verdict, the same in single-instance and batch mode.

    "OK" / "FAIL": proven optimal, objective equal to / different from the
    expected optimum. "-": proven optimal, no known optimum. "TLIM": timed out.
    "ERR": failed (nonzero exit, unreadable output). "UNPR": finished without
    proving optimality. Only "OK" and "-" count as solved, and only "FAIL" as a
    mismatch.
    """
    status: str = result["status"]
    if status == "TIMEOUT":
        return "TLIM"
    if status.startswith("ERROR"):
        return "ERR"
    if status != PROVEN_STATUS:
        return "UNPR"
    if result["obj"] is None:
        return "ERR"
    if result["expected"] is None:
        return "-"
    return "OK" if abs(result["obj"] - result["expected"]) < 0.01 else "FAIL"


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
        verdict = classify(result)
        if verdict == "OK":
            print(f"  MATCH (diff={abs(result['obj'] - expected):.6f})")
        elif verdict == "FAIL":
            print(f"  MISMATCH (diff={abs(result['obj'] - expected):.4f})")
        else:
            print(f"  NOT COMPARED ({verdict}: not proven optimal)")


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
        notes = ["unit_size buckets"] if unit_buckets else []
        settings_path = tmpdir / "pathwyse.set"
        if not write_settings(pw_bin, algo, settings_path, unit_buckets=unit_buckets):
            notes.append("no pathwyse.set: PathWyse defaults")
        note = "; ".join(notes)
        result = run_pathwyse(
            pw_bin, pw_instance, time_limit=time_limit, algo=algo, settings_path=settings_path
        )
    except (OSError, ValueError, AssertionError, IndexError, KeyError, OverflowError) as e:
        result = _no_result(f"ERROR: {e}")
    finally:
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
    result["expected"] = expected

    if verbose:
        print_result(result, expected)

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


def _is_data(line: str) -> bool:
    """True for a CSV line, header or row; False for a blank or comment line."""
    return bool(line.strip()) and not line.startswith("#")


def check_csv_header(csv_path: Path) -> None:
    """Exit (2), before anything runs, if csv_path has a header other than CSV_HEADER.

    benchmarks/pathwyse.csv (the paper's rows) has the older header without
    the algo and pathwyse columns, so new rows go to a file of their own.
    """
    if not csv_path.exists():
        return
    with open(csv_path) as check:
        for line in check:
            if _is_data(line):
                header = next(csv.reader([line]))
                if header != CSV_HEADER:
                    print(
                        f"Error: {csv_path} has columns {header}, expected {CSV_HEADER}; "
                        "write to another file",
                        file=sys.stderr,
                    )
                    sys.exit(2)
                return


def _kept_lines(csv_path: Path, replaced: set[tuple[str, str]]) -> tuple[list[str], list[str]]:
    """(comment lines above the header, lines below it) of csv_path, minus replaced rows.

    A row is replaced when its (instance, algo) is in replaced.
    """
    head: list[str] = []
    body: list[str] = []
    if not csv_path.exists():
        return head, body
    seen_header = False
    with open(csv_path, newline="") as f:
        for raw in f:
            line = raw if raw.endswith("\n") else raw + "\n"
            if not seen_header:
                seen_header = _is_data(line)
                if not seen_header:
                    head.append(line)
                continue
            if _is_data(line):
                fields = next(csv.reader([line]))
                if len(fields) > 2 and (fields[0], fields[2]) in replaced:
                    continue
            body.append(line)
    return head, body


def write_rows(csv_path: Path, rows: list[list[str]]) -> None:
    """Write rows (CSV_HEADER's columns) into csv_path, replacing earlier runs.

    A row replaces the file's rows of the same instance and algorithm; every
    other row, and every comment line, is kept as it was. The new file is
    written beside the old one and then moved over it, so a run that fails or
    is interrupted before this point, or while writing, leaves the existing
    rows untouched.
    """
    head, body = _kept_lines(csv_path, {(row[0], row[2]) for row in rows})
    tmp = csv_path.with_name(f".{csv_path.name}.{os.getpid()}.tmp")
    try:
        with open(tmp, "x", newline="") as f:
            f.writelines(head)
            writer = csv.writer(f)
            writer.writerow(CSV_HEADER)
            f.writelines(body)
            writer.writerows(rows)
        if csv_path.exists():
            shutil.copymode(csv_path, tmp)
        os.replace(tmp, csv_path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


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
    parser.add_argument(
        "--csv",
        help="Write results to this CSV file, replacing its rows of the same instances "
        "and algorithm (other rows are kept)",
    )
    parser.add_argument("--machine", help="Machine description (e.g. 'AMD Ryzen 9 3950X 32T')")
    parser.add_argument(
        "--require-pathwyse",
        metavar="COMMIT",
        help="Refuse to run unless the binary's provenance stamp proves it was built at "
        "COMMIT (run_pathwyse_reference.sh passes setup_pathwyse.sh's pin)",
    )
    args = parser.parse_args()

    # Every refusal before a run exits 2, so a caller can tell it from a mismatch (1).
    pw_bin = find_pathwyse_bin(args.pathwyse_bin)
    if not pw_bin.exists():
        print(f"Error: PathWyse binary not found at {pw_bin}", file=sys.stderr)
        print("Run ./benchmarks/setup_pathwyse.sh first", file=sys.stderr)
        sys.exit(2)

    if args.csv:
        check_csv_header(Path(args.csv))
    revision = pathwyse_revision(pw_bin)
    if args.require_pathwyse and revision != args.require_pathwyse:
        print(
            f"Error: {pw_bin} is PathWyse {revision}, not the pinned {args.require_pathwyse}: "
            "no provenance stamp proves it. Run ./benchmarks/setup_pathwyse.sh",
            file=sys.stderr,
        )
        sys.exit(2)

    instances = find_instances(Path(args.instance))

    # Table header for batch mode
    batch = len(instances) > 1
    if batch:
        print(
            f"Running PathWyse {args.algo} ({revision}) on {len(instances)} instances "
            f"(time limit: {args.time_limit}s)"
        )
        print(f"{'Instance':<30} {'PW Obj':>12} {'PW Time':>10} {'Expected':>12} {'Match':>7}")
        print("-" * 73)

    results = []
    verdicts: Counter[str] = Counter()

    for inst_path in instances:
        # --expected is only meaningful for a single instance; a batch uses optimal.csv
        expected = None if batch else args.expected

        r = run_single(
            pw_bin, inst_path, expected, args.time_limit, verbose=not batch, algo=args.algo
        )
        results.append(r)
        verdict = classify(r)
        verdicts[verdict] += 1

        if batch:
            obj_str = f"{r['obj']:.6g}" if r["obj"] is not None else "N/A"
            time_str = f"{r['time']:.2f}s" if r["time"] else "N/A"
            exp_str = f"{r['expected']:.6g}" if r["expected"] is not None else "N/A"
            print(f"{inst_path.stem:<30} {obj_str:>12} {time_str:>10} {exp_str:>12} {verdict:>7}")

    if batch:
        print("-" * 73)
        print(
            f"{verdicts['OK']} match, {verdicts['FAIL']} mismatch, {verdicts['-']} solved "
            f"without a known optimum, {verdicts['UNPR']} not proven optimal, "
            f"{verdicts['TLIM'] + verdicts['ERR']} timeout/error out of {len(instances)}"
        )

    if args.csv:
        timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
        machine = args.machine or ""
        write_rows(
            Path(args.csv), [csv_row(r, args.algo, timestamp, machine, revision) for r in results]
        )

    if verdicts["FAIL"] > 0:
        sys.exit(1)


def find_instances(target: Path) -> list[Path]:
    """The instance file, or a directory's .sppcc and .vrp files; exits (2) if none."""
    if target.is_dir():
        instances = sorted(target.glob("*.sppcc")) + sorted(target.glob("*.vrp"))
        if not instances:
            print(f"No .sppcc or .vrp files found in {target}", file=sys.stderr)
            sys.exit(2)
        return instances
    if target.is_file():
        return [target]
    print(f"Error: {target} not found", file=sys.stderr)
    sys.exit(2)


def csv_row(r: dict[str, Any], algo: str, timestamp: str, machine: str, revision: str) -> list[str]:
    """One result as a CSV_HEADER row."""
    return [
        r["instance"],
        r["set"],
        algo,
        f"{r['obj']:.12g}" if r["obj"] is not None else "",
        f"{r['time']:.3f}" if r["time"] is not None else "",
        r["status"],
        timestamp,
        machine,
        revision,
        r["note"],
    ]


if __name__ == "__main__":
    main()
