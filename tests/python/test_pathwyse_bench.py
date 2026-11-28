"""PathWyse benchmark scripts (#172): binary provenance, failed runs, row staging.

A stub shell script stands in for the PathWyse binary: it prints the file
$PW_STUB_OUT and exits with $PW_STUB_EXIT, or, with $PW_STUB_KILL_PARENT set,
kills the compare_pathwyse.py that runs it, as an interrupted run would end.
"""

import csv
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BENCH = REPO_ROOT / "benchmarks"
COMPARE = BENCH / "compare_pathwyse.py"
AGGREGATE = BENCH / "compare_pathwyse.sh"
REFERENCE = BENCH / "run_pathwyse_reference.sh"

OTHER_COMMIT = "0123456789abcdef0123456789abcdef01234567"
HEADER = "instance,set,algo,obj,time_s,status,timestamp,machine,pathwyse,note"

STUB = """#!/bin/sh
if [ -n "$PW_STUB_KILL_PARENT" ]; then kill -TERM "$PPID"; sleep 5; fi
cat "$PW_STUB_OUT"
exit "${PW_STUB_EXIT:-0}"
"""

# A 3-node SPPRCLIB-style instance; its optimal.csv says -42.
INSTANCE = """NAME : tiny
TYPE : SPPRC
DIMENSION : 3
CAPACITY : 10
EDGE_WEIGHT_SECTION
0 1 1
1 0 1
1 1 0
NODE_WEIGHT_SECTION
0 -20 -24
DEMAND_SECTION
1 0
2 1
3 1
DEPOT_SECTION
1
-1
EOF
"""


def pw_output(obj: float, status: str = "Optimal") -> str:
    """PathWyse 1.1.0's aligned output style."""
    return (
        f"Obj              : {obj}\n"
        "Tour             : 0 1 2 0\n"
        f"Solution Status  : {status}\n"
        "PWDefault Time   : 0.5\n"
    )


def pinned_commit() -> str:
    text = (BENCH / "setup_pathwyse.sh").read_text()
    match = re.search(r'^PATHWYSE_COMMIT="([0-9a-f]{40})"$', text, re.M)
    assert match is not None
    return match.group(1)


def write_stamp(checkout: Path, commit: str) -> None:
    """The stamp setup_pathwyse.sh writes: commit, bin/, then sha256sum lines relative to it."""
    bin_dir = checkout / "bin"
    files = [p.name for p in sorted(bin_dir.iterdir()) if p.name != "pathwyse.provenance"]
    lines = [f"commit {commit}", f"dir {bin_dir.resolve()}"]
    for rel in [*files, "../pathwyse.set"]:
        lines.append(f"{hashlib.sha256((bin_dir / rel).read_bytes()).hexdigest()}  {rel}")
    (bin_dir / "pathwyse.provenance").write_text("\n".join(lines) + "\n")


@pytest.fixture
def checkout(tmp_path: Path) -> Path:
    """A built PathWyse checkout: stub binary, shared library, pathwyse.set; no stamp."""
    root = tmp_path / "pathwyse"
    (root / "bin").mkdir(parents=True)
    binary = root / "bin" / "pathwyse"
    binary.write_text(STUB)
    binary.chmod(0o755)
    (root / "bin" / "libpathwyse_core.so").write_bytes(b"solver code")
    (root / "pathwyse.set").write_text("verbosity = 1\n")
    return root


@pytest.fixture
def instances(tmp_path: Path) -> Path:
    directory = tmp_path / "instances"
    directory.mkdir()
    (directory / "tiny.sppcc").write_text(INSTANCE)
    (directory / "optimal.csv").write_text("instance,optimal\ntiny,-42\ntiny2,-42\n")
    return directory


def add_second_instance(instances: Path) -> None:
    (instances / "tiny2.sppcc").write_text(INSTANCE.replace("NAME : tiny", "NAME : tiny2"))


def stub_env(tmp_path: Path, checkout: Path, output: str, exit_code: int = 0) -> dict[str, str]:
    out = tmp_path / "stub_output.txt"
    out.write_text(output)
    env = {k: v for k, v in os.environ.items() if not k.startswith("PW_STUB_")}
    env.update(
        PATHWYSE_BIN=str(checkout / "bin" / "pathwyse"),
        PW_STUB_OUT=str(out),
        PW_STUB_EXIT=str(exit_code),
    )
    return env


def run(cmd: list[str], env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=60, check=False)


def compare(target: Path, env: dict[str, str], *args: str) -> subprocess.CompletedProcess[str]:
    return run([sys.executable, str(COMPARE), str(target), "--time-limit", "5", *args], env)


def rows(csv_path: Path) -> list[dict[str, str]]:
    with open(csv_path) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


# --- provenance ---


def test_stamp_proves_the_commit(tmp_path: Path, checkout: Path, instances: Path) -> None:
    write_stamp(checkout, OTHER_COMMIT)
    out = tmp_path / "out.csv"
    env = stub_env(tmp_path, checkout, pw_output(-42))
    proc = compare(instances / "tiny.sppcc", env, "--csv", str(out))
    assert proc.returncode == 0, proc.stderr
    assert [r["pathwyse"] for r in rows(out)] == [OTHER_COMMIT]


@pytest.mark.parametrize("change", ["binary", "library", "settings", "no stamp"])
def test_changed_build_is_unknown(
    tmp_path: Path, checkout: Path, instances: Path, change: str
) -> None:
    """A stamp that no longer matches what runs proves nothing (the stale-binary case)."""
    write_stamp(checkout, OTHER_COMMIT)
    if change == "binary":
        (checkout / "bin" / "pathwyse").write_text(STUB + "# rebuilt\n")
    elif change == "library":
        (checkout / "bin" / "libpathwyse_core.so").write_bytes(b"rebuilt solver code")
    elif change == "settings":
        (checkout / "pathwyse.set").write_text("verbosity = 2\n")
    else:
        (checkout / "bin" / "pathwyse.provenance").unlink()
    out = tmp_path / "out.csv"
    env = stub_env(tmp_path, checkout, pw_output(-42))
    proc = compare(instances / "tiny.sppcc", env, "--csv", str(out))
    assert proc.returncode == 0, proc.stderr
    assert [r["pathwyse"] for r in rows(out)] == ["unknown"]


def test_binary_copied_out_of_its_build_is_unknown(
    tmp_path: Path, checkout: Path, instances: Path
) -> None:
    """The copy would still load the library of the build directory, which may have changed."""
    write_stamp(checkout, OTHER_COMMIT)
    copy = tmp_path / "copy"
    shutil.copytree(checkout, copy)
    out = tmp_path / "out.csv"
    env = stub_env(tmp_path, checkout, pw_output(-42))
    env["PATHWYSE_BIN"] = str(copy / "bin" / "pathwyse")
    proc = compare(instances / "tiny.sppcc", env, "--csv", str(out))
    assert proc.returncode == 0, proc.stderr
    assert [r["pathwyse"] for r in rows(out)] == ["unknown"]


def test_require_pathwyse_refuses_another_build(
    tmp_path: Path, checkout: Path, instances: Path
) -> None:
    write_stamp(checkout, OTHER_COMMIT)
    out = tmp_path / "out.csv"
    env = stub_env(tmp_path, checkout, pw_output(-42))
    proc = compare(instances, env, "--csv", str(out), "--require-pathwyse", pinned_commit())
    assert proc.returncode == 2
    assert "not the pinned" in proc.stderr
    assert f"proves it was built at {OTHER_COMMIT}" in proc.stderr
    assert not out.exists()


# --- failed and unproven runs ---


def test_nonzero_exit_is_an_error_row(tmp_path: Path, checkout: Path, instances: Path) -> None:
    """Output printed before a nonzero exit is not kept, so it is not counted."""
    out = tmp_path / "out.csv"
    env = stub_env(tmp_path, checkout, pw_output(-42), exit_code=2)
    proc = compare(instances / "tiny.sppcc", env, "--csv", str(out))
    assert proc.returncode == 0, proc.stderr
    assert "NOT COMPARED (ERR" in proc.stdout
    [row] = rows(out)
    assert (row["status"], row["obj"], row["time_s"]) == ("ERROR: exit 2", "", "")


def test_nonzero_exit_is_not_a_match_in_a_batch(
    tmp_path: Path, checkout: Path, instances: Path
) -> None:
    add_second_instance(instances)
    proc = compare(instances, stub_env(tmp_path, checkout, pw_output(-42), exit_code=2))
    assert proc.returncode == 0, proc.stderr
    assert "0 match, 0 mismatch" in proc.stdout
    assert "2 timeout/error out of 2" in proc.stdout


def test_unproven_status_is_not_a_match(tmp_path: Path, checkout: Path, instances: Path) -> None:
    add_second_instance(instances)
    proc = compare(instances, stub_env(tmp_path, checkout, pw_output(-42, status="Feasible")))
    assert proc.returncode == 0, proc.stderr
    assert "0 match, 0 mismatch, 0 solved without a known optimum, 2 not proven" in proc.stdout


@pytest.mark.parametrize("batch", [False, True])
def test_mismatch_exits_1(tmp_path: Path, checkout: Path, instances: Path, batch: bool) -> None:
    if batch:
        add_second_instance(instances)
        target = instances
    else:
        target = instances / "tiny.sppcc"
    proc = compare(target, stub_env(tmp_path, checkout, pw_output(-41)))
    assert proc.returncode == 1, proc.stdout
    assert ("FAIL" if batch else "MISMATCH") in proc.stdout


def test_match_exits_0(tmp_path: Path, checkout: Path, instances: Path) -> None:
    proc = compare(instances / "tiny.sppcc", stub_env(tmp_path, checkout, pw_output(-42)))
    assert proc.returncode == 0, proc.stdout
    assert "MATCH" in proc.stdout


# --- compare_pathwyse.sh aggregation ---

CPTP_HEADER = "instance,obj,gap_pct,time_s\n"


def aggregate(tmp_path: Path, pw_rows: str, cptp_rows: str) -> subprocess.CompletedProcess[str]:
    pw, cptp = tmp_path / "pw.csv", tmp_path / "cptp.csv"
    pw.write_text(HEADER + "\n" + pw_rows)
    cptp.write_text(CPTP_HEADER + cptp_rows)
    cmd = ["bash", str(AGGREGATE), "--pathwyse", str(pw), "--cptp", str(cptp)]
    return run(cmd, dict(os.environ))


def test_error_row_is_neither_match_nor_speedup(tmp_path: Path) -> None:
    """An ERROR row that still carries an objective and time (as rows did before #172)."""
    proc = aggregate(
        tmp_path,
        "a,s,PWDefault,-42,0.5,ERROR: exit 2,t,m,c,\nb,s,PWDefault,-7,1.0,Optimal,t,m,c,\n",
        "a,-42,0,2.0\nb,-7,0,2.0\n",
    )
    assert proc.returncode == 0, proc.stdout
    assert "Matches: 1  Mismatches: 0  Not both proven: 1" in proc.stdout
    assert "Proven optimal: PathWyse 1, cptp 2" in proc.stdout
    assert "0.50x  (over 1 matching instances)" in proc.stdout


def test_unproven_rows_are_not_compared(tmp_path: Path) -> None:
    proc = aggregate(
        tmp_path,
        "a,s,PWDefault,-40,9.0,Feasible,t,m,c,\nb,s,PWDefault,-7,1.0,Optimal,t,m,c,\n",
        "a,-42,0,2.0\nb,-9,1.5,3600\n",
    )
    assert proc.returncode == 0, proc.stdout
    assert "Matches: 0  Mismatches: 0  Not both proven: 2" in proc.stdout
    assert "Geometric mean" not in proc.stdout


def test_cptp_gap_within_highs_default_is_proven(tmp_path: Path) -> None:
    """gap_pct is a percent (#187): 0.005% is within HiGHS's 0.01%, 0.02% is not."""
    proc = aggregate(
        tmp_path,
        "a,s,PWDefault,-42,0.5,Optimal,t,m,c,\nb,s,PWDefault,-7,1.0,Optimal,t,m,c,\n",
        "a,-42,0.005,2.0\nb,-7,0.02,3600\n",
    )
    assert proc.returncode == 0, proc.stdout
    assert re.search(r"^a .* OK$", proc.stdout, re.M)
    assert re.search(r"^b .* TLIM$", proc.stdout, re.M)


def test_rows_without_an_objective_are_errors(tmp_path: Path) -> None:
    proc = aggregate(
        tmp_path,
        "a,s,PWDefault,,0.5,Optimal,t,m,c,\nb,s,PWDefault,-7,1.0,Optimal,t,m,c,\n",
        "a,-42,0,2.0\nb,,0,2.0\n",
    )
    assert proc.returncode == 0, proc.stdout
    assert re.search(r"^a .* ERR$", proc.stdout, re.M)
    assert re.search(r"^b .* ERR$", proc.stdout, re.M)
    assert "Proven optimal: PathWyse 1, cptp 1" in proc.stdout


def test_cptp_file_without_gap_is_refused(tmp_path: Path) -> None:
    pw, cptp = tmp_path / "pw.csv", tmp_path / "cptp.csv"
    pw.write_text(HEADER + "\na,s,PWDefault,-42,0.5,Optimal,t,m,c,\n")
    cptp.write_text("instance,obj,time_s\na,-42,2.0\n")
    proc = run(
        ["bash", str(AGGREGATE), "--pathwyse", str(pw), "--cptp", str(cptp)], dict(os.environ)
    )
    assert proc.returncode == 1
    assert "no gap_pct column" in proc.stderr


def test_aggregate_mismatch_exits_1(tmp_path: Path) -> None:
    proc = aggregate(tmp_path, "a,s,PWDefault,-41,0.5,Optimal,t,m,c,\n", "a,-42,0,2.0\n")
    assert proc.returncode == 1
    assert "Mismatches: 1" in proc.stdout


# --- run_pathwyse_reference.sh staging ---

OLD_ROWS = (
    "# kept comment\n"
    f"{HEADER}\n"
    "tiny,s,PWDefault,-42,0.100,Optimal,old,m,c,\n"
    "tiny,s,PWBucket,-42,0.200,Optimal,old,m,c,\n"
    "other,s,PWDefault,-1,0.300,Optimal,old,m,c,\n"
)


def reference(
    tmp_path: Path, instances: Path, env: dict[str, str]
) -> tuple[subprocess.CompletedProcess[str], Path]:
    out = tmp_path / "pathwyse-PWDefault.csv"
    out.write_text(OLD_ROWS)
    cmd = ["bash", str(REFERENCE), "--time-limit", "5", "--out", str(out), str(instances)]
    return run(cmd, env), out


def test_reference_refuses_an_unpinned_binary(
    tmp_path: Path, checkout: Path, instances: Path
) -> None:
    write_stamp(checkout, OTHER_COMMIT)
    proc, out = reference(tmp_path, instances, stub_env(tmp_path, checkout, pw_output(-42)))
    assert proc.returncode == 2, proc.stdout + proc.stderr
    assert out.read_text() == OLD_ROWS


def test_reference_keeps_rows_of_an_interrupted_rerun(
    tmp_path: Path, checkout: Path, instances: Path
) -> None:
    write_stamp(checkout, pinned_commit())
    env = stub_env(tmp_path, checkout, pw_output(-42))
    env["PW_STUB_KILL_PARENT"] = "1"
    proc, out = reference(tmp_path, instances, env)
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert out.read_text() == OLD_ROWS


def test_reference_replaces_only_this_algorithms_rows(
    tmp_path: Path, checkout: Path, instances: Path
) -> None:
    write_stamp(checkout, pinned_commit())
    proc, out = reference(tmp_path, instances, stub_env(tmp_path, checkout, pw_output(-42)))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert out.read_text().startswith("# kept comment\n")
    got = [(r["instance"], r["algo"], r["timestamp"] == "old", r["pathwyse"]) for r in rows(out)]
    assert got == [
        ("tiny", "PWBucket", True, "c"),
        ("other", "PWDefault", True, "c"),
        ("tiny", "PWDefault", False, pinned_commit()),
    ]
