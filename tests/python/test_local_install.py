"""Local editable-install smoke test for the Python package.

This is opt-in because it builds the C++ extension in a fresh virtualenv.
Run with:
    RCSPP_RUN_LOCAL_INSTALL_TEST=1 pytest tests/python/test_local_install.py -v
"""

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "tests" / "data"


@pytest.mark.skipif(
    os.environ.get("RCSPP_RUN_LOCAL_INSTALL_TEST") != "1",
    reason="Set RCSPP_RUN_LOCAL_INSTALL_TEST=1 to run local install smoke test.",
)
def test_local_editable_install_and_use(tmp_path: Path):
    """Create a venv, pip install -e ., and run a tiny API smoke check."""
    venv_dir = tmp_path / "venv"
    subprocess.run([sys.executable, "-m", "venv", str(venv_dir)], check=True, timeout=120)

    if os.name == "nt":
        py_exe = venv_dir / "Scripts" / "python.exe"
    else:
        py_exe = venv_dir / "bin" / "python"

    pip_env = dict(os.environ)
    pip_env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"

    subprocess.run(
        [str(py_exe), "-m", "pip", "install", "-U", "pip"],
        check=True,
        timeout=240,
        env=pip_env,
    )
    subprocess.run(
        [str(py_exe), "-m", "pip", "install", "-e", str(REPO_ROOT)],
        check=True,
        timeout=2400,
        env=pip_env,
    )

    smoke = textwrap.dedent(
        f"""
        from rcspp_bac import load, Model, has_highs
        prob = load(r"{DATA_DIR / 'tiny4.txt'}")
        model = Model()
        model.set_problem(prob)
        print("has_highs", has_highs)
        print("ok")
        """
    )
    result = subprocess.run(
        [str(py_exe), "-c", smoke],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"
    assert "ok" in result.stdout
