"""Sphinx configuration for the cptp documentation site.

Build it from the repository root (see docs/development.md):

    uv run --only-group docs sphinx-build -W -b html docs build/docs/html

The C++ reference needs `doxygen` on PATH; this file runs it on the
user-facing headers before Sphinx reads the pages.
"""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOXYGEN_OUT = ROOT / "build" / "docs" / "doxygen"

# The user-facing C++ headers; the rest of src/ is internal.
CPP_HEADERS = [
    "src/model/model.h",
    "src/model/resource.h",
    "src/core/problem.h",
    "src/core/solution.h",
    "src/core/io.h",
]

DOXYGEN_CONFIG = f"""
PROJECT_NAME     = cptp
INPUT            = {" ".join(f'"{ROOT / h}"' for h in CPP_HEADERS)}
STRIP_FROM_PATH  = "{ROOT / "src"}"
OUTPUT_DIRECTORY = "{DOXYGEN_OUT}"
EXTRACT_ALL      = YES
GENERATE_XML     = YES
GENERATE_HTML    = NO
GENERATE_LATEX   = NO
QUIET            = YES
WARN_IF_UNDOCUMENTED = NO
"""

DOXYGEN_OUT.mkdir(parents=True, exist_ok=True)
subprocess.run(["doxygen", "-"], input=DOXYGEN_CONFIG, text=True, check=True, cwd=ROOT)

project = "cptp"
author = "Simon Spoorendonk"
copyright = "2026, Simon Spoorendonk"

extensions = [
    "myst_parser",
    "sphinx.ext.napoleon",
    "breathe",
    "autoapi.extension",
    "sphinx_llms_txt",
]

html_theme = "furo"
html_title = "cptp"
html_baseurl = "https://spoorendonk.dk/cptp/"

myst_enable_extensions = ["dollarmath", "amsmath"]
myst_heading_anchors = 3

breathe_projects = {"cptp": str(DOXYGEN_OUT / "xml")}
breathe_default_project = "cptp"
# Most members have no doc comment; list them anyway.
breathe_default_members = ("members", "undoc-members")

# Read statically from the committed stub and sources: the docs build does
# not compile the extension.
autoapi_dirs = [str(ROOT / "python" / "cptp")]
autoapi_file_patterns = ["*.pyi", "*.py"]
autoapi_ignore = ["*__main__*"]
autoapi_add_toctree_entry = False

# The root page's first paragraph is the README's badge row.
llms_txt_summary = (
    "cptp is a branch-and-cut solver for the Capacitated Profitable Tour Problem, "
    "with its open s-t path variant and several capacity rows: a C++23 library, "
    "a command-line solver and Python bindings, on HiGHS."
)
# llms-full.txt concatenates the page sources, where the README's MyST
# include and the Breathe directives stay unexpanded: append what they render.
llms_txt_code_files = [f"+:../{f}" for f in ["README.md", *CPP_HEADERS]]
