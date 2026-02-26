# Getting Started

## Requirements

- GCC 14+ (C++23)
- CMake 3.25+
- TBB (Threading Building Blocks)

HiGHS and Catch2 are fetched automatically via CMake.

## Installation

=== "Python (pip)"

    ```bash
    pip install rcspp-bac
    ```

=== "C++ (from source)"

    ```bash
    apt install libtbb-dev
    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j$(nproc)
    ```

=== "Python bindings (from source)"

    ```bash
    cmake -B build -DCMAKE_BUILD_TYPE=Release -DRCSPP_BUILD_PYTHON=ON
    cmake --build build -j$(nproc)
    ```

    Or via pip:

    ```bash
    pip install .
    ```

## CLI Usage

```bash
./build/rcspp-solve <instance> [--source <node>] [--target <node>] [--<highs_option> <value> ...]
```

Accepts TSPLIB (`.vrp`, `.sppcc`) and numeric (`.txt`) instance formats. All options beyond `--source`/`--target` are forwarded to HiGHS.

### Examples

```bash
# Tour (closed loop from depot)
./build/rcspp-solve benchmarks/instances/spprclib/B-n45-k6-54.sppcc --time_limit 120

# s-t path (source/target read from file)
./build/rcspp-solve tests/data/tiny4_path.txt

# Override source/target via CLI
./build/rcspp-solve tests/data/tiny4.txt --source 0 --target 3

# Suppress HiGHS log output
./build/rcspp-solve tests/data/tiny4.txt --output_flag false
```

## Running Tests

```bash
./build/rcspp_tests
```

## Solver Options

All HiGHS options can be passed via CLI flags or the options map in the API:

| Option | Default | Description |
|---|---|---|
| `time_limit` | 600 | Time limit in seconds |
| `threads` | auto | Number of threads |
| `output_flag` | true | Print solver log |
| `max_cuts_per_separator` | 3 | Top-k cuts per separator per round |
| `separation_tol` | 0.1 | Fractional violation threshold |
| `enable_rglm` | false | Enable RGLM separator |
| `submip_separation` | true | SEC separation at sub-MIP root |
