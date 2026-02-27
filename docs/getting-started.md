# Getting Started

## Requirements

- GCC 14+ (C++23)
- CMake 3.25+
- TBB

HiGHS is required and fetched automatically.

## Build

```bash
apt install libtbb-dev
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

## Python Package

```bash
pip install .
```

## CLI Usage

```bash
./build/rcspp-solve <instance> [--source <node>] [--target <node>] [--<highs_option> <value> ...]
```

### Examples

```bash
./build/rcspp-solve tests/data/tiny4.txt --time_limit 30 --output_flag false
./build/rcspp-solve tests/data/tiny4_path.txt
./build/rcspp-solve tests/data/tiny4.txt --source 0 --target 3
```

## Tests

```bash
./build/rcspp_tests
pytest tests/python/test_solver.py
```
