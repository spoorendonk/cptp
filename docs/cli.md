# CLI Usage

## Command

```bash
./build/rcspp-solve <instance> [--source <node>] [--target <node>] [--<highs_option> <value> ...]
```

## Examples

```bash
./build/rcspp-solve tests/data/tiny4.txt --time_limit 30 --output_flag false
./build/rcspp-solve tests/data/tiny4_path.txt
./build/rcspp-solve tests/data/tiny4.txt --source 0 --target 3
```

## Parameters

| Parameter | Type | Description |
|---|---|---|
| `--source` | int | Override source node |
| `--target` | int | Override target node |
| `--<highs_option> <value>` | mixed | Forwarded to HiGHS |

TODO: expand this into a full parameter/options table (including defaults and solver-specific options).

## Instances

See [Instance Formats](instance-formats.md) for supported file formats.
