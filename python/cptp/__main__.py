"""CLI entry point: python -m cptp <instance> [options]."""

import argparse
import sys

from cptp._cptp import Model, load as _load


def _rebuild_with_endpoints(problem, source, target):
    """Rebuild a Problem with overridden source/target."""
    import numpy as np
    from cptp._cptp import Problem

    return Problem(
        num_nodes=problem.num_nodes,
        edges=np.array(problem.graph_edges()),
        edge_costs=np.array(problem.edge_costs),
        profits=np.array(problem.profits),
        demands=np.array(problem.demands),
        capacity=problem.capacity,
        source=source,
        target=target,
        name=problem.name,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="cptp",
        description="Branch-and-cut solver for the Capacitated Profitable Tour Problem",
    )
    parser.add_argument("instance", help="Path to instance file (.txt, .sppcc, or .vrp)")
    parser.add_argument("--source", type=int, default=None, help="Override source node")
    parser.add_argument("--target", type=int, default=None, help="Override target node")
    parser.add_argument("--time_limit", type=float, default=600.0, help="Time limit in seconds (default: 600)")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("--verbose", action="store_true", help="Show HiGHS output")

    args, extra = parser.parse_known_args(argv)

    try:
        problem = _load(args.instance)
    except Exception as e:
        print(f"Error loading instance: {e}", file=sys.stderr)
        return 1

    # Apply source/target overrides
    if args.source is not None or args.target is not None:
        src = args.source if args.source is not None else problem.source
        tgt = args.target if args.target is not None else problem.target
        problem = _rebuild_with_endpoints(problem, src, tgt)

    mode = "tour" if problem.is_tour else f"path {problem.source} -> {problem.target}"
    print(f"Instance: {problem.name} ({problem.num_nodes} nodes, {problem.num_edges} edges, {mode})")

    model = Model()
    model.set_problem(problem)

    options = [
        ("time_limit", str(args.time_limit)),
        ("threads", str(args.threads)),
        ("output_flag", "true" if args.verbose else "false"),
    ]
    # Forward extra --key value pairs to HiGHS
    i = 0
    while i < len(extra):
        if extra[i].startswith("--") and i + 1 < len(extra):
            options.append((extra[i][2:], extra[i + 1]))
            i += 2
        else:
            print(f"Unknown argument: {extra[i]}", file=sys.stderr)
            return 1

    result = model.solve(options)

    # Print results
    if result.has_solution():
        path_str = " -> ".join(str(n) for n in result.tour)
        label = "Tour" if problem.is_tour else "Path"
        print(f"\n{label}: {path_str}")
        print(f"Objective: {result.objective}  Bound: {result.bound}  "
              f"Gap: {result.gap * 100:.2f}%  Time: {result.time_seconds:.2f}s  Nodes: {result.nodes}")

        if result.total_cuts > 0:
            print(f"User cuts: {result.total_cuts} ({result.separation_rounds} rounds)")
            for name, stats in result.separator_stats.items():
                print(f"  {name:<10} {stats.cuts_added:>6} cuts {stats.rounds_called:>6} rounds {stats.time_seconds:>8.3f}s")

        return 0
    else:
        print(f"\nNo feasible solution found (status: {result.status.name}, "
              f"time: {result.time_seconds:.2f}s, nodes: {result.nodes})")
        return 1


if __name__ == "__main__":
    sys.exit(main())
