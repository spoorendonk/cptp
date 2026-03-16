# cptp — Capacitated Profitable Tour Problem Solver

## Quick Reference

```bash
# build
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)

# debug build
cmake -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build -j$(nproc)

# python editable install
pip install -e .

# format
clang-format -i src/**/*.{h,cpp}

# test
./build/cptp_tests && pytest tests/python/test_solver.py

# all tests (including extra: max flow, oracle, separators)
./build/cptp_tests && ./build/cptp_tests_extra && pytest tests/python/test_solver.py

# benchmarks
benchmarks/run_benchmarks.sh
```

The build patches HiGHS with custom user-cut and propagator callbacks from `third_party/highs_patch/`.

## Git

- Never commit directly to `main`. Always feature branches.
- Linear history only (squash-merge or rebase-merge).

## Workflow: Plan → Grind

Every task has two phases. Do not skip planning.

### Plan (default)

When given a task, **plan first**: investigate the code, propose an approach,
discuss with the user. Wait for approval before implementing (e.g. "grind", "go", "do it").

### Grind (on approval)

Execute autonomously. Build, test, fix, repeat until green.
Self-review, then fullgate: branch, PR, sync main, push.
Progress lives in files and git — not in your context window.

Only pause and ask a human when:
- A fix requires changing the public API or architecture
- You discover a bug in unrelated code you shouldn't touch
- You're stuck after multiple failed attempts

### Fullgate

Also runs standalone when user says **"fullgate"**:
branch → PR → sync (merge main **into** feature branch) → tests → docs →
push → review → build → test → push fixes → squash-merge → delete branch

### Claiming Work (GitHub)

- `gh issue edit <N> --add-label agent-wip` when starting on an issue or PR
- Check for `agent-wip` label before picking up work
- Remove label and close/merge when done

### Teams

For independent sub-tasks, launch a team. Each teammate works in its own
worktree. Lead integrates: merge, resolve conflicts, build/test the result.
