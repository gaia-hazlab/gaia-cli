"""Line coverage of ``src/gaia_cli`` by the test suite, using only the standard library.

    pixi run coverage                    # report; fail below the task's --fail-under
    pixi run -e dev python scripts/cover.py --fail-under 80 tests/test_io.py

Lines are recorded with ``sys.monitoring`` (Python 3.12+). Each line is reported
once and then disabled, so the tests run at nearly full speed.
"""

import argparse
import sys
import types
from collections import defaultdict
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src" / "gaia_cli"
mon = sys.monitoring
hits: dict[str, set[int]] = defaultdict(set)


def on_line(code: types.CodeType, line: int):
    if code.co_filename.startswith(str(SRC)):
        hits[code.co_filename].add(line)
    return mon.DISABLE


def executable_lines(path: Path) -> set[int]:
    lines, stack = set(), [compile(path.read_text(), str(path), "exec")]
    while stack:
        code = stack.pop()
        lines |= {line for _, _, line in code.co_lines() if line}
        stack += [c for c in code.co_consts if isinstance(c, types.CodeType)]
    return lines


def ranges(lines: list[int]) -> str:
    out, start = [], None
    for prev, line in zip([None, *lines], [*lines, None]):
        if start is None:
            start = line
        elif line != prev + 1:
            out.append(f"{start}" if start == prev else f"{start}-{prev}")
            start = line
    return ", ".join(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fail-under", type=float, default=0, help="minimum total %%")
    args, pytest_args = parser.parse_known_args()

    mon.use_tool_id(mon.COVERAGE_ID, "gaia-cover")
    mon.register_callback(mon.COVERAGE_ID, mon.events.LINE, on_line)
    mon.set_events(mon.COVERAGE_ID, mon.events.LINE)
    status = pytest.main(["-q", "-p", "no:cacheprovider", *pytest_args])
    mon.set_events(mon.COVERAGE_ID, 0)
    mon.free_tool_id(mon.COVERAGE_ID)

    total = covered = 0
    print(f"\n{'file':<28} {'lines':>5} {'cover':>6}  missing")
    for path in sorted(SRC.rglob("*.py")):
        lines = executable_lines(path)
        missing = sorted(lines - hits[str(path)])
        total, covered = total + len(lines), covered + len(lines) - len(missing)
        pct = 100 * (1 - len(missing) / len(lines)) if lines else 100
        name = str(path.relative_to(SRC))
        print(f"{name:<28} {len(lines):>5} {pct:>5.1f}%  {ranges(missing)}")
    pct = 100 * covered / total
    print(f"{'TOTAL':<28} {total:>5} {pct:>5.1f}%")
    if status:
        sys.exit(status)
    if pct < args.fail_under:
        sys.exit(f"Coverage {pct:.1f}% is below --fail-under {args.fail_under:g}%")


if __name__ == "__main__":
    main()
