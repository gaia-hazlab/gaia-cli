"""Times, durations and time windows, as written in a stage file and in
``gaia stage -s``/``-e``. See dev/stage-file.md#windows.
"""

from __future__ import annotations

import dataclasses
import datetime as dt
import difflib
import re

import pandas as pd

WINDOW_KEYS = ("start", "end", "before", "after", "select", "count")
NAMED_WINDOWS = ("pre", "post")
SELECT = ("all", "first", "last", "nearest")
DURATION = re.compile(r"(\d+(?:\.\d+)?)(s|min|h|d|w)")
SECONDS = {"s": 1, "min": 60, "h": 3600, "d": 86400, "w": 604800}
DATE_ONLY = re.compile(r"\d{4}-\d{2}-\d{2}")

# Timed datasets are read this far beyond each end of their window, so that
# periods straddling the window's edges (e.g. PRISM's 12 UTC day) are complete
# when `gaia align` aggregates one product into another's time steps
PAD = pd.Timedelta(days=1)


@dataclasses.dataclass
class Window:
    """A half-open UTC interval ``[start, end)``; ``None`` is an open end."""

    start: pd.Timestamp | None
    end: pd.Timestamp | None
    name: str | None = None
    select: str | None = None
    count: int | None = None

    def __str__(self) -> str:
        return "/".join(
            ".." if t is None else t.isoformat() for t in (self.start, self.end)
        )

    def padded(self) -> Window:
        """The interval actually read: ``PAD`` beyond each closed end."""
        return dataclasses.replace(
            self,
            start=None if self.start is None else self.start - PAD,
            end=None if self.end is None else self.end + PAD,
        )


def unknown(keys, allowed, what: str) -> list[str]:
    out = []
    for key in keys:
        if key not in allowed:
            close = difflib.get_close_matches(str(key), allowed, n=1)
            hint = f"; did you mean '{close[0]}'?" if close else ""
            out.append(f"unknown {what} '{key}'{hint}")
    return out


def parse_time(value, *, end: bool = False):
    """An RFC 3339 timestamp or date as naive UTC.

    A date-only *end* includes that whole day, so it becomes the next midnight.
    """
    date_only = (isinstance(value, dt.date) and not isinstance(value, dt.datetime)) or (
        isinstance(value, str) and DATE_ONLY.fullmatch(value) is not None
    )
    if not isinstance(value, (str, dt.date)):
        raise TypeError(f"{value!r} is not a date or timestamp")
    try:
        t = pd.Timestamp(value)
    except ValueError:
        raise ValueError(f"{value!r} is not a date or timestamp") from None
    if t.tzinfo is not None:
        t = t.tz_convert("UTC").tz_localize(None)
    return t + pd.Timedelta(days=1) if end and date_only else t


def parse_duration(value):
    """A non-negative duration such as ``30d``, ``17min`` or ``PT17M``."""
    s = str(value).strip()
    if m := DURATION.fullmatch(s):
        return pd.Timedelta(seconds=float(m[1]) * SECONDS[m[2]])
    if s.startswith("P"):
        try:
            return pd.Timedelta(s)
        except ValueError:
            pass
    raise ValueError(
        f"{value!r} is not a duration; use e.g. 30d, 17min or 6h (units s, min, h, d, w)"
    )


def parse_window(
    spec, time: pd.Timestamp | None, name: str | None = None
) -> list[Window]:
    """Resolve a window spec (see dev/stage-file.md) into one or more windows."""
    if isinstance(spec, str):
        parts = spec.split("/")
        if len(parts) != 2:
            raise ValueError(f"{spec!r} is not an interval; use 'start/end'")
        start, end = (None if p == ".." else p for p in parts)
        spec = {"start": start, "end": end}
    if not isinstance(spec, dict):
        raise TypeError(f"{spec!r} is not a window")
    if name is None and spec and set(spec) <= set(NAMED_WINDOWS):
        return [w for k, sub in spec.items() for w in parse_window(sub, time, k)]
    if problems := unknown(spec, WINDOW_KEYS, "window key"):
        raise ValueError("; ".join(problems))
    if "start" in spec and "before" in spec or "end" in spec and "after" in spec:
        raise ValueError("give either start or before, and either end or after")
    before, after = (
        parse_duration(spec[k]) if k in spec else None for k in ("before", "after")
    )
    relative = before is not None or after is not None
    if (relative or spec.get("select") == "nearest") and time is None:
        raise ValueError("relative windows need a top-level 'time'")
    anchor = time if relative else None
    start = (
        parse_time(spec["start"])
        if spec.get("start") is not None
        else (anchor - before if before is not None and anchor is not None else anchor)
    )
    end = (
        parse_time(spec["end"], end=True)
        if spec.get("end") is not None
        else (anchor + after if after is not None and anchor is not None else anchor)
    )
    if start is None and end is None:
        raise ValueError("window has neither a start nor an end")
    if start is not None and end is not None and start >= end:
        raise ValueError(f"window starts at {start} but ends at {end}")
    select, count = spec.get("select"), spec.get("count")
    if select is not None and select not in SELECT:
        raise ValueError(f"select must be one of {', '.join(SELECT)}, not {select!r}")
    if count is not None:
        if select is None:
            raise ValueError("count needs select")
        if isinstance(count, bool) or not isinstance(count, int) or count < 1:
            raise ValueError(f"count must be a positive integer, not {count!r}")
    return [Window(start, end, name, select, count)]
