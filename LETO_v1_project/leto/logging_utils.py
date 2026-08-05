"""Shared console, ArcGIS, and timing helpers."""
from __future__ import annotations

import time
from contextlib import contextmanager
from dataclasses import dataclass

try:
    import arcpy
except Exception:  # Allows lightweight non-ArcPy tests.
    arcpy = None


def message(text: object = "") -> None:
    value = str(text)
    print(value, flush=True)
    if arcpy is not None:
        try:
            arcpy.AddMessage(value)
        except Exception:
            pass


@dataclass
class StageTiming:
    name: str
    seconds: float


@contextmanager
def timed_stage(name: str, timings: list[StageTiming]):
    message("=" * 72)
    message(name.upper())
    message("=" * 72)
    start = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        timings.append(StageTiming(name=name, seconds=elapsed))
        message(f"{name} completed in {elapsed:,.1f} seconds")
