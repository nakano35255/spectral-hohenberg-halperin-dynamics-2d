#!/usr/bin/env python3
"""Kugui wrapper for the order-parameter-relaxation input generator."""

import runpy
from pathlib import Path


generator = Path(__file__).resolve().parents[1] / "ohtaka_jobs" / "prepare_ohtaka_inputs.py"
runpy.run_path(str(generator), run_name="__main__")
