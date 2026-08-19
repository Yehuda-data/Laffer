"""Parity checks for the browser-only GitHub Pages model engine."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess

import pytest

from backend.dashboard_model_service import baseline, capital_curve, labor_curve, presets
from backend.schemas import ModelSpecification


NODE = shutil.which("node")
PROBE = Path(__file__).with_name("static_model_probe.cjs")


@pytest.mark.skipif(NODE is None, reason="Node.js is required for browser-model parity tests")
def test_static_browser_engine_matches_python() -> None:
    specifications = [item["specification"] for item in presets()["presets"]]
    completed = subprocess.run(
        [NODE, str(PROBE)],
        input=json.dumps(specifications),
        text=True,
        capture_output=True,
        check=True,
    )
    browser_results = json.loads(completed.stdout)

    for raw_specification, browser in zip(specifications, browser_results, strict=True):
        specification = ModelSpecification.model_validate(raw_specification)
        python_baseline = baseline(specification)
        python_labor = labor_curve(specification)
        python_capital = capital_curve(specification)

        assert browser["baseline"]["valid"] == python_baseline["validity"]["valid"]
        for key in ("n", "y", "k_y", "c_y", "T_total_y"):
            assert browser["baseline"][key] == pytest.approx(
                python_baseline["equilibrium"][key], rel=1e-9, abs=1e-10
            )

        assert browser["labor"]["peak_tax"] == python_labor["summary"]["peak_tax"]
        assert browser["labor"]["peak_revenue"] == pytest.approx(
            python_labor["summary"]["peak_revenue"], rel=1e-9, abs=1e-9
        )
        assert browser["labor"]["valid_points"] == python_labor["validity"]["valid_points"]

        assert browser["capital"]["peak_tax"] == python_capital["summary"]["peak_tax"]
        assert browser["capital"]["peak_revenue"] == pytest.approx(
            python_capital["summary"]["peak_revenue"], rel=1e-9, abs=1e-9
        )
        assert browser["capital"]["valid_points"] == python_capital["validity"]["valid_points"]
