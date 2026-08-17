from __future__ import annotations

import math

import numpy as np
import pytest

from backend.dashboard_model_service import (
    _engine_solution,
    baseline,
    calibrate,
    capital_curve,
    labor_curve,
)
from backend.schemas import (
    CalibrationMethod,
    ExternalBalanceConvention,
    FiscalClosure,
    KappaMode,
    ModelSpecification,
)
from laffer_model import steady
from laffer_model_g import steady_g


def test_exact_us_benchmark(preset_specs):
    result = labor_curve(preset_specs[2])
    assert result["validity"]["valid"]
    assert result["summary"]["peak_tax"] == pytest.approx(0.633, abs=5e-4)
    assert result["summary"]["peak_revenue"] == pytest.approx(129.606122, abs=1e-5)


def test_s_laffer_baseline_equals_existing_engine(preset_specs):
    spec = preset_specs[0]
    context = calibrate(spec)
    dashboard = baseline(spec)["equilibrium"]
    direct = steady(
        context.b_bar, context.g_bar, context.tb_bar, context.theta, context.delta,
        context.kappa, context.waste_bar, spec.tau_n, spec.tau_c, spec.tau_k,
        spec.psi, spec.phi, spec.eta, "CFE", spec.R, spec.gamma,
    )
    assert dashboard["n"] == pytest.approx(direct.n_bar, abs=1e-12)
    assert dashboard["y"] == pytest.approx(direct.y_bar, abs=1e-12)
    assert dashboard["T_total"] == pytest.approx(direct.taxrev_bar, abs=1e-12)
    assert dashboard["g"] == pytest.approx(context.g_bar, abs=1e-12)


def test_g_laffer_baseline_equals_existing_engine(preset_specs):
    s_spec = preset_specs[1]
    g_spec = s_spec.model_copy(update={"closure": FiscalClosure.G_LAFFER})
    context = calibrate(g_spec)
    dashboard = baseline(g_spec)["equilibrium"]
    direct = steady_g(
        context.b_bar, context.s_bar, context.tb_bar, context.theta, context.delta,
        context.kappa, context.waste_bar, g_spec.tau_n, g_spec.tau_c, g_spec.tau_k,
        g_spec.psi, g_spec.phi, g_spec.eta, "CFE", g_spec.R, g_spec.gamma,
    )
    assert dashboard["n"] == pytest.approx(direct.n_bar, abs=1e-12)
    assert dashboard["g"] == pytest.approx(direct.govcons_bar, abs=1e-12)
    assert dashboard["s"] == pytest.approx(context.s_bar, abs=1e-12)
    assert dashboard["T_total"] == pytest.approx(direct.taxrev_bar, abs=1e-12)


def test_model_implied_theta_and_delta_formulas(preset_specs):
    spec = preset_specs[1]
    context = calibrate(spec)
    expected_delta = spec.x_y / spec.k_y - (spec.psi - 1.0)
    expected_theta = spec.k_y * ((spec.R - 1.0) / (1.0 - spec.tau_k) + expected_delta)
    assert context.delta == pytest.approx(expected_delta, abs=1e-14)
    assert context.theta == pytest.approx(expected_theta, abs=1e-14)
    payload = baseline(spec)["parameters"]
    assert payload["theta"]["status"] == "MODEL-IMPLIED"
    assert payload["delta"]["status"] == "MODEL-IMPLIED"
    assert payload["kappa"]["status"] == "CALIBRATED"


def test_trade_balance_and_net_imports_are_opposites(preset_specs):
    original = preset_specs[0]
    m_spec = original.model_copy(update={
        "external_balance_convention": ExternalBalanceConvention.NET_IMPORTS,
        "external_balance_y": 0.031,
    })
    tb_spec = original.model_copy(update={
        "external_balance_convention": ExternalBalanceConvention.TRADE_BALANCE,
        "external_balance_y": -0.031,
    })
    m_context = calibrate(m_spec)
    tb_context = calibrate(tb_spec)
    assert m_context.tb_bar == pytest.approx(tb_context.tb_bar, abs=1e-14)
    assert baseline(m_spec)["equilibrium"]["m_y"] == pytest.approx(0.031, abs=1e-12)


def test_invalid_equilibrium_is_not_clipped(preset_specs):
    bad = preset_specs[0].model_copy(update={"theta": 1.2, "delta": -0.01})
    result = baseline(bad)
    assert not result["validity"]["valid"]
    levels = {item["level"] for item in result["diagnostics"]}
    assert "INVALID EQUILIBRIUM" in levels
    theta = result.get("parameters", {}).get("theta", {}).get("value")
    assert theta in (1.2, None)  # Never silently repaired into (0, 1).


def test_level_invariant_over_curve(preset_specs):
    coarse_s = preset_specs[0].model_copy(update={"grid_step": 0.05})
    s_curve = labor_curve(coarse_s)["curve"]
    g_levels = [row["g"] for row in s_curve if row["g"] is not None]
    assert max(g_levels) - min(g_levels) < 1e-12

    coarse_g = preset_specs[1].model_copy(update={"closure": FiscalClosure.G_LAFFER, "grid_step": 0.05})
    g_curve = labor_curve(coarse_g)["curve"]
    s_levels = [row["s"] for row in g_curve if row["s"] is not None]
    assert max(s_levels) - min(s_levels) < 1e-12


def test_capital_g_laffer_is_explicitly_unsupported(preset_specs):
    spec = preset_specs[1].model_copy(update={"closure": FiscalClosure.G_LAFFER})
    with pytest.raises(NotImplementedError, match="not yet implemented"):
        capital_curve(spec)


def test_high_ky_is_displayed_and_warned(preset_specs):
    result = baseline(preset_specs[0])
    assert result["equilibrium"]["k_y"] == pytest.approx(7.2320755025, abs=1e-9)
    assert any(item["code"] == "HIGH_KY" for item in result["diagnostics"])
