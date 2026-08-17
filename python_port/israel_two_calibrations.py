from __future__ import annotations

import argparse
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from laffer_model import Steady, steady

ROOT = Path(__file__).resolve().parent
INPUT_FILE = ROOT / "input" / "Israel_laffer_input_2005_2019.xlsx"
OUTPUT_ROOT = ROOT / "output" / "israel_two_calibrations"
FIGURES_DIR = OUTPUT_ROOT / "figures"
TABLES_DIR = OUTPUT_ROOT / "tables"
DIAGNOSTICS_DIR = OUTPUT_ROOT / "diagnostics"
DOCS_DIR = OUTPUT_ROOT / "docs"

TAUC = 0.18
TAUN = 0.28
TAUK = 0.30
PSI = 1.035
GAMMA = 1.0
OTHER_WASTE = 0.0
THETA_V1 = 0.33
DELTA_V1 = 0.02
N_TARGET_V1 = 0.25
KY_V2 = 1.6
LABOR_NORMALIZATION_DENOMINATOR = 100.0
GRID = np.round(np.arange(0.0, 0.991, 0.001), 10)
EPS_TAX = 0.01

PREFERENCES = [
    {"key": "benchmark", "label": "Benchmark", "phi": 1.0, "eta": 2.0, "color": "#1f5f99", "style": "-"},
    {"key": "alt1", "label": "Alternative 1", "phi": 3.0, "eta": 1.0, "color": "#b52a2a", "style": "--"},
    {"key": "alt2", "label": "Alternative 2", "phi": 3.0, "eta": 2.0, "color": "#238b45", "style": "-."},
]

KY_V2_SOURCE = "User-supplied external calibration input"


@dataclass(frozen=True)
class IsraelInputs:
    raw: pd.DataFrame
    sample_start: int
    sample_end: int
    c_private_y: float
    x_private_y: float
    g_consumption_y: float
    g_investment_y: float
    g_total_y: float
    net_exports_y_raw_mean: float
    net_exports_y: float
    net_imports_y: float
    debt_y: float
    transfers_y: float
    weekly_hours: float
    participation_rate: float
    labor_normalization_denominator: float
    n_observed: float
    R: float
    psi: float
    ky_observed: float | None


@dataclass(frozen=True)
class Calibration:
    preference: dict[str, Any]
    theta: float
    delta: float
    kappa: float
    n_target: float
    ky: float
    yn: float
    baseline: Steady
    b_bar: float
    g_bar: float
    tb_bar: float
    labor_curve: dict[str, Any]
    capital_curve: dict[str, Any]
    labor_self_financing: float
    capital_self_financing: float


def load_inputs() -> IsraelInputs:
    raw = pd.read_excel(INPUT_FILE)
    expected = {
        "Unnamed: 0",
        "Consumption_private_to_GDP",
        "Invetment_private_to_GDP",
        "Consumption_government_to_GDP",
        "Investment_government_to_GDP",
        "Net_exprots_to_GDP",
        "Governent_debt_to_GDP",
        "Transfers_to_GDP",
        "Weekly_hours_worked_per_worker",
        "Participation rate",
        "R",
    }
    missing = expected.difference(raw.columns)
    if missing:
        raise ValueError(f"Missing required workbook columns: {sorted(missing)}")
    if raw[list(expected)].isna().any().any():
        raise ValueError("Israel input contains missing values in required columns")

    mean = lambda col: float(pd.to_numeric(raw[col], errors="raise").mean())
    net_exports_raw = mean("Net_exprots_to_GDP")
    net_exports = net_exports_raw / 100.0 if abs(net_exports_raw) > 0.2 else net_exports_raw
    weekly_hours = mean("Weekly_hours_worked_per_worker")
    participation = mean("Participation rate")
    n_observed = weekly_hours * participation / LABOR_NORMALIZATION_DENOMINATOR
    g_consumption = mean("Consumption_government_to_GDP")
    g_investment = mean("Investment_government_to_GDP")

    return IsraelInputs(
        raw=raw,
        sample_start=int(raw["Unnamed: 0"].min()),
        sample_end=int(raw["Unnamed: 0"].max()),
        c_private_y=mean("Consumption_private_to_GDP"),
        x_private_y=mean("Invetment_private_to_GDP"),
        g_consumption_y=g_consumption,
        g_investment_y=g_investment,
        g_total_y=g_consumption + g_investment,
        net_exports_y_raw_mean=net_exports_raw,
        net_exports_y=net_exports,
        net_imports_y=-net_exports,
        debt_y=mean("Governent_debt_to_GDP"),
        transfers_y=mean("Transfers_to_GDP"),
        weekly_hours=weekly_hours,
        participation_rate=participation,
        labor_normalization_denominator=LABOR_NORMALIZATION_DENOMINATOR,
        n_observed=n_observed,
        R=mean("R"),
        psi=PSI,
        ky_observed=KY_V2,
    )


def cfe_kappa(theta: float, taun: float, tauc: float, phi: float, eta: float, cy: float, n: float) -> float:
    alpha = (1.0 - theta) * ((1.0 - taun) / (1.0 + tauc)) * (phi / (1.0 + phi))
    denominator = eta * cy / alpha + 1.0 - eta
    if denominator <= 0.0:
        raise ValueError(f"Invalid CFE kappa denominator: {denominator}")
    return (1.0 / denominator) * n ** (-(1.0 + 1.0 / phi))


def is_valid(solution: Steady) -> bool:
    return (
        np.isfinite(solution.n_bar)
        and 0.0 <= solution.n_bar <= 1.0
        and np.isfinite(solution.cy_bar)
        and 0.0 <= solution.cy_bar <= 1.0
        and np.isfinite(solution.ky_bar)
        and solution.ky_bar >= 0.0
        and np.isfinite(solution.taxrev_bar)
    )


def solve_curve(cal: dict[str, Any], tax_kind: str) -> dict[str, Any]:
    revenue = np.full(GRID.size, np.nan)
    for i, rate in enumerate(GRID):
        taun = float(rate) if tax_kind == "labor" else TAUN
        tauk = float(rate) if tax_kind == "capital" else TAUK
        result = steady(
            b_bar=cal["b_bar"],
            govcons_bar=cal["g_bar"],
            tb_bar=cal["tb_bar"],
            theta=cal["theta"],
            delta=cal["delta"],
            kappa=cal["kappa"],
            othergovwaste_bar=OTHER_WASTE,
            taun_bar=taun,
            tauc_bar=TAUC,
            tauk_bar=tauk,
            psi=PSI,
            FRISCH=cal["phi"],
            eta=cal["eta"],
            utility="CFE",
            R_bar=cal["R"],
            gamma_bar=GAMMA,
        )
        if is_valid(result):
            revenue[i] = result.taxrev_bar

    normalized = revenue / cal["baseline"].taxrev_bar * 100.0
    peak_index = int(np.nanargmax(normalized))
    return {
        "grid": GRID.copy(),
        "revenue": revenue,
        "normalized": normalized,
        "peak_rate": float(GRID[peak_index]),
        "peak_normalized": float(normalized[peak_index]),
        "max_additional_pct": float(normalized[peak_index] - 100.0),
        "valid_points": int(np.isfinite(normalized).sum()),
    }


def self_financing(cal: dict[str, Any], tax_kind: str) -> float:
    base = cal["baseline"]
    if tax_kind == "labor":
        mechanical_loss = EPS_TAX * (1.0 - cal["theta"]) * base.y_bar
        taun, tauk = TAUN - EPS_TAX, TAUK
    else:
        mechanical_loss = EPS_TAX * (cal["theta"] - cal["delta"] * base.ky_bar) * base.y_bar
        taun, tauk = TAUN, TAUK - EPS_TAX

    lower_tax = steady(
        b_bar=cal["b_bar"],
        govcons_bar=cal["g_bar"],
        tb_bar=cal["tb_bar"],
        theta=cal["theta"],
        delta=cal["delta"],
        kappa=cal["kappa"],
        othergovwaste_bar=OTHER_WASTE,
        taun_bar=taun,
        tauc_bar=TAUC,
        tauk_bar=tauk,
        psi=PSI,
        FRISCH=cal["phi"],
        eta=cal["eta"],
        utility="CFE",
        R_bar=cal["R"],
        gamma_bar=GAMMA,
    )
    dynamic_loss = base.taxrev_bar - lower_tax.taxrev_bar
    return float(100.0 * (1.0 - dynamic_loss / mechanical_loss))


def calibrate_variant1(inputs: IsraelInputs, preference: dict[str, Any]) -> Calibration:
    theta, delta, n_target = THETA_V1, DELTA_V1, N_TARGET_V1
    ky = ((inputs.R - 1.0) / ((1.0 - TAUK) * theta) + delta / theta) ** (-1.0)
    yn = (GAMMA * ky**theta) ** (1.0 / (1.0 - theta))
    xy = (inputs.psi - 1.0 + delta) * ky
    cy = 1.0 - xy - inputs.g_total_y - inputs.net_exports_y - OTHER_WASTE
    kappa = cfe_kappa(theta, TAUN, TAUC, preference["phi"], preference["eta"], cy, n_target)

    y = yn * n_target
    b_bar = inputs.debt_y * y
    g_bar = inputs.g_total_y * y
    tb_bar = inputs.net_exports_y * y
    baseline = steady(
        b_bar=b_bar,
        govcons_bar=g_bar,
        tb_bar=tb_bar,
        theta=theta,
        delta=delta,
        kappa=kappa,
        othergovwaste_bar=OTHER_WASTE,
        taun_bar=TAUN,
        tauc_bar=TAUC,
        tauk_bar=TAUK,
        psi=inputs.psi,
        FRISCH=preference["phi"],
        eta=preference["eta"],
        utility="CFE",
        R_bar=inputs.R,
        gamma_bar=GAMMA,
    )
    if not is_valid(baseline):
        raise RuntimeError(f"Invalid Variant 1 baseline for {preference['label']}")
    if not np.isclose(baseline.n_bar, n_target, atol=1e-9, rtol=0.0):
        raise AssertionError(f"Variant 1 failed labor target: {baseline.n_bar} != {n_target}")

    common = {
        "b_bar": b_bar,
        "g_bar": g_bar,
        "tb_bar": tb_bar,
        "theta": theta,
        "delta": delta,
        "kappa": kappa,
        "phi": preference["phi"],
        "eta": preference["eta"],
        "R": inputs.R,
        "baseline": baseline,
    }
    labor_curve = solve_curve(common, "labor")
    capital_curve = solve_curve(common, "capital")

    return Calibration(
        preference=preference,
        theta=theta,
        delta=delta,
        kappa=kappa,
        n_target=n_target,
        ky=ky,
        yn=yn,
        baseline=baseline,
        b_bar=b_bar,
        g_bar=g_bar,
        tb_bar=tb_bar,
        labor_curve=labor_curve,
        capital_curve=capital_curve,
        labor_self_financing=self_financing(common, "labor"),
        capital_self_financing=self_financing(common, "capital"),
    )


def calibrate_variant2(inputs: IsraelInputs, preference: dict[str, Any]) -> Calibration:
    if inputs.ky_observed is None or inputs.ky_observed <= 0.0:
        raise ValueError("Variant 2 requires a positive externally supplied Israeli k/y")

    ky = inputs.ky_observed
    delta = inputs.x_private_y / ky - (inputs.psi - 1.0)
    theta = ky * ((inputs.R - 1.0) / (1.0 - TAUK) + delta)
    n_target = inputs.n_observed
    if not (0.0 < delta < 1.0):
        raise ValueError(f"Variant 2 implied invalid depreciation delta={delta}")
    if not (0.0 < theta < 1.0):
        raise ValueError(f"Variant 2 implied invalid capital share theta={theta}")

    yn = (GAMMA * ky**theta) ** (1.0 / (1.0 - theta))
    cy = 1.0 - inputs.x_private_y - inputs.g_total_y - inputs.net_exports_y - OTHER_WASTE
    kappa = cfe_kappa(theta, TAUN, TAUC, preference["phi"], preference["eta"], cy, n_target)
    y = yn * n_target
    b_bar = inputs.debt_y * y
    g_bar = inputs.g_total_y * y
    tb_bar = inputs.net_exports_y * y
    baseline = steady(
        b_bar=b_bar,
        govcons_bar=g_bar,
        tb_bar=tb_bar,
        theta=theta,
        delta=delta,
        kappa=kappa,
        othergovwaste_bar=OTHER_WASTE,
        taun_bar=TAUN,
        tauc_bar=TAUC,
        tauk_bar=TAUK,
        psi=inputs.psi,
        FRISCH=preference["phi"],
        eta=preference["eta"],
        utility="CFE",
        R_bar=inputs.R,
        gamma_bar=GAMMA,
    )
    if not is_valid(baseline):
        raise RuntimeError(f"Invalid Variant 2 baseline for {preference['label']}")
    if not np.isclose(baseline.n_bar, n_target, atol=1e-9, rtol=0.0):
        raise AssertionError(f"Variant 2 failed labor target: {baseline.n_bar} != {n_target}")
    if not np.isclose(baseline.ky_bar, ky, atol=1e-10, rtol=0.0):
        raise AssertionError(f"Variant 2 failed supplied k/y target: {baseline.ky_bar} != {ky}")

    common = {
        "b_bar": b_bar,
        "g_bar": g_bar,
        "tb_bar": tb_bar,
        "theta": theta,
        "delta": delta,
        "kappa": kappa,
        "phi": preference["phi"],
        "eta": preference["eta"],
        "R": inputs.R,
        "baseline": baseline,
    }
    labor_curve = solve_curve(common, "labor")
    capital_curve = solve_curve(common, "capital")
    return Calibration(
        preference=preference,
        theta=theta,
        delta=delta,
        kappa=kappa,
        n_target=n_target,
        ky=ky,
        yn=yn,
        baseline=baseline,
        b_bar=b_bar,
        g_bar=g_bar,
        tb_bar=tb_bar,
        labor_curve=labor_curve,
        capital_curve=capital_curve,
        labor_self_financing=self_financing(common, "labor"),
        capital_self_financing=self_financing(common, "capital"),
    )


def fmt(value: Any, digits: int = 8) -> Any:
    if value is None:
        return "N/A"
    if isinstance(value, (float, np.floating)):
        if not np.isfinite(value):
            return "N/A"
        return f"{float(value):.{digits}f}"
    return value


def model_metrics(calibration: Calibration) -> dict[str, float]:
    s = calibration.baseline
    y = s.y_bar
    wage = (1.0 - calibration.theta) * y / s.n_bar
    rental_rate = calibration.theta / s.ky_bar
    consumption_base = s.cy_bar * y
    labor_base = (1.0 - calibration.theta) * y
    capital_base = (calibration.theta - calibration.delta * s.ky_bar) * y
    return {
        "y": y,
        "wage": wage,
        "rental_rate": rental_rate,
        "consumption_base": consumption_base,
        "labor_base": labor_base,
        "capital_base": capital_base,
        "consumption_base_y": s.cy_bar,
        "labor_base_y": 1.0 - calibration.theta,
        "capital_base_y": calibration.theta - calibration.delta * s.ky_bar,
        "taxrev_y": s.taxrev_bar / y,
        "constaxrev_y": s.constaxrev_bar / y,
        "labtaxrev_y": s.labtaxrev_bar / y,
        "captaxrev_y": s.captaxrev_bar / y,
        "g_y": calibration.g_bar / y,
        "b_y": calibration.b_bar / y,
        "tb_y": calibration.tb_bar / y,
    }


def calibration_values(cal: Calibration, inputs: IsraelInputs) -> dict[str, Any]:
    s = cal.baseline
    m = model_metrics(cal)
    return {
        "tauc": TAUC,
        "taun": TAUN,
        "tauk": TAUK,
        "theta": cal.theta,
        "delta": cal.delta,
        "eta": cal.preference["eta"],
        "phi": cal.preference["phi"],
        "kappa": cal.kappa,
        "R": inputs.R,
        "Rminus1": inputs.R - 1.0,
        "psi": inputs.psi,
        "psiminus1": inputs.psi - 1.0,
        "debt_y": inputs.debt_y,
        "gcons_y": inputs.g_consumption_y,
        "ginv_y": inputs.g_investment_y,
        "gtotal_y": inputs.g_total_y,
        "netimports_y": inputs.net_imports_y,
        "transfers_data_y": inputs.transfers_y,
        "cobs_y": inputs.c_private_y,
        "xobs_y": inputs.x_private_y,
        "ky_observed": inputs.ky_observed,
        "hours": inputs.weekly_hours,
        "participation": inputs.participation_rate,
        "normalization": inputs.labor_normalization_denominator,
        "n_observed": inputs.n_observed,
        "n_target": cal.n_target,
        "n_model": s.n_bar,
        "ky_model": s.ky_bar,
        "xy_model": s.xy_bar,
        "cy_model": s.cy_bar,
        "gy_model": m["g_y"],
        "sy_model": s.sy_bar,
        "wage": m["wage"],
        "rental": m["rental_rate"],
        "consbase_level": m["consumption_base"],
        "labbase_level": m["labor_base"],
        "capbase_level": m["capital_base"],
        "consbase_y": m["consumption_base_y"],
        "labbase_y": m["labor_base_y"],
        "capbase_y": m["capital_base_y"],
        "constaxrev_level": s.constaxrev_bar,
        "labtaxrev_level": s.labtaxrev_bar,
        "captaxrev_level": s.captaxrev_bar,
        "taxrev_level": s.taxrev_bar,
        "constaxrev_y": m["constaxrev_y"],
        "labtaxrev_y": m["labtaxrev_y"],
        "captaxrev_y": m["captaxrev_y"],
        "taxrev_y": m["taxrev_y"],
        "labor_baseline": TAUN,
        "labor_peak": cal.labor_curve["peak_rate"],
        "labor_peak_norm": cal.labor_curve["peak_normalized"],
        "labor_additional": cal.labor_curve["max_additional_pct"],
        "labor_self_fin": cal.labor_self_financing,
        "capital_baseline": TAUK,
        "capital_peak": cal.capital_curve["peak_rate"],
        "capital_peak_norm": cal.capital_curve["peak_normalized"],
        "capital_additional": cal.capital_curve["max_additional_pct"],
        "capital_self_fin": cal.capital_self_financing,
    }


def comparison_rows(inputs: IsraelInputs, cal: Calibration) -> list[tuple[str, str, Any]]:
    s = cal.baseline
    m = model_metrics(cal)
    rows: list[tuple[str, str, Any]] = [
        ("Consumption tax rate tau_c [imposed]", "tauc", TAUC),
        ("Labor tax rate tau_n [imposed]", "taun", TAUN),
        ("Capital tax rate tau_k [imposed]", "tauk", TAUK),
        ("Capital share theta [V1 imposed; V2 calibrated if identified]", "theta", cal.theta),
        ("Depreciation delta [V1 imposed; V2 calibrated if identified]", "delta", cal.delta),
        ("Eta [imposed preference specification]", "eta", cal.preference["eta"]),
        ("Phi / Frisch elasticity [imposed preference specification]", "phi", cal.preference["phi"]),
        ("Kappa [calibrated to labor target]", "kappa", cal.kappa),
        ("Gross real return R [observed workbook mean]", "R", inputs.R),
        ("Real interest rate R-1 [observed/transformed]", "Rminus1", inputs.R - 1.0),
        ("Gross balanced-growth factor psi [externally imposed; absent from workbook]", "psi", inputs.psi),
        ("Growth rate psi-1 [externally imposed]", "psiminus1", inputs.psi - 1.0),
        ("Debt/GDP [observed]", "debt_y", inputs.debt_y),
        ("Government consumption/GDP [observed]", "gcons_y", inputs.g_consumption_y),
        ("Government investment/GDP [observed]", "ginv_y", inputs.g_investment_y),
        ("Total government spending/GDP [observed/transformed]", "gtotal_y", inputs.g_total_y),
        ("Net imports/GDP [observed/transformed; negative net exports]", "netimports_y", inputs.net_imports_y),
        ("Transfers/GDP [observed data comparison target]", "transfers_data_y", inputs.transfers_y),
        ("Private consumption/GDP [observed]", "cobs_y", inputs.c_private_y),
        ("Private investment/GDP [observed]", "xobs_y", inputs.x_private_y),
        ("Capital/GDP k/y [user-supplied external input]", "ky_observed", inputs.ky_observed),
        ("Raw weekly hours [observed]", "hours", inputs.weekly_hours),
        ("Participation rate [observed]", "participation", inputs.participation_rate),
        ("Labor normalization denominator [imposed]", "normalization", inputs.labor_normalization_denominator),
        ("Normalized Israeli labor n [observed/transformed]", "n_observed", inputs.n_observed),
        ("Labor calibration target n [V1 imposed; V2 observed]", "n_target", cal.n_target),
        ("Equilibrium labor n [model-implied]", "n_model", s.n_bar),
        ("Capital/output k/y [model-implied in V1]", "ky_model", s.ky_bar),
        ("Investment/output x/y [model-implied]", "xy_model", s.xy_bar),
        ("Consumption/output c/y [model-implied]", "cy_model", s.cy_bar),
        ("Government spending/output g/y [model baseline]", "gy_model", m["g_y"]),
        ("Transfers/output s/y [model-implied fiscal residual]", "sy_model", s.sy_bar),
        ("Wage w [model-implied]", "wage", m["wage"]),
        ("Rental rate d [model-implied]", "rental", m["rental_rate"]),
        ("Consumption tax base [model units]", "consbase_level", m["consumption_base"]),
        ("Labor tax base [model units]", "labbase_level", m["labor_base"]),
        ("Net capital-income tax base [model units]", "capbase_level", m["capital_base"]),
        ("Consumption tax base/GDP [model-implied]", "consbase_y", m["consumption_base_y"]),
        ("Labor tax base/GDP [model-implied]", "labbase_y", m["labor_base_y"]),
        ("Net capital-income tax base/GDP [model-implied]", "capbase_y", m["capital_base_y"]),
        ("Consumption tax revenue [model units]", "constaxrev_level", s.constaxrev_bar),
        ("Labor tax revenue [model units]", "labtaxrev_level", s.labtaxrev_bar),
        ("Capital tax revenue [model units]", "captaxrev_level", s.captaxrev_bar),
        ("Total tax revenue [model units]", "taxrev_level", s.taxrev_bar),
        ("Consumption tax revenue/GDP [model-implied]", "constaxrev_y", m["constaxrev_y"]),
        ("Labor tax revenue/GDP [model-implied]", "labtaxrev_y", m["labtaxrev_y"]),
        ("Capital tax revenue/GDP [model-implied]", "captaxrev_y", m["captaxrev_y"]),
        ("Total tax revenue/GDP [model-implied]", "taxrev_y", m["taxrev_y"]),
        ("Baseline labor tax rate", "labor_baseline", TAUN),
        ("Revenue-maximizing labor tax rate [model-implied]", "labor_peak", cal.labor_curve["peak_rate"]),
        ("Maximum normalized labor-Laffer revenue [baseline=100]", "labor_peak_norm", cal.labor_curve["peak_normalized"]),
        ("Maximum additional labor-Laffer revenue [% baseline]", "labor_additional", cal.labor_curve["max_additional_pct"]),
        ("Labor-tax self-financing rate [%]", "labor_self_fin", cal.labor_self_financing),
        ("Baseline capital tax rate", "capital_baseline", TAUK),
        ("Revenue-maximizing capital tax rate [model-implied]", "capital_peak", cal.capital_curve["peak_rate"]),
        ("Maximum normalized capital-Laffer revenue [baseline=100]", "capital_peak_norm", cal.capital_curve["peak_normalized"]),
        ("Maximum additional capital-Laffer revenue [% baseline]", "capital_additional", cal.capital_curve["max_additional_pct"]),
        ("Capital-tax self-financing rate [%]", "capital_self_fin", cal.capital_self_financing),
    ]
    return rows


def build_comparison(inputs: IsraelInputs, cal1: Calibration, cal2: Calibration) -> pd.DataFrame:
    variant2_values = calibration_values(cal2, inputs)
    records = []
    for label, key, variant1 in comparison_rows(inputs, cal1):
        variant2 = variant2_values[key]
        records.append({"Variable": label, "Variant 1": fmt(variant1), "Variant 2": fmt(variant2)})
    return pd.DataFrame(records, columns=["Variable", "Variant 1", "Variant 2"])


def observed_model_table(inputs: IsraelInputs, cal: Calibration) -> pd.DataFrame:
    rows = [
        ("n", inputs.n_observed, cal.baseline.n_bar),
        ("c/y", inputs.c_private_y, cal.baseline.cy_bar),
        ("x/y (private investment)", inputs.x_private_y, cal.baseline.xy_bar),
        ("k/y", inputs.ky_observed, cal.baseline.ky_bar),
        ("s/y", inputs.transfers_y, cal.baseline.sy_bar),
    ]
    out = []
    for variable, observed, model in rows:
        difference: Any
        if isinstance(observed, (float, np.floating)) and isinstance(model, (float, np.floating)):
            difference = float(model - observed)
        else:
            difference = "N/A"
        out.append({"Variable": variable, "Observed Israel": fmt(observed), "Model-implied": fmt(model), "Difference": fmt(difference)})
    return pd.DataFrame(out, columns=["Variable", "Observed Israel", "Model-implied", "Difference"])


def transformed_inputs_table(inputs: IsraelInputs) -> pd.DataFrame:
    rows = [
        ("Sample period", "Unnamed: 0", "year", f"{inputs.sample_start}-{inputs.sample_end}", "No transformation", f"{inputs.sample_start}-{inputs.sample_end}", "observed"),
        ("Private consumption/GDP", "Consumption_private_to_GDP", "ratio", inputs.c_private_y, "Arithmetic mean", inputs.c_private_y, "observed"),
        ("Private investment/GDP", "Invetment_private_to_GDP", "ratio", inputs.x_private_y, "Arithmetic mean", inputs.x_private_y, "observed"),
        ("Government consumption/GDP", "Consumption_government_to_GDP", "ratio", inputs.g_consumption_y, "Arithmetic mean", inputs.g_consumption_y, "observed"),
        ("Government investment/GDP", "Investment_government_to_GDP", "ratio", inputs.g_investment_y, "Arithmetic mean", inputs.g_investment_y, "observed"),
        ("Total government spending/GDP", "two government columns", "ratio", "N/A", "mean consumption + mean investment", inputs.g_total_y, "observed/transformed"),
        ("Net exports/GDP", "Net_exprots_to_GDP", "percentage points in workbook", inputs.net_exports_y_raw_mean, "divide by 100 because |mean| > 0.2", inputs.net_exports_y, "observed/transformed"),
        ("Net imports/GDP", "Net_exprots_to_GDP", "ratio", inputs.net_exports_y, "multiply net exports by -1", inputs.net_imports_y, "observed/transformed"),
        ("Debt/GDP", "Governent_debt_to_GDP", "ratio", inputs.debt_y, "Arithmetic mean", inputs.debt_y, "observed"),
        ("Transfers/GDP", "Transfers_to_GDP", "ratio", inputs.transfers_y, "Arithmetic mean; comparison target only", inputs.transfers_y, "observed"),
        ("Weekly hours", "Weekly_hours_worked_per_worker", "hours/week", inputs.weekly_hours, "Arithmetic mean", inputs.weekly_hours, "observed"),
        ("Participation rate", "Participation rate", "fraction", inputs.participation_rate, "Arithmetic mean", inputs.participation_rate, "observed"),
        ("Normalized labor n", "hours and participation", "fraction of 100-hour endowment", "N/A", "mean hours * mean participation / 100", inputs.n_observed, "observed/transformed"),
        ("Gross real return R", "R", "gross annual factor", inputs.R, "Arithmetic mean", inputs.R, "observed"),
        ("Balanced-growth factor psi", "not present", "gross annual factor", "N/A", "existing Israel assumption", inputs.psi, "externally imposed"),
        ("Capital/output k/y", "not present in workbook", "ratio", "N/A", "User-supplied external calibration input", inputs.ky_observed, "externally supplied"),
    ]
    return pd.DataFrame(rows, columns=["Variable", "Source column", "Raw unit", "Raw/mean value", "Transformation", "Model value", "Status"])


def identification_table(inputs: IsraelInputs) -> pd.DataFrame:
    rows = [
        ("Observed Israeli private investment/GDP x/y", True, inputs.x_private_y, "Available from workbook"),
        ("Israeli capital/output k/y", True, inputs.ky_observed, KY_V2_SOURCE),
        ("Gross growth factor psi", True, inputs.psi, "Externally imposed existing Israel assumption; not in workbook"),
        ("Gross real return R", True, inputs.R, "Observed workbook mean"),
        ("Capital tax tau_k", True, TAUK, "Externally imposed"),
        ("Depreciation delta", True, inputs.x_private_y / inputs.ky_observed - (inputs.psi - 1.0), "Identified from x/y, supplied k/y, and psi"),
        ("Capital share theta", True, inputs.ky_observed * ((inputs.R - 1.0) / (1.0 - TAUK) + inputs.x_private_y / inputs.ky_observed - (inputs.psi - 1.0)), "Identified from supplied k/y, R, tau_k, and delta"),
        ("Kappa for each preference", True, "See comparison tables", "Calibrated to normalized observed Israeli labor"),
        ("Variant 2 status", True, "IDENTIFIED", "All Variant 2 baselines and curves produced"),
    ]
    return pd.DataFrame(rows, columns=["Identification item", "Available/identified", "Value", "Diagnostic"])


def validation_table(inputs: IsraelInputs, calibration_groups: list[tuple[str, list[Calibration]]]) -> pd.DataFrame:
    rows = []
    for variant, calibrations in calibration_groups:
        for cal in calibrations:
            s = cal.baseline
            m = model_metrics(cal)
            tax_sum = s.constaxrev_bar + s.labtaxrev_bar + s.captaxrev_bar
            fiscal_rhs = m["taxrev_y"] - m["b_y"] * (inputs.R - inputs.psi) - m["g_y"] - OTHER_WASTE / s.y_bar
            checks = [
                ("Equation (15): k/y", s.ky_bar, ((inputs.R - 1.0) / ((1.0 - TAUK) * cal.theta) + cal.delta / cal.theta) ** -1.0, 1e-10),
                ("Capital accumulation: x/y", s.xy_bar, (inputs.psi - 1.0 + cal.delta) * s.ky_bar, 1e-10),
                ("Labor target", s.n_bar, cal.n_target, 1e-9),
                ("Tax revenue decomposition", s.taxrev_bar, tax_sum, 1e-10),
                ("s-Laffer government budget", s.sy_bar, fiscal_rhs, 1e-10),
                ("Labor curve baseline normalization", float(np.interp(TAUN, cal.labor_curve["grid"], cal.labor_curve["normalized"])), 100.0, 1e-8),
                ("Capital curve baseline normalization", float(np.interp(TAUK, cal.capital_curve["grid"], cal.capital_curve["normalized"])), 100.0, 1e-8),
            ]
            if variant == "Variant 2":
                checks.extend([
                    ("User-supplied k/y target", s.ky_bar, inputs.ky_observed, 1e-10),
                    ("Observed private investment/GDP target", s.xy_bar, inputs.x_private_y, 1e-10),
                ])
            for name, actual, expected, tolerance in checks:
                difference = float(actual - expected)
                rows.append({
                    "Variant": variant,
                    "Preference": cal.preference["label"],
                    "Validation": name,
                    "Actual": actual,
                    "Expected": expected,
                    "Difference": difference,
                    "Tolerance": tolerance,
                    "Passed": abs(difference) <= tolerance,
                })
            rows.append({
                "Variant": variant,
                "Preference": cal.preference["label"],
                "Validation": "Labor curve recomputed points",
                "Actual": cal.labor_curve["valid_points"],
                "Expected": len(GRID),
                "Difference": cal.labor_curve["valid_points"] - len(GRID),
                "Tolerance": 0,
                "Passed": cal.labor_curve["valid_points"] > 0,
            })
            rows.append({
                "Variant": variant,
                "Preference": cal.preference["label"],
                "Validation": "Capital curve recomputed points",
                "Actual": cal.capital_curve["valid_points"],
                "Expected": len(GRID),
                "Difference": cal.capital_curve["valid_points"] - len(GRID),
                "Tolerance": 0,
                "Passed": cal.capital_curve["valid_points"] > 0,
            })
    return pd.DataFrame(rows)


def render_table(df: pd.DataFrame, path: Path, title: str, font_size: float = 8.0) -> None:
    display = df.copy().fillna("N/A").astype(str)
    rows, cols = display.shape
    width = max(12.0, min(28.0, 3.2 * cols))
    height = max(3.0, 0.34 * (rows + 2))
    fig, ax = plt.subplots(figsize=(width, height), dpi=180)
    ax.axis("off")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    table = ax.table(cellText=display.values, colLabels=display.columns, loc="center", cellLoc="left", colLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(font_size)
    table.scale(1.0, 1.25)
    for (row, _col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor("#23374e")
            cell.set_text_props(color="white", weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#f2f5f8")
    table.auto_set_column_width(col=list(range(cols)))
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def save_table(df: pd.DataFrame, directory: Path, stem: str, title: str, font_size: float = 8.0) -> None:
    df.to_csv(directory / f"{stem}.csv", index=False, encoding="utf-8-sig")
    df.to_excel(directory / f"{stem}.xlsx", index=False)
    render_table(df, directory / f"{stem}.png", title, font_size=font_size)


def plot_curves(calibrations: list[Calibration], tax_kind: str, variant: str, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 6.2), dpi=180)
    for cal in calibrations:
        curve = cal.labor_curve if tax_kind == "labor" else cal.capital_curve
        pref = cal.preference
        ax.plot(curve["grid"], curve["normalized"], color=pref["color"], ls=pref["style"], lw=2.3,
                label=rf"{pref['label']}: $\varphi={pref['phi']:.0f},\ \eta={pref['eta']:.0f}$")
        ax.plot(curve["peak_rate"], curve["peak_normalized"], "o", color=pref["color"], ms=7, mfc="none")
    baseline_tax = TAUN if tax_kind == "labor" else TAUK
    ax.axvline(baseline_tax, color="black", ls="--", lw=1.5, label=f"Baseline tax={baseline_tax:.0%}")
    label = "Labor" if tax_kind == "labor" else "Capital"
    ax.set_title(
        f"{variant} — Israel {label}-Tax Laffer Curve\n"
        rf"$\theta={calibrations[0].theta:.4f}$, $\delta={calibrations[0].delta:.4f}$, "
        rf"$n^{{target}}={calibrations[0].n_target:.4f}$; "
        rf"$\tau^c=0.18$, $\tau^n=0.28$, $\tau^k=0.30$",
        fontsize=12,
    )
    ax.set_xlabel(f"Steady-State {label} Tax Rate")
    ax.set_ylabel(f"Total Tax Revenue ({variant} baseline = 100)")
    ax.set_xlim(0.0, 0.99)
    ax.grid(alpha=0.22)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def baseline_summary(calibration_groups: list[tuple[str, list[Calibration]]]) -> pd.DataFrame:
    rows = []
    for variant, calibrations in calibration_groups:
        for cal in calibrations:
            s = cal.baseline
            m = model_metrics(cal)
            rows.append({
                "Variant": variant,
                "Preference": cal.preference["label"],
                "phi": cal.preference["phi"],
                "eta": cal.preference["eta"],
                "theta": cal.theta,
                "delta": cal.delta,
                "kappa [calibrated]": cal.kappa,
                "n target": cal.n_target,
                "n model": s.n_bar,
                "k/y model": s.ky_bar,
                "x/y model": s.xy_bar,
                "c/y model": s.cy_bar,
                "s/y model": s.sy_bar,
                "tax revenue/GDP": m["taxrev_y"],
                "labor peak tax": cal.labor_curve["peak_rate"],
                "capital peak tax": cal.capital_curve["peak_rate"],
                "labor self-financing %": cal.labor_self_financing,
                "capital self-financing %": cal.capital_self_financing,
            })
    return pd.DataFrame(rows)


def write_readme(inputs: IsraelInputs, calibration_groups: list[tuple[str, list[Calibration]]]) -> None:
    variant2 = calibration_groups[1][1][0]
    lines = [
        "# Israel Two-Calibrations Workflow",
        "",
        "## Status",
        "",
        "- **Variant 1 completed successfully** for all three CFE specifications.",
        "- **Variant 2 completed successfully** for all three CFE specifications using user-supplied `k/y=1.6`.",
        "- `k/y=1.6` is an external calibration input supplied by the user; it is not described as an observation from the workbook.",
        "",
        "## Source and sample",
        "",
        f"- Source: `{INPUT_FILE.relative_to(ROOT)}`.",
        f"- Sample: {inputs.sample_start}–{inputs.sample_end}; arithmetic means use all {len(inputs.raw)} complete annual observations.",
        "- Raw annual values are preserved in `tables/raw_annual_israel_inputs.*`.",
        "- Transformations are listed in `tables/transformed_input_values.*`.",
        "",
        "## Variant 1 — Baseline-style Israel calibration",
        "",
        "- `theta=0.33` is externally imposed.",
        "- `delta=0.02` is externally imposed.",
        "- `n_target=0.25` is externally imposed as the U.S.-style normalized labor target.",
        "- Taxes are externally imposed at `tau_c=0.18`, `tau_n=0.28`, and `tau_k=0.30`.",
        "- Israeli workbook means provide `R`, government spending/GDP, debt/GDP, and trade balance/GDP.",
        "- `psi=1.035` is the existing Israel assumption; the workbook contains no growth column.",
        "- For each `(phi, eta)`, `kappa` is calibrated using the original common-parameter CFE condition and model-consistent baseline `c/y` so equilibrium labor equals 0.25.",
        "",
        "## Variant 2 — Israel-specific calibration",
        "",
        "The intended identifying equations are:",
        "",
        "- `delta = (x/y)/(k/y) - (psi-1)`.",
        "- `theta = (k/y) * ((R-1)/(1-tau_k) + delta)`.",
        "- `kappa` then matches observed normalized Israeli labor for each preference case.",
        "",
        f"The workbook contains private investment/GDP but no capital/GDP field. With externally supplied `k/y={inputs.ky_observed:.1f}`, Variant 2 implies `delta={variant2.delta:.8f}` and `theta={variant2.theta:.8f}`. This assumption should be replaced by a documented compatible capital-stock/GDP series if one becomes available.",
        "",
        "## Labor normalization",
        "",
        f"- Mean raw weekly hours: `{inputs.weekly_hours:.8f}`.",
        f"- Mean participation rate: `{inputs.participation_rate:.8f}`.",
        f"- Denominator: `{inputs.labor_normalization_denominator:.1f}` weekly hours.",
        f"- Existing transformation: `n = hours * participation / 100 = {inputs.n_observed:.8f}`.",
        "- The denominator is an explicit modeling assumption; the original MATLAB data only labels `n` as hours worked per person and does not document its denominator.",
        "",
        "## Treatment of R and psi",
        "",
        f"- `R={inputs.R:.8f}` is the arithmetic mean of the workbook's gross annual real-return column; `R-1={inputs.R-1:.8f}`.",
        f"- `psi={inputs.psi:.8f}` is externally imposed; `psi-1={inputs.psi-1:.8f}`.",
        "- In the original MATLAB baseline, `R=1.04` and `psi=1.02` are independently imposed. The static implementation does not calculate or use `beta`.",
        "- Current Israel inputs have `R < psi`; reviewers should assess the economic interpretation, especially debt-service and any implied discount factor.",
        "",
        "## Fiscal closure",
        "",
        "Both requested variants use the original `s`-Laffer closure. In each successfully calibrated curve, government spending `g_bar`, government debt `b_bar`, trade balance `tb_bar`, and other government waste are held fixed in levels. Each tax-grid point recomputes the full steady state, and transfers `s/y` are the government-budget residual.",
        "",
        "## Deviations from original MATLAB baseline",
        "",
        "- Variant 1 uses Israel macro inputs, `theta=0.33`, `delta=0.02`, taxes 18%/28%/30%, and `psi=1.035`, rather than the MATLAB defaults.",
        "- Variant 2 additionally uses externally supplied `k/y=1.6`, calibrates `delta` and `theta`, and targets normalized observed Israeli labor.",
        "- The Python CFE solver uses an algebraically transformed residual to remain real while preserving admissible roots.",
        "- Israel labor normalization is explicit and differs from the undocumented original denominator.",
        "- Government consumption plus government investment is treated as total fixed government spending, matching the existing project logic.",
        "- Net exports enter the solver as `tb`; reported net imports have the opposite sign.",
        "",
        "## Outputs",
        "",
        "- `figures/`: separate Variant 1 and Variant 2 labor and capital Laffer figures.",
        "- `tables/`: raw/transformed inputs, comprehensive comparison tables, and baseline summary in CSV/XLSX/PNG.",
        "- `diagnostics/`: observed-vs-model tables, identification check, and equation/closure validation.",
        "",
        "## Reproduction",
        "",
        "Run `python_port/israel_two_calibrations.py` from the project environment. Add `--overwrite` only when intentionally replacing the existing task output directory.",
    ]
    (DOCS_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run both Israel calibration variants")
    parser.add_argument("--overwrite", action="store_true", help="replace the existing task output directory")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if OUTPUT_ROOT.exists() and any(OUTPUT_ROOT.iterdir()):
        if not args.overwrite:
            raise FileExistsError(f"Refusing to overwrite existing output directory without --overwrite: {OUTPUT_ROOT}")
        for existing_file in OUTPUT_ROOT.rglob("*"):
            if existing_file.is_file():
                existing_file.unlink()
    for directory in (FIGURES_DIR, TABLES_DIR, DIAGNOSTICS_DIR, DOCS_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    inputs = load_inputs()
    variant1 = [calibrate_variant1(inputs, preference) for preference in PREFERENCES]
    variant2 = [calibrate_variant2(inputs, preference) for preference in PREFERENCES]
    calibration_groups = [("Variant 1", variant1), ("Variant 2", variant2)]

    shutil.copy2(INPUT_FILE, DOCS_DIR / INPUT_FILE.name)
    save_table(inputs.raw, TABLES_DIR, "raw_annual_israel_inputs", "Raw Israel Input Data (2005–2019)", font_size=6.5)
    save_table(transformed_inputs_table(inputs), TABLES_DIR, "transformed_input_values", "Observed and Transformed Israel Inputs", font_size=7.0)
    save_table(baseline_summary(calibration_groups), TABLES_DIR, "all_variants_baseline_summary", "Both Variants: Baseline and Laffer Summary", font_size=7.0)

    for cal1, cal2 in zip(variant1, variant2, strict=True):
        key = cal1.preference["key"]
        comparison = build_comparison(inputs, cal1, cal2)
        save_table(comparison, TABLES_DIR, f"comparison_{key}", f"{cal1.preference['label']}: Variant 1 vs Variant 2", font_size=7.0)
        save_table(observed_model_table(inputs, cal1), DIAGNOSTICS_DIR,
                   f"variant1_{key}_observed_vs_model", f"Variant 1 {cal1.preference['label']}: Observed vs Model")
        save_table(observed_model_table(inputs, cal2), DIAGNOSTICS_DIR,
                   f"variant2_{key}_observed_vs_model", f"Variant 2 {cal2.preference['label']}: Observed vs Model")

    identification = identification_table(inputs)
    save_table(identification, DIAGNOSTICS_DIR, "variant2_identification_check", "Variant 2 Identification Check", font_size=7.5)
    validation = validation_table(inputs, calibration_groups)
    save_table(validation, DIAGNOSTICS_DIR, "all_variants_model_validation", "Both Variants: Equation and Closure Validation", font_size=7.0)
    if not validation["Passed"].all():
        failed = validation.loc[~validation["Passed"]]
        raise AssertionError(f"Validation failures:\n{failed.to_string(index=False)}")

    plot_curves(variant1, "labor", "Variant 1", FIGURES_DIR / "variant1_labor_laffer.png")
    plot_curves(variant1, "capital", "Variant 1", FIGURES_DIR / "variant1_capital_laffer.png")
    plot_curves(variant2, "labor", "Variant 2", FIGURES_DIR / "variant2_labor_laffer.png")
    plot_curves(variant2, "capital", "Variant 2", FIGURES_DIR / "variant2_capital_laffer.png")
    write_readme(inputs, calibration_groups)

    print("Both variants completed for all preference specifications.")
    print(f"Variant 2 external k/y input: {inputs.ky_observed}")
    print(baseline_summary(calibration_groups).to_string(index=False))
    print(f"Outputs: {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
