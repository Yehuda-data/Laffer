"""Dashboard adapter around the validated Trabandt–Uhlig Python engines.

This module contains presentation/calibration glue only. Every equilibrium is
computed by the existing ``laffer_model.steady`` or ``laffer_model_g.steady_g``
function; neither engine is copied or modified here.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
import math
import sys
from typing import Any, Iterable

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PYTHON_PORT = ROOT / "python_port"
if str(PYTHON_PORT) not in sys.path:
    sys.path.insert(0, str(PYTHON_PORT))

from israel_two_calibrations import (  # noqa: E402
    DELTA_V1,
    GAMMA,
    KY_V2,
    N_TARGET_V1,
    OTHER_WASTE,
    PREFERENCES,
    PSI,
    TAUC,
    TAUK,
    TAUN,
    THETA_V1,
    calibrate_variant1,
    calibrate_variant2,
    cfe_kappa,
    load_inputs,
)
from laffer_data import IDX_USA  # noqa: E402
from laffer_model import params, steady  # noqa: E402
from laffer_model_g import steady_g  # noqa: E402

from .schemas import (  # noqa: E402
    CalibrationMethod,
    CompareRequest,
    ExternalBalanceConvention,
    FiscalClosure,
    KappaMode,
    ModelSpecification,
    SensitivityRequest,
)

TOL = 1e-8
BRANCH_JUMP_THRESHOLD = 0.05
KEY_RATES = (0.80, 0.90, 0.95, 0.99)
STATUS_INPUT = "INPUT"
STATUS_CALIBRATED = "CALIBRATED"
STATUS_IMPLIED = "MODEL-IMPLIED"
STATUS_OUTPUT = "EQUILIBRIUM OUTPUT"


@dataclass(frozen=True)
class CalibrationContext:
    specification: ModelSpecification
    theta: float
    delta: float
    kappa: float
    kappa_status: str
    n_anchor: float
    k_y_anchor: float
    x_y_anchor: float
    c_y_anchor: float
    y_anchor: float
    b_bar: float
    g_bar: float | None
    s_bar: float | None
    tb_bar: float
    waste_bar: float


def _finite(value: Any) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def _safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if _finite(denominator) and denominator != 0 else math.nan


def json_safe(value: Any) -> Any:
    """Recursively replace non-finite numbers with JSON null."""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def _trade_balance_y(spec: ModelSpecification) -> float:
    if spec.external_balance_convention == ExternalBalanceConvention.NET_IMPORTS:
        return -spec.external_balance_y
    return spec.external_balance_y


def _structural_values(spec: ModelSpecification) -> tuple[float, float, str, str]:
    if spec.calibration == CalibrationMethod.MODEL_IMPLIED:
        assert spec.k_y is not None and spec.x_y is not None
        delta = spec.x_y / spec.k_y - (spec.psi - 1.0)
        theta = spec.k_y * ((spec.R - 1.0) / (1.0 - spec.tau_k) + delta)
        return theta, delta, STATUS_IMPLIED, STATUS_IMPLIED
    assert spec.theta is not None and spec.delta is not None
    return spec.theta, spec.delta, STATUS_INPUT, STATUS_INPUT


def _anchor_ratios(
    spec: ModelSpecification, theta: float, delta: float
) -> tuple[float, float, float, float]:
    """Compute ratios needed to translate user ratios into fixed levels."""
    k_y = ((spec.R - 1.0) / ((1.0 - spec.tau_k) * theta) + delta / theta) ** -1.0
    x_y = (spec.psi - 1.0 + delta) * k_y
    tb_y = _trade_balance_y(spec)
    if spec.closure == FiscalClosure.S_LAFFER:
        assert spec.g_y is not None
        c_y = 1.0 - x_y - spec.g_y - tb_y - spec.other_waste_y
    else:
        assert spec.s_y is not None
        income_tax_y = spec.tau_n * (1.0 - theta) + spec.tau_k * (theta - delta * k_y)
        c_y = (
            1.0
            - x_y
            - income_tax_y
            + spec.debt_y * (spec.R - spec.psi)
            + spec.s_y
            - tb_y
        ) / (1.0 + spec.tau_c)
    return float(k_y), float(x_y), float(c_y), float(tb_y)


def _labor_from_kappa(
    spec: ModelSpecification, theta: float, c_y: float, kappa: float
) -> float:
    alpha = (
        (1.0 - theta)
        * ((1.0 - spec.tau_n) / (1.0 + spec.tau_c))
        * (spec.phi / (1.0 + spec.phi))
    )
    base = kappa * (spec.eta * c_y / alpha + 1.0 - spec.eta)
    if not _finite(base) or base <= 0.0:
        return math.nan
    return float(base ** (-spec.phi / (spec.phi + 1.0)))


def calibrate(spec: ModelSpecification) -> CalibrationContext:
    """Translate an explicit ratio specification into engine-level inputs."""
    theta, delta, _theta_status, _delta_status = _structural_values(spec)
    k_y, x_y, c_y, tb_y = _anchor_ratios(spec, theta, delta)

    if spec.kappa_mode == KappaMode.LABOR_TARGET:
        assert spec.n_target is not None
        n_anchor = spec.n_target
        kappa = cfe_kappa(
            theta, spec.tau_n, spec.tau_c, spec.phi, spec.eta, c_y, n_anchor
        )
        kappa_status = STATUS_CALIBRATED
    else:
        assert spec.kappa is not None
        kappa = spec.kappa
        n_anchor = _labor_from_kappa(spec, theta, c_y, kappa)
        kappa_status = STATUS_INPUT

    y_per_n = (spec.gamma * k_y**theta) ** (1.0 / (1.0 - theta))
    y_anchor = float(y_per_n * n_anchor)
    b_bar = spec.debt_y * y_anchor
    tb_bar = tb_y * y_anchor
    waste_bar = spec.other_waste_y * y_anchor
    g_bar = spec.g_y * y_anchor if spec.closure == FiscalClosure.S_LAFFER else None
    s_bar = spec.s_y * y_anchor if spec.closure == FiscalClosure.G_LAFFER else None
    return CalibrationContext(
        specification=spec,
        theta=float(theta),
        delta=float(delta),
        kappa=float(kappa),
        kappa_status=kappa_status,
        n_anchor=float(n_anchor),
        k_y_anchor=k_y,
        x_y_anchor=x_y,
        c_y_anchor=c_y,
        y_anchor=y_anchor,
        b_bar=float(b_bar),
        g_bar=None if g_bar is None else float(g_bar),
        s_bar=None if s_bar is None else float(s_bar),
        tb_bar=float(tb_bar),
        waste_bar=float(waste_bar),
    )


def _engine_solution(
    context: CalibrationContext,
    *,
    tau_n: float | None = None,
    tau_k: float | None = None,
) -> tuple[Any, float, float]:
    spec = context.specification
    labor_tax = spec.tau_n if tau_n is None else tau_n
    capital_tax = spec.tau_k if tau_k is None else tau_k
    common = dict(
        b_bar=context.b_bar,
        tb_bar=context.tb_bar,
        theta=context.theta,
        delta=context.delta,
        kappa=context.kappa,
        othergovwaste_bar=context.waste_bar,
        taun_bar=labor_tax,
        tauc_bar=spec.tau_c,
        tauk_bar=capital_tax,
        psi=spec.psi,
        FRISCH=spec.phi,
        eta=spec.eta,
        utility="CFE",
        R_bar=spec.R,
        gamma_bar=spec.gamma,
    )
    with np.errstate(all="ignore"):
        if spec.closure == FiscalClosure.S_LAFFER:
            assert context.g_bar is not None
            solution = steady(govcons_bar=context.g_bar, **common)
            g_bar = context.g_bar
            s_bar = float(solution.sy_bar * solution.y_bar)
        else:
            assert context.s_bar is not None
            solution = steady_g(s_bar=context.s_bar, **common)
            g_bar = float(solution.govcons_bar)
            s_bar = context.s_bar
    return solution, float(g_bar), float(s_bar)


def _point(
    context: CalibrationContext,
    *,
    tau_n: float | None = None,
    tau_k: float | None = None,
) -> dict[str, Any]:
    spec = context.specification
    actual_tau_n = spec.tau_n if tau_n is None else float(tau_n)
    actual_tau_k = spec.tau_k if tau_k is None else float(tau_k)
    try:
        solution, g, s = _engine_solution(context, tau_n=actual_tau_n, tau_k=actual_tau_k)
        y = float(solution.y_bar)
        n = float(solution.n_bar)
        k_y = float(solution.ky_bar)
        x_y = float(solution.xy_bar)
        c_y = float(solution.cy_bar)
        k = k_y * y
        x = x_y * y
        c = c_y * y
        w = (1.0 - context.theta) * y / n
        d = context.theta / k_y
        labor_base = (1.0 - context.theta) * y
        capital_base = (context.theta - context.delta * k_y) * y
        consumption_base = c
        t_n = float(solution.labtaxrev_bar)
        t_k = float(solution.captaxrev_bar)
        t_c = float(solution.constaxrev_bar)
        total = float(solution.taxrev_bar)
        resource_residual = y - c - x - g - context.tb_bar - context.waste_bar
        budget_residual = (
            total
            - context.b_bar * (spec.R - spec.psi)
            - g
            - s
            - context.waste_bar
        )
        values = (y, n, k_y, x_y, c_y, k, x, c, g, s, total)
        finite = all(_finite(value) for value in values)
        reasons: list[str] = []
        if not finite:
            reasons.append("solver failure or non-finite equilibrium value")
        if finite and not (0.0 < n <= 1.0):
            reasons.append("labor is outside 0 < n <= 1")
        if finite and c < 0.0:
            reasons.append("consumption is negative")
        if finite and g < 0.0:
            reasons.append("government spending is negative")
        if finite and k_y < 0.0:
            reasons.append("capital/output is negative")
        if _finite(resource_residual) and abs(resource_residual) > TOL:
            reasons.append("resource residual exceeds tolerance")
        if _finite(budget_residual) and abs(budget_residual) > TOL:
            reasons.append("government-budget residual exceeds tolerance")
        valid = not reasons
        return {
            "tau_n": actual_tau_n,
            "tau_k": actual_tau_k,
            "point_valid": valid,
            "valid": valid,
            "invalid_reasons": reasons,
            "n": n,
            "y": y,
            "k": k,
            "k_y": k_y,
            "x": x,
            "x_y": x_y,
            "c": c,
            "c_y": c_y,
            "g": g,
            "g_y": _safe_ratio(g, y),
            "s": s,
            "s_y": _safe_ratio(s, y),
            "m_y": -_safe_ratio(context.tb_bar, y),
            "tb_y": _safe_ratio(context.tb_bar, y),
            "w": w,
            "d": d,
            "labor_tax_base": labor_base,
            "capital_tax_base": capital_base,
            "consumption_tax_base": consumption_base,
            "labor_tax_base_y": _safe_ratio(labor_base, y),
            "capital_tax_base_y": _safe_ratio(capital_base, y),
            "consumption_tax_base_y": _safe_ratio(consumption_base, y),
            "T_n": t_n,
            "T_k": t_k,
            "T_c": t_c,
            "T_total": total,
            "T_n_y": _safe_ratio(t_n, y),
            "T_k_y": _safe_ratio(t_k, y),
            "T_c_y": _safe_ratio(t_c, y),
            "T_total_y": _safe_ratio(total, y),
            "resource_residual": resource_residual,
            "government_budget_residual": budget_residual,
        }
    except Exception as error:  # Invalid calibrations must be reported, not clipped.
        return {
            "tau_n": actual_tau_n,
            "tau_k": actual_tau_k,
            "point_valid": False,
            "valid": False,
            "invalid_reasons": [f"solver failure: {type(error).__name__}: {error}"],
        }


def _diagnostic(level: str, code: str, message: str) -> dict[str, str]:
    return {"level": level, "code": code, "message": message}


def baseline_diagnostics(context: CalibrationContext, point: dict[str, Any]) -> list[dict[str, str]]:
    spec = context.specification
    diagnostics: list[dict[str, str]] = []
    k_y = point.get("k_y", context.k_y_anchor)
    x_y = point.get("x_y", context.x_y_anchor)
    c_y = point.get("c_y")
    n = point.get("n")
    g = point.get("g")
    c = point.get("c")

    if _finite(k_y) and k_y > 4.0:
        diagnostics.append(_diagnostic(
            "WARNING", "HIGH_KY",
            f"k/y = {k_y:.4g} is unusually high. It follows from the supplied theta, delta, R, and tau_k through equation (15).",
        ))
    if _finite(x_y) and x_y > 0.30:
        diagnostics.append(_diagnostic("WARNING", "HIGH_XY", f"x/y = {x_y:.4g} exceeds 0.30."))
    if context.delta <= 0.0:
        diagnostics.append(_diagnostic("WARNING", "NONPOSITIVE_DELTA", f"delta = {context.delta:.4g} is non-positive."))
    elif context.delta > 0.20:
        diagnostics.append(_diagnostic("WARNING", "HIGH_DELTA", f"delta = {context.delta:.4g} exceeds 0.20."))
    if not (0.0 < context.theta < 1.0):
        diagnostics.append(_diagnostic("INVALID EQUILIBRIUM", "INVALID_THETA", f"theta = {context.theta:.4g} is outside (0, 1)."))
    if _finite(n) and not (0.0 < n <= 1.0):
        diagnostics.append(_diagnostic("INVALID EQUILIBRIUM", "INVALID_LABOR", f"n = {n:.4g} is outside (0, 1]."))
    if _finite(c) and c < 0.0:
        diagnostics.append(_diagnostic("INVALID EQUILIBRIUM", "NEGATIVE_CONSUMPTION", f"c = {c:.4g} is negative."))
    if _finite(g) and g < 0.0:
        diagnostics.append(_diagnostic("INVALID EQUILIBRIUM", "NEGATIVE_GOVERNMENT", f"g = {g:.4g} is negative."))
    if _finite(c_y) and abs(c_y) < 0.01:
        diagnostics.append(_diagnostic("WARNING", "NEAR_ZERO_CY", f"c/y = {c_y:.4g} is close to zero."))
    if spec.R < spec.psi:
        diagnostics.append(_diagnostic("WARNING", "R_BELOW_PSI", f"R = {spec.R:.4g} is below psi = {spec.psi:.4g}."))
    beta = spec.psi**spec.eta / spec.R if spec.R != 0 else math.inf
    if _finite(beta) and beta > 1.0:
        diagnostics.append(_diagnostic("WARNING", "BETA_ABOVE_ONE", f"Implied beta = psi^eta/R = {beta:.4g} exceeds 1."))
    if not point.get("valid", False):
        diagnostics.append(_diagnostic(
            "INVALID EQUILIBRIUM", "SOLVER_OR_ADMISSIBILITY",
            "; ".join(point.get("invalid_reasons", ["Equilibrium is invalid."])),
        ))
    if not diagnostics:
        diagnostics.append(_diagnostic("INFO", "ADMISSIBLE", "No configured diagnostic threshold was triggered."))
    return diagnostics


def _parameter_payload(context: CalibrationContext) -> dict[str, dict[str, Any]]:
    spec = context.specification
    theta_status = STATUS_IMPLIED if spec.calibration == CalibrationMethod.MODEL_IMPLIED else STATUS_INPUT
    delta_status = theta_status
    return {
        "tau_c": {"value": spec.tau_c, "status": STATUS_INPUT},
        "tau_n": {"value": spec.tau_n, "status": STATUS_INPUT},
        "tau_k": {"value": spec.tau_k, "status": STATUS_INPUT},
        "theta": {"value": context.theta, "status": theta_status},
        "delta": {"value": context.delta, "status": delta_status},
        "eta": {"value": spec.eta, "status": STATUS_INPUT},
        "phi": {"value": spec.phi, "status": STATUS_INPUT},
        "kappa": {"value": context.kappa, "status": context.kappa_status},
        "R": {"value": spec.R, "status": STATUS_INPUT},
        "psi": {"value": spec.psi, "status": STATUS_INPUT},
        "gamma": {"value": spec.gamma, "status": STATUS_INPUT},
    }


def baseline(spec: ModelSpecification) -> dict[str, Any]:
    try:
        context = calibrate(spec)
        point = _point(context)
        diagnostics = baseline_diagnostics(context, point)
        output_status = {key: STATUS_OUTPUT for key in point if key not in {"invalid_reasons", "valid", "point_valid"}}
        return json_safe({
            "inputs": spec.model_dump(mode="json"),
            "parameters": _parameter_payload(context),
            "equilibrium": point,
            "statuses": output_status,
            "diagnostics": diagnostics,
            "validity": {"valid": bool(point.get("valid")), "reasons": point.get("invalid_reasons", [])},
        })
    except Exception as error:
        return json_safe({
            "inputs": spec.model_dump(mode="json"),
            "parameters": {},
            "equilibrium": {},
            "statuses": {},
            "diagnostics": [_diagnostic("INVALID EQUILIBRIUM", "CALIBRATION_FAILURE", f"{type(error).__name__}: {error}")],
            "validity": {"valid": False, "reasons": [str(error)]},
        })


def _grid(spec: ModelSpecification, baseline_rate: float) -> np.ndarray:
    count = int(round((spec.grid_max - spec.grid_min) / spec.grid_step))
    values = spec.grid_min + np.arange(count + 1) * spec.grid_step
    values = values[(values >= spec.grid_min - 1e-12) & (values <= spec.grid_max + 1e-12)]
    if not np.any(np.isclose(values, baseline_rate, atol=1e-12, rtol=0.0)):
        values = np.append(values, baseline_rate)
    return np.unique(np.round(values, 10))


def _first_rate(points: Iterable[dict[str, Any]], predicate) -> float | None:
    for point in points:
        try:
            if predicate(point):
                return float(point["tau_n"])
        except (KeyError, TypeError, ValueError):
            continue
    return None


def labor_curve(spec: ModelSpecification) -> dict[str, Any]:
    try:
        context = calibrate(spec)
        base = _point(context)
        points = [_point(context, tau_n=float(rate)) for rate in _grid(spec, spec.tau_n)]
        baseline_index = next(i for i, point in enumerate(points) if np.isclose(point["tau_n"], spec.tau_n))

        if spec.closure == FiscalClosure.G_LAFFER:
            connected = [False] * len(points)
            if points[baseline_index].get("point_valid", False):
                connected[baseline_index] = True
                for index in range(baseline_index + 1, len(points)):
                    connected[index] = connected[index - 1] and bool(points[index].get("point_valid"))
                for index in range(baseline_index - 1, -1, -1):
                    connected[index] = connected[index + 1] and bool(points[index].get("point_valid"))
            for point, is_connected in zip(points, connected, strict=True):
                point["valid"] = is_connected
                if point.get("point_valid") and not is_connected:
                    point.setdefault("invalid_reasons", []).append("valid root is not connected to the baseline branch")

        base_total = base.get("T_total", math.nan)
        for key in ("n", "y", "k", "c", "w"):
            denominator = base.get(key, math.nan)
            for point in points:
                point[f"{key}_index"] = _safe_ratio(point.get(key, math.nan), denominator) * 100.0
        for point in points:
            for key in ("T_n", "T_k", "T_c", "T_total"):
                point[f"{key}_index"] = _safe_ratio(point.get(key, math.nan), base_total) * 100.0
            point["T_n_own_index"] = _safe_ratio(point.get("T_n", math.nan), base.get("T_n", math.nan)) * 100.0

        jumps: list[dict[str, float]] = []
        solver_failures: list[float] = []
        for previous, current in zip(points, points[1:]):
            if not _finite(current.get("n")):
                solver_failures.append(float(current["tau_n"]))
            if _finite(previous.get("n")) and _finite(current.get("n")):
                change = abs(float(current["n"]) - float(previous["n"]))
                if change > BRANCH_JUMP_THRESHOLD:
                    jumps.append({"tau_n": float(current["tau_n"]), "delta_n": change})

        valid_points = [point for point in points if point.get("valid") and _finite(point.get("T_total_index"))]
        peak = max(valid_points, key=lambda point: point["T_total_index"]) if valid_points else None
        valid_rates = [point["tau_n"] for point in valid_points]
        diagnostics = baseline_diagnostics(context, base)
        if jumps:
            diagnostics.append(_diagnostic("WARNING", "POSSIBLE_BRANCH_SWITCH", "Possible branch/root switch: adjacent labor jump exceeds 0.05."))
        disconnected = sum(bool(point.get("point_valid")) and not bool(point.get("valid")) for point in points)
        if disconnected:
            diagnostics.append(_diagnostic("WARNING", "DISCONNECTED_ROOTS", f"{disconnected} raw valid roots are outside the baseline-connected branch."))

        summary = {
            "baseline_tax": spec.tau_n,
            "peak_tax": None if peak is None else peak["tau_n"],
            "peak_revenue": None if peak is None else peak["T_total_index"],
            "valid_tau_min": min(valid_rates) if valid_rates else None,
            "valid_tau_max": max(valid_rates) if valid_rates else None,
            "first_g_nonpositive": _first_rate(points, lambda p: _finite(p.get("g")) and p["g"] <= 0),
            "first_n_at_least_one": _first_rate(points, lambda p: _finite(p.get("n")) and p["n"] >= 1),
            "first_c_nonpositive": _first_rate(points, lambda p: _finite(p.get("c")) and p["c"] <= 0),
            "solver_failure_rates": solver_failures,
            "possible_branch_switches": jumps,
        }
        return json_safe({
            "inputs": spec.model_dump(mode="json"),
            "parameters": _parameter_payload(context),
            "baseline": base,
            "curve": points,
            "summary": summary,
            "diagnostics": diagnostics,
            "validity": {"valid": bool(valid_points), "valid_points": len(valid_points), "total_points": len(points)},
        })
    except Exception as error:
        return json_safe({
            "inputs": spec.model_dump(mode="json"), "parameters": {}, "baseline": {}, "curve": [], "summary": {},
            "diagnostics": [_diagnostic("INVALID EQUILIBRIUM", "CURVE_FAILURE", f"{type(error).__name__}: {error}")],
            "validity": {"valid": False, "valid_points": 0, "total_points": 0},
        })


def capital_curve(spec: ModelSpecification) -> dict[str, Any]:
    if spec.closure == FiscalClosure.G_LAFFER:
        raise NotImplementedError("Capital-tax g-Laffer not yet implemented.")
    context = calibrate(spec)
    base = _point(context)
    points = [_point(context, tau_k=float(rate)) for rate in _grid(spec, spec.tau_k)]
    base_total = base.get("T_total", math.nan)
    for point in points:
        for key in ("T_n", "T_k", "T_c", "T_total"):
            point[f"{key}_index"] = _safe_ratio(point.get(key, math.nan), base_total) * 100.0
        for key in ("n", "y", "k", "c", "w"):
            point[f"{key}_index"] = _safe_ratio(point.get(key, math.nan), base.get(key, math.nan)) * 100.0
        point["T_n_own_index"] = _safe_ratio(point.get("T_n", math.nan), base.get("T_n", math.nan)) * 100.0
    valid_points = [point for point in points if point.get("valid") and _finite(point.get("T_total_index"))]
    peak = max(valid_points, key=lambda point: point["T_total_index"]) if valid_points else None
    return json_safe({
        "inputs": spec.model_dump(mode="json"),
        "parameters": _parameter_payload(context),
        "baseline": base,
        "curve": points,
        "summary": {
            "baseline_tax": spec.tau_k,
            "peak_tax": None if peak is None else peak["tau_k"],
            "peak_revenue": None if peak is None else peak["T_total_index"],
        },
        "diagnostics": baseline_diagnostics(context, base),
        "validity": {"valid": bool(valid_points), "valid_points": len(valid_points), "total_points": len(points)},
    })


def sensitivity(request: SensitivityRequest) -> dict[str, Any]:
    allowed = {"theta", "delta", "R", "psi", "eta", "phi", "kappa", "k_y", "x_y", "g_y", "s_y", "m_y"}
    if request.parameter not in allowed:
        raise ValueError(f"unsupported sensitivity parameter: {request.parameter}")
    if request.specification.calibration == CalibrationMethod.MODEL_IMPLIED and request.parameter in {"theta", "delta"}:
        raise ValueError("theta and delta are model-implied and cannot be edited directly")
    if request.parameter == "g_y" and request.specification.closure != FiscalClosure.S_LAFFER:
        raise ValueError("g_y sensitivity is available only for s-Laffer")
    if request.parameter == "s_y" and request.specification.closure != FiscalClosure.G_LAFFER:
        raise ValueError("s_y sensitivity is available only for g-Laffer")

    scenarios = []
    for value in np.linspace(request.minimum, request.maximum, request.scenarios):
        updates: dict[str, Any] = {}
        parameter = request.parameter
        if parameter == "m_y":
            updates["external_balance_convention"] = ExternalBalanceConvention.NET_IMPORTS
            updates["external_balance_y"] = float(value)
        else:
            updates[parameter] = float(value)
        if parameter == "kappa":
            updates["kappa_mode"] = KappaMode.KAPPA
            updates["kappa"] = float(value)
        scenario_spec = request.specification.model_copy(update=updates)
        result = labor_curve(scenario_spec)
        scenarios.append({"label": f"{parameter} = {value:.6g}", "value": float(value), "result": result})
    return json_safe({"parameter": request.parameter, "scenarios": scenarios})


def compare(request: CompareRequest) -> dict[str, Any]:
    def scenario(spec: ModelSpecification) -> dict[str, Any]:
        base = baseline(spec)
        curve = labor_curve(spec)
        key_rows = []
        for rate in KEY_RATES:
            rows = curve.get("curve", [])
            row = min(rows, key=lambda item: abs(item["tau_n"] - rate)) if rows else {}
            key_rows.append({
                "tau_n": rate,
                "n": row.get("n") if row.get("valid") else None,
                "T_total_index": row.get("T_total_index") if row.get("valid") else None,
                "valid": bool(row.get("valid")),
                "reasons": row.get("invalid_reasons", []),
            })
        return {"baseline": base, "curve": curve, "key_rates": key_rows}
    return json_safe({"scenario_a": scenario(request.scenario_a), "scenario_b": scenario(request.scenario_b)})


def _spec_from_calibration(name: str, calibration, inputs, method: CalibrationMethod) -> ModelSpecification:
    return ModelSpecification(
        name=name,
        closure=FiscalClosure.S_LAFFER,
        calibration=method,
        external_balance_convention=ExternalBalanceConvention.NET_IMPORTS,
        kappa_mode=KappaMode.LABOR_TARGET,
        tau_c=TAUC,
        tau_n=TAUN,
        tau_k=TAUK,
        eta=float(calibration.preference["eta"]),
        phi=float(calibration.preference["phi"]),
        theta=float(calibration.theta),
        delta=float(calibration.delta),
        kappa=float(calibration.kappa),
        n_target=float(calibration.n_target),
        k_y=float(calibration.ky),
        x_y=float(inputs.x_private_y),
        R=float(inputs.R),
        psi=float(inputs.psi),
        gamma=GAMMA,
        debt_y=float(inputs.debt_y),
        g_y=float(inputs.g_total_y),
        s_y=float(calibration.baseline.sy_bar),
        external_balance_y=float(inputs.net_imports_y),
        other_waste_y=OTHER_WASTE,
    )


@lru_cache(maxsize=1)
def presets() -> dict[str, Any]:
    inputs = load_inputs()
    preference = PREFERENCES[0]
    input_calibration = calibrate_variant1(inputs, preference)
    implied_calibration = calibrate_variant2(inputs, preference)

    p = params(FRISCH=1.0, eta=2.0, utility="CFE", R_bar=1.04, startdate=1995, enddate=2007, use_common_parameters=True)
    i = IDX_USA
    us_baseline = steady(
        b_bar=float(p.b_bar[i]), govcons_bar=float(p.govcons_bar[i]), tb_bar=float(p.tb_bar[i]),
        theta=float(p.theta[i]), delta=float(p.delta[i]), kappa=float(p.kappa[i]),
        othergovwaste_bar=float(p.othergovwaste_bar[i]), taun_bar=float(p.taun_target[i]),
        tauc_bar=float(p.tauc_target[i]), tauk_bar=float(p.tauk_target[i]), psi=float(p.psi),
        FRISCH=float(p.FRISCH), eta=float(p.eta), utility=p.utility, R_bar=float(p.R_bar), gamma_bar=float(p.gamma_bar),
    )
    us = ModelSpecification(
        name="Paper / US benchmark", closure=FiscalClosure.S_LAFFER,
        calibration=CalibrationMethod.EXTERNAL,
        external_balance_convention=ExternalBalanceConvention.TRADE_BALANCE,
        kappa_mode=KappaMode.KAPPA,
        tau_c=float(p.tauc_target[i]), tau_n=float(p.taun_target[i]), tau_k=float(p.tauk_target[i]),
        eta=float(p.eta), phi=float(p.FRISCH), theta=float(p.theta[i]), delta=float(p.delta[i]),
        kappa=float(p.kappa[i]), n_target=float(p.n_target[i]), k_y=float(us_baseline.ky_bar),
        x_y=float(us_baseline.xy_bar), R=float(p.R_bar), psi=float(p.psi), gamma=float(p.gamma_bar),
        debt_y=float(p.b_bar[i] / us_baseline.y_bar), g_y=float(p.govcons_bar[i] / us_baseline.y_bar),
        s_y=float(us_baseline.sy_bar), external_balance_y=float(p.tb_bar[i] / us_baseline.y_bar),
        other_waste_y=float(p.othergovwaste_bar[i] / us_baseline.y_bar),
    )
    values = [
        _spec_from_calibration("Israel — Input assumptions", input_calibration, inputs, CalibrationMethod.EXTERNAL),
        _spec_from_calibration("Israel — Model-implied", implied_calibration, inputs, CalibrationMethod.MODEL_IMPLIED),
        us,
    ]
    return json_safe({"presets": [{"id": f"preset_{index + 1}", "label": spec.name, "specification": spec.model_dump(mode="json")} for index, spec in enumerate(values)]})


def equations(closure: FiscalClosure) -> dict[str, Any]:
    common = [
        {"name": "Capital/output", "latex": r"k/y=[(R-1)/(theta(1-tau_k))+delta/theta]^{-1}", "source": "laffer_model.steady / equation (15)"},
        {"name": "Capital accumulation", "latex": r"x/y=(psi-1+delta)k/y", "source": "laffer_model.steady"},
        {"name": "Production", "latex": r"y/n=[gamma(k/y)^theta]^{1/(1-theta)}", "source": "laffer_model.steady"},
        {"name": "Tax revenue", "latex": r"T/y=tau_c(c/y)+tau_n(1-theta)+tau_k[theta-delta(k/y)]", "source": "laffer_model.steady and laffer_model_g.steady_g"},
        {"name": "CFE labor", "latex": r"n^{-(1+1/phi)}=kappa[eta(c/y)/alpha+1-eta]", "source": "existing CFE residual solvers"},
    ]
    if closure == FiscalClosure.S_LAFFER:
        common.append({"name": "s-Laffer closure", "latex": r"s/y=T/y-[b/y](R-psi)-g/y-q/y; g fixed in levels", "source": "laffer_model.steady"})
    else:
        common.extend([
            {"name": "g-Laffer consumption", "latex": r"c/y=(1+tau_c)^{-1}{1-x/y-tau_n(1-theta)-tau_k[theta-delta(k/y)]+[b(R-psi)+s-tb]/y}", "source": "laffer_model_g._cy_g, equation (19)"},
            {"name": "g-Laffer closure", "latex": r"g=T-b(R-psi)-s-q; s fixed in levels", "source": "laffer_model_g.steady_g / Proposition 3"},
        ])
    return {"closure": closure.value, "equations": common}
