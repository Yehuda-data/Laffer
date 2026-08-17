"""Explicit request schemas for the research dashboard API."""

from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, model_validator


class FiscalClosure(str, Enum):
    S_LAFFER = "s_laffer"
    G_LAFFER = "g_laffer"


class CalibrationMethod(str, Enum):
    EXTERNAL = "external"
    MODEL_IMPLIED = "model_implied"


class ExternalBalanceConvention(str, Enum):
    NET_IMPORTS = "net_imports"
    TRADE_BALANCE = "trade_balance"


class KappaMode(str, Enum):
    KAPPA = "kappa"
    LABOR_TARGET = "labor_target"


class ModelSpecification(BaseModel):
    """A complete, explicit model specification.

    Economic admissibility is deliberately not enforced as schema bounds. The
    service returns unusual calibrations unchanged with diagnostics.
    """

    model_config = ConfigDict(extra="forbid")

    name: str = "Custom"
    closure: FiscalClosure = FiscalClosure.S_LAFFER
    calibration: CalibrationMethod = CalibrationMethod.EXTERNAL
    external_balance_convention: ExternalBalanceConvention = ExternalBalanceConvention.NET_IMPORTS
    kappa_mode: KappaMode = KappaMode.LABOR_TARGET

    tau_c: float = 0.18
    tau_n: float = 0.28
    tau_k: float = 0.30
    eta: float = 2.0
    phi: float = 1.0

    theta: float | None = 0.33
    delta: float | None = 0.02
    kappa: float | None = None
    n_target: float | None = 0.25
    k_y: float | None = 1.6
    x_y: float | None = 0.12

    R: float = 1.04
    psi: float = 1.02
    gamma: float = 1.0
    debt_y: float = 0.60
    g_y: float | None = 0.20
    s_y: float | None = 0.10
    external_balance_y: float = 0.0
    other_waste_y: float = 0.0

    grid_min: float = 0.0
    grid_max: float = 0.99
    grid_step: float = 0.001

    @model_validator(mode="after")
    def require_method_fields(self) -> "ModelSpecification":
        if self.calibration == CalibrationMethod.EXTERNAL:
            if self.theta is None or self.delta is None:
                raise ValueError("external calibration requires theta and delta")
        else:
            if self.k_y is None or self.x_y is None:
                raise ValueError("model-implied calibration requires k_y and x_y")
        if self.kappa_mode == KappaMode.KAPPA and self.kappa is None:
            raise ValueError("kappa mode requires kappa")
        if self.kappa_mode == KappaMode.LABOR_TARGET and self.n_target is None:
            raise ValueError("labor-target mode requires n_target")
        if self.closure == FiscalClosure.S_LAFFER and self.g_y is None:
            raise ValueError("s-Laffer requires baseline g_y")
        if self.closure == FiscalClosure.G_LAFFER and self.s_y is None:
            raise ValueError("g-Laffer requires baseline s_y")
        if self.grid_step <= 0 or self.grid_max < self.grid_min:
            raise ValueError("invalid tax grid")
        return self


class SensitivityRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    specification: ModelSpecification
    parameter: str
    minimum: float
    maximum: float
    scenarios: int = Field(default=3, ge=2, le=12)


class CompareRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    scenario_a: ModelSpecification
    scenario_b: ModelSpecification


class ApiEnvelope(BaseModel):
    """Documentation model; endpoints return dictionaries with this shape."""

    inputs: dict[str, Any] = Field(default_factory=dict)
    parameters: dict[str, Any] = Field(default_factory=dict)
    equilibrium: dict[str, Any] = Field(default_factory=dict)
    diagnostics: list[dict[str, Any]] = Field(default_factory=list)
    validity: dict[str, Any] = Field(default_factory=dict)
