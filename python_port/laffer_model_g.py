"""Parallel g-Laffer steady state for Trabandt & Uhlig (2011).

This module leaves ``laffer_model.py`` unchanged.  It mirrors its notation and
numerical style, but fixes transfers ``s_bar`` and solves government spending
``govcons_bar`` residually, as in equation (19) and Proposition 3.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import fsolve

# Initial guesses and tolerances deliberately match laffer_model.py / steady.m.
_GUESS_VEC = [0.25, 0.05, 0.10, 0.15, 0.20, 0.30, 0.35, 0.40, 0.45,
              0.50, 0.55, 0.60, 0.65]


@dataclass
class SteadyG:
    y_bar: float
    sy_bar: float
    ky_bar: float
    xy_bar: float
    cy_bar: float
    n_bar: float
    govcons_bar: float
    othergovwaste_bar: float
    taxrev_bar: float
    constaxrev_bar: float
    labtaxrev_bar: float
    captaxrev_bar: float


def _cy_g(n, xy_bar, yn_bar, b_bar, s_bar, tb_bar, taun_bar, tauc_bar,
          tauk_bar, theta, delta, ky_bar, R_bar, psi):
    """Equation (19), in the source code's net-export sign convention.

    Paper net imports m equal ``-tb_bar``.  Waste cancels when the government
    budget is substituted into feasibility because it is both a government
    outlay and a resource use.
    """
    income_tax_y = (taun_bar * (1.0 - theta)
                    + tauk_bar * (theta - delta * ky_bar))
    tax_unaffected_income = b_bar * (R_bar - psi) + s_bar - tb_bar
    return ((1.0 - xy_bar - income_tax_y)
            + tax_unaffected_income / (yn_bar * n)) / (1.0 + tauc_bar)


def _resid_nbar_cfe_g(n_bar, theta, taun_bar, tauc_bar, tauk_bar, FRISCH,
                       xy_bar, yn_bar, ky_bar, delta, b_bar, s_bar, tb_bar,
                       kappa, eta, R_bar, psi):
    """Residual of equations (16), (17), and (19) under the g closure."""
    n = np.asarray(n_bar, dtype=float)
    n = np.where(n <= 0.0, 1e-12, n)
    alpha = ((1.0 - theta) * ((1.0 - taun_bar) / (1.0 + tauc_bar))
             * (FRISCH / (1.0 + FRISCH)))
    cy = _cy_g(n, xy_bar, yn_bar, b_bar, s_bar, tb_bar, taun_bar,
               tauc_bar, tauk_bar, theta, delta, ky_bar, R_bar, psi)
    rhs = kappa * (eta * cy / alpha + 1.0 - eta)
    return n ** (-(1.0 + 1.0 / FRISCH)) - rhs


def _solve_nbar_cfe_g(theta, taun_bar, tauc_bar, tauk_bar, FRISCH, xy_bar,
                       yn_bar, ky_bar, delta, b_bar, s_bar, tb_bar, kappa,
                       eta, R_bar, psi):
    args = (theta, taun_bar, tauc_bar, tauk_bar, FRISCH, xy_bar, yn_bar,
            ky_bar, delta, b_bar, s_bar, tb_bar, kappa, eta, R_bar, psi)
    for guess in _GUESS_VEC:
        with np.errstate(all="ignore"):
            x, _info, ier, _msg = fsolve(
                _resid_nbar_cfe_g, [guess], args=args, full_output=True,
                xtol=1e-12,
            )
            if ier != 1:
                continue
            n = float(x[0])
            if not np.isfinite(n) or n <= 0.0:
                continue
            if abs(float(_resid_nbar_cfe_g(n, *args))) < 1e-9:
                return n
    return np.nan


def steady_g(b_bar, s_bar, tb_bar, theta, delta, kappa,
             othergovwaste_bar, taun_bar, tauc_bar, tauk_bar, psi, FRISCH,
             eta, utility, R_bar, gamma_bar) -> SteadyG:
    """Solve a balanced-growth g-Laffer equilibrium with fixed transfers.

    The unchanged equations determine ``k/y``, ``y/n``, and ``x/y``.  For CFE
    preferences, equations (16), (17), and (19) determine labor.  Government
    spending is then the exact government-budget residual:

    ``g = T - b*(R-psi) - s - waste``.
    """
    ky_bar = ((R_bar - 1.0) / ((1.0 - tauk_bar) * theta)
              + delta / theta) ** -1.0
    yn_bar = (gamma_bar * ky_bar ** theta) ** (1.0 / (1.0 - theta))
    xy_bar = (psi - 1.0 + delta) * ky_bar

    if utility == "CFE":
        n_bar = _solve_nbar_cfe_g(
            theta, taun_bar, tauc_bar, tauk_bar, FRISCH, xy_bar, yn_bar,
            ky_bar, delta, b_bar, s_bar, tb_bar, kappa, eta, R_bar, psi,
        )
    else:
        raise ValueError("steady_g currently implements the requested CFE cases only")

    y_bar = yn_bar * n_bar
    cy_bar = float(_cy_g(
        n_bar, xy_bar, yn_bar, b_bar, s_bar, tb_bar, taun_bar, tauc_bar,
        tauk_bar, theta, delta, ky_bar, R_bar, psi,
    ))
    taxrevy_bar = (tauc_bar * cy_bar + taun_bar * (1.0 - theta)
                   + tauk_bar * (theta - delta * ky_bar))
    taxrev_bar = taxrevy_bar * y_bar
    govcons_bar = (taxrev_bar - b_bar * (R_bar - psi) - s_bar
                   - othergovwaste_bar)
    sy_bar = s_bar / y_bar

    return SteadyG(
        y_bar=y_bar,
        sy_bar=sy_bar,
        ky_bar=ky_bar,
        xy_bar=xy_bar,
        cy_bar=cy_bar,
        n_bar=n_bar,
        govcons_bar=govcons_bar,
        othergovwaste_bar=othergovwaste_bar,
        taxrev_bar=taxrev_bar,
        constaxrev_bar=tauc_bar * cy_bar * y_bar,
        labtaxrev_bar=taun_bar * (1.0 - theta) * y_bar,
        captaxrev_bar=tauk_bar * (theta - delta * ky_bar) * y_bar,
    )
