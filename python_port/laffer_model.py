"""Python port of the Trabandt & Uhlig (2011) baseline model.

Ports LafferCodeWeb/1_BaselineCode/{params.m, get_further_params.m, steady.m,
steady_fsolve_nbar_CFE.m}.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import fsolve

from laffer_data import COUNTRY, COUNTRY_SELECTION, IDX_USA, data

# Initial guesses for n_bar, exactly as in steady.m
_GUESS_VEC = [0.25, 0.05, 0.10, 0.15, 0.20, 0.30, 0.35, 0.40, 0.45,
              0.50, 0.55, 0.60, 0.65]


@dataclass
class Params:
    R_bar: float
    gamma_bar: float
    psi: float
    b_bar: np.ndarray
    govcons_bar: np.ndarray
    tb_bar: np.ndarray
    theta: np.ndarray
    delta: np.ndarray
    kappa: np.ndarray
    othergovwaste_bar: np.ndarray
    country: list
    taun_target: np.ndarray
    tauc_target: np.ndarray
    tauk_target: np.ndarray
    by_target: np.ndarray
    govconsy_target: np.ndarray
    tby_target: np.ndarray
    n_target: np.ndarray
    ky_target: np.ndarray
    cy_target: np.ndarray
    xy_target: np.ndarray
    sy_target: np.ndarray
    othergovwastey_target: float
    FRISCH: float
    eta: float
    utility: str
    use_common_parameters: bool
    constaxrevy_target: np.ndarray
    labtaxrevy_target: np.ndarray
    captaxrevy_target: np.ndarray
    startdate: int
    enddate: int


@dataclass
class Steady:
    y_bar: float
    sy_bar: float
    ky_bar: float
    xy_bar: float
    cy_bar: float
    n_bar: float
    othergovwaste_bar: float
    taxrev_bar: float
    constaxrev_bar: float
    labtaxrev_bar: float
    captaxrev_bar: float


def _resid_nbar_cfe(n_bar, theta, taun_bar, tauc_bar, FRISCH, xy_bar,
                    abbrev1_bar, yn_bar, kappa, eta):
    """Residual of the labour-leisure condition under CFE preferences.

    steady_fsolve_nbar_CFE.m writes it as
        n - (kappa*(eta*cy/alpha + 1 - eta))**(-FRISCH/(FRISCH+1))
    which becomes complex whenever the base is negative.  Raising both sides to
    the power -(FRISCH+1)/FRISCH (a monotone transform for n > 0) gives an
    equivalent equation that stays real, with exactly the same roots.
    """
    n = np.asarray(n_bar, dtype=float)
    n = np.where(n <= 0.0, 1e-12, n)  
    alpha = (1.0 - theta) * ((1.0 - taun_bar) / (1.0 + tauc_bar)) * (FRISCH / (1.0 + FRISCH))
    cy = 1.0 - xy_bar - abbrev1_bar / yn_bar / n
    rhs = kappa * (eta * cy / alpha + 1.0 - eta)
    return n ** (-(1.0 + 1.0 / FRISCH)) - rhs


def _solve_nbar_cfe(theta, taun_bar, tauc_bar, FRISCH, xy_bar, abbrev1_bar,
                    yn_bar, kappa, eta):
    args = (theta, taun_bar, tauc_bar, FRISCH, xy_bar, abbrev1_bar, yn_bar,
            kappa, eta)
    for guess in _GUESS_VEC:
        with np.errstate(all="ignore"):
            x, _info, ier, _msg = fsolve(_resid_nbar_cfe, [guess], args=args,
                                         full_output=True, xtol=1e-12)
            if ier != 1:
                continue
            n = float(x[0])
            if not np.isfinite(n) or n <= 0.0:
                continue
            if abs(float(_resid_nbar_cfe(n, *args))) < 1e-9:
                return n
    return np.nan


def steady(b_bar, govcons_bar, tb_bar, theta, delta, kappa, othergovwaste_bar,
           taun_bar, tauc_bar, tauk_bar, psi, FRISCH, eta, utility, R_bar,
           gamma_bar) -> Steady:
    """Port of steady.m."""
    ky_bar = ((R_bar - 1.0) / ((1.0 - tauk_bar) * theta) + delta / theta) ** -1.0
    yn_bar = (gamma_bar * ky_bar ** theta) ** (1.0 / (1.0 - theta))
    xy_bar = (psi - 1.0 + delta) * ky_bar
    abbrev1_bar = govcons_bar + tb_bar + othergovwaste_bar

    if utility == "CFE":
        n_bar = _solve_nbar_cfe(theta, taun_bar, tauc_bar, FRISCH, xy_bar,
                                abbrev1_bar, yn_bar, kappa, eta)
    elif utility == "C-D":
        alpha2 = (1.0 - theta) * ((1.0 - taun_bar) / (1.0 + tauc_bar))
        n_bar = ((1.0 + (1.0 - kappa) / kappa / alpha2 * abbrev1_bar / yn_bar)
                 / (1.0 + (1.0 - kappa) / kappa / alpha2 * (1.0 - xy_bar)))
    else:
        raise ValueError(f"unknown utility '{utility}'")

    y_bar = yn_bar * n_bar
    cy_bar = 1.0 - xy_bar - abbrev1_bar / yn_bar / n_bar
    taxrevy_bar = (tauc_bar * cy_bar + taun_bar * (1.0 - theta)
                   + tauk_bar * (theta - delta * ky_bar))
    sy_bar = (taxrevy_bar - b_bar / y_bar * (R_bar - psi)
              - govcons_bar / y_bar - othergovwaste_bar / y_bar)

    return Steady(
        y_bar=y_bar,
        sy_bar=sy_bar,
        ky_bar=ky_bar,
        xy_bar=xy_bar,
        cy_bar=cy_bar,
        n_bar=n_bar,
        othergovwaste_bar=othergovwaste_bar,
        taxrev_bar=taxrevy_bar * y_bar,
        constaxrev_bar=tauc_bar * cy_bar * y_bar,
        labtaxrev_bar=taun_bar * (1.0 - theta) * y_bar,
        captaxrev_bar=tauk_bar * (theta - delta * ky_bar) * y_bar,
    )


def get_further_params(by_t, govconsy_t, tby_t, n_t, taun_t, tauc_t, tauk_t,
                       ky_t, cy_t, xy_t, othergovwastey_t, theta0, delta0,
                       R_bar, gamma_bar, psi, use_common, utility, FRISCH, eta):
    """Port of get_further_params.m."""
    ncty = by_t.shape[0]
    b_bar = np.zeros(ncty)
    govcons_bar = np.zeros(ncty)
    tb_bar = np.zeros(ncty)
    theta_v = np.zeros(ncty)
    delta_v = np.zeros(ncty)
    kappa_v = np.zeros(ncty)
    waste_bar = np.zeros(ncty)

    kappa = np.nan  # inherited from the US iteration

    for i in COUNTRY_SELECTION:
        taun_bar, tauc_bar, tauk_bar = taun_t[i], tauc_t[i], tauk_t[i]

        if use_common:
            theta, delta = theta0, delta0
            ky_bar = ((R_bar - 1.0) / ((1.0 - tauk_bar) * theta) + delta / theta) ** -1.0
            yn_bar = (gamma_bar * ky_bar ** theta) ** (1.0 / (1.0 - theta))
            xy_bar = (psi - 1.0 + delta) * ky_bar
            cy_bar = 1.0 - xy_bar - govconsy_t[i] - tby_t[i] - othergovwastey_t

            if utility == "CFE":
                alpha = ((1.0 - theta) * ((1.0 - taun_bar) / (1.0 + tauc_bar))
                         * (FRISCH / (1.0 + FRISCH)))
                if i == IDX_USA:
                    n_bar = n_t[i]
                    kappa = 1.0 / (eta * cy_bar / alpha + 1.0 - eta) * n_bar ** (-(1.0 + 1.0 / FRISCH))
                else:
                    n_bar = (kappa * (eta * cy_bar / alpha + 1.0 - eta)) ** (-FRISCH / (FRISCH + 1.0))
            elif utility == "C-D":
                alpha2 = (1.0 - theta) * ((1.0 - taun_bar) / (1.0 + tauc_bar))
                if i == IDX_USA:
                    n_bar = n_t[i]
                    kappa = 1.0 / (1.0 + alpha2 * (1.0 - n_bar) / n_bar / cy_bar)
                else:
                    n_bar = 1.0 / (1.0 + (1.0 - kappa) / kappa / alpha2 * cy_bar)
            else:
                raise ValueError(f"unknown utility '{utility}'")

            y_bar = yn_bar * n_bar
            waste_y = othergovwastey_t
        else:
            delta = xy_t[i] / ky_t[i] - psi + 1.0
            theta = ((R_bar - 1.0) / (1.0 - tauk_bar) + delta) * ky_t[i]
            yn_bar = (gamma_bar * ky_t[i] ** theta) ** (1.0 / (1.0 - theta))
            waste_y = 1.0 - xy_t[i] - govconsy_t[i] - tby_t[i] - cy_t[i]

            if utility == "CFE":
                alpha = ((1.0 - theta) * ((1.0 - taun_bar) / (1.0 + tauc_bar))
                         * (FRISCH / (1.0 + FRISCH)))
                kappa = 1.0 / (eta * cy_t[i] / alpha + 1.0 - eta) * n_t[i] ** (-(1.0 + 1.0 / FRISCH))
            elif utility == "C-D":
                alpha2 = (1.0 - theta) * ((1.0 - taun_bar) / (1.0 + tauc_bar))
                kappa = 1.0 / (1.0 + alpha2 * (1.0 - n_t[i]) / n_t[i] / cy_t[i])
            else:
                raise ValueError(f"unknown utility '{utility}'")

            y_bar = yn_bar * n_t[i]

        b_bar[i] = by_t[i] * y_bar
        govcons_bar[i] = govconsy_t[i] * y_bar
        tb_bar[i] = tby_t[i] * y_bar
        theta_v[i] = theta
        delta_v[i] = delta
        kappa_v[i] = kappa
        waste_bar[i] = waste_y * y_bar

    return b_bar, govcons_bar, tb_bar, theta_v, delta_v, kappa_v, waste_bar


def params(FRISCH=None, eta=None, utility=None, R_bar=None, startdate=None,
           enddate=None, use_common_parameters=None) -> Params:
    """Port of params.m.  Passing None reproduces MATLAB's `[]` defaults."""
    if startdate is None:
        startdate = 1995
    if enddate is None:
        enddate = 2007
    if R_bar is None:
        R_bar = 1.04
    if eta is None:
        eta = 2.0
    if FRISCH is None:
        FRISCH = 1.0
    if utility is None:
        utility = "CFE"
    if use_common_parameters is None:
        use_common_parameters = True

    d = data(startdate, enddate)
    psi = 1.02
    gamma_bar = 1.0
    theta0, delta0 = 0.380, 0.070

    mean = lambda a: np.nanmean(a, axis=1)

    by_t = mean(d["by"])
    govconsy_t = mean(d["govconsy"] + d["xy_gov"])
    tby_t = mean(d["tby"])
    n_t = mean(d["n"])
    taun_t = mean(d["taun"])
    tauc_t = mean(d["tauc"])
    tauk_t = mean(d["tauk"])
    ky_t = mean(d["ky"])
    cy_t = mean(d["cy"])
    xy_t = mean(d["xy_priv"])
    sy_t = mean(d["sy"])
    constaxrevy_t = mean(d["constaxrevy"])
    labtaxrevy_t = mean(d["labtaxrevy"])
    captaxrevy_t = mean(d["captaxrevy"])
    othergovwastey_t = 0.0

    b_bar, govcons_bar, tb_bar, theta_v, delta_v, kappa_v, waste_bar = get_further_params(
        by_t, govconsy_t, tby_t, n_t, taun_t, tauc_t, tauk_t, ky_t, cy_t, xy_t,
        othergovwastey_t, theta0, delta0, R_bar, gamma_bar, psi,
        use_common_parameters, utility, FRISCH, eta)

    return Params(
        R_bar=R_bar, gamma_bar=gamma_bar, psi=psi, b_bar=b_bar,
        govcons_bar=govcons_bar, tb_bar=tb_bar, theta=theta_v, delta=delta_v,
        kappa=kappa_v, othergovwaste_bar=waste_bar, country=list(COUNTRY),
        taun_target=taun_t, tauc_target=tauc_t, tauk_target=tauk_t,
        by_target=by_t, govconsy_target=govconsy_t, tby_target=tby_t,
        n_target=n_t, ky_target=ky_t, cy_target=cy_t, xy_target=xy_t,
        sy_target=sy_t, othergovwastey_target=othergovwastey_t, FRISCH=FRISCH,
        eta=eta, utility=utility, use_common_parameters=use_common_parameters,
        constaxrevy_target=constaxrevy_t, labtaxrevy_target=labtaxrevy_t,
        captaxrevy_target=captaxrevy_t, startdate=startdate, enddate=enddate,
    )
