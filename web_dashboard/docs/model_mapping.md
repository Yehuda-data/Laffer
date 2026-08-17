# Dashboard-to-model mapping

## Status labels

- **INPUT** — supplied by the user or selected preset.
- **CALIBRATED** — chosen to reproduce the labor target.
- **MODEL-IMPLIED** — inferred from model identities and supplied targets.
- **EQUILIBRIUM OUTPUT** — returned by `steady()` or `steady_g()` or derived directly from those returned levels.

## Controls and parameters

| Dashboard variable | Python/API variable | Engine/source | Status |
|---|---|---|---|
| Fiscal closure | `closure` | selects `laffer_model.steady` or `laffer_model_g.steady_g` | INPUT |
| Calibration | `calibration` | dashboard adapter selection | INPUT |
| Consumption tax | `tau_c` → `tauc_bar` | both steady-state engines | INPUT |
| Labor tax | `tau_n` → `taun_bar` | both steady-state engines | INPUT |
| Capital tax | `tau_k` → `tauk_bar` | both steady-state engines | INPUT |
| Intertemporal preference | `eta` | CFE labor equation | INPUT |
| Frisch elasticity | `phi` → `FRISCH` | CFE labor equation | INPUT |
| Capital share | `theta` | production and equation (15) | INPUT or MODEL-IMPLIED |
| Depreciation | `delta` | equation (15), capital accumulation and capital-tax base | INPUT or MODEL-IMPLIED |
| Labor weight | `kappa` | CFE labor equation | INPUT or CALIBRATED |
| Labor target | `n_target` | argument to existing `cfe_kappa()` calibration equation | INPUT |
| Capital/output target | `k_y` | model-implied calibration | INPUT |
| Investment/output target | `x_y` | model-implied calibration | INPUT |
| Gross real return | `R` → `R_bar` | equation (15), government budget | INPUT |
| Balanced-growth factor | `psi` | accumulation and government budget | INPUT |
| Production scale | `gamma` → `gamma_bar` | production / `y/n` | INPUT |
| Debt/GDP | `debt_y` | converted to fixed `b_bar` at baseline | INPUT |
| Government spending/GDP | `g_y` | converted to `govcons_bar`; active only for s-Laffer | INPUT |
| Transfers/GDP | `s_y` | converted to `s_bar`; active only for g-Laffer | INPUT |
| Other waste/GDP | `other_waste_y` | converted to `othergovwaste_bar` | INPUT |
| Net imports/GDP | `external_balance_y` under `net_imports` | converted to `tb/y = -m/y` | INPUT |
| Trade balance/GDP | `external_balance_y` under `trade_balance` | passed as `tb/y` | INPUT |

## Model-implied calibration

| Dashboard output | Formula/source | Status |
|---|---|---|
| `delta` | `(x/y)/(k/y) - (psi - 1)`; existing Variant 2 logic | MODEL-IMPLIED |
| `theta` | `(k/y) * ((R - 1)/(1 - tau_k) + delta)`; existing Variant 2 logic | MODEL-IMPLIED |
| `kappa` | `israel_two_calibrations.cfe_kappa()` evaluated at `n_target` | CALIBRATED |

The adapter does not alter an externally entered `theta`, `delta`, `R`, or tax rate to produce a more conventional `k/y`.

## Baseline equilibrium

| Dashboard variable | Python expression/source | Status |
|---|---|---|
| Labor `n` | `solution.n_bar` | EQUILIBRIUM OUTPUT |
| Output `y` | `solution.y_bar` | EQUILIBRIUM OUTPUT |
| Capital/output `k/y` | `solution.ky_bar` | EQUILIBRIUM OUTPUT |
| Investment/output `x/y` | `solution.xy_bar` | EQUILIBRIUM OUTPUT |
| Consumption/output `c/y` | `solution.cy_bar` | EQUILIBRIUM OUTPUT |
| Capital `k` | `solution.ky_bar * solution.y_bar` | EQUILIBRIUM OUTPUT |
| Investment `x` | `solution.xy_bar * solution.y_bar` | EQUILIBRIUM OUTPUT |
| Consumption `c` | `solution.cy_bar * solution.y_bar` | EQUILIBRIUM OUTPUT |
| Government spending `g` | fixed input level for s; `SteadyG.govcons_bar` for g | EQUILIBRIUM OUTPUT |
| Transfers `s` | `Steady.sy_bar*y` for s; fixed input level for g | EQUILIBRIUM OUTPUT |
| Net imports/output `m/y` | `-tb_bar/y` | EQUILIBRIUM OUTPUT |
| Wage `w` | `(1-theta)y/n` | EQUILIBRIUM OUTPUT |
| Rental rate `d` | `theta/(k/y)` | EQUILIBRIUM OUTPUT |

## Tax bases and revenues

| Dashboard variable | Python expression/source | Status |
|---|---|---|
| Labor tax base | `(1-theta)y` | EQUILIBRIUM OUTPUT |
| Capital tax base | `[theta-delta(k/y)]y` | EQUILIBRIUM OUTPUT |
| Consumption tax base | `c` | EQUILIBRIUM OUTPUT |
| `T_n` | `solution.labtaxrev_bar` | EQUILIBRIUM OUTPUT |
| `T_k` | `solution.captaxrev_bar` | EQUILIBRIUM OUTPUT |
| `T_c` | `solution.constaxrev_bar` | EQUILIBRIUM OUTPUT |
| `T_total` | `solution.taxrev_bar` | EQUILIBRIUM OUTPUT |
| Revenue/GDP fields | corresponding returned revenue divided by `y` | EQUILIBRIUM OUTPUT |

## Equations used by both closures

Capital/output, from equation (15):

```text
k/y = [(R-1)/(theta(1-tau_k)) + delta/theta]^(-1)
```

Balanced-growth capital accumulation:

```text
x/y = (psi - 1 + delta) k/y
```

Production / output per labor:

```text
y/n = [gamma (k/y)^theta]^[1/(1-theta)]
```

Tax revenue:

```text
T/y = tau_c(c/y) + tau_n(1-theta) + tau_k[theta-delta(k/y)]
```

CFE labor equation:

```text
n^[-(1+1/phi)] = kappa[eta(c/y)/alpha + 1-eta]
alpha = (1-theta)[(1-tau_n)/(1+tau_c)] phi/(1+phi)
```

## Closure mapping

### s-Laffer

`dashboard_model_service` passes a fixed `govcons_bar` to `laffer_model.steady()` at every tax-grid point. Transfers are returned as the fiscal residual:

```text
s/y = T/y - (b/y)(R-psi) - g/y - waste/y
```

Consequently, `g` is constant in levels while `g/y` can move.

### g-Laffer

`dashboard_model_service` passes a fixed `s_bar` to `laffer_model_g.steady_g()` at every labor-tax point. Consumption follows the existing equation (19) implementation in `_cy_g()`, and government spending is:

```text
g = T - b(R-psi) - s - waste
```

Consequently, `s` is constant in levels while `s/y` can move. Capital-tax g-Laffer is not mapped because no validated implementation currently exists.

## Curves and normalization

| Display | Backend field | Normalization |
|---|---|---|
| Main Laffer curve | `T_total_index` | total revenue / baseline total revenue × 100 |
| Revenue decomposition | `T_n_index`, `T_k_index`, `T_c_index`, `T_total_index` | each component / baseline total revenue × 100 |
| Macro response | `n_index`, `y_index`, `k_index`, `c_index` | each variable / its own baseline × 100 |
| Fiscal response | `g_y`, `s_y`, `T_total_y` | ratios to current output, not indices |

A point with `valid=false` is rendered as `null`, with `connectgaps=false`; the frontend performs no interpolation or economic calculation.
