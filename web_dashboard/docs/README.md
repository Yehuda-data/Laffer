# Trabandt–Uhlig Laffer Curve Lab

A local, research-oriented dashboard for the Trabandt and Uhlig (2011) steady-state Laffer-curve model. The dashboard is a new adapter and user interface; it does not replace or modify the existing economic engines.

## Architecture

```text
Existing Python engines
  python_port/laffer_model.py       (s-Laffer)
  python_port/laffer_model_g.py     (g-Laffer)
            ↓
web_dashboard/backend/dashboard_model_service.py
  ratio-to-level calibration, serialization, diagnostics
            ↓
FastAPI / Pydantic REST API
            ↓
HTML + CSS + vanilla JavaScript + local Plotly.js
```

Economic equilibrium calculations remain in `steady()` and `steady_g()`. JavaScript only collects inputs, calls the API, renders returned values and creates browser-side downloads.

## Install on Windows

Open a PowerShell terminal in `סקירת ספרות/uhlig2011` and use the existing project virtual environment:

```powershell
.\.venv\Scripts\python.exe -m pip install -r .\web_dashboard\requirements.txt
```

Do not use global `pip`. The dashboard requirements reuse `../requirements.txt` and add FastAPI, Uvicorn, Pydantic, Plotly, pytest and httpx.

## Start the dashboard

From `סקירת ספרות/uhlig2011/web_dashboard`:

```powershell
..\.venv\Scripts\python.exe -m uvicorn backend.app:app --host 127.0.0.1 --port 8000
```

Open <http://127.0.0.1:8000/>. Interactive API documentation is available at <http://127.0.0.1:8000/docs>.

Plotly JavaScript is served from the installed local Python package at `/assets/plotly.min.js`; no internet connection or CDN is required.

## API

- `GET /api/presets`
- `POST /api/baseline`
- `POST /api/laffer/labor`
- `POST /api/laffer/capital`
- `POST /api/sensitivity`
- `POST /api/compare`
- `GET /api/equations/{s_laffer|g_laffer}`
- `GET /api/health`

Each model request contains an explicit complete specification. Responses include submitted inputs, parameter values and statuses, equilibrium/curve data, diagnostics and validity metadata. Non-finite numbers become JSON `null`; economically invalid values are not clipped or repaired.

## Fiscal closures

### s-Laffer

- Input: baseline `g/y`.
- Fixed throughout a curve: government spending `g` in **levels**.
- Endogenous: transfers `s`, obtained from the government-budget residual.
- Labor- and capital-tax curves are supported.

### g-Laffer

- Input: baseline `s/y`.
- Fixed throughout a curve: transfers `s` in **levels**.
- Endogenous: government spending `g`.
- Labor is solved using equation (19), as implemented in `laffer_model_g.py`.
- Only the baseline-connected valid branch is plotted. Raw isolated roots remain visible through validity metadata.

The UI disables the irrelevant fiscal input. It never fixes both fiscal instruments independently.

## Calibration methods

### External parameters

`theta` and `delta` are user inputs. The user can enter `kappa` or provide `n_target`; in the latter case the existing CFE calibration equation is used and `kappa` is labelled `CALIBRATED`.

### Model-implied parameters

The user enters `k/y` and `x/y`. The service applies:

```text
delta = (x/y)/(k/y) - (psi - 1)
theta = (k/y) * ((R - 1)/(1 - tau_k) + delta)
```

It then calibrates `kappa` to the labor target. `theta` and `delta` are read-only in the UI and labelled `MODEL-IMPLIED`.

## Presets

- **Israel — Input assumptions** uses `load_inputs()` and `calibrate_variant1()` from `israel_two_calibrations.py`.
- **Israel — Model-implied** uses `load_inputs()` and `calibrate_variant2()`.
- **Paper / US benchmark** uses `params()` and the original US row from `laffer_data.py`.
- **Custom** retains the current editable form.

The current Israel workflow uses `tau_c=0.18`, `tau_n=0.28`, and `tau_k=0.30`. Older provisional values in `israel_laffer_2005_2019.py` are intentionally not used for these presets.

## Diagnostics

Diagnostics are classified as `INFO`, `WARNING`, or `INVALID EQUILIBRIUM`. Configured checks include unusual `k/y` and `x/y`, depreciation and capital-share ranges, labor admissibility, negative consumption/government spending, near-zero consumption, solver failures, `R < psi`, and `beta = psi^eta/R > 1`.

For g-Laffer curves the response also reports:

- baseline-connected valid tax range;
- first tax with `g <= 0`;
- first tax with `n >= 1`;
- first tax with `c <= 0`;
- solver-failure rates;
- disconnected raw roots;
- adjacent labor jumps larger than `0.05`, labelled as a possible branch/root switch.

Invalid regions are returned as chart gaps. The service does not interpolate across them.

## Exports

The browser can download baseline, full curve, sensitivity and comparison CSV files. Curve exports retain raw equilibrium values, validity flags and invalidity reasons. Plotly exports the active Laffer chart as PNG.

## Tests

From `web_dashboard`:

```powershell
..\.venv\Scripts\python.exe -m pytest tests -q
```

The suite checks the exact US result (`peak tau_n=0.633`, normalized maximum revenue `129.606122`), direct equality with existing s/g engines, model-implied formulas, trade-sign conversion, preservation of invalid calibrations, closure level invariants, unsupported capital g-Laffer behavior, API contracts and offline static assets.

## Known limitations

1. Capital-tax g-Laffer is not implemented in the existing model and is therefore deliberately unavailable. The API returns HTTP 501: `Capital-tax g-Laffer not yet implemented.`
2. `steady_g()` currently supports CFE preferences only. The dashboard exposes CFE with `phi` mapped to the existing `FRISCH` argument.
3. The existing solvers return the first converged root from their fixed initial-guess sequence. The dashboard diagnoses suspicious jumps and filters to the baseline-connected branch, but does not invent a global root-selection algorithm.
4. Sensitivity runs can be computationally expensive because every displayed point is recomputed by the Python engine.
5. The Israel model-implied preset uses the existing externally supplied `k/y=1.6`; the workbook itself does not contain a capital-stock series.
