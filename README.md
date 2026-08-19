# Laffer Curve Israel Dashboard

A local research dashboard for exploring Trabandt–Uhlig (2011) steady-state Laffer curves with Israel and paper/US calibrations.

The repository contains the dashboard plus the minimal model and data files required to run it. Generated research outputs and unrelated source material are intentionally excluded.

## Features

- FastAPI backend with an offline HTML/CSS/JavaScript frontend
- Israel input-assumption and model-implied calibrations
- Paper/US benchmark calibration
- Labor- and capital-tax Laffer curves
- s-Laffer and g-Laffer fiscal closures
- Sensitivity analysis, diagnostics, scenario comparison, and CSV/PNG export

## Setup on Windows

From the repository root:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r .\web_dashboard\requirements.txt
```

## Run

```powershell
Set-Location .\web_dashboard
..\.venv\Scripts\python.exe -m uvicorn backend.app:app --host 127.0.0.1 --port 8000
```

Open <http://127.0.0.1:8000/>. API documentation is available at <http://127.0.0.1:8000/docs>.

## GitHub Pages

The same interface can run entirely in the browser without a Python server. On
GitHub Pages, `static-model.js` performs the steady-state calculations locally
and no model inputs are sent to a backend. The Pages workflow publishes
`web_dashboard/frontend` after every push to `main`.

Enable **Settings → Pages → Source: GitHub Actions** in the GitHub repository.
The public URL then has the form `https://<account>.github.io/<repository>/`.

## Test

From the repository root:

```powershell
.\.venv\Scripts\python.exe -m pytest .\web_dashboard\tests -q
```

See [`web_dashboard/docs/README.md`](web_dashboard/docs/README.md) for model architecture, calibration details, API routes, diagnostics, and known limitations.
