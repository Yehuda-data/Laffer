"""FastAPI application for the Trabandt–Uhlig Laffer Curve Lab."""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from .dashboard_model_service import (
    baseline,
    capital_curve,
    compare,
    equations,
    labor_curve,
    presets,
    sensitivity,
)
from .schemas import CompareRequest, FiscalClosure, ModelSpecification, SensitivityRequest

APP_ROOT = Path(__file__).resolve().parents[1]
FRONTEND = APP_ROOT / "frontend"

app = FastAPI(
    title="Trabandt–Uhlig Laffer Curve Lab API",
    version="1.0.0",
    description="A thin research-dashboard adapter over the existing Python model engines.",
)


@app.get("/api/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/api/presets")
def get_presets():
    return presets()


@app.post("/api/baseline")
def post_baseline(specification: ModelSpecification):
    return baseline(specification)


@app.post("/api/laffer/labor")
def post_labor_curve(specification: ModelSpecification):
    return labor_curve(specification)


@app.post("/api/laffer/capital")
def post_capital_curve(specification: ModelSpecification):
    try:
        return capital_curve(specification)
    except NotImplementedError as error:
        raise HTTPException(status_code=501, detail=str(error)) from error


@app.post("/api/sensitivity")
def post_sensitivity(request: SensitivityRequest):
    try:
        return sensitivity(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@app.post("/api/compare")
def post_compare(request: CompareRequest):
    return compare(request)


@app.get("/api/equations/{closure}")
def get_equations(closure: FiscalClosure):
    return equations(closure)


@app.get("/assets/plotly.min.js", include_in_schema=False)
def plotly_javascript():
    try:
        import plotly
    except ImportError as error:
        raise HTTPException(status_code=503, detail="The local plotly package is not installed.") from error
    asset = Path(plotly.__file__).resolve().parent / "package_data" / "plotly.min.js"
    if not asset.exists():
        raise HTTPException(status_code=503, detail="Local Plotly JavaScript asset was not found.")
    return FileResponse(asset, media_type="application/javascript")


app.mount("/", StaticFiles(directory=FRONTEND, html=True), name="frontend")
