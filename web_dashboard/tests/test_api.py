from __future__ import annotations

from fastapi.testclient import TestClient

from backend.app import app

client = TestClient(app)


def test_health_and_static_frontend():
    assert client.get("/api/health").json() == {"status": "ok"}
    index = client.get("/")
    assert index.status_code == 200
    assert "Trabandt–Uhlig Laffer Curve Lab" in index.text
    assert client.get("/styles.css").status_code == 200
    assert client.get("/app.js").status_code == 200


def test_local_plotly_asset():
    response = client.get("/assets/plotly.min.js")
    assert response.status_code == 200
    assert "plotly" in response.text[:1000].lower()


def test_presets_and_baseline_endpoint():
    presets = client.get("/api/presets")
    assert presets.status_code == 200
    items = presets.json()["presets"]
    assert len(items) == 3
    response = client.post("/api/baseline", json=items[1]["specification"])
    assert response.status_code == 200
    payload = response.json()
    assert payload["validity"]["valid"]
    assert {"inputs", "parameters", "equilibrium", "diagnostics", "validity"}.issubset(payload)


def test_labor_endpoint_and_equations():
    spec = client.get("/api/presets").json()["presets"][2]["specification"]
    spec["grid_step"] = 0.05
    response = client.post("/api/laffer/labor", json=spec)
    assert response.status_code == 200
    payload = response.json()
    assert payload["curve"]
    assert "peak_tax" in payload["summary"]
    equations = client.get("/api/equations/s_laffer")
    assert equations.status_code == 200
    assert equations.json()["equations"]


def test_capital_g_endpoint_returns_501():
    spec = client.get("/api/presets").json()["presets"][1]["specification"]
    spec["closure"] = "g_laffer"
    response = client.post("/api/laffer/capital", json=spec)
    assert response.status_code == 501
    assert response.json()["detail"] == "Capital-tax g-Laffer not yet implemented."


def test_schema_rejects_fixing_wrong_closure_without_required_field():
    spec = client.get("/api/presets").json()["presets"][0]["specification"]
    spec["closure"] = "g_laffer"
    spec["s_y"] = None
    response = client.post("/api/baseline", json=spec)
    assert response.status_code == 422
