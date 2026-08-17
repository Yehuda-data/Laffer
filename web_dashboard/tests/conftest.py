from __future__ import annotations

from pathlib import Path
import sys

import pytest

DASHBOARD = Path(__file__).resolve().parents[1]
if str(DASHBOARD) not in sys.path:
    sys.path.insert(0, str(DASHBOARD))

from backend.dashboard_model_service import presets
from backend.schemas import ModelSpecification


@pytest.fixture(scope="session")
def preset_specs() -> list[ModelSpecification]:
    return [ModelSpecification.model_validate(item["specification"]) for item in presets()["presets"]]
