"""Python port of LafferCodeWeb/1_BaselineCode/data.m

The numeric matrices are not re-typed by hand: they are parsed directly out of
the original MATLAB source file, so the data are byte-for-byte the researchers'
data.  Rows are countries (16), columns are years 1995..2007.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

MATLAB_DIR = Path(__file__).resolve().parent.parent / "LafferCodeWeb" / "1_BaselineCode"
DATA_M = MATLAB_DIR / "data.m"

FIRST_YEAR = 1995
LAST_YEAR = 2007

COUNTRY = [
    "GER", "FRA", "ITA", "GBR", "AUT", "BEL", "DNK", "FIN",
    "GRE", "IRL", "NET", "PRT", "ESP", "SWE", "USA", "EU-14",
]

# MATLAB 1-based country indices used throughout the original code.
IDX_USA = COUNTRY.index("USA")      # MATLAB 15 -> python 14
IDX_EU14 = COUNTRY.index("EU-14")   # MATLAB 16 -> python 15

# Order in which get_further_params.m loops over countries.  USA must come
# first because every other country inherits the US value of kappa.
COUNTRY_SELECTION = [IDX_USA, IDX_EU14] + list(range(0, 14))

_MATRIX_NAMES = [
    "ky", "sy", "by", "xy_priv", "xy_gov", "cy", "govconsy", "tby", "n",
    "tauk", "taun", "tauc", "captaxrevy", "labtaxrevy", "constaxrevy",
    "gov_interest_y",
]


def _strip_comments(text: str) -> str:
    return "\n".join(line.split("%")[0] for line in text.splitlines())


def _parse_matlab_matrices(path: Path) -> dict[str, np.ndarray]:
    src = _strip_comments(path.read_text(encoding="latin-1"))
    out: dict[str, np.ndarray] = {}
    for name in _MATRIX_NAMES:
        m = re.search(rf"^\s*{re.escape(name)}\s*=\s*\[(.*?)\]\s*;", src,
                      re.DOTALL | re.MULTILINE)
        if m is None:
            raise ValueError(f"matrix '{name}' not found in {path}")
        rows = []
        for line in m.group(1).splitlines():
            tokens = line.split()
            if not tokens:
                continue
            rows.append([np.nan if t.lower() == "nan" else float(t) for t in tokens])
        arr = np.array(rows, dtype=float)
        if arr.shape != (len(COUNTRY), LAST_YEAR - FIRST_YEAR + 1):
            raise ValueError(f"matrix '{name}' has shape {arr.shape}, expected "
                             f"{(len(COUNTRY), LAST_YEAR - FIRST_YEAR + 1)}")
        out[name] = arr
    return out


_CACHE: dict[str, np.ndarray] | None = None


def raw_data() -> dict[str, np.ndarray]:
    global _CACHE
    if _CACHE is None:
        _CACHE = _parse_matlab_matrices(DATA_M)
    return _CACHE


def data(startdate: int = FIRST_YEAR, enddate: int = LAST_YEAR) -> dict[str, np.ndarray]:
    """Equivalent of data.m: return the sub-sample [startdate, enddate]."""
    d = raw_data()
    lo = startdate - FIRST_YEAR
    hi = enddate - FIRST_YEAR + 1
    return {k: v[:, lo:hi] for k, v in d.items()}
