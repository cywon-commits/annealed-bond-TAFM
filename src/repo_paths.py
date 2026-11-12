"""Repository-relative paths shared by the analysis and figure scripts."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGDATA = DATA / "figdata"
OUT = Path(os.environ.get("OUT_DIR", ROOT / "output"))
FIGS = OUT / "figs"
FIGS.mkdir(parents=True, exist_ok=True)
