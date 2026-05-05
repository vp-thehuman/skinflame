"""Console entry points."""
from __future__ import annotations

import os
import sys
from pathlib import Path


def run_demo():  # pragma: no cover
    """`skinflame-demo` -> launch the bundled Streamlit demo."""
    here = Path(__file__).resolve().parent
    # demo/streamlit_app.py lives outside the package src tree, so prefer the
    # installed copy if present, else fall back to the repo path.
    candidates = [
        here / "demo" / "streamlit_app.py",
        here.parent.parent / "demo" / "streamlit_app.py",
    ]
    app = next((c for c in candidates if c.exists()), None)
    if app is None:
        sys.exit(
            "Could not locate streamlit_app.py. If you installed from source, run:\n"
            "  streamlit run demo/streamlit_app.py"
        )
    os.execvp("streamlit", ["streamlit", "run", str(app)])
