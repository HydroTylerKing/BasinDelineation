#!/bin/bash
# Basin Delineator - macOS / Linux launcher. Double-click (macOS) or run ./start.command
cd "$(dirname "$0")"
if [ ! -x .venv/bin/python ]; then
  echo "First run: setting up Python environment ..."
  python3 -m venv .venv || { echo "Install Python 3.10+ first."; exit 1; }
  .venv/bin/python -m pip install --upgrade pip
  .venv/bin/python -m pip install -r requirements.txt
fi
exec .venv/bin/python basin_app.py "$@"
