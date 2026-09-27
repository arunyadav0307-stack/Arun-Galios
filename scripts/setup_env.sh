#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# P0 environment bootstrap for the Arun-Galios execution plan.
#
# Installs the pinned dependency set from requirements.txt. Idempotent: safe
# to re-run; it is a no-op when the pinned versions are already importable.
#
# WHY THIS SCRIPT EXISTS (P0 finding, see docs/P0_REPORT.md §B.3):
#   * Debian 12 enforces PEP 668 (externally-managed-environment), so a bare
#     `pip install -r requirements.txt` is REFUSED with
#     "error: externally-managed-environment".
#   * The sandbox image does not ship sympy / numpy / matplotlib / pytest;
#     they were installed at runtime into /usr/local/lib/python3.11/
#     dist-packages, which is OUTSIDE the persisted workspace root
#     (/home/user). A fresh sandbox therefore starts WITHOUT them.
#   Consequence: environment reproducibility is NOT guaranteed by cloning the
#   repository alone. Run this script first.
#
# Usage:  bash scripts/setup_env.sh
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REQ="$REPO_ROOT/requirements.txt"

echo "== Arun-Galios environment bootstrap =="
echo "repo    : $REPO_ROOT"
echo "python  : $(python3 -V 2>&1)"
echo "pep668  : detecting whether pip is externally managed"

# --- Detect PEP 668 -------------------------------------------------------
PEP668=0
if python3 -c "import sysconfig,os;sys.exit(0 if os.path.exists(sysconfig.get_path('stdlib')+'/EXTERNALLY-MANAGED') else 1)"; then
  PEP668=1
fi

# Prefer a venv when the venv module is usable; fall back to
# --break-system-packages (how the P0 baseline was actually installed).
if [ "$PEP668" -eq 1 ]; then
  echo "          PEP 668 active (EXTERNALLY-MANAGED present)"
  if python3 -m venv --help >/dev/null 2>&1; then
    VENV="$REPO_ROOT/.venv"
    echo "          creating/using venv at $VENV"
    python3 -m venv "$VENV" || { echo "venv creation failed; falling back"; VENV=""; }
    if [ -n "${VENV:-}" ]; then
      # shellcheck disable=SC1091
      source "$VENV/bin/activate"
      python3 -m pip install --quiet --upgrade pip
      python3 -m pip install --quiet -r "$REQ"
      echo "OK: installed into $VENV"
      echo "    activate with:  source $VENV/bin/activate"
      python3 -m pip list --format=columns 2>/dev/null | grep -Ei 'sympy|numpy|matplotlib|pymupdf|pytest|mpmath' || true
      exit 0
    fi
  fi
  echo "          falling back to --break-system-packages (matches the P0 baseline)"
  python3 -m pip install --quiet --break-system-packages -r "$REQ"
else
  python3 -m pip install --quiet -r "$REQ"
fi

echo "OK: dependencies installed."
echo
echo "== Verified import check =="
python3 - <<'PY'
import importlib, importlib.metadata as md
req = ["sympy", "numpy", "matplotlib", "pymupdf", "pytest", "mpmath"]
bad = []
for p in req:
    try:
        importlib.import_module(p)
        print(f"  {p:12s} {md.version(p):10s} import OK")
    except Exception as exc:
        bad.append(p)
        print(f"  {p:12s} {'':10s} FAILED: {type(exc).__name__}")
raise SystemExit(1 if bad else 0)
PY
