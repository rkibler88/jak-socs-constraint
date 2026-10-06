#!/usr/bin/env bash
#
# JAK-SOCS evolutionary constraint and specificity analysis - full pipeline.
#
# Usage:
#   ./run_pipeline.sh                  # full run, writes to results/
#   ./run_pipeline.sh --from 05        # resume from stage 05
#   ./run_pipeline.sh --only 08        # run a single stage
#   ./run_pipeline.sh --quick          # reduced permutation/replicate counts (~5 min)
#
# Stages 01-02 query UniProt and take ~50 min in total (14,338 sequence
# fetches); everything downstream is minutes. Stage 07 downloads two mmCIF
# files from RCSB. All other stages are offline.
#
# Every stage is idempotent and writes only into results/ and data/.
# Random seed is fixed at 42 inside each script.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

WORKDIR="work/analysis"     # stages read/write here using short internal names
OUTDIR="$WORKDIR"          # published tree is produced by stage 12
SRC="src"
THREADS="${THREADS:-8}"

N_PERM=10000        # label permutations (specificity, motifs, interface)
N_BOOT=1000         # sequence bootstraps (constraint CIs)
N_SCRAMBLE=200      # species-scrambled replicates (coupling null)

FROM="01"
ONLY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from)  FROM="$2"; shift 2 ;;
    --only)  ONLY="$2"; shift 2 ;;
    --quick) N_PERM=1000; N_BOOT=200; N_SCRAMBLE=20; shift ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

mkdir -p "$WORKDIR" results/tables results/figures results/logs data/sequences data/alignments work

# Seed the working directory from the shipped published tree so that any single
# stage can be re-run offline without repeating stages 01-02.
if [[ "$FROM" != "01" || -n "$ONLY" ]]; then
  echo "== seeding $WORKDIR from published tree =="
  python "$SRC/12_package.py" --unpack --workdir "$WORKDIR" --root .
fi

# --------------------------------------------------------------- preflight

require() {
  command -v "$1" >/dev/null 2>&1 || { echo "MISSING required tool: $1" >&2; exit 1; }
}
echo "== preflight =="
require python
require mafft
require mmseqs
python - <<'PY'
import sys
missing = []
for m in ("numpy", "scipy", "pandas", "matplotlib", "Bio", "statsmodels"):
    try:
        __import__(m)
    except ImportError:
        missing.append(m)
if missing:
    sys.exit("MISSING required python packages: " + ", ".join(missing))
print("  python packages OK")
PY
echo "  mafft   $(mafft --version 2>&1 | head -1)"
echo "  mmseqs  $(mmseqs version 2>&1 | head -1)"
echo "  threads $THREADS ; permutations $N_PERM ; bootstraps $N_BOOT ; scrambles $N_SCRAMBLE"

# ------------------------------------------------------------------ driver

run_stage() {
  local id="$1"; shift
  local name="$1"; shift
  if [[ -n "$ONLY" && "$ONLY" != "$id" ]]; then return 0; fi
  if [[ -z "$ONLY" && "$id" < "$FROM" ]]; then
    echo "-- [$id] $name (skipped)"
    return 0
  fi
  echo ""
  echo "== [$id] $name =="
  local log="results/logs/stage${id}.log"
  if "$@" 2>&1 | tee "$log"; then
    echo "-- [$id] done (log: $log)"
  else
      echo "!! [$id] FAILED - see $log" >&2
      exit 1
  fi
}

PY="python"

run_stage 01 "depth-maximisation scoping survey (UniProt, ~5 min)" \
  $PY "$SRC/01_depth_scoping.py" --outdir "$OUTDIR"

run_stage 02 "paralogue-resolved retrieval and assignment (~45 min)" \
  $PY "$SRC/02_retrieval.py" --scoping-raw "$OUTDIR/depth_scoping_raw.csv" \
      --outdir "$OUTDIR" --threads "$THREADS"

run_stage 03 "per-family alignment and coordinate mapping" \
  $PY "$SRC/03_align.py" --outdir "$OUTDIR" --threads "$THREADS"

run_stage 04 "domain extraction and column provenance" \
  $PY "$SRC/04_provenance.py" --outdir "$OUTDIR"

run_stage 05 "per-residue constraint" \
  $PY "$SRC/05_constraint.py" --outdir "$OUTDIR" --n-boot "$N_BOOT"

run_stage 06 "specificity-determining positions" \
  $PY "$SRC/06_specificity.py" --outdir "$OUTDIR" --n-perm "$N_PERM"

run_stage 07 "structural interfaces from 4GL9 and 6C7Y" \
  $PY "$SRC/07_structures.py" --outdir "$OUTDIR" --cifdir work/structures

run_stage 08 "interface constraint test against matched controls" \
  $PY "$SRC/08_interface_test.py" --outdir "$OUTDIR" --n-perm "$N_PERM"

run_stage 09 "declared motif tests" \
  $PY "$SRC/09_motif_tests.py" --outdir "$OUTDIR" --n-perm "$N_PERM"

run_stage 10 "coupling depth gate and species-scrambled null" \
  $PY "$SRC/10_coupling_power.py" --outdir "$OUTDIR" --n-rep "$N_SCRAMBLE"

run_stage 11 "cross-paralogue consistency and meta-analysis" \
  $PY "$SRC/11_consistency.py" --outdir "$OUTDIR"

run_stage 12 "package outputs into the published tree" \
  $PY "$SRC/12_package.py" --pack --workdir "$WORKDIR" --root .

echo ""
echo "== pipeline complete =="
echo "tables   results/tables   (manifest: results/tables/MANIFEST.csv)"
echo "data     data/sequences, data/alignments, data/reference"
echo "figures  results/figures"
echo "logs     results/logs"
