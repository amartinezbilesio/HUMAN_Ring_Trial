#!/bin/bash
#SBATCH --job-name=consensus_build
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=32000
#SBATCH --partition=slow-mc1
#SBATCH --time=48:00:00
#SBATCH --output=consensus_build_%j.log


date
PTH=$(pwd)
cd "$PTH"

# Try several candidate Rscript locations (check the BioC install too)
export PATH=/shared/bioinf/R/bin:$PATH

# Candidate paths (adjust if your site uses a different tree)
R_CANDIDATES=(
  "/shared/bioinf/R/bin/Rscript-4.5-BioC3.22"
  "/shared/bioinf/R/bin/R-4.5-BioC3.22"
  "/shared/bioinf/R/bin/Rscript"
  "/usr/bin/Rscript"
)

RSCRIPT=""
for p in "${R_CANDIDATES[@]}"; do
  if [ -x "$p" ]; then
    RSCRIPT="$p"
    break
  fi
done

# Try module system if available
if [ -z "$RSCRIPT" ]; then
  if command -v module >/dev/null 2>&1; then
    module load R-4.5-BioC3.22 2>/dev/null || true
  fi
  # fallback to whatever is on PATH
  RSCRIPT=$(command -v Rscript || true)
fi

if [ -z "$RSCRIPT" ]; then
  echo "ERROR: Could not find Rscript. Check the install path or load the R module." >&2
  echo "Tried candidates: ${R_CANDIDATES[*]}" >&2
  echo "You can inspect available R installations with: ls -l /shared/bioinf/R || module avail R" >&2
  exit 1
fi

echo "Using Rscript: $RSCRIPT"
echo "Rscript version: $($RSCRIPT --version 2>&1 | head -n1)"

# Run the consensus builder. Adjust --cores and --chunk_size as needed.
"$RSCRIPT" 5_downstream_analysis/scripts/build_consensus_no_bin.R \
  --input_dir 5_downstream_analysis/object \
  --out_prefix 5_downstream_analysis/object/consensus_all \
  --cores 10 \
  --chunk_size 5000

date
