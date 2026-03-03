#!/bin/bash
#SBATCH --job-name=ring_trial_5a_adduct
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

date
PTH=$(pwd)
cd "$PTH"

mkdir -p logs

echo "[START] $(date '+%Y-%m-%d %H:%M:%S')"
echo "JobID: ${SLURM_JOB_ID:-NA}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-NA}"
echo "Node: ${SLURMD_NODENAME:-NA}"
echo "PWD: $(pwd)"

/shared/bioinf/R/bin/R-4.5-BioC3.22 -e "quarto::quarto_render('5_downstream_analysis/5a_adduct_method_development.qmd')"

echo "[END] $(date '+%Y-%m-%d %H:%M:%S')"
date
