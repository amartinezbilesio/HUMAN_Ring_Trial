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

# Ensure R from the requested installation is used
export PATH=/shared/bioinf/R/bin:$PATH
RSCRIPT=${RSCRIPT:-/shared/bioinf/R/bin/Rscript}

echo "Using Rscript: $(which $RSCRIPT)"

# Run the consensus builder. Adjust --cores and --chunk_size as needed.
"$RSCRIPT" 5_downstream_analysis/scripts/build_consensus_no_bin.R \
  --input_dir 5_downstream_analysis/object \
  --out_prefix 5_downstream_analysis/object/consensus_all \
  --cores 10 \
  --chunk_size 5000

date
