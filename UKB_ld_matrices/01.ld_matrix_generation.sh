#!/bin/bash
#SBATCH --output=/work/users/m/e/melange/bcb784_proj/results/slurm/ld_matrix_generation.R_%j.out
#SBATCH --error=/work/users/m/e/melange/bcb784_proj/results/slurm/ld_matrix_generation.R_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=1:00:00

# ======== Edit your region here (hg19) ========
NAME=APOE
CHR=19
START=44409053
END=46409053
# ==============================================

LD_DIR=/work/users/m/e/melange/bcb784_proj/data/ld/raw_ld     # UKB downloads
OUT_DIR=/work/users/m/e/melange/bcb784_proj/data/ld/ld_matrices     # LD matrix output

module load r/4.4.0
module load anaconda
conda activate base
export RETICULATE_PYTHON=$(which python)

Rscript 01.ld_matrix_generation.R "$CHR" "$START" "$END" "$LD_DIR" "${OUT_DIR}/${NAME}_ld.rds"
