#!/bin/bash
#SBATCH --job-name=coloc_fads1
#SBATCH --output=/hickory/users/a/b/abskin/coloc_fads1_%j.out
#SBATCH --error=/hickory/users/a/b/abskin/coloc_fads1_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=4:00:00
# =============================================================================
# FADS1: coloc.abf -> coloc + SuSiE (UK Biobank LD) -> plots.   Submit:  sbatch run_coloc_fads1.sl
# UK Biobank LD comes from UKB_ld_matrices/01.ld_matrix_generation.R (downloads
# the 3 Mb UKB window the first time, ~GBs, so the first run is slower).
# =============================================================================
GENE=fads1
cd /hickory/users/a/b/abskin          # folder with the .R scripts and UKB_ld_matrices/
export COLOC_DIR=/hickory/users/a/b/abskin/subsetted_data   # your data folder (used by run_coloc_susie.R and plot_locus.R)
module load r
module load anaconda                   # Python + numpy/scipy for reading the UKB .npz files
conda activate base
export RETICULATE_PYTHON=$(which python)

python -c "import numpy, scipy" || { echo "numpy/scipy not available in $(which python)"; exit 1; }
Rscript -e 'for (p in c("coloc","susieR","data.table","ggplot2","patchwork","reticulate")) if (!requireNamespace(p, quietly = TRUE)) stop("R package missing: ", p)' || exit 1

echo "===== 1. coloc.abf ====="
Rscript run_coloc_Abbey.R $GENE
echo "===== 2. coloc + SuSiE (UK Biobank LD) ====="
Rscript run_coloc_susie.R $GENE
echo "===== 3. plots ====="
Rscript plot_locus.R $GENE ukb
