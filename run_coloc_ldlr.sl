#!/bin/bash
#SBATCH --job-name=coloc_ldlr
#SBATCH --output=/hickory/users/a/b/abskin/coloc_ldlr_%j.out
#SBATCH --error=/hickory/users/a/b/abskin/coloc_ldlr_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
# =============================================================================
# LDLR: coloc.abf -> coloc + SuSiE -> plots.   Submit with:  sbatch run_ldlr.sl
# Uses the group's analysis (github.com/beamimc/784). Before the first run:
#   sbatch get_LD_ref.sh   (downloads the 1000G EUR LD panel; once for all genes)
# =============================================================================
GENE=ldlr
cd /hickory/users/a/b/abskin          # folder holding the .R scripts
module load r
module load plink/1.90b6.21   # must be PLINK 1.9: the default (plink/2.00) is named plink2 and lacks --r square
which plink || { echo "plink not found"; exit 1; }
Rscript -e 'for (p in c("coloc","susieR","ieugwasr","data.table","ggplot2","patchwork")) if (!requireNamespace(p, quietly = TRUE)) stop("R package missing: ", p)' || exit 1

echo "===== 1. coloc.abf ====="
Rscript run_coloc_Abbey.R $GENE
echo "===== 2. coloc + SuSiE ====="
Rscript run_coloc_susie.R $GENE
echo "===== 3. plots ====="
Rscript plot_locus.R $GENE
