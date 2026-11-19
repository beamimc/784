#!/bin/bash
#SBATCH --output=/work/users/b/e/beacm/784/784/subsetted_data_%j.out
#SBATCH --error=/work/users/b/e/beacm/784/784/subsetted_data_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00

module load samtools

INDIR="/work/users/b/e/beacm/784/784"
OUTDIR="/work/users/b/e/beacm/784/784"

# GWAS data subset:
## This code is for two separate genes:
### Note: Gene 1. SORT1 (chr1, minus strand), hg19 TSS = 109940540, +- 1MB:
{ zcat "$INDIR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" | head -n1; tabix "$INDIR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" 1:108940540-110940540; } > "$OUTDIR/gwas_ldl_sort1_region.tsv"
### Note: Gene 2. HMGCR (chr5, plus strand), hg19 TSS = 74632993, +- 1MB:
{ zcat "$INDIR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" | head -n1; tabix "$INDIR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" 5:73632993-75632993; } > "$OUTDIR/gwas_ldl_hmgcr_region.tsv"

# eQTL data subset:
## Note: Gene 1. SORT1, chr1:
tar -xzOf "$INDIR/Liver_eQTL_Meta_Marginal_20240530.tar.gz" "./chr1_marginal_summary_results.tsv" | awk -F'\t' 'NR==1 || ($4>=108940540 && $4<=110940540)' > "$OUTDIR/eqtl_liver_sort1_region.tsv"
## Note: Gene 2. HMGCR, chr5:
tar -xzOf "$INDIR/Liver_eQTL_Meta_Marginal_20240530.tar.gz" "./chr5_marginal_summary_results.tsv" | awk -F'\t' 'NR==1 || ($4>=73632993 && $4<=75632993)' > "$OUTDIR/eqtl_liver_hmgcr_region.tsv"

wc -l "$OUTDIR"/*.tsv