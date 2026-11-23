#!/bin/bash
#SBATCH --job-name=subset_regions
#SBATCH --output=/hickory/users/a/b/abskin/subsetted_data_%j.out
#SBATCH --error=/hickory/users/a/b/abskin/subsetted_data_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00

module load samtools

INDIR="/hickory/users/a/b/abskin"
OUTDIR="/hickory/users/a/b/abskin/subsetted_data"
mkdir -p "$OUTDIR"

GWAS_FILE="$INDIR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"
EQTL_FILE="$INDIR/Liver_eQTL_Meta_Marginal_20240530.tar.gz"

# GWAS data subset (+-1 Mb around the hg19 TSS). The header comes from
# `zcat | head -n1` because the GLGC header line doesn't start with "#",
# so `tabix -H` prints nothing.
## Gene 1. FADS1 (chr11, minus strand), hg19 TSS = 61596790:
{ zcat "$GWAS_FILE" | head -n1; tabix "$GWAS_FILE" 11:60596790-62596790; } > "$OUTDIR/gwas_ldl_fads1_region.tsv"
## Gene 2. LDLR (chr19, plus strand), hg19 TSS = 11200038:
{ zcat "$GWAS_FILE" | head -n1; tabix "$GWAS_FILE" 19:10200038-12200038; } > "$OUTDIR/gwas_ldl_ldlr_region.tsv"

# eQTL data subset (same windows):
## Gene 1. FADS1, chr11:
tar -xzOf "$EQTL_FILE" "./chr11_marginal_summary_results.tsv" | awk -F'\t' 'NR==1 || ($4>=60596790 && $4<=62596790)' > "$OUTDIR/eqtl_liver_fads1_region.tsv"
## Gene 2. LDLR, chr19:
tar -xzOf "$EQTL_FILE" "./chr19_marginal_summary_results.tsv" | awk -F'\t' 'NR==1 || ($4>=10200038 && $4<=12200038)' > "$OUTDIR/eqtl_liver_ldlr_region.tsv"

wc -l "$OUTDIR"/*.tsv
