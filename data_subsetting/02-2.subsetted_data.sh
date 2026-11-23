#!/bin/bash
#SBATCH --job-name=subset_genes
#SBATCH --output=/hickory/users/a/b/abskin/subsetted_genes_%j.out
#SBATCH --error=/hickory/users/a/b/abskin/subsetted_genes_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:10:00

DIR="/hickory/users/a/b/abskin/subsetted_data"

# eQTL data subset by gene:
## Gene 1. FADS1:
awk -F'\t' 'NR==1 || $13=="FADS1"' "$DIR/eqtl_liver_fads1_region.tsv" > "$DIR/eqtl_liver_fads1_only.tsv"
## Gene 2. LDLR:
awk -F'\t' 'NR==1 || $13=="LDLR"' "$DIR/eqtl_liver_ldlr_region.tsv" > "$DIR/eqtl_liver_ldlr_only.tsv"

wc -l "$DIR"/eqtl_liver_*_only.tsv
