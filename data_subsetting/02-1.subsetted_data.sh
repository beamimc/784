#!/bin/bash
#SBATCH --output=/work/users/b/e/beacm/784/784/02.subsetted_data_%j.out
#SBATCH --error=/work/users/b/e/beacm/784/784/02.subsetted_data_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:10:00

DIR="/work/users/b/e/beacm/784/784"

# eQTL data subset by gene:
## Note: Gene 1. SORT1:
awk -F'\t' 'NR==1 || $13=="SORT1"' "$DIR/eqtl_liver_sort1_region.tsv" > "$DIR/eqtl_liver_sort1_only.tsv"
## Note: Gene 2. HMGCR:
awk -F'\t' 'NR==1 || $13=="HMGCR"' "$DIR/eqtl_liver_hmgcr_region.tsv" > "$DIR/eqtl_liver_hmgcr_only.tsv"

wc -l "$DIR"/eqtl_liver_*_only.tsv