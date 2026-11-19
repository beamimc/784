#!/bin/bash
#SBATCH --output=/work/users/m/e/melange/bcb784_proj/results/slurm/02.subsetted_data_%j.out
#SBATCH --error=/work/users/m/e/melange/bcb784_proj/results/slurm/02.subsetted_data_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:10:00

DIR="/work/users/m/e/melange/bcb784_proj/data/subsetted_data"

# eQTL data subset by gene:
## Note: Gene 1. Select only SNPs corresponding to gene. Update "$13=="APOB"" to gene name and "eqtl_liver_apob_region.tsv" to file name. You'll also need to update output file name:
awk -F'\t' 'NR==1 || $13=="APOB"' "$DIR/eqtl_liver_apob_region.tsv" > "$DIR/eqtl_liver_apob_only.tsv"
## Note: Gene 2. Select only SNPs corresponding to gene. Update "$13=="APOE"" to gene name and "eqtl_liver_apoe_region.tsv" to file name. You'll also need to update output file name:
awk -F'\t' 'NR==1 || $13=="APOE"' "$DIR/eqtl_liver_apoe_region.tsv" > "$DIR/eqtl_liver_apoe_only.tsv"

wc -l eqtl_liver_*_only.tsv
