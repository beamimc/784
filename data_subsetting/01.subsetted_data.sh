#!/bin/bash
#SBATCH --output=/work/users/m/e/melange/bcb784_proj/results/slurm/subsetted_data_%j.out
#SBATCH --error=/work/users/m/e/melange/bcb784_proj/results/slurm/subsetted_data_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00

module load samtools

INDIR="/work/users/m/e/melange/bcb784_proj/data/original_data"
OUTDIR="/work/users/m/e/melange/bcb784_proj/data/subsetted_data"

# GWAS data subset:
## This code is for two separate genes:
### Note: Gene 1. Change "2:20266945-22266945" to your gene location with +- 1MB of the TSS. The format is [chr]:[start_position-end_positon]. You'll also need to change the output file name:
{ zcat "$INDIR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" | head -n1; tabix "$INDIR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" 2:20266945-22266945; } > "$OUTDIR/gwas_ldl_apob_region.tsv"
### Note: Gene 2. Change "19:44409053-46409053" to your gene location with +- 1MB of the TSS. The format is [chr]:[start_position-end_positon]. You'll also need to change the output file name:
{ zcat "$INDIR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" | head -n1; tabix "$INDIR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" 19:44409053-46409053; } > "$OUTDIR/gwas_ldl_apoe_region.tsv"

# eQTL data subset:
## Note: Gene 1. You'll need to change the chromosome number depending on the chromosome your gene is on for file "./chr2_marginal_summary_results.tsv". $4>=20266945 && $4<=22266945 are the start and stop position of the gene (use the same start and stop position you used for subsetting the GWAS data). You'll also need to change the output file name:
tar -xzOf "$INDIR/Liver_eQTL_Meta_Marginal_20240530.tar.gz" "./chr2_marginal_summary_results.tsv" | awk -F'\t' 'NR==1 || ($4>=20266945 && $4<=22266945)' > "$OUTDIR/eqtl_liver_apob_region.tsv"
## Note: Gene 2. You'll need to change the chromosome number depending on the chromosome your gene is on for file "./chr19_marginal_summary_results.tsv". $4>=44409053 && $4<=46409053 are the start and stop position of the gene (use the same start and stop position you used for subsetting the GWAS data). You'll also need to change the output file name:
tar -xzOf "$INDIR/Liver_eQTL_Meta_Marginal_20240530.tar.gz" "./chr19_marginal_summary_results.tsv" | awk -F'\t' 'NR==1 || ($4>=44409053 && $4<=46409053)' > "$OUTDIR/eqtl_liver_apoe_region.tsv"

wc -l "$OUTDIR"/*.tsv
