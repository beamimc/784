#!/bin/bash
#SBATCH --output=/work/users/b/e/beacm/784/784/03.ld_ref_%j.out
#SBATCH --error=/work/users/b/e/beacm/784/784/03.ld_ref_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00

# Downloads the MRCIEU / OpenGWAS 1000 Genomes LD reference panel (hg19, plink format)
# and keeps the EUR files. Used by run_coloc_susie.R via ieugwasr::ld_matrix_local().
# Needs internet.

DIR="/work/users/b/e/beacm/784/784"
mkdir -p "$DIR/1kg"
cd "$DIR/1kg"

[ -s 1kg.v3.tgz ] || wget -c http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz

# Extract only the EUR files (archive also has AFR, AMR, EAS, SAS)
tar -xzvf 1kg.v3.tgz --wildcards '*EUR.*'

# Move them to 1kg/ in case the archive puts them in a subfolder
find . -mindepth 2 -name 'EUR.*' -exec mv {} . \;

ls -lh EUR.*
wc -l EUR.fam
