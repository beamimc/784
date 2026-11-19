#!/bin/bash
#SBATCH --output=/work/users/m/e/melange/bcb784_proj/results/slurm/download_%j.out
#SBATCH --error=/work/users/m/e/melange/bcb784_proj/results/slurm/download_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=16:00:00

set -euo pipefail

OUTDIR="/work/users/m/e/melange/bcb784_proj/data"

cd "$OUTDIR"

wget -c --tries=10 "https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"

wget -c --tries=10 "https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.tbi"

wget -c --tries=10 "ftp://mohlkeanon:anon@rc-ns-ftp.its.unc.edu/Liver_eQTL_Meta_Marginal_20240530.tar.gz"

ls -lh
