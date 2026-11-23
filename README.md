# 784: GWAS–eQTL colocalization and disease gene prioritization
 
Group project: Abbey Skinner, Heming Pei, Bea Campillo, Marisa Lange.
 
We test whether known LDL-cholesterol genes are picked up by colocalizing the **LDL GWAS** (GLGC 2021, European ancestry) with **liver eQTL** data (Broadaway et al. 2024 meta-analysis, marginal results), using **coloc.abf** (one causal variant per trait) and **coloc + SuSiE** (several causal variants, with UK Biobank LD). Each person runs two genes.
 
Everything is **hg19 / GRCh37**.
 
## What's in this repo
 
| Path | What it does |
|---|---|
| `download_data/download.sh` | Downloads the GWAS (+ index) and eQTL raw files |
| `data_subsetting/01*.subsetted_data.sh` | Cuts the GWAS and eQTL files to ±1 Mb around each gene's TSS (one script per person) |
| `data_subsetting/02*.subsetted_data.sh` | Keeps only the eQTL rows for the gene itself |
| `run_coloc_Abbey.R` | coloc.abf for one or more genes (same analysis as `run_coloc_Bea.R`, but takes genes on the command line) |
| `run_coloc_susie.R` | coloc + SuSiE; LD from UK Biobank by default (`LD_SOURCE <- "ukb"`) |
| `UKB_ld_matrices/01.ld_matrix_generation.R` | Downloads/extracts the UK Biobank LD matrix (used by `run_coloc_susie.R`) |
| `plot_locus.R` | Regional Manhattan (GWAS over eQTL, colored by LD) + coloc.abf vs SuSiE plot |
| `run_coloc_fads1.sl`, `run_coloc_ldlr.sl` | SLURM job per gene: coloc.abf → coloc + SuSiE → plots. **Copy one of these for your genes.** |
| `03.get_LD_ref.sh` | 1000 Genomes LD panel (only for the older `LD_SOURCE <- "1kg"` option) |
| `results/` | Results for finished genes |
 
---
 
# How to run your genes (step by step)
 
Do one gene per job. In everything below, replace `X/Y/ONYEN` with your own path (e.g. `/work/users/h/p/hpei`) and `<gene>` with your lowercase gene name (e.g. `pcsk9`).
 
## 0. Before you start: your genes' coordinates
 
For each of your 2 genes you need the **chromosome and hg19 TSS**:
 
1. Go to **https://grch37.ensembl.org**. Use this GRCh37 site, *not* the default ensembl.org, which is GRCh38.
2. Search the gene and open its page. Near the top it says `Location: Chromosome <chr>: <start>-<end> <forward|reverse> strand`.
3. TSS = **start** if forward (+) strand, **end** if reverse (−) strand.
4. Window = **TSS − 1,000,000 to TSS + 1,000,000**.
   Example, FADS1: chr11, reverse strand, TSS 61,596,790, so the window is `11:60596790-62596790`.
| | Gene 1 | Gene 2 |
|---|---|---|
| Gene symbol (e.g. `FADS1`) | | |
| lowercase name for files (e.g. `fads1`) | | |
| Chromosome | | |
| hg19 TSS | | |
| Window start / end | | |
 
You'll need about **15 GB** free: raw data ~5 GB, plus UK Biobank LD ~2–4 GB per gene.
 
## 1. Get the code and make folders (once)
 
```bash
ssh ONYEN@longleaf.unc.edu
cd /work/users/X/Y/ONYEN
git clone https://github.com/beamimc/784.git
mkdir -p 784_data/original_data 784_data/subsetted_data 784_logs
```
 
Keep data in `784_data/`, **outside** the repo folder, so it never ends up on GitHub.
 
## 2. Install R packages and check Python (once)
 
```bash
module load r
Rscript -e 'install.packages(c("coloc","susieR","data.table","ggplot2","patchwork","reticulate"), repos = "https://cloud.r-project.org", lib = Sys.getenv("R_LIBS_USER"))'
Rscript -e 'for (p in c("coloc","susieR","data.table","ggplot2","patchwork","reticulate")) cat(p, requireNamespace(p, quietly = TRUE), "\n")'
module load anaconda
conda activate base            # if this errors: source activate base
python -c "import numpy, scipy; print('numpy/scipy OK')"
```
 
Every package line should say `TRUE`. If the install says the library "is not writable", run `mkdir -p $(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')` and install again.
 
## 3. Get the raw data (once)
 
You need these three files in `784_data/original_data/`:
- `LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz`
- `LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.tbi`
- `Liver_eQTL_Meta_Marginal_20240530.tar.gz`
**Easiest:** ask a teammate who has them and copy (`cp`) them from their folder.
**Or download them:** copy `download_data/download.sh`, change `OUTDIR` and the two `#SBATCH --output/--error` paths to your folders, then `sbatch` it. This can take hours.
 
## 4. Subset the data to your genes
 
Use Abbey's subsetting scripts as the template. They contain the header fix listed in Troubleshooting.
 
```bash
cd /work/users/X/Y/ONYEN/784/data_subsetting
cp 01-2.subsetted_data.sh 01-3.subsetted_data.sh
cp 02-2.subsetted_data.sh 02-3.subsetted_data.sh
```
 
**Edit `01-3.subsetted_data.sh`** (`nano 01-3.subsetted_data.sh`):
- `#SBATCH --output` / `--error`: point them to your `784_logs` folder.
- `INDIR` = `.../784_data/original_data`; `OUTDIR` = `.../784_data/subsetted_data`.
- For each gene:
  - GWAS line: change the region (`11:60596790-62596790`) and the output name (`gwas_ldl_fads1_region.tsv`).
  - eQTL line: change the chromosome file (`./chr11_marginal_summary_results.tsv`), both positions (`$4>=60596790 && $4<=62596790`), and the output name.
**Edit `02-3.subsetted_data.sh`:**
- the log paths, and `DIR` = `.../784_data/subsetted_data`;
- for each gene: the gene symbol (`$13=="FADS1"`) and the input/output file names.
Use the **same lowercase name** in all file names: `gwas_ldl_<gene>_region.tsv` and `eqtl_liver_<gene>_only.tsv`. The later scripts find the files by these names.
 
Run the two scripts in order:
```bash
sbatch 01-3.subsetted_data.sh
squeue -u ONYEN                 # wait until it's gone
sbatch 02-3.subsetted_data.sh
```
 
**Check before going on:**
```bash
cd /work/users/X/Y/ONYEN/784_data/subsetted_data
head -2 gwas_ldl_*_region.tsv | cut -c1-120     # first line must be column names: rsID  CHROM  POS_b37 ...
wc -l eqtl_liver_*_only.tsv                      # must be more than 1 line
```
 
## 5. Point the coloc.abf script at your data (one line)
 
`run_coloc_Abbey.R` has Abbey's data folder written in. In **your copy** (don't push this change), set line 13 to your folder:
 
```bash
cd /work/users/X/Y/ONYEN/784
nano run_coloc_Abbey.R
#   DIR   <- "/work/users/X/Y/ONYEN/784_data/subsetted_data"
```
 
(`run_coloc_susie.R` and `plot_locus.R` don't need editing; they read your folder from `COLOC_DIR` in the job file.)
 
## 6. Make one job file per gene
 
Copy Abbey's FADS1 job file once for each of your genes:
 
```bash
cd /work/users/X/Y/ONYEN/784
cp run_coloc_fads1.sl run_coloc_<gene1>.sl
cp run_coloc_fads1.sl run_coloc_<gene2>.sl
```
 
In each copy, change these 6 lines:
 
| Line in the file | Change to |
|---|---|
| `#SBATCH --job-name=coloc_fads1` | `coloc_<gene>` |
| `#SBATCH --output=/hickory/users/a/b/abskin/coloc_fads1_%j.out` | `/work/users/X/Y/ONYEN/784_logs/coloc_<gene>_%j.out` |
| `#SBATCH --error=/hickory/users/a/b/abskin/coloc_fads1_%j.err` | `/work/users/X/Y/ONYEN/784_logs/coloc_<gene>_%j.err` |
| `GENE=fads1` | `GENE=<gene>` (lowercase, same as your file names) |
| `cd /hickory/users/a/b/abskin` | `cd /work/users/X/Y/ONYEN/784` (your repo clone) |
| `export COLOC_DIR=/hickory/users/a/b/abskin/subsetted_data` | `export COLOC_DIR=/work/users/X/Y/ONYEN/784_data/subsetted_data` |
 
The `784_logs` folder must exist (step 1 made it). If it doesn't, SLURM fails without writing a log.
 
## 7. Run
 
```bash
sbatch run_coloc_<gene1>.sl
sbatch run_coloc_<gene2>.sl
squeue -u ONYEN          # PD = waiting, R = running, gone = finished
```
 
The first run downloads that gene's UK Biobank LD window (several GB) into `subsetted_data/ukb_ld/`. Expect a few minutes to about an hour per gene.
 
## 8. Check it worked
 
```bash
cd /work/users/X/Y/ONYEN
grep -iE "error|skipped|did NOT|retrying|UKB LD file|LD consistency|Removed" 784_logs/coloc_*.out 784_logs/coloc_*.err
ls 784_data/subsetted_data/coloc_susie_results_ukb/ 784_data/subsetted_data/plots/
```
 
| Log message | Meaning |
|---|---|
| `UKB LD file: ... SNPs found: N of M` | N should be close to M |
| `LD consistency s = ...` | How well the LD fits the data: **< 0.05 good**, > 0.1 means mismatch, so be cautious with SuSiE. FADS1/LDLR were 0.005–0.023. |
| `retrying with coverage 0.1` | SuSiE found no signal at the usual 95% level (common for a weak eQTL) |
| `did NOT converge` | SuSiE results for that gene are unreliable |
| `Error` / `skipped:` | Something failed. See Troubleshooting. |
 
**Outputs** (in `784_data/subsetted_data/`):
 
| File | What it is |
|---|---|
| `coloc_summary_<gene>.tsv` | coloc.abf PP.H0–H4 |
| `coloc_<gene>_credible_set.tsv` | coloc.abf 95% credible set (only meaningful if H4 is high) |
| `coloc_<gene>_sensitivity.pdf`, `coloc_<gene>_datasets.pdf` | coloc.abf prior-sensitivity plot; GWAS and eQTL Manhattan plots |
| `coloc_susie_results_ukb/coloc_susie_<gene>_summary.tsv` | coloc + SuSiE PP.H0–H4, one row per GWAS signal × eQTL signal |
| `coloc_susie_results_ukb/susie_<gene>_{gwas,eqtl}_credible_sets.tsv` | SuSiE signals per trait (lead SNP, PIP) |
| `coloc_susie_results_ukb/ld_check_<gene>.pdf` | LD-fit diagnostic plots |
| `plots/<gene>_ukb_regional_manhattan.png` | GWAS over eQTL, colored by LD; purple diamond = top shared SNP, gold triangles = SuSiE leads |
| `plots/<gene>_ukb_pp_comparison.png` | coloc.abf vs coloc + SuSiE posteriors |
 
No `coloc_susie_<gene>_summary.tsv` means one trait (usually the eQTL) had no SuSiE signal. That's a result ("not testable"), not an error.
 
## 9. Reading the results
 
| | Meaning |
|---|---|
| H0 | no association with either trait |
| H1 | LDL signal only, no eQTL |
| H2 | eQTL only |
| H3 | both, **different** causal variants |
| H4 | both, **shared** causal variant (colocalized) |
 
- Group rule: **H4 > 0.5** = colocalized (0.8 is the stricter common choice; say which you use).
- coloc + SuSiE: each row is one GWAS signal × one eQTL signal; a high H4 in any row means colocalization.
- Examples:
  - **FADS1** colocalizes: coloc.abf H4 = 0.76; coloc + SuSiE finds 2 shared signals (H4 = 0.85 and 1.00).
  - **LDLR** does not: H1 = 0.70, a strong LDL signal but no liver eQTL.
## 10. Upload to GitHub
 
Through the website: **Add file → Upload files**.
- `data_subsetting/`: your `01-3` and `02-3` scripts.
- Top level: your `run_coloc_<gene>.sl` files.
- `results/`: `coloc_summary_<gene>.tsv`, `coloc_susie_<gene>_summary.tsv`, `susie_<gene>_*_credible_sets.tsv`, `<gene>_ukb_regional_manhattan.png`, `<gene>_ukb_pp_comparison.png`.
**Don't upload** data (`*_region.tsv`, `*_only.tsv`, `ukb_ld/`, raw downloads) or your edited `run_coloc_Abbey.R`.
 
---
 
## Troubleshooting (problems we already hit)
 
| Symptom | Cause | Fix |
|---|---|---|
| `No position column ... Columns present: V1, V2, ...` | GWAS region file has no header (`tabix -H` prints nothing for this file) | Use the `{ zcat ... \| head -n1; tabix ...; }` lines from `01-2.subsetted_data.sh` |
| `there is no package called 'X'` | R package not installed | Install it as in step 2 (after `module load r`) |
| `conda activate` error | conda not set up in batch jobs | Replace that line in your job file with `source activate base` |
| `numpy/scipy not available` | Python environment missing them | `module load anaconda; conda activate base; conda install numpy scipy` |
| `Download failed for chr..._....npz` | UK Biobank window has a different file name (e.g. near a chromosome end) | Run the `aws s3 ls` command printed in the error to find the right name |
| coloc.abf reads Abbey's files or "File not found" | `DIR` in `run_coloc_Abbey.R` not changed | Step 5 |
| eQTL `_only` file has 1 line | Gene symbol not found in the eQTL file | `cut -f13 eqtl_liver_<gene>_region.tsv \| sort -u \| grep -i <gene>` |
| Job stopped at the time limit | Usually SuSiE not converging | Check `LD consistency s`; raise `--time` |
 
The 1000 Genomes LD option still exists (`LD_SOURCE <- "1kg"` in `run_coloc_susie.R`, plus `03.get_LD_ref.sh` and `module load plink/1.90b6.21`; the default `plink` module is PLINK 2 and won't work). We use UK Biobank LD because it fits the GWAS much better (s ≈ 0.01 vs ≈ 0.2).
 
