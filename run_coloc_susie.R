# =============================================================================
# Colocalisation of a GWAS and an eQTL dataset with SuSiE, gene by gene
# (allows more than one causal variant per trait)
# Follows the coloc SuSiE vignette:
#   https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html
# LD: reference panel in PLINK format (e.g. 1000 Genomes EUR, hg19), via ieugwasr::ld_matrix_local()
#
# Run:   module load plink
#        Rscript run_coloc_susie.R                 # uses the genes listed in SETTINGS
#        Rscript run_coloc_susie.R sort1 hmgcr     # or give genes on the command line
#
# Per gene it expects two region TSVs (names set by GWAS_FILE / EQTL_FILE below).
# Outputs go to OUT_DIR:
#   ld_check_<gene>.pdf, coloc_susie_<gene>_sensitivity.pdf, coloc_susie_<gene>_summary.tsv,
#   susie_<gene>_gwas_credible_sets.tsv, susie_<gene>_eqtl_credible_sets.tsv, coloc_susie_summary.tsv
#
# QC before SuSiE: drop strand-ambiguous (A/T, C/G) SNPs, SNPs with MAF < MAF_MIN,
# SNPs with N < N_FRAC_MIN x max N, and SNPs whose z-scores don't fit the LD panel (kriging_rss).
# (install.packages("Rfast") makes SuSiE much faster on large LD matrices)
# =============================================================================

library(coloc)
library(susieR)
library(ieugwasr)
library(data.table)

# =============================================================================
# SETTINGS - the only part you should need to edit
# =============================================================================

# --- Where things are -------------------------------------------------------
DIR      <- "/work/users/b/e/beacm/784/784"          # working directory with the input TSVs
genes    <- c("sort1", "hmgcr")                      # overridden by command-line arguments, if given
LD_BFILE <- file.path(DIR, "1kg", "EUR")             # PLINK prefix of the LD reference (.bed/.bim/.fam)
OUT_DIR  <- file.path(DIR, "coloc_susie_results")    # created if it doesn't exist

# --- Input file names: {gene} is replaced by each gene name -----------------
GWAS_FILE <- "gwas_ldl_{gene}_region.tsv"
EQTL_FILE <- "eqtl_liver_{gene}_only.tsv"

# --- Column names in the input files (right-hand side = your column) --------
GWAS_COLS <- c(rsid = "rsID", pos = "POS_b37", ref = "REF", alt = "ALT",   # alt = effect allele
               af = "POOLED_ALT_AF", beta = "EFFECT_SIZE", se = "SE", n = "N")
EQTL_COLS <- c(snp = "Variant", pos = "POS", nea = "NEA", ea = "EA",
               eaf = "EAF", beta = "Beta", se = "SE", n = "N")

# --- Trait details ----------------------------------------------------------
GWAS_LABEL <- "LDL GWAS"
EQTL_LABEL <- "Liver eQTL"
GWAS_SDY   <- 1      # SD of the GWAS trait; 1 for an inverse-normal transformed trait, NULL = estimate
EQTL_SDY   <- NULL   # NULL = coloc estimates it from MAF and N

# --- QC and SuSiE settings (usually leave as they are) ----------------------
MAF_MIN            <- 0.01   # a small LD panel can't estimate LD for rare variants
N_FRAC_MIN         <- 0.9    # keep SNPs with N >= 90% of max N (lower to 0.8 if too few remain)
DROP_PALINDROMIC   <- TRUE   # A/T and C/G SNPs can't be strand-aligned reliably
REMOVE_LD_OUTLIERS <- TRUE   # drop SNPs whose z-scores don't fit the LD (kriging_rss)
MAX_OUTLIER_FRAC   <- 0.05   # if more than this share is flagged, the panel is the problem: drop nothing
SUSIE_MAXIT        <- 1000   # iterations before giving up
H4_RULE            <- "H4 > 0.5"   # decision rule for the sensitivity plots (0.8 is also common)
MIN_SNPS           <- 50     # skip a gene if fewer SNPs than this survive QC

# =============================================================================
# Nothing below should need editing
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) genes <- tolower(args)

PLINK <- Sys.which("plink")
if (PLINK == "") PLINK <- genetics.binaRies::get_plink_binary()
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_file <- function(...) file.path(OUT_DIR, paste0(...))

# Read only the needed columns and rename them to the script's internal names
read_cols <- function(file, cols, prefix, keep_plain = c("snp", "rsid", "pos")) {
  if (!file.exists(file)) stop("File not found: ", file)
  d <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop(basename(file), ": missing column(s) ", paste(miss, collapse = ", "))
  d <- d[, cols]
  names(d) <- ifelse(names(cols) %in% keep_plain, names(cols), paste0(prefix, names(cols)))
  d
}

# Lead SNP and size of each SuSiE credible set
cs_table <- function(S) {
  cs <- S$sets$cs
  if (length(cs) == 0) return(data.frame())
  data.frame(cs       = names(cs),
             n_snps   = sapply(cs, length),
             lead_snp = sapply(cs, function(i) names(S$pip)[i][which.max(S$pip[i])]),
             lead_pip = sapply(cs, function(i) max(S$pip[i])))
}

# Keep a subset of SNPs in a coloc dataset (vectors + LD)
subset_ds <- function(D, keep) {
  for (f in c("beta", "varbeta", "snp", "position", "MAF")) D[[f]] <- D[[f]][keep]
  D$LD <- D$LD[keep, keep]
  D
}

# How well do this trait's z-scores fit the LD panel? (susieR summary-stat diagnostics)
#   s near 0 = consistent; s > ~0.1 = LD mismatch between panel and GWAS/eQTL sample
ld_check <- function(D, label) {
  z  <- D$beta / sqrt(D$varbeta)
  s  <- estimate_s_rss(z, D$LD, n = D$N)
  kr <- kriging_rss(z, D$LD, n = D$N, s = s)
  cd <- kr$conditional_dist
  bad <- D$snp[cd$logLR > 2 & abs(cd$z) > 2]
  cat(sprintf("%s: LD consistency s = %.3f | outlier SNPs = %d | max |z| = %.1f\n",
              label, s, length(bad), max(abs(z))))
  print(kr$plot + ggplot2::ggtitle(sprintf("%s: observed vs expected z (s = %.3f)", label, s)))
  list(s = s, outliers = bad, lead = D$snp[which.max(abs(z))])
}

# SuSiE with an iteration cap + the vignette's low-coverage retry
fit_susie <- function(D, label) {
  S <- runsusie(D, maxit = SUSIE_MAXIT, repeat_until_convergence = FALSE)
  if (length(S$sets$cs) == 0) {
    message(label, ": no credible sets at coverage 0.95; retrying with coverage 0.1 (interpret with caution)")
    S <- runsusie(D, coverage = 0.1, maxit = SUSIE_MAXIT, repeat_until_convergence = FALSE)
  }
  if (!isTRUE(S$converged))
    message(label, ": SuSiE did NOT converge in ", SUSIE_MAXIT,
            " iterations - credible sets may be unreliable (likely LD mismatch)")
  S
}

# -----------------------------------------------------------------------------
# Full analysis for one gene; returns its summary table (or NULL)
# -----------------------------------------------------------------------------
run_gene <- function(g) {
  cat("\n==================== ", toupper(g), " ====================\n")
  
  # ---- 1. Read and harmonise GWAS + eQTL ----
  eqtl <- read_cols(file.path(DIR, gsub("{gene}", g, EQTL_FILE, fixed = TRUE)), EQTL_COLS, "eqtl_")
  gwas <- read_cols(file.path(DIR, gsub("{gene}", g, GWAS_FILE, fixed = TRUE)), GWAS_COLS, "gwas_")
  for (a in c("eqtl_nea", "eqtl_ea")) eqtl[[a]] <- toupper(eqtl[[a]])
  for (a in c("gwas_ref", "gwas_alt")) gwas[[a]] <- toupper(gwas[[a]])
  
  m <- merge(eqtl, gwas, by = "pos")
  same <- m$eqtl_nea == m$gwas_ref & m$eqtl_ea == m$gwas_alt
  flip <- m$eqtl_nea == m$gwas_alt & m$eqtl_ea == m$gwas_ref
  m <- m[same | flip, ]
  flip <- flip[same | flip]
  m$gwas_beta <- ifelse(flip, -m$gwas_beta, m$gwas_beta)   # align GWAS to eQTL effect allele
  
  m$gwas_maf <- pmin(m$gwas_af, 1 - m$gwas_af)
  m$eqtl_maf <- pmin(m$eqtl_eaf, 1 - m$eqtl_eaf)
  m <- m[complete.cases(m[, c("gwas_beta", "gwas_se", "eqtl_beta", "eqtl_se",
                              "gwas_maf", "eqtl_maf", "gwas_n", "eqtl_n")]), ]
  m <- m[m$gwas_se > 0 & m$eqtl_se > 0 & m$gwas_maf > 0 & m$eqtl_maf > 0, ]
  m <- m[!duplicated(m$snp), ]
  m <- m[grepl("^rs", m$rsid) & !duplicated(m$rsid), ]      # LD panel is indexed by rsID
  cat("Matched GWAS-eQTL SNPs with rsIDs:", nrow(m), "\n")
  
  # ---- 1b. QC: strand-ambiguous, rare, low-N SNPs ----
  if (DROP_PALINDROMIC) {
    pal <- paste0(m$eqtl_ea, m$eqtl_nea) %in% c("AT", "TA", "CG", "GC")
    m <- m[!pal, ]
    cat("Dropped strand-ambiguous (A/T, C/G) SNPs:", sum(pal), "\n")
  }
  n_before <- nrow(m)
  m <- m[m$gwas_maf >= MAF_MIN & m$eqtl_maf >= MAF_MIN, ]
  cat("Dropped SNPs with MAF <", MAF_MIN, ":", n_before - nrow(m), "| remaining:", nrow(m), "\n")
  
  # per-SNP N varies in a meta-analysis, but SuSiE gets a single N,
  # so low-N SNPs look inconsistent with the LD
  n_before <- nrow(m)
  m <- m[m$gwas_n >= N_FRAC_MIN * max(m$gwas_n) & m$eqtl_n >= N_FRAC_MIN * max(m$eqtl_n), ]
  cat("Dropped SNPs with N <", N_FRAC_MIN, "x max N:", n_before - nrow(m), "| remaining:", nrow(m), "\n")
  if (nrow(m) < MIN_SNPS) stop("only ", nrow(m), " SNPs left after QC (MIN_SNPS = ", MIN_SNPS, ")")
  
  # ---- 2. Signed LD matrix (names are rsID_A1_A2, A1 = counted allele) ----
  ld <- ld_matrix_local(m$rsid, bfile = LD_BFILE, plink_bin = PLINK, with_alleles = TRUE)
  
  ld_info <- do.call(rbind, strsplit(colnames(ld), "_"))
  ld_info <- data.frame(rsid = ld_info[, 1], ld_a1 = ld_info[, 2], ld_a2 = ld_info[, 3],
                        stringsAsFactors = FALSE)
  m <- m[match(ld_info$rsid, m$rsid), ]                     # same order as the LD matrix
  
  # Align both traits' betas to the LD panel's A1 (betas are currently for the eQTL EA)
  a1_is_ea  <- m$eqtl_ea == ld_info$ld_a1 & m$eqtl_nea == ld_info$ld_a2
  a1_is_nea <- m$eqtl_ea == ld_info$ld_a2 & m$eqtl_nea == ld_info$ld_a1
  m$gwas_beta[a1_is_nea] <- -m$gwas_beta[a1_is_nea]
  m$eqtl_beta[a1_is_nea] <- -m$eqtl_beta[a1_is_nea]
  
  keep <- a1_is_ea | a1_is_nea                               # drop allele mismatches
  m  <- m[keep, ]
  ld <- ld[keep, keep]
  dimnames(ld) <- list(m$snp, m$snp)
  cat("SNPs in LD panel (allele-matched):", nrow(m),
      " (flipped to LD A1:", sum(a1_is_nea), ")\n")
  
  # ---- 3. Build coloc datasets (with LD) ----
  D1 <- list(beta = m$gwas_beta, varbeta = m$gwas_se^2, snp = m$snp, position = m$pos,
             type = "quant", N = median(m$gwas_n), MAF = m$gwas_maf, LD = ld)
  D2 <- list(beta = m$eqtl_beta, varbeta = m$eqtl_se^2, snp = m$snp, position = m$pos,
             type = "quant", N = max(m$eqtl_n), MAF = m$eqtl_maf, LD = ld)
  if (!is.null(GWAS_SDY)) D1$sdY <- GWAS_SDY
  if (!is.null(EQTL_SDY)) D2$sdY <- EQTL_SDY
  
  check_dataset(D1, req = "LD")
  check_dataset(D2, req = "LD")
  
  # ---- 3b. Check z-scores against the LD panel; optionally drop outlier SNPs ----
  pdf(out_file("ld_check_", g, ".pdf"), width = 7, height = 6)
  chk1 <- ld_check(D1, GWAS_LABEL)
  chk2 <- ld_check(D2, EQTL_LABEL)
  dev.off()
  
  drop <- union(chk1$outliers, chk2$outliers)
  leads_flagged <- intersect(drop, c(chk1$lead, chk2$lead))
  if (length(leads_flagged) > 0)
    message("Lead SNP(s) flagged as LD outliers and kept anyway: ", paste(leads_flagged, collapse = ", "))
  drop <- setdiff(drop, c(chk1$lead, chk2$lead))             # never drop the lead SNPs
  
  if (REMOVE_LD_OUTLIERS && length(drop) > 0) {
    if (length(drop) > MAX_OUTLIER_FRAC * length(D1$snp)) {
      message(length(drop), " SNPs flagged (> ", 100 * MAX_OUTLIER_FRAC, "%): the LD panel itself doesn't ",
              "fit these summary stats, so nothing was removed. Treat SuSiE results with caution.")
    } else {
      keep <- !(D1$snp %in% drop)
      D1 <- subset_ds(D1, keep)
      D2 <- subset_ds(D2, keep)
      cat("Removed", length(drop), "LD-outlier SNPs; remaining:", length(D1$snp), "\n")
    }
  }
  
  # ---- 4. Fine-map each trait with SuSiE ----
  S1 <- fit_susie(D1, GWAS_LABEL)
  S2 <- fit_susie(D2, EQTL_LABEL)
  
  cat("\n", GWAS_LABEL, " credible sets:\n", sep = "");  print(cs_table(S1))
  cat("\n", EQTL_LABEL, " credible sets:\n", sep = "");  print(cs_table(S2))
  fwrite(cs_table(S1), out_file("susie_", g, "_gwas_credible_sets.tsv"), sep = "\t")
  fwrite(cs_table(S2), out_file("susie_", g, "_eqtl_credible_sets.tsv"), sep = "\t")
  
  # ---- 5. Colocalise every pair of credible sets ----
  susie.res <- coloc.susie(S1, S2)
  if (is.null(susie.res$summary) || nrow(susie.res$summary) == 0) {
    cat("\nNo credible-set pairs to test for", toupper(g), "\n")
    return(NULL)
  }
  cat("\ncoloc.susie results (one row per GWAS signal x eQTL signal):\n")
  print(susie.res$summary)
  
  # ---- 6. Sensitivity plot for each pair ----
  pdf(out_file("coloc_susie_", g, "_sensitivity.pdf"), width = 11, height = 7)
  for (r in seq_len(nrow(susie.res$summary))) {
    tryCatch(
      sensitivity(susie.res, H4_RULE, row = r, dataset1 = D1, dataset2 = D2),
      error = function(e) {
        message("Row ", r, ": Manhattan panels failed (", conditionMessage(e), "); p12 panel only")
        sensitivity(susie.res, H4_RULE, row = r, plot.manhattans = FALSE)
      })
  }
  dev.off()
  
  # ---- 7. Save ----
  out <- data.frame(gene = toupper(g), susie.res$summary,
                    gwas_converged = isTRUE(S1$converged), eqtl_converged = isTRUE(S2$converged),
                    gwas_ld_s = chk1$s, eqtl_ld_s = chk2$s)
  fwrite(out, out_file("coloc_susie_", g, "_summary.tsv"), sep = "\t")
  out
}

# -----------------------------------------------------------------------------
# Run all genes; a failure in one gene is reported and the rest still run
# -----------------------------------------------------------------------------
all_summaries <- list()
for (g in genes) {
  res <- tryCatch(run_gene(g), error = function(e) {
    while (dev.cur() > 1) dev.off()                         # close any PDF left open
    message("\n!! ", toupper(g), " skipped: ", conditionMessage(e))
    NULL
  })
  if (!is.null(res)) all_summaries[[g]] <- res
}

if (length(all_summaries) > 0) {
  summary_tab <- rbindlist(all_summaries, fill = TRUE)
  cat("\n==================== SUMMARY ====================\n")
  print(summary_tab)
  fwrite(summary_tab, out_file("coloc_susie_summary.tsv"), sep = "\t")
} else {
  message("No genes produced coloc.susie results.")
}
