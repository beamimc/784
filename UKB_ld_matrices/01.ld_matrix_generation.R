#!/usr/bin/env Rscript
# =============================================================================
# 01.ld_matrix_generation.R — Extract a UK Biobank LD matrix for a genomic window
#
# Source: Broad/PolyFun precomputed UKB LD (Weissbrod et al. 2020)
#   ~337k unrelated British-ancestry UKB participants, hg19 (b37)
#   3 Mb windows in 1 Mb steps: chr{C}_{S}_{S+3000000}.{npz,gz}
#   https://github.com/omerwe/polyfun
#
# What it does:
#   1. Picks the UKB 3 Mb file that contains your window (best centered)
#   2. Downloads the .gz (SNP list) and .npz (sparse LD) if not already present
#   3. Loads the sparse matrix via Python/SciPy (reticulate), symmetrizes it
#   4. Subsets to SNPs inside [start, end] (and optionally a SNP list)
#      BEFORE densifying, to keep memory low
#   5. Returns signed LD (r) + the matching SNP table (same order)
#
# Requirements (one-time):
#   install.packages(c("reticulate", "data.table"))
#   reticulate::py_install(c("numpy", "scipy"))
#
# Use in R:
#   source("01.ld_matrix_generation.R")
#   ld <- get_ukb_ld(chr = 1, start = 108274968, end = 110274968)
#   ld$R      # SNP x SNP signed correlation matrix
#   ld$snps   # rsid, chromosome, position, allele1, allele2 (row order of R)
#
# Use from the command line (saves an .rds with list(file, R, snps)):
#   Rscript 01.ld_matrix_generation.R <chr> <start> <end> [ld_dir] [out.rds | outdir/] [snp_keep.txt]
#   - 5th arg as a file path:  saves exactly there
#   - 5th arg as a folder:     saves <outdir>/ld_chr<C>_<start>_<end>.rds
#
# NOT done here (do these afterwards):
#   - Allele harmonization of GWAS/eQTL z-scores to ld$snps alleles
#   - z/R consistency checks (susieR::estimate_s_rss, kriging_rss)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(reticulate)
})

UKB_LD_URL <- "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/"

# -----------------------------------------------------------------------------
# Choose the UKB window file that contains [start, end]
#   File k spans [k*1e6 + 1, k*1e6 + 1 + 3e6]
#   Valid k: window start <= start AND window end >= end
#   Among valid k, choose the one that best centers the region
# -----------------------------------------------------------------------------
ukb_window <- function(chr, start, end) {
  chr   <- sub("^chr", "", as.character(chr))
  start <- as.numeric(start); end <- as.numeric(end)
  if (end <= start) stop("end must be greater than start")

  k_max <- floor((start - 1) / 1e6)
  k_min <- max(ceiling((end - 3e6 - 1) / 1e6), 0)
  if (k_min > k_max)
    stop(sprintf("Region %s:%d-%d (%.2f Mb) does not fit in one 3 Mb UKB window; narrow it.",
                 chr, start, end, (end - start) / 1e6))

  # Best-centered k (same integer rounding as the bash script)
  k <- floor((start + end - 3000002 + 1e6) / 2e6)
  k <- min(max(k, k_min), k_max)

  s <- k * 1e6 + 1
  sprintf("chr%s_%.0f_%.0f", chr, s, s + 3e6)
}

# -----------------------------------------------------------------------------
# Download the .gz and .npz for a window prefix (skips files already present)
# -----------------------------------------------------------------------------
download_ukb_ld <- function(prefix, ld_dir = "ukb_ld") {
  dir.create(ld_dir, showWarnings = FALSE, recursive = TRUE)
  old_to <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(old_to))

  for (ext in c(".gz", ".npz")) {
    f <- file.path(ld_dir, paste0(prefix, ext))
    if (file.exists(f) && file.size(f) > 0) next
    message("Downloading ", prefix, ext, " ...")
    ok <- tryCatch({
      download.file(paste0(UKB_LD_URL, prefix, ext), f, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok || !file.exists(f) || file.size(f) == 0) {
      if (file.exists(f)) file.remove(f)
      stop("Download failed for ", prefix, ext,
           ". The file name may differ (e.g., at chromosome ends). List the bucket:\n",
           "  aws s3 ls --no-sign-request s3://broad-alkesgroup-ukbb-ld/UKBB_LD/ | grep chr",
           sub("_.*", "", sub("^chr", "", prefix)), "_")
    }
  }
  invisible(file.path(ld_dir, paste0(prefix, c(".gz", ".npz"))))
}

# -----------------------------------------------------------------------------
# Python helper: load .npz, symmetrize, subset rows/cols, densify
# -----------------------------------------------------------------------------
.ukb_py_helper <- function() {
  # Start Python explicitly and fail with a clear message if it is unusable
  ok <- tryCatch(py_available(initialize = TRUE), error = function(e) FALSE)
  if (!isTRUE(ok))
    stop("Python could not be started by reticulate.\n",
         "  RETICULATE_PYTHON = '", Sys.getenv("RETICULATE_PYTHON"), "'\n",
         "  Activate an environment with numpy + scipy and set RETICULATE_PYTHON=$(which python).")
  for (m in c("numpy", "scipy"))
    if (!py_module_available(m))
      stop("Python module '", m, "' not found in ", py_config()$python,
           "\n  Install it there (e.g., conda install numpy scipy) or point RETICULATE_PYTHON to an env that has it.")

  main <- import_main(convert = TRUE)
  if (py_has_attr(main, "ukb_subset")) return(invisible(main))
  py_run_string("
import numpy as np
import scipy.sparse as sp

def ukb_subset(npz_file, idx, out_bin):
    M = sp.load_npz(npz_file)
    n = M.shape[0]
    M = (M + M.T).tocsr()
    idx = np.atleast_1d(np.asarray(idx, dtype=np.int64))
    A = M[idx][:, idx].toarray().astype('<f8')
    A.T.tofile(out_bin)          # column-major doubles, read with readBin() in R
    return int(n)
")
  invisible(main)
}

# -----------------------------------------------------------------------------
# Main: LD matrix for SNPs in [start, end]
#   snp_keep: optional character vector of rsIDs to retain (e.g., harmonized
#             SNPs shared by GWAS and eQTL)
# -----------------------------------------------------------------------------
get_ukb_ld <- function(chr, start, end, ld_dir = "ukb_ld", snp_keep = NULL,
                       download = TRUE) {
  prefix <- ukb_window(chr, start, end)
  if (download) download_ukb_ld(prefix, ld_dir)

  gz_file  <- file.path(ld_dir, paste0(prefix, ".gz"))
  npz_file <- file.path(ld_dir, paste0(prefix, ".npz"))
  if (!file.exists(gz_file) || !file.exists(npz_file))
    stop("Missing ", prefix, ".gz/.npz in ", ld_dir)

  # SNP table (one row per SNP in the 3 Mb file; row i = matrix row i)
  snps <- fread(gz_file)
  if (!"position" %in% names(snps))
    stop("Unexpected SNP file columns: ", paste(names(snps), collapse = ", "))

  # SNPs inside the requested window (+ optional SNP list)
  keep <- which(snps$position >= start & snps$position <= end)
  if (!is.null(snp_keep)) keep <- keep[snps$rsid[keep] %in% snp_keep]
  if (length(keep) == 0) stop("No UKB SNPs left in ", chr, ":", start, "-", end)

  # Sparse LD (stored as one triangle): symmetrize, subset, THEN densify.
  # Done in Python so only the subset ever becomes a dense matrix.
  main <- .ukb_py_helper()
  # Matrix is passed back via a temporary binary file, which avoids
  # reticulate/numpy version-conversion issues and is fast for large matrices
  tmp <- tempfile(fileext = ".bin"); on.exit(unlink(tmp), add = TRUE)
  n_total <- main$ukb_subset(npz_file, as.list(as.integer(keep - 1L)), tmp)  # 0-based
  if (as.integer(n_total) != nrow(snps))
    stop("Matrix dimension does not match SNP table for ", prefix)
  m <- length(keep)
  R <- matrix(readBin(tmp, "double", n = m * m, size = 8, endian = "little"),
              nrow = m, ncol = m)

  # Diagonal: should be 1 (or 2 if the triangle included the diagonal)
  d <- diag(R)
  if (!all(abs(d - 1) < 1e-4 | abs(d - 2) < 1e-4))
    warning("Unexpected diagonal values (range ", paste(signif(range(d), 4), collapse = "-"),
            "); check the LD file.")
  diag(R) <- 1

  ids <- if ("rsid" %in% names(snps)) snps$rsid[keep] else
    paste(snps$chromosome[keep], snps$position[keep], sep = ":")
  dimnames(R) <- list(ids, ids)

  message(sprintf("%s: %d SNPs in %s:%.0f-%.0f", prefix, length(keep),
                  sub("^chr", "", chr), start, end))
  list(file = prefix, R = R, snps = snps[keep])
}

# -----------------------------------------------------------------------------
# Command-line entry point
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    cat("Usage: Rscript 01.ld_matrix_generation.R <chr> <start> <end> [ld_dir] [out.rds] [snp_keep.txt]\n")
    quit(status = 1)
  }
  chr    <- args[1]
  start  <- as.numeric(args[2])
  end    <- as.numeric(args[3])
  ld_dir <- if (length(args) >= 4) args[4] else "ukb_ld"
  default_name <- sprintf("ld_chr%s_%.0f_%.0f.rds", sub("^chr", "", chr), start, end)
  out    <- if (length(args) >= 5) args[5] else default_name
  # 5th argument may be an output folder (existing, or ending in "/") or a file path
  if (dir.exists(out) || grepl("/$", out)) out <- file.path(sub("/$", "", out), default_name)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  snp_keep <- if (length(args) >= 6) readLines(args[6]) else NULL

  ld <- get_ukb_ld(chr, start, end, ld_dir = ld_dir, snp_keep = snp_keep)
  saveRDS(ld, out)
  message("Saved ", out, " (", nrow(ld$R), " x ", ncol(ld$R), ")")
}
