# =============================================================================
# plot_locus.R -- extra figures for one gene, made from the outputs of the
# group's scripts (coloc_abf.R and run_coloc_susie.R). Run it after both.
#
#   Rscript plot_locus.R fads1          # uses the UK Biobank-LD SuSiE results (default)
#   Rscript plot_locus.R fads1 1kg      # uses the 1000 Genomes-LD SuSiE results
#
# Makes, in <DIR>/plots/:
#   <gene>_<ld>_regional_manhattan.png  LDL GWAS (top) and liver eQTL (bottom), same x-axis,
#                                       coloured by LD (r^2) with the lead GWAS SNP;
#                                       diamond = top coloc.abf shared SNP,
#                                       triangles = SuSiE credible-set lead SNPs
#   <gene>_<ld>_pp_comparison.png       PP.H0-H4: coloc.abf vs each coloc.SuSiE pair
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

DIR      <- Sys.getenv("COLOC_DIR", "/hickory/users/a/b/abskin/subsetted_data")  # same DIR as run_coloc_susie.R (set COLOC_DIR to override)
LD_BFILE <- file.path(DIR, "1kg", "EUR")                  # 1000G panel, only a fallback for LD colouring
PLOT_DIR <- file.path(DIR, "plots"); dir.create(PLOT_DIR, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
g <- tolower(args[1])
if (is.na(g)) stop("Usage: Rscript plot_locus.R <gene> [ukb|1kg]")
LD_SOURCE <- if (length(args) >= 2) tolower(args[2]) else "ukb"
SUS_DIR  <- file.path(DIR, paste0("coloc_susie_results_", LD_SOURCE))   # same as run_coloc_susie.R
G <- toupper(g)
logp <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)

# ---- coloc.abf output (per-SNP z-scores for both traits) ----
res <- fread(file.path(DIR, paste0("coloc_", g, "_snp_results.tsv")))
res[, `:=`(logp_g = logp(z.df1), logp_e = logp(z.df2))]
gw  <- fread(file.path(DIR, paste0("gwas_ldl_", g, "_region.tsv")), select = c("rsID", "POS_b37"))
gw  <- gw[!duplicated(POS_b37)]
res[, rsid := gw$rsID[match(position, gw$POS_b37)]]
res[, label := fifelse(grepl("^rs", rsid), rsid, snp)]
CHR <- sub("_.*", "", res$snp[1])

pp_abf <- fread(file.path(DIR, paste0("coloc_summary_", g, ".tsv")))
pp_abf <- unlist(pp_abf[gene == G, paste0("PP.H", 0:4, ".abf"), with = FALSE])

# ---- coloc.SuSiE output (if it ran) ----
f_sus <- file.path(SUS_DIR, paste0("coloc_susie_", g, "_summary.tsv"))
sus   <- if (file.exists(f_sus)) fread(f_sus) else NULL
cs_lead <- function(which) {
  f <- file.path(SUS_DIR, paste0("susie_", g, "_", which, "_credible_sets.tsv"))
  if (file.exists(f) && file.size(f) > 1) fread(f)$lead_snp else character(0)
}
lead_gwas_cs <- cs_lead("gwas"); lead_eqtl_cs <- cs_lead("eqtl")

# ---- LD (r^2) with the lead GWAS SNP ----
# Preferred: the LD matrix run_coloc_susie.R used (saved as ld_<gene>.rds), so the
# colours match the LD that SuSiE saw. Fallback: 1000G panel via plink.
lead <- res[which.max(logp_g)]
res[, r2 := NA_real_]
plink <- Sys.which("plink")
f_ld <- file.path(SUS_DIR, paste0("ld_", g, ".rds"))
if (file.exists(f_ld)) {
  R  <- readRDS(f_ld)
  inR <- res$snp %in% colnames(R)
  lead <- res[inR][which.max(logp_g)]                    # strongest GWAS SNP that is in the matrix
  res[inR, r2 := R[lead$snp, snp]^2]
  message("LD colouring from ", basename(f_ld), " (", LD_SOURCE, "): r^2 for ", sum(inR), " of ", nrow(res),
          " SNPs; grey = not in the LD matrix (removed by SuSiE QC or not in the panel)")
} else if (nzchar(plink) && file.exists(paste0(LD_BFILE, ".bim")) && grepl("^rs", lead$rsid)) {
  tmp <- tempfile()
  system2(plink, c("--bfile", LD_BFILE, "--r2", "--ld-snp", lead$rsid,
                   "--ld-window-kb", "5000", "--ld-window", "999999", "--ld-window-r2", "0",
                   "--out", tmp), stdout = FALSE, stderr = FALSE)
  if (file.exists(paste0(tmp, ".ld"))) {
    ld <- fread(paste0(tmp, ".ld"))
    res[, r2 := ld$R2[match(position, ld$BP_B)]]
    message("LD colouring: r^2 found for ", sum(!is.na(res$r2)), " of ", nrow(res), " SNPs")
  } else message("plink ran but made no .ld file (is ", lead$rsid, " in the panel?) -- see ", tmp, ".log")
} else message("No LD colouring (needs plink + ", LD_BFILE, ".bim + an rsID for the lead SNP)")

ld_cols <- c("<0.2" = "#2c3e91", "0.2-0.4" = "#56b4e9", "0.4-0.6" = "#009e73",
             "0.6-0.8" = "#e69f00", ">0.8" = "#d55e00")
res[, ld := cut(r2, c(-Inf, .2, .4, .6, .8, Inf), labels = names(ld_cols))]
top <- res[which.max(SNP.PP.H4)]

panel <- function(y, ttl, cs_snps) {
  ggplot(res, aes(position / 1e6, .data[[y]])) +
    geom_hline(yintercept = -log10(5e-8), colour = "red", linetype = 2, linewidth = 0.3) +
    geom_point(aes(colour = ld), size = 1.2, alpha = 0.85) +
    geom_point(data = res[snp %in% cs_snps], shape = 24, size = 3, fill = "gold", colour = "black") +
    geom_point(data = top, shape = 23, size = 3.5, fill = "purple", colour = "black") +
    geom_text(data = top, aes(label = label), vjust = -1.1, size = 3) +
    scale_colour_manual(values = ld_cols, na.value = "grey75", drop = FALSE,
                        name = paste0("r^2 with\n", lead$label)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = sprintf("chr%s position (Mb, hg19)", CHR), y = expression(-log[10](p)), title = ttl) +
    theme_bw(11)
}
p <- (panel("logp_g", "LDL-C GWAS (GLGC 2021, EUR)", lead_gwas_cs) /
      panel("logp_e", sprintf("%s expression in liver (Broadaway 2024)", G), lead_eqtl_cs)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = sprintf("%s locus", G),
    subtitle = sprintf("coloc.abf PP.H4 = %.2f | purple diamond = top shared SNP\ngold triangles = SuSiE credible-set leads | red line = p 5e-8",
                       pp_abf[5]))
ggsave(file.path(PLOT_DIR, paste0(g, "_", LD_SOURCE, "_regional_manhattan.png")), p, width = 9, height = 8, dpi = 200)

# ---- coloc.abf vs coloc.SuSiE posteriors ----
pp <- data.table(method = "coloc.abf", hyp = paste0("H", 0:4), pp = as.numeric(pp_abf))
if (!is.null(sus) && nrow(sus)) for (i in seq_len(nrow(sus)))
  pp <- rbind(pp, data.table(method = sprintf("coloc.SuSiE: GWAS CS%s x eQTL CS%s", sus$idx1[i], sus$idx2[i]),
                             hyp = paste0("H", 0:4),
                             pp = as.numeric(unlist(sus[i, paste0("PP.H", 0:4, ".abf"), with = FALSE]))))
pp[, method := factor(method, levels = rev(unique(method)))]
p2 <- ggplot(pp, aes(pp, method, fill = hyp)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c(H0 = "#d9d9d9", H1 = "#9ecae1", H2 = "#a1d99b", H3 = "#fdae6b", H4 = "#cb181d"),
    labels = c("H0 none", "H1 LDL only", "H2 eQTL only", "H3 different variants", "H4 shared variant"), name = NULL) +
  scale_x_continuous(limits = c(0, 1.0001), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Posterior probability", y = NULL, title = sprintf("%s: coloc.abf vs coloc.SuSiE (%s LD)", G, toupper(LD_SOURCE))) +
  theme_bw(11) + theme(legend.position = "bottom")
ggsave(file.path(PLOT_DIR, paste0(g, "_", LD_SOURCE, "_pp_comparison.png")), p2,
       width = 8.5, height = 2 + 0.45 * uniqueN(pp$method), dpi = 200)
cat("Plots written to", PLOT_DIR, "\n")
