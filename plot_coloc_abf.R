# =============================================================================
# plot_coloc_abf.R -- coloc.abf figure for one gene, made only from the
# outputs of run_coloc_Bea.R / run_coloc_Abbey.R (no SuSiE, no LD needed).
#
#   Rscript plot_coloc_abf.R sort1
#
# Reads from <DIR>: coloc_<gene>_snp_results.tsv, coloc_summary_<gene>.tsv,
#                   gwas_ldl_<gene>_region.tsv (only for rsIDs in the labels)
# Writes <DIR>/plots/<gene>_coloc_abf.png with two panels:
#   A. LDL-C GWAS vs liver eQTL -log10(p), one point per SNP. If one variant
#      drives both traits, the strongest SNPs sit together in the top-right
#      corner; if different variants do, the points form an "L" along the axes.
#   B. SNP.PP.H4 along the region: the probability that each SNP is the shared
#      causal variant, assuming H4 is true.
# Purple = 95% credible set for the shared variant; diamond = top SNP
# (the same SNP as the purple diamond in plot_locus.R).
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

DIR <- Sys.getenv("COLOC_DIR", "/work/users/b/e/beacm/784/784") # set COLOC_DIR to override
PLOT_DIR <- file.path(DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE)

g <- tolower(commandArgs(trailingOnly = TRUE)[1])
if (is.na(g)) stop("Usage: Rscript plot_coloc_abf.R <gene>")
G <- toupper(g)

COL_CS <- "#7B2CBF" # credible-set SNPs (purple, as in plot_locus.R)
COL_OTHER <- "grey72" # all other SNPs

# -log10(two-sided p) from z, without underflow at huge z
logp <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)

# ---- Read coloc.abf outputs ----
res <- fread(file.path(DIR, paste0("coloc_", g, "_snp_results.tsv")))
pp <- fread(file.path(DIR, paste0("coloc_summary_", g, ".tsv")))
pp <- unlist(pp[gene == G, paste0("PP.H", 0:4, ".abf"), with = FALSE])
res[, `:=`(logp_g = logp(z.df1), logp_e = logp(z.df2))]
CHR <- sub("_.*", "", res$snp[1])

# rsIDs for labels (fall back to chr_pos_ref_alt)
f_gwas <- file.path(DIR, paste0("gwas_ldl_", g, "_region.tsv"))
res[, label := snp]
if (file.exists(f_gwas)) {
  gw <- fread(f_gwas, select = c("rsID", "POS_b37"))
  gw <- gw[!duplicated(POS_b37)]
  rs <- gw$rsID[match(res$position, gw$POS_b37)]
  res[grepl("^rs", rs), label := rs[grepl("^rs", rs)]]
}

# ---- 95% credible set for the shared variant (same rule as run_coloc_Bea.R) ----
setorder(res, -SNP.PP.H4)
res[, in_cs := seq_len(.N) <= which(cumsum(SNP.PP.H4) > 0.95)[1]]
n_cs <- sum(res$in_cs)
top <- res[1]
res <- res[order(in_cs)] # draw credible-set SNPs on top

# ---- Panel A: GWAS vs eQTL significance ----
pa <- ggplot(res, aes(logp_g, logp_e)) +
  geom_point(aes(colour = in_cs), size = 1.4, alpha = 0.8) +
  geom_point(data = top, shape = 23, size = 4, fill = COL_CS, colour = "black") +
  geom_text(data = top, aes(label = label), vjust = -1.2, hjust = "inward", size = 3.8) +
  scale_colour_manual(values = c(`FALSE` = COL_OTHER, `TRUE` = COL_CS), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    x = expression("LDL-C GWAS " * -log[10](p)),
    y = bquote("Liver" ~ italic(.(G)) ~ "eQTL" ~ -log[10](p)),
    title = "A. Do the same SNPs drive both traits?"
  )

# ---- Panel B: SNP.PP.H4 along the region ----
# Zoomed to the credible set (+-100 kb, or 10% of its span if wider), and the
# y-axis fitted to the data so that SNP.PP.H4 spread over many SNPs is visible
cs_range <- range(res[in_cs == TRUE]$position)
pad <- max(1e5, 0.1 * diff(cs_range))
x_win <- (cs_range + c(-pad, pad)) / 1e6
y_max <- max(res$SNP.PP.H4) * 1.18

pb <- ggplot(res, aes(position / 1e6, SNP.PP.H4)) +
  geom_segment(aes(xend = position / 1e6, yend = 0, colour = in_cs), linewidth = 0.5) +
  geom_point(data = res[in_cs == TRUE], colour = COL_CS, size = 1.6) +
  geom_point(data = top, shape = 23, size = 4, fill = COL_CS, colour = "black") +
  geom_text(
    data = top, aes(label = sprintf("%s (%.2f)", label, SNP.PP.H4)),
    vjust = -1.2, hjust = "inward", size = 3.8
  ) +
  scale_colour_manual(values = c(`FALSE` = COL_OTHER, `TRUE` = COL_CS), guide = "none") +
  scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
  coord_cartesian(xlim = x_win) +
  labs(
    x = sprintf("chr%s position (Mb, hg19), zoomed to the credible set", CHR),
    y = "SNP.PP.H4",
    title = "B. Which SNP is the shared variant (if H4 is true)?"
  )

# ---- Combine and save ----
note <- if (pp[5] < 0.5) {
  "\nH4 < 0.5: not colocalised by the H4 > 0.5 rule, so read panel B (which assumes H4) with caution"
} else {
  ""
}
p <- (pa | pb) +
  plot_layout(widths = c(1, 1.6)) +
  plot_annotation(
    title = sprintf("%s: coloc.abf (one causal variant per trait)", G),
    subtitle = sprintf(
      "PP.H4 (shared variant) = %.2f | PP.H3 (different variants) = %.2f | purple = 95%% credible set (%d SNP%s) | diamond = top SNP%s",
      pp[5], pp[4], n_cs, if (n_cs == 1) "" else "s", note
    )
  ) &
  theme_bw(13) &
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 13)
  )

out <- file.path(PLOT_DIR, paste0(g, "_coloc_abf.png"))
ggsave(out, p, width = 13, height = 5.5, dpi = 200)
cat("Wrote", out, "\n")
