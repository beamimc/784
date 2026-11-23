# Colocalisation of LDL GWAS and liver eQTL (coloc.abf), one or more genes
# Follows the coloc enumeration vignette (coloc.abf):
#   https://chr1swallace.github.io/coloc/articles/a03_enumeration.html
#
# Same analysis as the group's run_coloc_Bea.R (github.com/beamimc/784).
# Only changes: DIR points to Abbey's data, genes can be given on the command
# line, and the summary file is named per gene so separate jobs don't overwrite it.
#
# Run:  Rscript run_coloc_Abbey.R fads1        (or: Rscript run_coloc_Abbey.R fads1 ldlr)

library(coloc)

DIR   <- "/hickory/users/a/b/abskin/subsetted_data"
genes <- c("fads1", "ldlr")
args  <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) genes <- tolower(args)

all_summaries <- list()

for (g in genes) {
  cat("\n==================== ", toupper(g), " ====================\n")
  
  # ---- 1. Read files ----
  eqtl <- read.delim(file.path(DIR, paste0("eqtl_liver_", g, "_only.tsv")),
                     stringsAsFactors = FALSE)
  gwas <- read.delim(file.path(DIR, paste0("gwas_ldl_", g, "_region.tsv")),
                     stringsAsFactors = FALSE)
  
  # Keep needed columns with clear names (both files have N and SE)
  eqtl <- data.frame(snp = eqtl$Variant, pos = eqtl$POS,
                     eqtl_nea = eqtl$NEA, eqtl_ea = eqtl$EA, eqtl_eaf = eqtl$EAF,
                     eqtl_beta = eqtl$Beta, eqtl_se = eqtl$SE, eqtl_n = eqtl$N)
  gwas <- data.frame(rsid = gwas$rsID, pos = gwas$POS_b37,
                     gwas_ref = gwas$REF, gwas_alt = gwas$ALT, gwas_af = gwas$POOLED_ALT_AF,
                     gwas_beta = gwas$EFFECT_SIZE, gwas_se = gwas$SE, gwas_n = gwas$N)
  
  # ---- 2. Match SNPs by position + alleles, align GWAS to eQTL effect allele ----
  m <- merge(eqtl, gwas, by = "pos")
  same <- m$eqtl_nea == m$gwas_ref & m$eqtl_ea == m$gwas_alt
  flip <- m$eqtl_nea == m$gwas_alt & m$eqtl_ea == m$gwas_ref
  m <- m[same | flip, ]
  flip <- flip[same | flip]
  
  # GWAS effect is for ALT; eQTL effect is for EA. Flip GWAS sign where ALT != EA.
  m$gwas_beta <- ifelse(flip, -m$gwas_beta, m$gwas_beta)
  
  # Clean up: complete rows, positive SE, non-zero MAF, one row per SNP
  m$gwas_maf <- pmin(m$gwas_af, 1 - m$gwas_af)
  m$eqtl_maf <- pmin(m$eqtl_eaf, 1 - m$eqtl_eaf)
  m <- m[complete.cases(m[, c("gwas_beta", "gwas_se", "eqtl_beta", "eqtl_se",
                              "gwas_maf", "eqtl_maf")]), ]
  m <- m[m$gwas_se > 0 & m$eqtl_se > 0 & m$gwas_maf > 0 & m$eqtl_maf > 0, ]
  m <- m[!duplicated(m$snp), ]
  m <- m[order(m$pos), ]
  
  cat("eQTL SNPs:", nrow(eqtl), " GWAS SNPs:", nrow(gwas),
      " Matched:", nrow(m), " (allele-flipped:", sum(flip), ")\n")
  
  # ---- 3. Build coloc datasets ----
  # GWAS: LDL is inverse-normal transformed, so sdY = 1
  D1 <- list(beta = m$gwas_beta, varbeta = m$gwas_se^2,
             snp = m$snp, position = m$pos,
             type = "quant", sdY = 1, MAF = m$gwas_maf,
             pvalues = 2 * pnorm(-abs(m$gwas_beta / m$gwas_se)))
  # eQTL: coloc estimates sdY from N and MAF
  D2 <- list(beta = m$eqtl_beta, varbeta = m$eqtl_se^2,
             snp = m$snp, position = m$pos,
             type = "quant", N = max(m$eqtl_n), MAF = m$eqtl_maf,
             pvalues = 2 * pnorm(-abs(m$eqtl_beta / m$eqtl_se)))
  
  check_dataset(D1)
  check_dataset(D2)
  
  # ---- 4. Plot each dataset ----
  pdf(file.path(DIR, paste0("coloc_", g, "_datasets.pdf")), width = 8, height = 8)
  par(mfrow = c(2, 1))
  plot_dataset(D1, main = paste(toupper(g), "- LDL GWAS"))
  plot_dataset(D2, main = paste(toupper(g), "- liver eQTL"))
  dev.off()
  
  # ---- 5. Run coloc.abf (default priors) ----
  res <- coloc.abf(dataset1 = D1, dataset2 = D2)
  print(res)
  
  # ---- 5b. Sensitivity plot (a04_sensitivity vignette) ----
  # If the Manhattan panels fail, fall back to the p12 panel only
  pdf(file.path(DIR, paste0("coloc_", g, "_sensitivity.pdf")), width = 11, height = 7)
  tryCatch(
    sensitivity(res, rule = "H4 > 0.5"),                                   # exactly as in vignette
    error = function(e1) tryCatch(
      sensitivity(res, rule = "H4 > 0.5", dataset1 = D1, dataset2 = D2),   # datasets passed explicitly
      error = function(e2) {
        message("Manhattan panels failed (", conditionMessage(e2), "); plotting p12 panel only")
        sensitivity(res, rule = "H4 > 0.5", plot.manhattans = FALSE)
      }))
  dev.off()
  
  # ---- 6. SNPs contributing to H4 ----
  cat("\nSNPs with SNP.PP.H4 > 0.01:\n")
  top <- subset(res$results, SNP.PP.H4 > 0.01)
  top <- top[order(top$SNP.PP.H4, decreasing = TRUE), ]
  print(top[, c("snp", "position", "SNP.PP.H4")])
  
  # ---- 7. 95% credible set (meaningful if PP.H4 is high) ----
  o  <- order(res$results$SNP.PP.H4, decreasing = TRUE)
  cs <- cumsum(res$results$SNP.PP.H4[o])
  w  <- which(cs > 0.95)[1]
  cred <- res$results[o, ][1:w, c("snp", "position", "SNP.PP.H4")]
  cat("\n95% credible set (", nrow(cred), "SNPs):\n")
  print(cred)
  
  # ---- 8. Save ----
  write.table(res$results, file.path(DIR, paste0("coloc_", g, "_snp_results.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cred, file.path(DIR, paste0("coloc_", g, "_credible_set.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  all_summaries[[g]] <- c(gene = toupper(g), res$summary)
}

summary_tab <- do.call(rbind, all_summaries)
cat("\n==================== SUMMARY ====================\n")
print(summary_tab)
write.table(summary_tab, file.path(DIR, paste0("coloc_summary_", paste(genes, collapse = "_"), ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

