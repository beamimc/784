# Simple Manhattan plot of the GLGC 2021 LDL-C GWAS (European ancestry).
# Positions are GRCh37. Run from the project root:
#   Rscript scripts/manhattan_simple.R

library(data.table)
library(ggplot2)

gwas_path <- "LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"

# 1. Read only the three columns a Manhattan plot needs. We use
#    pvalue_neg_log10 because the raw pvalue underflows to 0 at strong loci.
gwas <- fread(gwas_path, select = c("CHROM", "POS_b37", "pvalue_neg_log10"))
setnames(gwas, c("chr", "pos", "logp"))

# 2. Thin the data: keep every variant with P < 0.01 and a 1% random sample
#    of the rest. Most of the 47M variants sit near the bottom of the plot,
#    so it looks the same but draws much faster.
set.seed(784)
gwas <- gwas[logp > 2 | runif(.N) < 0.01]

# 3. Plot. One facet per chromosome puts them side by side without computing
#    cumulative positions; alternating colors separate neighboring chromosomes.
p <- ggplot(gwas, aes(pos, logp, colour = factor(chr %% 2))) +
  geom_point(size = 0.4) +
  geom_hline(
    yintercept = -log10(5e-8), linetype = "dashed", colour = "firebrick"
  ) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x", switch = "x") +
  scale_colour_manual(values = c("grey30", "steelblue"), guide = "none") +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_sqrt(
    breaks = c(10, 100, 1000, 2500, 5000, 7500),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(P))),
    title = "LDL-C GWAS, European ancestry (GLGC 2021)",
    subtitle = "Square-root y-axis. Dashed line: P = 5e-8"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 10)
  )

ggsave("LDL_EUR_manhattan_simple.png", p, width = 14, height = 5, dpi = 300)
