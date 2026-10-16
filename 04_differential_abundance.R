library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)

COL_HYPO <- "#2C7FB8"; COL_HYPER <- "#D7301F"; COL_NS <- "grey65"
theme_pub <- theme_bw(base_size=11) + theme(
  panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth=0.25, color="grey85"),
  panel.border = element_rect(color="black", linewidth=0.5), plot.title = element_text(size=11, face="bold"),
  axis.text = element_text(color="black"))

mc <- fread("final_analytic_dataset.csv")
mc[, phenotype := factor(phenotype, levels=c("Hypoinflammatory","Hyperinflammatory"))]
mito_genes <- readLines("mito_full_panel.txt")

zmat <- scale(mc[, ..mito_genes])
res <- rbindlist(lapply(mito_genes, function(g) {
  x_hyper <- zmat[mc$phenotype=="Hyperinflammatory", g]
  x_hypo  <- zmat[mc$phenotype=="Hypoinflammatory", g]
  wt <- wilcox.test(x_hyper, x_hypo)
  data.table(protein=g, mean_diff=mean(x_hyper,na.rm=TRUE)-mean(x_hypo,na.rm=TRUE), p=wt$p.value)
}))
res[, q := p.adjust(p, method="BH")]
setorder(res, q)
res[, col := fifelse(q>=0.05, COL_NS, fifelse(mean_diff>0, COL_HYPER, COL_HYPO))]
cat("N FDR-significant (q<0.05):", sum(res$q<0.05), "/", nrow(res), "\n")
fwrite(res[, .(Protein=protein, `Standardized mean difference`=mean_diff,
               `p-value`=p, `FDR q-value`=q)],
       "TableS3_differential_abundance.csv")

res[, lbl := ""]
top_proteins <- res[order(q)][1:15, protein]
res[protein %in% top_proteins, lbl := protein]
pS3a <- ggplot(res, aes(x=mean_diff, y=-log10(q), color=col)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50", linewidth=0.4) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50", linewidth=0.4) +
  geom_point(size=2, alpha=0.85) + scale_color_identity() +
  ggrepel::geom_text_repel(aes(label=lbl), size=2.6, max.overlaps=30, min.segment.length=0, segment.size=0.3, color="black") +
  labs(x="Standardized mean difference (hyperinflammatory - hypoinflammatory)",
       y=expression(-log[10]~"(FDR q-value)"), title="A") + theme_pub

topN <- res[order(-abs(mean_diff))][1:20]
topN[, protein := factor(protein, levels=rev(protein[order(mean_diff)]))]
pS3b <- ggplot(topN, aes(x=mean_diff, y=protein, color=col)) +
  geom_vline(xintercept=0, color="grey50", linewidth=0.4) +
  geom_segment(aes(x=0, xend=mean_diff, y=protein, yend=protein), linewidth=0.6) +
  geom_point(size=2.6) + scale_color_identity() +
  labs(x="Standardized mean difference", y=NULL, title="B") + theme_pub

figS3 <- pS3a | pS3b
ggsave("FigureS3_differential_abundance.png", figS3, width=11, height=5.2, dpi=350, bg="white")
