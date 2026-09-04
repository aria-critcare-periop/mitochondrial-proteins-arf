library(data.table); library(ggplot2); library(patchwork)

COL_HYPO <- "#2C7FB8"; COL_HYPER <- "#D7301F"
COL_MITO_LOW <- "#1B9E77"; COL_MITO_HIGH <- "#7570B3"
theme_pub <- theme_classic(base_size=12) + theme(
  plot.title = element_text(size=12, face="bold"),
  legend.position = "top", legend.title = element_blank(), legend.text = element_text(size=9.5, color="black"),
  axis.title = element_text(size=10.5, color="black"), axis.text = element_text(size=9.5, color="black"),
  axis.ticks = element_line(color="black"))

## ---- Panel A: standardized biomarker profile ----
target_proteins <- c("IL6","CXCL8","TNFRSF1A","PROC","SERPINE1","ICAM1","AGER","SFTPD","VWF")
lbl_map <- c(IL6="IL-6", CXCL8="IL-8", TNFRSF1A="sTNFR1", PROC="Protein C", SERPINE1="PAI-1",
             ICAM1="sICAM-1", AGER="sRAGE", SFTPD="SP-D", VWF="vWF")
m4 <- fread("final_merged_v4_restricted.csv")
zcols <- paste0("z_", target_proteins)
zdt <- as.data.table(scale(m4[, ..zcols])); setnames(zdt, zcols, target_proteins)
zdt[, phenotype := m4$phenotype]
long <- melt(zdt, id.vars="phenotype", variable.name="protein", value.name="z")
long[, protein_lbl := lbl_map[as.character(protein)]]
summ <- long[, .(mean_z=mean(z,na.rm=TRUE)), by=.(phenotype,protein_lbl)]
ordr <- summ[phenotype=="Hyperinflammatory"][order(-mean_z), protein_lbl]
summ[, protein_lbl := factor(protein_lbl, levels=ordr)]
pA <- ggplot(summ, aes(x=protein_lbl, y=mean_z, color=phenotype, group=phenotype)) +
  geom_hline(yintercept=0, color="grey50", linewidth=0.3) +
  geom_line(linewidth=0.9) + geom_point(size=2.2) +
  scale_color_manual(values=c("Hypoinflammatory"=COL_HYPO,"Hyperinflammatory"=COL_HYPER)) +
  scale_y_continuous(limits=c(-1,1)) +
  labs(x=NULL, y="Standardized value (z-score)", title="A") +
  theme_pub + theme(axis.text.x=element_text(angle=45,hjust=1,color="black"))

## ---- Data ----
mc <- fread("final_analytic_dataset.csv")
mc[, phenotype := factor(phenotype, levels=c("Hypoinflammatory","Hyperinflammatory"))]
mito_genes <- readLines("mito_full_panel.txt")
zmito <- scale(mc[, ..mito_genes])
zmito_imp <- zmito
zmito_imp[is.na(zmito_imp)] <- 0
pca <- prcomp(zmito_imp, center=FALSE, scale.=FALSE)

## ---- Panel B: PCA scatter ----
pcdt <- as.data.table(pca$x[,1:2]); pcdt[, mito_named := mc$mito_named]
ve <- summary(pca)$importance[2,1:2]*100
pB <- ggplot(pcdt, aes(x=PC1, y=PC2, color=mito_named)) +
  geom_point(size=1.8, alpha=0.85) +
  scale_color_manual(values=c("Hypomitochondrial"=COL_MITO_LOW,"Hypermitochondrial"=COL_MITO_HIGH)) +
  labs(x=sprintf("PC1 (%.1f%%)", ve[1]), y=sprintf("PC2 (%.1f%%)", ve[2]), title="B") +
  theme_pub

## ---- Panel C: mito score by phenotype ----
wt <- wilcox.test(mito_score_full ~ phenotype, data=mc)
pC <- ggplot(mc, aes(x=phenotype, y=mito_score_full, fill=phenotype)) +
  geom_boxplot(outlier.shape=NA, alpha=0.75) + geom_jitter(width=0.15, alpha=0.35, size=1.3) +
  scale_fill_manual(values=c("Hypoinflammatory"=COL_HYPO,"Hyperinflammatory"=COL_HYPER), guide="none") +
  labs(x=NULL, y="Continuous mitochondrial score", title="C") +
  annotate("text", x=1.5, y=max(mc$mito_score_full,na.rm=TRUE)*1.08,
           label=ifelse(wt$p.value<0.001,"p<0.001",sprintf("p=%.3f", wt$p.value)), size=3.2) +
  theme_pub + theme(axis.text.x=element_text(size=8.5, angle=20, hjust=1))

## ---- Panel D: crude vs mutually adjusted odds ratios ----
mc[, mito_named := factor(mito_named, levels=c("Hypomitochondrial","Hypermitochondrial"))]
mod_cont <- glm(death ~ phenotype + scale(mito_score_full), data=mc, family=binomial)
co1 <- exp(cbind(OR=coef(mod_cont), confint(mod_cont)))[-1,]
p1 <- summary(mod_cont)$coefficients[-1,4]
mod_cat <- glm(death ~ phenotype + mito_named, data=mc, family=binomial)
co2 <- exp(cbind(OR=coef(mod_cat), confint(mod_cat)))[-1,]
p2 <- summary(mod_cat)$coefficients[-1,4]

extract1 <- function(mod, term_name) {
  co <- exp(cbind(OR=coef(mod), confint(mod)))
  p <- summary(mod)$coefficients[,4]
  list(OR=unname(co[term_name,"OR"]), lo=unname(co[term_name,2]), hi=unname(co[term_name,3]), p=unname(p[term_name]))
}
e_c1 <- extract1(glm(death ~ phenotype, data=mc, family=binomial), "phenotypeHyperinflammatory")
e_c2 <- extract1(glm(death ~ scale(mito_score_full), data=mc, family=binomial), "scale(mito_score_full)")
e_c3 <- extract1(glm(death ~ mito_named, data=mc, family=binomial), "mito_namedHypermitochondrial")

fdt <- data.table(
  term=rep(c("Hyperinflammatory\n(vs Hypoinflammatory)","Mitochondrial score, continuous\n(per 1-SD)",
         "Hypermitochondrial cluster\n(vs Hypomitochondrial)"), each=2),
  model=rep(c("Crude","Mutually\nadjusted"), 3),
  OR=c(e_c1$OR, co1[1,1], e_c2$OR, co1[2,1], e_c3$OR, co2[2,1]),
  lo=c(e_c1$lo, co1[1,2], e_c2$lo, co1[2,2], e_c3$lo, co2[2,2]),
  hi=c(e_c1$hi, co1[1,3], e_c2$hi, co1[2,3], e_c3$hi, co2[2,3]),
  p=c(e_c1$p, p1[1], e_c2$p, p1[2], e_c3$p, p2[2]),
  col=rep(c(COL_HYPER, "#7570B3", "grey45"), each=2))
fdt[, term := factor(term, levels=rev(unique(term)))]
fdt[, model := factor(model, levels=c("Crude","Mutually\nadjusted"))]
fdt[, plab := ifelse(p<0.001, "p<0.001", sprintf("p=%.2f", p))]

pD <- ggplot(fdt, aes(x=OR, y=term, group=model)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50", linewidth=0.4) +
  geom_errorbarh(aes(xmin=lo, xmax=hi, alpha=model), height=0.15, linewidth=0.55, color="black",
                 position=position_dodge(width=0.55)) +
  geom_point(aes(color=col, shape=model, alpha=model), size=3, position=position_dodge(width=0.55)) +
  scale_color_identity() +
  scale_shape_manual(values=c("Crude"=1,"Mutually\nadjusted"=16), name=NULL) +
  scale_alpha_manual(values=c("Crude"=0.55,"Mutually\nadjusted"=1), guide="none") +
  geom_text(aes(label=sprintf("%.2f [%.2f-%.2f], %s", OR, lo, hi, plab)), position=position_dodge(width=0.55),
            hjust=-0.06, vjust=-0.7, size=2.3) +
  scale_x_log10(breaks=c(0.5,1,2,5,10,20), limits=c(0.3,65)) +
  labs(x="Odds ratio for 28-day mortality; crude vs mutually adjusted", y=NULL, title="D") + theme_pub +
  theme(axis.text.y=element_text(size=8.3), legend.position="bottom", legend.text=element_text(size=8))

## ---- Panel E: crossed 2x2 mortality ----
cross_tab <- mc[, .(N=.N, Death=sum(death), death_pct=mean(death)*100), by=.(phenotype, mito_named)]
cross_tab[, label := paste0(as.character(phenotype), "\n", as.character(mito_named))]
cross_tab[, label := factor(label, levels=unique(label[order(phenotype,mito_named)]))]
cross_tab[, ntxt := sprintf("%d/%d (%.0f%%)", Death, N, death_pct)]
pE <- ggplot(cross_tab, aes(x=label, y=death_pct, fill=phenotype, alpha=mito_named)) +
  geom_col(width=0.65, color="black", linewidth=0.35) +
  geom_text(aes(label=ntxt), vjust=-0.3, size=2.5, lineheight=0.8, color="black", alpha=1) +
  scale_fill_manual(values=c("Hypoinflammatory"=COL_HYPO,"Hyperinflammatory"=COL_HYPER), guide="none") +
  scale_alpha_manual(values=c("Hypomitochondrial"=0.45,"Hypermitochondrial"=1), guide="none") +
  scale_y_continuous(limits=c(0,45), expand=expansion(mult=c(0,0.12))) +
  labs(x=NULL, y="28-day mortality (%)", title="E") +
  theme_pub + theme(axis.text.x=element_text(size=6.3, color="black", angle=30, hjust=1, vjust=1, lineheight=0.8))

## ---- Assemble ----
top <- pA + pB + plot_layout(widths=c(1.15,1))
bottom <- pC + pD + pE + plot_layout(widths=c(0.7,1.2,1.2))
fig1 <- top / bottom + plot_layout(heights=c(1,0.95))

ggsave("Figure1_5panel.png", fig1, width=12.5, height=9.2, dpi=350, bg="white")
ggsave("Figure1_5panel.tiff", fig1, width=12.5, height=9.2, dpi=350, bg="white", compression="lzw")
