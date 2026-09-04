# Inputs: data/Human.MitoCarta3.0.xls (Rath et al. 2021, https://www.broadinstitute.org/mitocarta)
library(data.table); library(readxl); library(cluster); library(ggplot2); library(patchwork)

COL_HYPER <- "#D7301F"; COL_MITO_LOW <- "#1B9E77"; COL_MITO_HIGH <- "#7570B3"
theme_pub <- theme_bw(base_size=11) + theme(
  panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth=0.25, color="grey85"),
  panel.border = element_rect(color="black", linewidth=0.5), plot.title = element_text(size=11, face="bold"),
  legend.position = "top", legend.title = element_blank(), axis.text = element_text(color="black"))

## ---- 59-protein mitochondrial panel: Olink assay list x MitoCarta3.0 ----
npx <- fread("data/MGH_COVID_OLINK_NPX.txt", sep=";")
mc3 <- as.data.table(suppressWarnings(read_excel("data/Human.MitoCarta3.0.xls", sheet="A Human MitoCarta3.0")))
mito_genes <- sort(intersect(unique(npx$Assay), unique(mc3$Symbol)))
writeLines(mito_genes, "mito_full_panel.txt")
cat("N mitochondrial-associated proteins:", length(mito_genes), "\n")

npx_mito <- npx[Timepoint=="D0" & Assay %in% mito_genes]
npx_mito[, npx_clean := ifelse(QC_Warning=="PASS" & Assay_Warning=="PASS", NPX, NA_real_)]
wide_mito <- dcast(npx_mito, subject_id ~ Assay, value.var="npx_clean", fun.aggregate=mean, na.rm=TRUE)

m4 <- fread("final_merged_v4_restricted.csv")
m4[, subject_id := as.character(subject_id)]
wide_mito[, subject_id := as.character(subject_id)]
mc <- merge(m4, wide_mito, by="subject_id")
mc[, phenotype := factor(phenotype, levels=c("Hypoinflammatory","Hyperinflammatory"))]
cat("N =", nrow(mc), "; with >=1 missing protein:", sum(!complete.cases(mc[, ..mito_genes])), "\n")

## continuous score: available-protein average (no imputation)
zmito <- scale(mc[, ..mito_genes])
mc[, mito_score_full := rowMeans(zmito, na.rm=TRUE)]

## PCA/k-means: zero-imputed (standardized-scale mean) to retain the full cohort
zmito_imp <- zmito
zmito_imp[is.na(zmito_imp)] <- 0
pca <- prcomp(zmito_imp, center=FALSE, scale.=FALSE)
pc3 <- pca$x[, 1:3]

set.seed(42)
sil <- sapply(2:8, function(k) { km <- kmeans(pc3, centers=k, nstart=100); mean(silhouette(km$cluster, dist(pc3))[,3]) })
km2 <- kmeans(pc3, centers=2, nstart=100)
mc[, mito_cluster := factor(km2$cluster)]
cl_means <- mc[, .(mean_score=mean(mito_score_full)), by=mito_cluster][order(mean_score)]
high_cl <- cl_means$mito_cluster[2]
mc[, mito_named := factor(fifelse(mito_cluster==high_cl,"Hypermitochondrial","Hypomitochondrial"),
                           levels=c("Hypomitochondrial","Hypermitochondrial"))]
mc[, Acuity_max := as.numeric(Acuity_max)]
mc[, death := as.numeric(Acuity_max==1)]
fwrite(mc, "final_analytic_dataset.csv")
cat("N per mito cluster:\n"); print(table(mc$mito_named))

print(cor.test(mc$mito_score_full, mc$CPROB_hyper, method="spearman"))
print(mc[, .(median=median(mito_score_full), q1=quantile(mito_score_full,.25), q3=quantile(mito_score_full,.75)), by=.(Acuity_0_cat==1)])
print(wilcox.test(mito_score_full ~ (Acuity_0_cat==1), data=mc))
print(mc[, .(median=median(mito_score_full), q1=quantile(mito_score_full,.25), q3=quantile(mito_score_full,.75)), by=death])
print(wilcox.test(mito_score_full ~ death, data=mc))

## ---- Figure S3: PCA diagnostics ----
ve <- summary(pca)$importance[2,1:10]*100
scree_dt <- data.table(PC=factor(1:10), var=ve)
pS3a <- ggplot(scree_dt, aes(x=PC, y=var)) +
  geom_col(fill="grey60", color="black", linewidth=0.3, width=0.7) +
  geom_text(aes(label=sprintf("%.1f%%", var)), vjust=-0.4, size=2.8) +
  scale_y_continuous(limits=c(0,50), expand=expansion(mult=c(0,0.12))) +
  labs(x="Principal component", y="Variance explained (%)", title="A") + theme_pub
sil_dt <- data.table(k=2:8, sil=sil)
pS3b <- ggplot(sil_dt, aes(x=k, y=sil)) +
  geom_line(linewidth=0.6, color="grey40") + geom_point(size=2.5, color="black") +
  geom_point(data=sil_dt[k==2], size=4, shape=21, fill=COL_MITO_HIGH, color="black", stroke=1) +
  scale_x_continuous(breaks=2:8) + scale_y_continuous(limits=c(0,0.5)) +
  labs(x="Number of clusters (k)", y="Average silhouette width", title="B") + theme_pub
figS3 <- pS3a | pS3b
ggsave("FigureS3_mito_PCA_diagnostics.png", figS3, width=8.5, height=4, dpi=350, bg="white")

## ---- Table S2: crude vs mutually adjusted odds ratios ----
extract1 <- function(mod, term) {
  co <- exp(cbind(OR=coef(mod), confint(mod))); p <- summary(mod)$coefficients[,4]
  sprintf("%.2f [%.2f-%.2f], p%s", co[term,1], co[term,2], co[term,3],
          ifelse(p[term]<0.001,"<0.001",sprintf("=%.2f",p[term])))
}
m_c1 <- glm(death ~ phenotype, data=mc, family=binomial)
m_c2 <- glm(death ~ scale(mito_score_full), data=mc, family=binomial)
m_c3 <- glm(death ~ mito_named, data=mc, family=binomial)
m_a1 <- glm(death ~ phenotype + scale(mito_score_full), data=mc, family=binomial)
m_a2 <- glm(death ~ phenotype + mito_named, data=mc, family=binomial)

tabS2 <- data.table(
  Term=c("Hyperinflammatory (vs Hypoinflammatory)","Mitochondrial score, continuous (per 1-SD)","Hypermitochondrial cluster (vs Hypomitochondrial)"),
  Crude=c(extract1(m_c1,"phenotypeHyperinflammatory"), extract1(m_c2,"scale(mito_score_full)"), extract1(m_c3,"mito_namedHypermitochondrial")),
  Adjusted=c(extract1(m_a1,"phenotypeHyperinflammatory"), extract1(m_a1,"scale(mito_score_full)"), extract1(m_a2,"mito_namedHypermitochondrial"))
)
fwrite(tabS2, "TableS2_crude_vs_adjusted_OR.csv")
print(tabS2)

mod_int <- glm(death ~ phenotype * scale(mito_score_full), data=mc, family=binomial)
cat("LRT interaction p =", anova(m_a1, mod_int, test="Chisq")[2,"Pr(>Chi)"], "\n")

## ---- Table S4: confusion matrix ----
conf <- table(mc$phenotype, mc$mito_named)
Po <- (conf[1,1]+conf[2,2])/sum(conf)
Pe <- (sum(conf[1,])*sum(conf[,1]) + sum(conf[2,])*sum(conf[,2])) / sum(conf)^2
kappa <- (Po-Pe)/(1-Pe)
ft <- fisher.test(conf)
tabS4 <- as.data.table(conf); setnames(tabS4, c("Phenotype","Mito_cluster","N"))
fwrite(tabS4, "TableS4_confusion_matrix.csv")
cat(sprintf("Observed agreement %.1f%% (%d/%d); kappa=%.2f; Fisher OR=%.2f [%.2f-%.2f], p=%.3g\n",
    100*Po, conf[1,1]+conf[2,2], sum(conf), kappa, ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
