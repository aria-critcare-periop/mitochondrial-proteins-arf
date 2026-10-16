# Input: data/Human.MitoCarta3.0.xls (Rath et al. 2021)
library(data.table); library(readxl); library(cluster)

COL_HYPER <- "#D7301F"; COL_MITO_LOW <- "#1B9E77"; COL_MITO_HIGH <- "#7570B3"

## ---- 59-protein mitochondrial panel: Olink assay list x MitoCarta3.0 ----
npx <- fread("data/MGH_COVID_OLINK_NPX.txt", sep=";")
mc3 <- as.data.table(suppressWarnings(read_excel("data/Human.MitoCarta3.0.xls", sheet="A Human MitoCarta3.0")))
mito_genes <- sort(intersect(unique(npx$Assay), unique(mc3$Symbol)))
stopifnot(length(mito_genes)==59, all(mito_genes %in% unique(npx$Assay)))
writeLines(mito_genes, "mito_full_panel.txt")
cat("N mitochondrial-associated proteins:", length(mito_genes), "\n")

npx_mito <- npx[Timepoint=="D0" & Assay %in% mito_genes]
npx_mito[, npx_clean := ifelse(QC_Warning=="PASS" & Assay_Warning=="PASS", NPX, NA_real_)]
wide_mito <- dcast(npx_mito, subject_id ~ Assay, value.var="npx_clean", fun.aggregate=mean, na.rm=TRUE)

m4 <- fread("final_merged.csv")
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
print(mc[, .(median=median(mito_score_full), q1=quantile(mito_score_full,.25),
             q3=quantile(mito_score_full,.75)), by=phenotype])
print(mc[, .(median=median(mito_score_full), q1=quantile(mito_score_full,.25), q3=quantile(mito_score_full,.75)), by=.(Acuity_0_cat==1)])
print(wilcox.test(mito_score_full ~ (Acuity_0_cat==1), data=mc))
print(mc[, .(median=median(mito_score_full), q1=quantile(mito_score_full,.25), q3=quantile(mito_score_full,.75)), by=death])
print(wilcox.test(mito_score_full ~ death, data=mc))

## ---- Figure S4: PCA diagnostics ----
ve <- 100*pca$sdev^2/sum(pca$sdev^2)
png("FigureS4_PCA_diagnostics.png", width=2900, height=1250, res=300)
par(mfrow=c(1,2), mar=c(5,5,3,1), cex.axis=0.8, cex.lab=1)
plot(seq_along(ve), ve, type="h", lwd=3,
     col=ifelse(seq_along(ve)<=3,COL_MITO_HIGH,"gray60"),
     xlab="Principal component", ylab="Variance explained (%)",
     main="A  All 59 principal components", xaxt="n")
axis(1, at=c(1,3,10,20,30,40,50,59))
points(1:3, ve[1:3], pch=19, col=COL_MITO_HIGH)
legend("topright",sprintf("PC1-PC3: %.1f%%",sum(ve[1:3])), bty="n",
       text.col=COL_MITO_HIGH, cex=0.9)
plot(2:8, sil, type="b", pch=19, col=COL_MITO_HIGH, ylim=c(0,0.46),
     xlab="Number of clusters", ylab="Average silhouette width",
     main="B  Cluster separation")
dev.off()

## ---- Table S5: crude vs mutually adjusted odds ratios ----
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

tabS5 <- data.table(
  Term=c("Hyperinflammatory (vs Hypoinflammatory)","Mitochondrial score, continuous (per 1-SD)","Hypermitochondrial cluster (vs Hypomitochondrial)"),
  Crude=c(extract1(m_c1,"phenotypeHyperinflammatory"), extract1(m_c2,"scale(mito_score_full)"), extract1(m_c3,"mito_namedHypermitochondrial")),
  Adjusted=c(extract1(m_a1,"phenotypeHyperinflammatory"), extract1(m_a1,"scale(mito_score_full)"), extract1(m_a2,"mito_namedHypermitochondrial"))
)
fwrite(tabS5[, .(Term, `Crude odds ratio`=Crude,
                  `Mutually adjusted odds ratio`=Adjusted)],
       "TableS5_crude_vs_adjusted_OR.csv")
print(tabS5)

mod_int <- glm(death ~ phenotype * scale(mito_score_full), data=mc, family=binomial)
cat("LRT interaction p =", anova(m_a1, mod_int, test="Chisq")[2,"Pr(>Chi)"], "\n")

## ---- Table S4: crossed phenotypes and mortality ----
conf <- table(mc$phenotype, mc$mito_named)
Po <- (conf[1,1]+conf[2,2])/sum(conf)
Pe <- (sum(conf[1,])*sum(conf[,1]) + sum(conf[2,])*sum(conf[,2])) / sum(conf)^2
kappa <- (Po-Pe)/(1-Pe)
ft <- fisher.test(conf)
cell <- function(p, c) {
  x <- mc[phenotype==p & mito_named==c]
  n <- nrow(x); deaths <- sum(x$death)
  sprintf("%d (%.0f%%); deaths %d/%d (%.1f%%)", n,
          100*n/sum(mc$phenotype==p), deaths, n, 100*deaths/n)
}
tabS4 <- data.table(
  Phenotype=c("Hypoinflammatory","Hyperinflammatory","Total"),
  Hypomitochondrial=c(cell("Hypoinflammatory","Hypomitochondrial"),
                      cell("Hyperinflammatory","Hypomitochondrial"),"179"),
  Hypermitochondrial=c(cell("Hypoinflammatory","Hypermitochondrial"),
                       cell("Hyperinflammatory","Hypermitochondrial"),"106"),
  Total=c("185","100","285"))
fwrite(tabS4, "TableS4_crossed_phenotypes.csv")
cat(sprintf("Observed agreement %.1f%% (%d/%d); kappa=%.2f; Fisher OR=%.2f [%.2f-%.2f], p=%.3g\n",
    100*Po, conf[1,1]+conf[2,2], sum(conf), kappa, ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))

## ---- Figure S2: individual scores ----
png("FigureS2_individual_scores.png", width=2600, height=1250, res=300)
par(mfrow=c(1,2), mar=c(5,5,3,1), cex.axis=0.9, cex.lab=1.1)
boxplot(mito_score_full~phenotype, data=mc, outline=FALSE, col=c("#A7C8DE","#E9998E"), xaxt="n",
        ylab="Mean standardized protein abundance", xlab="Inflammatory phenotype",
        main="A  Individual score distribution")
axis(1, at=1:2, labels=c("Hypoinflammatory","Hyperinflammatory"), cex.axis=0.75)
points(jitter(as.numeric(mc$phenotype), amount=0.16), mc$mito_score_full,
       pch=16, cex=0.5, col=adjustcolor("black", alpha.f=0.35))
plot(mc$CPROB_hyper, mc$mito_score_full, pch=16, cex=0.6,
     col=adjustcolor("#424242", alpha.f=0.5),
     xlab="Posterior hyperinflammatory probability", ylab="Mitochondrial protein score",
     main="B  Score and inflammatory probability")
abline(lm(mito_score_full~CPROB_hyper, data=mc), col=COL_HYPER, lwd=2)
dev.off()

## ---- Figure S5: protein loadings ----
png("FigureS5_protein_loadings.png", width=2700, height=1800, res=300)
par(mfrow=c(1,3), mar=c(5,7,3,1), cex.axis=0.8, cex.lab=0.9)
for (j in 1:3) {
  selected <- order(abs(pca$rotation[,j]), decreasing=TRUE)[1:12]
  vals <- pca$rotation[selected,j]
  ord <- order(vals)
  barplot(vals[ord], horiz=TRUE, names.arg=mito_genes[selected][ord], las=1,
          col=ifelse(vals[ord]>=0,COL_MITO_HIGH,COL_MITO_LOW), border=NA,
          xlab="PCA loading", main=paste0("PC",j," (",sprintf("%.1f",ve[j]),"%)"))
  abline(v=0, col="grey40")
}
dev.off()

## ---- Table S6: sensitivity analyses ----
agreement <- function(a,b) max(mean(a==b), mean(a!=b))
set.seed(42)
km_pc2 <- kmeans(pca$x[,1:2], centers=2, nstart=100)$cluster
set.seed(42)
km_pc4 <- kmeans(pca$x[,1:4], centers=2, nstart=100)$cluster
z_median <- zmito
for (j in seq_len(ncol(z_median))) {
  z_median[is.na(z_median[,j]),j] <- median(z_median[,j], na.rm=TRUE)
}
set.seed(42)
km_median <- kmeans(prcomp(z_median, center=FALSE, scale.=FALSE)$x[,1:3],
                    centers=2, nstart=100)$cluster
complete <- complete.cases(zmito)
set.seed(42)
km_complete <- kmeans(prcomp(zmito[complete,], center=FALSE, scale.=FALSE)$x[,1:3],
                      centers=2, nstart=100)$cluster
median_score <- apply(zmito, 1, median, na.rm=TRUE)
valid <- mc[Acuity_0!=1]
valid[, vent := as.numeric(Acuity_0==2)]
fit_score <- glm(death~scale(mito_score_full)+Age_cat+COVID+vent+creat0c+neut0c,
                 data=valid, family=binomial)
fit_cluster <- glm(death~mito_named+Age_cat+COVID+vent+creat0c+neut0c,
                   data=valid, family=binomial)
wald_or <- function(fit, term) {
  x <- summary(fit)$coefficients[term,]
  sprintf("OR %.2f (95%% CI %.2f-%.2f), p=%.3f", exp(x[1]), exp(x[1]-1.96*x[2]),
          exp(x[1]+1.96*x[2]), x[4])
}
tabS6 <- data.table(
  Analysis=c("Two versus three PCs, cluster agreement", "Four versus three PCs, cluster agreement",
             "Median versus zero imputation, cluster agreement", "Complete-case PCA, cluster agreement",
             "Mean versus median protein score, Spearman rho", "Clinical-adjusted mitochondrial abundance score",
             "Clinical-adjusted mitochondrial cluster"),
  Result=c(sprintf("%.1f%%",100*agreement(km_pc2,km2$cluster)),
           sprintf("%.1f%%",100*agreement(km_pc4,km2$cluster)),
           sprintf("%.1f%%",100*agreement(km_median,km2$cluster)),
           sprintf("%.1f%%",100*agreement(km_complete,km2$cluster[complete])),
           sprintf("%.3f",cor(mc$mito_score_full,median_score,method="spearman")),
           wald_or(fit_score,"scale(mito_score_full)"),
           wald_or(fit_cluster,"mito_namedHypermitochondrial")))
fwrite(tabS6, "TableS6_sensitivity_analyses.csv")
