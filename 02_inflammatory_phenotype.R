# Run after Mplus has produced lca_v4_k1.out ... lca_v4_k5.out and class_v4_k2.dat
library(data.table)
library(MplusAutomation)
library(ggplot2)
library(patchwork)

target_proteins <- c("IL6","CXCL8","TNFRSF1A","PROC","SERPINE1","ICAM1","AGER","SFTPD","VWF")
zvars <- paste0("z_", target_proteins)
COL_HYPO <- "#2C7FB8"; COL_HYPER <- "#D7301F"
theme_pub <- theme_bw(base_size=11) + theme(
  panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth=0.25, color="grey85"),
  panel.border = element_rect(color="black", linewidth=0.5), plot.title = element_text(size=11, face="bold"),
  axis.text = element_text(color="black"))

## ---- class enumeration across k=1..5: Table S1 + Figure S1 ----
p4 <- readRDS("primary_sample_v4.rds")
models <- readModels(target=getwd(), filefilter="lca_v4_k")
enum <- data.table()
for (k in 1:5) {
  key <- grep(sprintf("lca_v4_k%d.out$", k), names(models), ignore.case=TRUE, value=TRUE)
  m <- models[[key]]; s <- m$summaries; cp <- m$class_counts$mostLikely
  enum <- rbind(enum, data.table(Classes=k, LL=s$LL, BIC=s$BIC, aBIC=s$aBIC,
    Entropy=ifelse(is.null(s$Entropy),NA,s$Entropy),
    VLMR_p=ifelse(is.null(s$T11_LMR_PValue),NA,s$T11_LMR_PValue),
    N_per_class=ifelse(is.null(cp), NA, paste(cp$count, collapse="/"))), fill=TRUE)
}
fwrite(enum, "TableS1_LCA_enumeration.csv")

pS1a <- ggplot(enum[!is.na(BIC)], aes(x=Classes, y=BIC)) +
  geom_line(color="grey40") + geom_point(size=2.5) +
  geom_point(data=enum[Classes==2], size=4, shape=21, fill=COL_HYPER, color="black", stroke=1) +
  scale_x_continuous(breaks=1:5) + labs(x="Number of classes (k)", y="BIC", title="A") + theme_pub
pS1b <- ggplot(enum[!is.na(Entropy)], aes(x=Classes, y=Entropy)) +
  geom_line(color="grey40") + geom_point(size=2.5) +
  geom_point(data=enum[Classes==2], size=4, shape=21, fill=COL_HYPER, color="black", stroke=1) +
  geom_hline(yintercept=0.8, linetype="dashed", color="grey50") +
  scale_x_continuous(breaks=1:5) + scale_y_continuous(limits=c(0.75,1)) +
  labs(x="Number of classes (k)", y="Entropy", title="B") + theme_pub
pS1c <- ggplot(enum[!is.na(VLMR_p)], aes(x=Classes, y=VLMR_p)) +
  geom_col(fill="grey60", color="black", linewidth=0.3, width=0.6) +
  geom_hline(yintercept=0.05, linetype="dashed", color=COL_HYPER) +
  geom_text(aes(label=sprintf("%.4f", VLMR_p)), vjust=-0.4, size=2.8) +
  scale_x_continuous(breaks=1:5) + scale_y_continuous(limits=c(0,0.6)) +
  labs(x="Number of classes (k, vs k-1)", y="VLMR p-value", title="C") + theme_pub
figS1 <- pS1a | pS1b | pS1c
ggsave("FigureS1_LCA_enumeration_diagnostics.png", figS1, width=11, height=3.8, dpi=350, bg="white")

## ---- phenotype from the retained k=2 model ----
cls <- fread("class_v4_k2.dat", header=FALSE)
setnames(cls, c(zvars, "acuity0c","creat0c","neut0c","CPROB1","CPROB2","C","mplus_id"))
m <- merge(p4, cls[,.(mplus_id,C,CPROB1,CPROB2)], by="mplus_id")
means_by_class <- m[, lapply(.SD, mean, na.rm=TRUE), by=C, .SDcols=zvars]
hyper_c <- means_by_class[which.max(rowMeans(means_by_class[,..zvars])), C]
m[, phenotype := factor(fifelse(C==hyper_c, "Hyperinflammatory","Hypoinflammatory"),
                         levels=c("Hypoinflammatory","Hyperinflammatory"))]
m[, CPROB_hyper := if (hyper_c==2) CPROB2 else CPROB1]
fwrite(m, "final_merged_v4_restricted.csv")
cat("N =", nrow(m), "\n"); print(m[, .N, by=phenotype])

ft_vent <- fisher.test(table(m$phenotype, m$Acuity_0_cat==1))
cat("Invasive ventilation by phenotype:\n"); print(table(m$phenotype, m$Acuity_0_cat==1))
cat("Fisher p =", ft_vent$p.value, "\n")

for (pr in target_proteins) {
  wt <- wilcox.test(as.formula(paste0(pr," ~ phenotype")), data=m)
  cat(pr, "p =", wt$p.value, "\n")
}

## ---- Table 1 ----
m[, Age_cat := as.numeric(Age_cat)]; m[, BMI_cat := as.numeric(BMI_cat)]
m[, creat_0_cat := as.numeric(creat_0_cat)]; m[, abs_neut_0_cat := as.numeric(abs_neut_0_cat)]
m[, Acuity_max := as.numeric(Acuity_max)]
hy <- m[phenotype=="Hypoinflammatory"]; he <- m[phenotype=="Hyperinflammatory"]

rows <- list()
addrow <- function(...) rows[[length(rows)+1]] <<- data.table(...)

cat_row <- function(var, lbl, levels_lbl) {
  p <- wilcox.test(as.formula(paste0(var," ~ phenotype")), data=m)$p.value
  addrow(Variable=lbl, Hypo="", Hyper="", p=p)
  codes <- sort(unique(m[[var]][!is.na(m[[var]])]))
  stopifnot(length(codes) == length(levels_lbl))
  for (lv in seq_along(levels_lbl)) {
    nH <- sum(hy[[var]]==codes[lv],na.rm=TRUE); nHe <- sum(he[[var]]==codes[lv],na.rm=TRUE)
    addrow(Variable=paste0("   ", levels_lbl[lv]),
           Hypo=sprintf("%d/%d (%.0f%%)", nH, nrow(hy), 100*nH/nrow(hy)),
           Hyper=sprintf("%d/%d (%.0f%%)", nHe, nrow(he), 100*nHe/nrow(he)), p=NA)
  }
}
cat_row("Age_cat", "Age category, n/N (%)", c("20-34 yr","35-49 yr","50-64 yr","65-79 yr","80+ yr"))
cat_row("BMI_cat", "Body mass index category, n/N (%)",
        c("<18.5 kg/m²","18.5-24.9 kg/m²","25.0-29.9 kg/m²","30.0-39.9 kg/m²","≥40 kg/m²","Unknown"))
comorbid_lbl <- c(HEART="Chronic heart disease, n/N (%)", LUNG="Chronic lung disease, n/N (%)",
                   KIDNEY="Chronic kidney disease, n/N (%)", HTN="Hypertension, n/N (%)", IMMUNO="Immunocompromised, n/N (%)")
for (v in names(comorbid_lbl)) {
  ft <- fisher.test(table(m$phenotype, m[[v]]))
  nH <- sum(hy[[v]]==1,na.rm=T); nHe <- sum(he[[v]]==1,na.rm=T)
  addrow(Variable=comorbid_lbl[[v]], Hypo=sprintf("%d/%d (%.0f%%)", nH, nrow(hy), 100*nH/nrow(hy)),
         Hyper=sprintf("%d/%d (%.0f%%)", nHe, nrow(he), 100*nHe/nrow(he)), p=ft$p.value)
}
cat_row("creat_0_cat", "Creatinine category, n/N (%)",
        c("0-0.79 mg/dL","0.80-1.19 mg/dL","1.20-1.79 mg/dL","1.80-2.99 mg/dL","≥3.00 mg/dL"))
cat_row("abs_neut_0_cat", "Neutrophil count category, n/N (%)",
        c("0-0.99 ×10⁹/L","1.00-3.99 ×10⁹/L","4.00-7.99 ×10⁹/L","8.00-11.99 ×10⁹/L","≥12.00 ×10⁹/L"))
ft <- fisher.test(table(m$phenotype, m$Acuity_max==1))
addrow(Variable="28-day mortality, n/N (%)",
       Hypo=sprintf("%d/%d (%.1f%%)", sum(hy$Acuity_max==1), nrow(hy), 100*mean(hy$Acuity_max==1)),
       Hyper=sprintf("%d/%d (%.1f%%)", sum(he$Acuity_max==1), nrow(he), 100*mean(he$Acuity_max==1)), p=ft$p.value)

tab1 <- rbindlist(rows)
tab1[, p_fmt := ifelse(is.na(p), "", ifelse(p<0.001, "<0.001", sprintf("%.3f", p)))]
fwrite(tab1[,.(Variable,Hypo,Hyper,p_fmt)], "Table1_baseline_characteristics.csv")
print(tab1[,.(Variable,Hypo,Hyper,p_fmt)])
