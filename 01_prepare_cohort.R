# Inputs from the public archive in data/
# (Zenodo, https://doi.org/10.5281/zenodo.7076472, CC-BY 4.0)
library(data.table)

npx <- fread("data/MGH_COVID_OLINK_NPX.txt", sep=";")
clin <- fread("data/MGH_COVID_Clinical_Info.txt", sep=";")

target_proteins <- c("IL6","CXCL8","TNFRSF1A","PROC","SERPINE1","ICAM1","AGER","SFTPD","VWF")
zvars <- paste0("z_", target_proteins)

npx_d0 <- npx[Timepoint=="D0" & Assay %in% target_proteins]
npx_d0[, npx_clean := ifelse(QC_Warning=="PASS" & Assay_Warning=="PASS", NPX, NA_real_)]
wide <- dcast(npx_d0, subject_id ~ Assay, value.var="npx_clean", fun.aggregate=mean, na.rm=TRUE)

clin[, subject_id := as.character(subject_id)]
wide[, subject_id := as.character(subject_id)]
primary <- merge(wide, clin, by="subject_id", all.x=TRUE)

primary[, Acuity_0 := as.numeric(Acuity_0)]
primary[, Acuity_max := as.numeric(Acuity_max)]
primary[, Acuity_0_cat := fcase(
  Acuity_0 %in% c(1,2), 1L,
  Acuity_0 == 3, 2L,
  Acuity_0 == 4, 3L,
  Acuity_0 == 5, 4L
)]
primary[, creat0c := as.numeric(creat_0_cat)]
primary[, neut0c := as.numeric(abs_neut_0_cat)]

p4 <- primary[Acuity_0 <= 3]
p4[, mplus_id := .I]
for (i in seq_along(target_proteins)) {
  p4[, (zvars[i]) := as.numeric(scale(get(target_proteins[i])))]
}

mplus_vars <- c("mplus_id", zvars, "Acuity_0_cat", "creat0c", "neut0c", "Acuity_max")
dat <- copy(p4[, ..mplus_vars])
for (v in c(zvars,"Acuity_0_cat","creat0c","neut0c","Acuity_max")) dat[is.na(get(v)), (v) := -999]
fwrite(dat, "mgh_lca_input.dat", sep=" ", col.names=FALSE, na="-999")

saveRDS(p4, "primary_sample.rds")
cat("N =", nrow(p4), "\n")
