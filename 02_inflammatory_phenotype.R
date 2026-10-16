# Run this script after Mplus
library(data.table)
library(MplusAutomation)

target_proteins <- c(
  "IL6", "CXCL8", "TNFRSF1A", "PROC", "SERPINE1",
  "ICAM1", "AGER", "SFTPD", "VWF"
)
zvars <- paste0("z_", target_proteins)

# Compare models with 1 to 5 classes
patients <- as.data.frame(readRDS("primary_sample.rds"))
models <- readModels(target = getwd(), filefilter = "lca_k")
model_results <- data.frame()

for (k in 1:5) {
  model_name <- grep(
    paste0("lca_k", k, ".out$"),
    names(models),
    ignore.case = TRUE,
    value = TRUE
  )
  
  model <- models[[model_name]]
  model_summary <- model$summaries
  class_counts <- model$class_counts$mostLikely
  
  entropy <- NA_real_
  if (!is.null(model_summary$Entropy)) {
    entropy <- model_summary$Entropy
  }
  
  vlmr_p <- NA_real_
  if (!is.null(model_summary$T11_LMR_PValue)) {
    vlmr_p <- model_summary$T11_LMR_PValue
  }
  
  n_per_class <- NA_character_
  if (!is.null(class_counts)) {
    n_per_class <- paste(class_counts$count, collapse = " / ")
  }
  
  model_results <- rbind(model_results, data.frame(
    Classes = k,
    LL = model_summary$LL,
    BIC = model_summary$BIC,
    aBIC = model_summary$aBIC,
    Entropy = entropy,
    VLMR_p = vlmr_p,
    N_per_class = n_per_class
  ))
}

# Table S2
table_s2 <- data.frame(
  Classes = model_results$Classes,
  check.names = FALSE
)
table_s2[["Log-likelihood"]] <- sprintf("%.1f", model_results$LL)
table_s2$BIC <- sprintf("%.1f", model_results$BIC)
table_s2$aBIC <- sprintf("%.1f", model_results$aBIC)
table_s2$Entropy <- sprintf("%.3f", model_results$Entropy)
table_s2[["VLMR p"]] <- "NA"
table_s2[["N per class"]] <- model_results$N_per_class

for (i in 1:nrow(model_results)) {
  p_value <- model_results$VLMR_p[i]
  
  if (!is.na(p_value)) {
    if (p_value < 0.01) {
      table_s2[["VLMR p"]][i] <- sprintf("%.4f", p_value)
    } else {
      table_s2[["VLMR p"]][i] <- sprintf("%.3f", p_value)
    }
  }
}

fwrite(table_s2, "TableS2_LCA_enumeration.csv")

# Plot model comparison
png("FigureS1_LCA_enumeration_diagnostics.png",
    width = 2700, height = 1100, res = 300)
par(mfrow = c(1, 3), mar = c(5, 5, 3, 1), cex.axis = 0.8, cex.lab = 0.95)

plot(model_results$Classes, model_results$BIC,
     type = "b", pch = 19, xaxt = "n",
     xlab = "Number of classes (k)", ylab = "BIC", main = "A")
axis(1, 1:5)

plot(model_results$Classes, model_results$Entropy,
     type = "b", pch = 19, xaxt = "n", xlim = c(1, 5), ylim = c(0.75, 1),
     xlab = "Number of classes (k)", ylab = "Entropy", main = "B")
axis(1, 1:5)
abline(h = 0.8, lty = 2, col = "gray50")

bar_positions <- barplot(
  model_results$VLMR_p[2:5], names.arg = 2:5,
  ylim = c(0, 0.43), col = "gray65", border = NA,
  xlab = "Number of classes (k vs k-1)", ylab = "VLMR p value", main = "C"
)
abline(h = 0.05, lty = 2, col = "#D7301F")

bar_labels <- c(
  sprintf("%.4f", model_results$VLMR_p[2]),
  sprintf("%.4f", round(model_results$VLMR_p[3:5], 3))
)
text(bar_positions, model_results$VLMR_p[2:5],
     labels = bar_labels, pos = 3, cex = 0.75)
dev.off()

classes <- as.data.frame(fread("class_k2.dat", header = FALSE))
names(classes) <- c(
  zvars, "acuity0c", "creat0c", "neut0c",
  "CPROB1", "CPROB2", "C", "mplus_id"
)

class_columns <- c("mplus_id", "C", "CPROB1", "CPROB2")
data <- merge(patients, classes[, class_columns], by = "mplus_id")

class_numbers <- unique(data$C)
class_scores <- numeric(length(class_numbers))

for (i in seq_along(class_numbers)) {
  class_data <- data[data$C == class_numbers[i], zvars, drop = FALSE]
  protein_means <- colMeans(class_data, na.rm = TRUE)
  class_scores[i] <- mean(protein_means)
}

hyper_class <- class_numbers[which.max(class_scores)]
data$phenotype <- ifelse(
  data$C == hyper_class, "Hyperinflammatory", "Hypoinflammatory"
)
data$phenotype <- factor(
  data$phenotype,
  levels = c("Hypoinflammatory", "Hyperinflammatory")
)

if (hyper_class == 2) {
  data$CPROB_hyper <- data$CPROB2
} else {
  data$CPROB_hyper <- data$CPROB1
}

fwrite(data, "final_merged.csv")
cat("N =", nrow(data), "\n")
print(table(data$phenotype))

# protein levels per phenotype
hypo <- subset(data, phenotype == "Hypoinflammatory")
hyper <- subset(data, phenotype == "Hyperinflammatory")

protein_labels <- c(
  "IL-6", "IL-8", "sTNFR1", "Protein C", "PAI-1",
  "sICAM-1", "sRAGE", "SP-D", "vWF"
)
table_s1 <- data.frame()

for (i in seq_along(zvars)) {
  hypo_mean <- mean(hypo[[zvars[i]]], na.rm = TRUE)
  hyper_mean <- mean(hyper[[zvars[i]]], na.rm = TRUE)
  
  table_s1 <- rbind(table_s1, data.frame(
    Indicator = protein_labels[i],
    Hypoinflammatory = sprintf("%.2f", hypo_mean),
    Hyperinflammatory = sprintf("%.2f", hyper_mean),
    Difference = sprintf("%+.2f SD", hyper_mean - hypo_mean)
  ))
}

indicator_row <- function(label, variable, code) {
  hypo_values <- hypo[[variable]]
  hyper_values <- hyper[[variable]]
  
  if (is.na(code)) {
    n_hypo <- sum(is.na(hypo_values))
    n_hyper <- sum(is.na(hyper_values))
  } else {
    n_hypo <- sum(hypo_values == code, na.rm = TRUE)
    n_hyper <- sum(hyper_values == code, na.rm = TRUE)
  }
  
  total_hypo <- length(hypo_values)
  total_hyper <- length(hyper_values)
  percent_hypo <- 100 * n_hypo / total_hypo
  percent_hyper <- 100 * n_hyper / total_hyper
  
  return(data.frame(
    Indicator = label,
    Hypoinflammatory = sprintf("%d/%d (%.1f%%)", n_hypo, total_hypo, percent_hypo),
    Hyperinflammatory = sprintf("%d/%d (%.1f%%)", n_hyper, total_hyper, percent_hyper),
    Difference = sprintf("%+.1f pp", percent_hyper - percent_hypo)
  ))
}

table_s1 <- rbind(table_s1, indicator_row("Invasive ventilation", "Acuity_0_cat", 1))
table_s1 <- rbind(table_s1, indicator_row("Supplemental oxygen", "Acuity_0_cat", 2))

creatinine_labels <- c("0–0.79 mg/dL", "0.80–1.19", "1.20–1.79", "1.80–2.99", "≥3.00")
neutrophil_labels <- c("0–0.99 ×10⁹/L", "1.00–3.99", "4.00–7.99", "8.00–11.99", "≥12.00")

for (i in 1:5) {
  label <- paste("Creatinine", creatinine_labels[i])
  table_s1 <- rbind(table_s1, indicator_row(label, "creat0c", i))
}
table_s1 <- rbind(table_s1, indicator_row("Creatinine missing", "creat0c", NA))

for (i in 1:5) {
  label <- paste("Neutrophils", neutrophil_labels[i])
  table_s1 <- rbind(table_s1, indicator_row(label, "neut0c", i))
}
table_s1 <- rbind(table_s1, indicator_row("Neutrophils missing", "neut0c", NA))
fwrite(table_s1, "TableS1_class_indicators.csv")

# Compare ventilation and protein levels between phenotypes
ventilation_table <- table(data$phenotype, data$Acuity_0_cat == 1)
ventilation_test <- fisher.test(ventilation_table)
cat("Invasive ventilation by phenotype:\n")
print(ventilation_table)
cat("Fisher p =", ventilation_test$p.value, "\n")

for (protein in target_proteins) {
  formula <- as.formula(paste(protein, "~ phenotype"))
  protein_test <- wilcox.test(formula, data = data)
  cat(protein, "p =", protein_test$p.value, "\n")
}

# Table 1: baseline characteristics
data$Age_cat <- as.numeric(data$Age_cat)
data$BMI_cat <- as.numeric(data$BMI_cat)
data$creat_0_cat <- as.numeric(data$creat_0_cat)
data$abs_neut_0_cat <- as.numeric(data$abs_neut_0_cat)
data$Acuity_max <- as.numeric(data$Acuity_max)

hypo <- subset(data, phenotype == "Hypoinflammatory")
hyper <- subset(data, phenotype == "Hyperinflammatory")

category_rows <- function(variable, label, category_labels) {
  formula <- as.formula(paste(variable, "~ phenotype"))
  test <- wilcox.test(formula, data = data)
  rows <- data.frame(Variable = label, Hypo = "", Hyper = "", p = test$p.value)
  
  values <- data[[variable]]
  codes <- sort(unique(values[!is.na(values)]))
  stopifnot(length(codes) == length(category_labels))
  
  for (i in seq_along(category_labels)) {
    n_hypo <- sum(hypo[[variable]] == codes[i], na.rm = TRUE)
    n_hyper <- sum(hyper[[variable]] == codes[i], na.rm = TRUE)
    percent_hypo <- 100 * n_hypo / nrow(hypo)
    percent_hyper <- 100 * n_hyper / nrow(hyper)
    
    rows <- rbind(rows, data.frame(
      Variable = paste0("   ", category_labels[i]),
      Hypo = sprintf("%d/%d (%.0f%%)", n_hypo, nrow(hypo), percent_hypo),
      Hyper = sprintf("%d/%d (%.0f%%)", n_hyper, nrow(hyper), percent_hyper),
      p = NA_real_
    ))
  }
  
  return(rows)
}

table1 <- category_rows(
  "Age_cat", "Age category, n/N (%)",
  c("20-34 yr", "35-49 yr", "50-64 yr", "65-79 yr", "80+ yr")
)

bmi_rows <- category_rows(
  "BMI_cat", "Body mass index category, n/N (%)",
  c("<18.5 kg/m²", "18.5-24.9 kg/m²", "25.0-29.9 kg/m²",
    "30.0-39.9 kg/m²", "≥40 kg/m²", "Unknown")
)
table1 <- rbind(table1, bmi_rows)

comorbidity_labels <- c(
  HEART = "Chronic heart disease, n/N (%)",
  LUNG = "Chronic lung disease, n/N (%)",
  KIDNEY = "Chronic kidney disease, n/N (%)",
  HTN = "Hypertension, n/N (%)",
  IMMUNO = "Immunocompromised, n/N (%)"
)

for (variable in names(comorbidity_labels)) {
  test <- fisher.test(table(data$phenotype, data[[variable]]))
  n_hypo <- sum(hypo[[variable]] == 1, na.rm = TRUE)
  n_hyper <- sum(hyper[[variable]] == 1, na.rm = TRUE)
  percent_hypo <- 100 * n_hypo / nrow(hypo)
  percent_hyper <- 100 * n_hyper / nrow(hyper)
  
  table1 <- rbind(table1, data.frame(
    Variable = comorbidity_labels[[variable]],
    Hypo = sprintf("%d/%d (%.0f%%)", n_hypo, nrow(hypo), percent_hypo),
    Hyper = sprintf("%d/%d (%.0f%%)", n_hyper, nrow(hyper), percent_hyper),
    p = test$p.value
  ))
}

covid_test <- fisher.test(table(data$phenotype, data$COVID))
table1 <- rbind(table1, data.frame(
  Variable = "Confirmed COVID-19, n/N (%)",
  Hypo = sprintf("%d/%d (%.1f%%)",
                 sum(hypo$COVID == 1), nrow(hypo), 100 * mean(hypo$COVID == 1)),
  Hyper = sprintf("%d/%d (%.1f%%)",
                  sum(hyper$COVID == 1), nrow(hyper), 100 * mean(hyper$COVID == 1)),
  p = covid_test$p.value
))

creatinine_rows <- category_rows(
  "creat_0_cat", "Creatinine category, n/N (%)",
  c("0-0.79 mg/dL", "0.80-1.19 mg/dL", "1.20-1.79 mg/dL",
    "1.80-2.99 mg/dL", "≥3.00 mg/dL")
)
table1 <- rbind(table1, creatinine_rows)

neutrophil_rows <- category_rows(
  "abs_neut_0_cat", "Neutrophil count category, n/N (%)",
  c("0-0.99 ×10⁹/L", "1.00-3.99 ×10⁹/L", "4.00-7.99 ×10⁹/L",
    "8.00-11.99 ×10⁹/L", "≥12.00 ×10⁹/L")
)
table1 <- rbind(table1, neutrophil_rows)

mortality_test <- fisher.test(table(data$phenotype, data$Acuity_max == 1))
table1 <- rbind(table1, data.frame(
  Variable = "28-day mortality, n/N (%)",
  Hypo = sprintf("%d/%d (%.1f%%)",
                 sum(hypo$Acuity_max == 1), nrow(hypo), 100 * mean(hypo$Acuity_max == 1)),
  Hyper = sprintf("%d/%d (%.1f%%)",
                  sum(hyper$Acuity_max == 1), nrow(hyper), 100 * mean(hyper$Acuity_max == 1)),
  p = mortality_test$p.value
))

table1$p_fmt <- ""
for (i in 1:nrow(table1)) {
  p_value <- table1$p[i]
  
  if (!is.na(p_value)) {
    if (p_value < 0.001) {
      table1$p_fmt[i] <- "<0.001"
    } else {
      table1$p_fmt[i] <- sprintf("%.3f", p_value)
    }
  }
}

table1_export <- table1[, c("Variable", "Hypo", "Hyper", "p_fmt")]
fwrite(table1_export, "Table1_baseline_characteristics.csv")
print(table1_export)

