#Packages
library(readxl)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(corrplot)
library(Hmisc)
library(survival)
library(ggpubr)
library(survminer)
library(maxstat)

#importazione Moffit dataset
dataset <- read.table("C:\\Users\\miche\\Desktop\\CNR lab\\Trascriptomic_R\\Moffit_database\\GSE71729_series_matrix.txt\\GSE71729_series_matrix.txt", header = TRUE, skip = 26, sep = "\t", fill = TRUE)
dataset <- dataset[-nrow(dataset), ]
numeric_matrix_Moffit <- apply(as.matrix(dataset[44:nrow(dataset), 2:ncol(dataset)]), 2, as.numeric)
rownames(numeric_matrix_Moffit) <- dataset[44:nrow(dataset), 1]

#creazione matrici SLIT/ROBO
SLIT_ROBO <- numeric_matrix_Moffit [c("SLIT1", "SLIT2", "SLIT3", "ROBO1", "ROBO2", "ROBO3", "ROBO4"),]
SLITROBO_normal_pancreas <- SLIT_ROBO [, grep("Normal.Pancreas", colnames(SLIT_ROBO), value = TRUE)]
SLITROBO_primary_pancreas <- SLIT_ROBO [, grep("Primary.Pancreas", colnames(SLIT_ROBO), value = TRUE)]
SLITROBO_metastasis_liver <- SLIT_ROBO [, grep("Met.Liver", colnames(SLIT_ROBO), value = TRUE)]
SLITROBO_normal_liver <- SLIT_ROBO [, grep("Normal.Liver", colnames(SLIT_ROBO), value = TRUE)]

#1.VIOLIN PLOT GENE EXPRESSION
valori_normal <- c(SLITROBO_normal_pancreas["ROBO4",])
valori_primary <- c(SLITROBO_primary_pancreas["ROBO4",])
valori_met_liver <- c(SLITROBO_metastasis_liver["ROBO4",])
valori_liver <- c(SLITROBO_normal_liver["ROBO4",])
valori <- c(valori_normal, valori_primary, valori_met_liver, valori_liver)
gruppi <- factor(c(rep("1.Normal Pancreas", length(valori_normal)),
                   rep("2.Primary Pancreas", length(valori_primary)),
                   rep("3.Liver Metastasis", length(valori_met_liver)),
                   rep("4.Normal Liver", length(valori_liver))))
dati <- data.frame(Gruppo = gruppi, Valore = valori)
dati$X <- factor(dati$Gruppo)
p <- ggplot(dati, aes(x = Gruppo, y = Valore)) +
  geom_violin(aes(fill = Gruppo), alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2) +
  stat_summary(fun = median, geom = "crossbar", color = "black", size = 0.3) +
  labs(title = "Violin Plot ROBO4",
       x = "",
       y = "log expressione value") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))

#Kruskal-Wallis test
num_gruppi <- length(unique(gruppi))
max_valore <- max(valori)  
p <- p + stat_compare_means(method = "kruskal.test", label.y = max_valore + 1.5)

#Pairwise test di Wilcoxon
comparisons <- list(c("1.Normal Pancreas", "2.Primary Pancreas"),
                    c("2.Primary Pancreas", "3.Liver Metastasis"),
                    c("3.Liver Metastasis", "4.Normal Liver"))

p <- p + stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif")
print(p)
ggsave("ROBO4_expression.jpg", plot = p, width = 10, height = 7, dpi = 300)

#1b.Violin Plot genes for each group
geni <- c("SLIT1", "SLIT2", "SLIT3", "ROBO1", "ROBO2", "ROBO3", "ROBO4")
valori <- unlist(lapply(geni, function(gene) c(SLITROBO_normal_pancreas[gene, ])))
gruppi <- factor(rep(geni, sapply(geni, function(gene) length(SLITROBO_normal_pancreas[gene, ]))),
                 levels = geni) 
dati <- data.frame(Gruppo = gruppi, Valore = valori)
p <- ggplot(dati, aes(x = Gruppo, y = Valore)) +
  geom_violin(aes(fill = Gruppo), alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2) +
  stat_summary(fun = median, geom = "crossbar", color = "black", size = 0.3) +
  labs(title = "Normal Pancreas",
       x = "",
       y = "log expression value") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))
ggsave("NormalPancreas_SLITROBO.jpg", plot = p, width = 10, height = 7, dpi = 300)

#2.GENE CORRELATION
cor_matrix <- cor(t(SLITROBO_metastasis_liver))
cor_test <- function(mat) {
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      test <- cor.test(mat[, i], mat[, j])
      p.mat[i, j] <- p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- NA  
  return(p.mat)
}
p_matrix <- cor_test(t(SLITROBO_metastasis_liver))
plot <- corrplot(cor_matrix, method="color", type="upper", order = "original" , 
                 p.mat = p_matrix, sig.level = c(0.0001, 0.001, 0.01, 0.05), 
                 insig = "label_sig", pch.col = "black", pch.cex = 2, 
                 tl.cex = 1, tl.col = "black", 
                 na.label = " ")

#3.KAPLAN-MEIER
GEO_ID <- dataset[1,2:ncol(dataset)]
ROBO4 <-  numeric_matrix_Moffit ["ROBO4",]
names(ROBO4) <- GEO_ID
dataset_survival <- as.matrix(read_xlsx("C:\\Users\\miche\\Desktop\\CNR lab\\Trascriptomic_R\\Moffit_database\\SLIT_ROBO\\Santiago Poster\\Tumor_survival.xlsx"))
row.names(dataset_survival) <- as.vector(dataset_survival[,1])
dataset_survival <- dataset_survival [,2:3]
ROBO4 <- t(ROBO4 [rownames(dataset_survival)])

data <- data.frame(
  time = as.numeric(dataset_survival [,"time"]),
  status = as.logical(dataset_survival [,"status"]),
  value = as.numeric(ROBO4)
)
maxstat_test <- maxstat.test(Surv(time, status) ~ value, data = data, 
                             smethod="LogRank", pmethod="exactGauss")
print(maxstat_test)
cutpoint <- maxstat_test$estimate
data$group <- ifelse(data$value <= cutpoint, "Low", "High")
fit <- survfit(Surv(time, status) ~ group, data = data)
cox_model <- coxph(Surv(time, status) ~ group, data = data)
cox_summary <- summary(cox_model)
hr <- cox_summary$coefficients[1, "exp(coef)"]
hr_confint <- confint(cox_model)
hr_lower <- hr_confint[1, 1]
hr_upper <- hr_confint[1, 2]
cat("Hazard Ratio (HR):", hr, "\n")
cat("95% Confidence Interval for HR:", hr_lower, "-", hr_upper, "\n")

km_plot <- ggsurvplot(
  fit, 
  data = data, 
  pval = TRUE, 
  conf.int = TRUE,
  risk.table = TRUE,
  legend.title = "Group",
  legend.labs = levels(data$group),
  xlab = "Time (days)",
  ylab = "Survival Probability",
  ggtheme = theme_bw(),
  title = "Kaplan-Meier ROBO4"
)
hr_text <- paste("HR:", round(hr, 2))
km_plot$plot <- km_plot$plot + 
  annotate("text", x = max(data$time) * 0.08, y = 0.3, label = hr_text, size = 5, color = "black")
print(km_plot)
ggsave("kaplan_meier_SLIT1.png", plot = km_plot$plot, bg = "white")

#4 CIBERSORT
Cibersort_all <- read.table("C:\\Users\\miche\\Desktop\\CNR lab\\Cibersort\\Moffit\\CIBERSORTx_Job11_Results.txt", 
                            header = TRUE, sep = "\t", fill = TRUE)
Cibersort <- Cibersort_all [, 2:23]
rownames(Cibersort) <- Cibersort_all [1:nrow(Cibersort_all), 1]
rows_primary_pancreas <- grep("Primary.Pancreas", rownames(Cibersort), value = TRUE)
Cibersort_Primary.Pancreas <-  t(Cibersort [rows_primary_pancreas,])
rows_normal_pancreas <- grep("Normal.Pancreas", rownames(Cibersort), value = TRUE)
Cibersort_Normal.Pancreas <-  t(Cibersort [rows_normal_pancreas,])

#4 CIBERSORT CORRELATION
cor_matrix <- cor(t(Cibersort_Primary.Pancreas))
diag(cor_matrix) <- NA
cor_test <- function(mat) {
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      test <- cor.test(mat[, i], mat[, j])
      p.mat[i, j] <- p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- NA  
  return(p.mat)
}
p_matrix <- cor_test(t(Cibersort_Primary.Pancreas))
png(filename = "Corr_Cibersort_PrimaryPancreas.jpeg", width = 1600, height = 1200, res = 300)
corrplot(cor_matrix, type = "upper", order = "original",
         p.mat = p_matrix, sig.level = c(0.0001, 0.001, 0.01, 0.05), 
         insig = "label_sig", pch.col = "black", pch.cex = 0.6, 
         tl.cex = 0.5, tl.col = "black", 
         na.label = " ")
dev.off()

cor_matrix <- cor(t(Cibersort_Primary.Pancreas [c("Macrophages.M0", "Monocytes", "T.cells.CD4.memory.resting", "NK.cells.activated", "B.cells.memory", "T.cells.regulatory..Tregs."),]))
diag(cor_matrix) <- NA
cor_test <- function(mat) {
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      test <- cor.test(mat[, i], mat[, j])
      p.mat[i, j] <- p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- NA  
  return(p.mat)
}
p_matrix <- cor_test(t(Cibersort_Primary.Pancreas [c("Macrophages.M0", "Monocytes", "T.cells.CD4.memory.resting", "NK.cells.activated", "B.cells.memory", "T.cells.regulatory..Tregs."),]))
png(filename = "Corr.sign.PrimaryPancreas.jpeg", width = 1600, height = 1200, res = 200)
corrplot(cor_matrix, type = "upper", order = "original",
         p.mat = p_matrix, sig.level = c(0.0001, 0.001, 0.01, 0.05), 
         insig = "label_sig", pch.col = "black", pch.cex = 0.8, 
         tl.cex = 0.8, tl.col = "black", 
         na.label = " ")
dev.off()


#4 Heatmap Correlation
correlations_SLIT1 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[1,]))
correlations_SLIT2 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[2,]))
correlations_SLIT3 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[3,]))
correlations_ROBO1 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[4,]))
correlations_ROBO2 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[5,]))
correlations_ROBO3 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[6,]))
correlations_ROBO4 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[7,]))

extract_pvalue <- function(row, target) {
  result <- cor.test(row, target)
  return(result$p.value)
}
pvalue_SLIT1 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[1,]))
pvalue_SLIT2 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[2,]))
pvalue_SLIT3 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[3,]))
pvalue_ROBO1 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[4,]))
pvalue_ROBO2 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[5,]))
pvalue_ROBO3 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[6,]))
pvalue_ROBO4 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[7,]))

dataframe_correlations_normal <- data.frame(SLIT1 = correlations_SLIT1, SLIT2 = correlations_SLIT2, SLIT3 = correlations_SLIT3, ROBO1= correlations_ROBO1, ROBO2 = correlations_ROBO2, ROBO3 = correlations_ROBO3, ROBO4 = correlations_ROBO4)

heatmap_result <- pheatmap(dataframe_correlations_normal, 
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE,
                           main = "Heatmap Normal Pancreas")
row_order <- heatmap_result$tree_row$order
col_order <- heatmap_result$tree_col$order
tiff("Heatmap_normal_pancreas.tiff", width = 1200, height = 800)
pheatmap(dataframe_correlations_normal, main = "Heatmap Normal Pancreas", fontsize = 20)
dev.off()

correlationsP_SLIT1 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[1,]))
correlationsP_SLIT2 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[2,]))
correlationsP_SLIT3 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[3,]))
correlationsP_ROBO1 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[4,]))
correlationsP_ROBO2 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[5,]))
correlationsP_ROBO3 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[6,]))
correlationsP_ROBO4 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[7,]))
dataframe_correlations_primary <- data.frame(SLIT1 = correlationsP_SLIT1, SLIT2 = correlationsP_SLIT2, SLIT3 = correlationsP_SLIT3, ROBO1= correlationsP_ROBO1, ROBO2 = correlationsP_ROBO2, ROBO3 = correlationsP_ROBO3, ROBO4 = correlationsP_ROBO4)
dataframe_correlations_primary_ordered <- dataframe_correlations_primary [row_order,col_order]
tiff("Heatmap_primary_pancreas.tiff", width = 1200, height = 800)
pheatmap(dataframe_correlations_primary_ordered, cluster_rows = FALSE, cluster_cols = FALSE, main = "Heatmap Primary Pancreas", fontsize = 20)
dev.off()

pvalueP_SLIT1 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[1,]))
pvalueP_SLIT2 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[2,]))
pvalueP_SLIT3 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[3,]))
pvalueP_ROBO1 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[4,]))
pvalueP_ROBO2 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[5,]))
pvalueP_ROBO3 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[6,]))
pvalueP_ROBO4 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[7,]))

get_pheatmap_colors <- function(num_colors) {
  default_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(num_colors)
  return(default_colors)
}
default_pheatmap_colors <- get_pheatmap_colors(100)

#4 Heatmap Correlation pesata
correlations_SLIT1 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[1,]))
correlations_SLIT2 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[2,]))
correlations_SLIT3 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[3,]))
correlations_ROBO1 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[4,]))
correlations_ROBO2 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[5,]))
correlations_ROBO3 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[6,]))
correlations_ROBO4 <- apply(Cibersort_Normal.Pancreas, 1, function(row) cor(row, SLITROBO_normal_pancreas[7,]))

extract_pvalue <- function(row, target) {
  result <- cor.test(row, target)
  return(result$p.value)
}
pvalue_SLIT1 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[1,]))
pvalue_SLIT2 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[2,]))
pvalue_SLIT3 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[3,]))
pvalue_ROBO1 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[4,]))
pvalue_ROBO2 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[5,]))
pvalue_ROBO3 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[6,]))
pvalue_ROBO4 <- apply(Cibersort_Normal.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_normal_pancreas[7,]))

dataframe_correlations_normal <- data.frame(SLIT1 = correlations_SLIT1, SLIT2 = correlations_SLIT2, SLIT3 = correlations_SLIT3, ROBO1= correlations_ROBO1, ROBO2 = correlations_ROBO2, ROBO3 = correlations_ROBO3, ROBO4 = correlations_ROBO4)
dataframe_pvalue_normal <- data.frame(SLIT1 = pvalue_SLIT1, SLIT2 = pvalue_SLIT2, SLIT3 = pvalue_SLIT3, ROBO1= pvalue_ROBO1, ROBO2 = pvalue_ROBO2, ROBO3 = pvalue_ROBO3, ROBO4 = pvalue_ROBO4)

matrix_correlations <- as.matrix(dataframe_correlations_normal)
matrix_pvalues <- as.matrix(dataframe_pvalue_normal)
log_pvalues <- -log10(matrix_pvalues)
product_matrix <- matrix_correlations * log_pvalues
dataframe_product <- as.data.frame(product_matrix)

heatmap_result <- pheatmap(dataframe_product, cluster_rows = TRUE, cluster_cols = TRUE, main = "Heatmap pesata Normal Pancreas")
row_order <- heatmap_result$tree_row$order
col_order <- heatmap_result$tree_col$order

correlationsP_SLIT1 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[1,]))
correlationsP_SLIT2 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[2,]))
correlationsP_SLIT3 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[3,]))
correlationsP_ROBO1 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[4,]))
correlationsP_ROBO2 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[5,]))
correlationsP_ROBO3 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[6,]))
correlationsP_ROBO4 <- apply(Cibersort_Primary.Pancreas, 1, function(row) cor(row, SLITROBO_primary_pancreas[7,]))

extract_pvalue <- function(row, target) {
  result <- cor.test(row, target)
  return(result$p.value)
}
pvalueP_SLIT1 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[1,]))
pvalueP_SLIT2 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[2,]))
pvalueP_SLIT3 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[3,]))
pvalueP_ROBO1 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[4,]))
pvalueP_ROBO2 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[5,]))
pvalueP_ROBO3 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[6,]))
pvalueP_ROBO4 <- apply(Cibersort_Primary.Pancreas, 1, function(row) extract_pvalue(row, SLITROBO_primary_pancreas[7,]))

dataframe_correlations_primary <- data.frame(SLIT1 = correlationsP_SLIT1, SLIT2 = correlationsP_SLIT2, SLIT3 = correlationsP_SLIT3, ROBO1= correlationsP_ROBO1, ROBO2 = correlationsP_ROBO2, ROBO3 = correlationsP_ROBO3, ROBO4 = correlationsP_ROBO4)
dataframe_pvalue_primary <- data.frame(SLIT1 = pvalueP_SLIT1, SLIT2 = pvalueP_SLIT2, SLIT3 = pvalueP_SLIT3, ROBO1= pvalueP_ROBO1, ROBO2 = pvalueP_ROBO2, ROBO3 = pvalueP_ROBO3, ROBO4 = pvalueP_ROBO4)

matrix_correlations <- as.matrix(dataframe_correlations_primary)
matrix_pvalues <- as.matrix(dataframe_pvalue_primary)
log_pvalues <- -log10(matrix_pvalues)
product_matrix <- matrix_correlations * log_pvalues
product_matrix_ordered <- product_matrix [row_order, col_order]
dataframe_product <- as.data.frame(product_matrix_ordered)

pheatmap(dataframe_product, cluster_rows = FALSE, cluster_cols = FALSE, main = "Heatmap pesata Primary Pancreas")
tiff("Heatmappesata_primary_pancreas.tiff", width = 1200, height = 800)
pheatmap(dataframe_product, main = "Heatmap pesata primary Pancreas", cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 25)
dev.off()

#4 CIBERSORT-genes correlations
correlations <- cbind(correlations_ROBO1, correlations_ROBO4, correlations_SLIT2, correlations_SLIT3, correlations_ROBO3, correlations_SLIT1, correlations_ROBO2)
correlations_ordered <- (correlations)[row_order,]
pvalues <- cbind(pvalue_ROBO1, pvalue_ROBO4, pvalue_SLIT2, pvalue_SLIT3, pvalue_ROBO3, pvalue_SLIT1, pvalue_ROBO2)
pvalues_ordered <- pvalues [row_order,]

correlations_sign <- correlations_ordered [c("Macrophages.M0", "Monocytes", "T.cells.CD4.memory.resting", "NK.cells.activated", "B.cells.memory", "T.cells.regulatory..Tregs."),]
pvalues_sign <- pvalues_ordered [c("Macrophages.M0", "Monocytes", "T.cells.CD4.memory.resting", "NK.cells.activated", "B.cells.memory", "T.cells.regulatory..Tregs."),]
tiff("Correlations_sign.tiff", width = 1600, height = 1200, res = 200)
corrplot(correlations_sign, order = "original",
         p.mat = pvalues_sign, sig.level = c(0.0001, 0.001, 0.01, 0.05), 
         insig = "label_sig", pch.col = "black", pch.cex = 1, 
         tl.cex = 1, tl.col = "black", tl.srt = 45, 
         na.label = " ",
         cl.pos = "r",
         xaxt = "n",
         col = default_pheatmap_colors)
dev.off()

tiff("Correlations_genes_cibersort_normal_pancreas.tiff", width = 1600, height = 1200, res = 200)
corrplot(correlations_ordered, order = "original",
         p.mat = pvalues_ordered, sig.level = c(0.0001, 0.001, 0.01, 0.05), 
         insig = "label_sig", pch.col = "black", pch.cex = 1, 
         tl.cex = 1, tl.col = "black", tl.srt = 45, 
         na.label = " ",
         cl.pos = "r",
         xaxt = "n",
         col = default_pheatmap_colors)
dev.off()

correlationsP <- cbind(correlationsP_ROBO1, correlationsP_ROBO4, correlationsP_SLIT2, correlationsP_SLIT3, correlationsP_ROBO3, correlationsP_SLIT1, correlationsP_ROBO2)
correlations_orderedP <- (correlationsP)[row_order,]
pvaluesP <- cbind(pvalueP_ROBO1, pvalueP_ROBO4, pvalueP_SLIT2, pvalueP_SLIT3, pvalueP_ROBO3, pvalueP_SLIT1, pvalueP_ROBO2)
pvalues_orderedP <- pvaluesP [row_order,]

tiff("Correlations_genes_cibersort_primary_pancreas.tiff", width = 1600, height = 1200, res = 200)
corrplot(correlations_orderedP, order = "original",
         p.mat = pvalues_orderedP, sig.level = c(0.0001, 0.001, 0.01, 0.05), 
         insig = "label_sig", pch.col = "black", pch.cex = 1, 
         tl.cex = 1, tl.col = "black", tl.srt = 45, 
         na.label = " ",
         cl.pos = "r",
         xaxt = "n",
         col = default_pheatmap_colors)
dev.off()

correlations_signP <- correlations_orderedP [c("Macrophages.M0", "Monocytes", "T.cells.CD4.memory.resting", "NK.cells.activated", "B.cells.memory", "T.cells.regulatory..Tregs."),]
pvalues_signP <- pvalues_orderedP [c("Macrophages.M0", "Monocytes", "T.cells.CD4.memory.resting", "NK.cells.activated", "B.cells.memory", "T.cells.regulatory..Tregs."),]
tiff("Correlations_signP.tiff", width = 1600, height = 1200, res = 200)
corrplot(correlations_signP, order = "original",
         p.mat = pvalues_signP, sig.level = c(0.0001, 0.001, 0.01, 0.05), 
         insig = "label_sig", pch.col = "black", pch.cex = 1, 
         tl.cex = 1, tl.col = "black", tl.srt = 45, 
         na.label = " ",
         cl.pos = "r",
         xaxt = "n",
         col = default_pheatmap_colors)
dev.off()