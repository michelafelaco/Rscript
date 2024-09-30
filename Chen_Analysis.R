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

#File richiesti
CHEN_linear_raw_data <- readRDS("C:/Users/miche/Desktop/CNR lab/Trascriptomic_R/Chen_database/Database/CHEN_linear_raw_data.rds")
write.table(CHEN_linear_raw_data, file = "C:\\Users\\miche\\Desktop\\CNR lab\\Trascriptomic_R\\Chen_database\\Chen_matrix.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
CHEN_SLITROBO <- CHEN_linear_raw_data [c("SLIT1", "SLIT2", "SLIT3", "ROBO1", "ROBO2", "ROBO3", "ROBO4"),]

CHEN_meta <- readRDS("C:/Users/miche/Desktop/CNR lab/Trascriptomic_R/Chen_database/Database/CHEN_meta.rds")
write.table(CHEN_meta, file = "C:\\Users\\miche\\Desktop\\CNR lab\\Trascriptomic_R\\Chen_database\\Chen_metadata.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

Cibersort_CHEN <- read_xlsx("C:\\Users\\miche\\Desktop\\CNR lab\\Trascriptomic_R\\Chen_database\\Cibersort\\CHENmicro_Cibersort.xlsx")
Cibersort_CHEN <- as.data.frame(Cibersort_CHEN)
rownames(Cibersort_CHEN) <- Cibersort_CHEN[, 1]
Cibersort_CHEN <- Cibersort_CHEN[, -1]
Cibersort_CHEN_filtered <- Cibersort_CHEN [Cibersort_CHEN$`P-value` <= 0.05, 1:22]
CHEN_SLITROBO_filtered <- CHEN_SLITROBO [, rownames(Cibersort_CHEN_filtered)]

#2.Gene Correlation
cor_matrix <- cor(t(CHEN_SLITROBO_filtered))
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
p_matrix <- cor_test(t(CHEN_SLITROBO_filtered))
plot <- corrplot(cor_matrix, method="color", type="upper", order = "original" , 
                 p.mat = p_matrix, sig.level = c(0.0001, 0.001, 0.01, 0.05), 
                 insig = "label_sig", pch.col = "black", pch.cex = 2, 
                 tl.cex = 1, tl.col = "black", 
                 na.label = " ")

#KM maxstat
#filtraggio matrici
rownames(CHEN_meta) <- CHEN_meta [,2]
CHEN_meta_filtered <- CHEN_meta [rownames(Cibersort_CHEN_filtered),]

dataframe <- data.frame(
  time = as.numeric(CHEN_meta_filtered [,"time"]),
  status = as.logical(CHEN_meta_filtered [,"status"]),
  value = as.numeric(CHEN_SLITROBO_filtered ["ROBO4",])
)
surv_object <- Surv(time = dataframe$time, event = dataframe$status)
maxstat_result <- maxstat.test(surv_object ~ dataframe$value,
                               data = dataframe, 
                               smethod = "LogRank",
                               pmethod = "exactGauss", 
                               minprop = 0.1,
                               maxprop = 0.9) 
optimal_cutpoint <- maxstat_result$estimate
print(paste("Optimal cutpoint:", optimal_cutpoint))
dataframe$group <- ifelse(dataframe$value > optimal_cutpoint, "High", "Low")
fit <- survfit(surv_object ~ group, data = dataframe)
logrank_test <- survdiff(surv_object ~ group, data = dataframe)
p_value <- 1 - pchisq(logrank_test$chisq, df = 1)
print(paste("Log-rank p-value:", p_value))

cox_model <- coxph(surv_object ~ group, data = dataframe)
hr <- exp(coef(cox_model))
hr_ci <- exp(confint(cox_model))
hr_label <- paste("HR:", round(hr, 2), 
                  "95% CI:", paste(round(hr_ci, 2), collapse = "-"))

plot_title <- paste("Survival Analysis for ROBO4\n", 
                    "Log-rank p-value:", round(p_value, 4), 
                    "-", hr_label)

ggsurvplot(fit, data = dataframe, 
           pval = paste("p =", round(p_value, 4)), 
           conf.int = TRUE,
           risk.table = TRUE, 
           xlab = "Time (months)", 
           ylab = "Survival Probability",
           ggtheme = theme_minimal(),
           title = plot_title)

#Cibersort boxplot
df <- as.data.frame(Cibersort_CHEN_filtered)
df_long <- melt(df)
ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  theme_minimal() +  # Usa un tema minimale
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),  
        plot.title = element_text(size = 16),    
        legend.position = "none") +  
  labs(title = "Boxplot Cibersort CHEN filtrato (n=52)", x = "Variabili", y = "Valori") +
  scale_y_continuous(limits = c(NA, 0.5)) +  
  scale_fill_manual(values = rainbow(22))

#4. Cibersort correlations
cor_matrix <- cor(Cibersort_CHEN_filtered)
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
p_matrix <- cor_test(Cibersort_CHEN_filtered)
png(filename = "Corr_Cibersort_CHEN.jpeg", width = 1600, height = 1200, res = 300)
corrplot(cor_matrix, type = "upper", order = "original",
         p.mat = p_matrix, sig.level = c(0.0001, 0.001, 0.01, 0.05), 
         insig = "label_sig", pch.col = "black", pch.cex = 0.6, 
         tl.cex = 0.5, tl.col = "black", 
         na.label = " ")
dev.off()

#heatmap
gene_names <- c("SLIT1", "SLIT2", "SLIT3", "ROBO1", "ROBO2", "ROBO3", "ROBO4")
correlations_list <- list()
pvalues_list <- list()

for (i in 1:length(gene_names)) {
  correlations <- apply(t(Cibersort_CHEN_filtered), 1, function(row) cor(row, CHEN_SLITROBO_filtered[i,]))
  correlations_list[[gene_names[i]]] <- correlations
  
  pvalues <- apply(t(Cibersort_CHEN_filtered), 1, function(row) {
    test <- cor.test(row, CHEN_SLITROBO_filtered[i,])
    return(test$p.value)
  })
  pvalues_list[[gene_names[i]]] <- pvalues
}
dataframe_correlations <- as.data.frame(correlations_list)
dataframe_pvalues <- as.data.frame(pvalues_list)
pvalue_to_stars <- function(p) {
  if (p < 0.0001) {
    return("****")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*") 
  } else {
    return("")
  }
}
stars_matrix <- apply(dataframe_pvalues, 2, function(col) sapply(col, pvalue_to_stars))
heatmap_result <- pheatmap(dataframe_correlations, 
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE,
                           main = "Correlations Cibersort-Genes",
                           display_numbers = stars_matrix,  
                           number_color = "black",
                           fontsize_number = 12)          