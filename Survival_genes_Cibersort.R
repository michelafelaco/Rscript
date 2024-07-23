#packages
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
PDAC <- numeric_matrix_Moffit [, grep("Primary.Pancreas", colnames(numeric_matrix_Moffit), value = TRUE)]
SLIT2 <-  PDAC ["SLIT2",]
SLIT3 <- PDAC ["SLIT3",]
ROBO1 <- PDAC ["ROBO1",]

#importazione survival
GEO_ID <- dataset[1,2:ncol(dataset)]
GEO_ID <- GEO_ID [,grep("Primary.Pancreas", colnames(GEO_ID), value = TRUE)]
dataset_survival <- as.matrix(read_xlsx("C:\\Users\\miche\\Desktop\\CNR lab\\Trascriptomic_R\\Moffit_database\\SLIT_ROBO\\Santiago Poster\\Tumor_survival.xlsx"))
row.names(dataset_survival) <- as.vector(dataset_survival[,1])
dataset_survival <- dataset_survival [,2:3]
colnames(PDAC) <- GEO_ID
PDAC <- PDAC [, rownames(dataset_survival)]
time <- as.numeric(dataset_survival[, 1])
status <- as.logical(dataset_survival[, 2])
names(SLIT2) <- GEO_ID
SLIT2 <- t(SLIT2 [rownames(dataset_survival)])
names(SLIT3) <- GEO_ID
SLIT3 <- (SLIT3 [rownames(dataset_survival)])
names(ROBO1) <- GEO_ID
ROBO1 <- (ROBO1 [rownames(dataset_survival)])

#importazione CIBERSORT
Cibersort_all <- read.table("C:\\Users\\miche\\Desktop\\CNR lab\\Cibersort\\Moffit\\CIBERSORTx_Job11_Results.txt", 
                            header = TRUE, sep = "\t", fill = TRUE)
Cibersort <- Cibersort_all [, 2:23]
rownames(Cibersort) <- Cibersort_all [1:nrow(Cibersort_all), 1]
rows_PDAC <- grep("Primary.Pancreas", rownames(Cibersort), value = TRUE)
Cibersort_PDAC <-  t(Cibersort [rows_PDAC,])
colnames(Cibersort_PDAC) <- GEO_ID
Cibersort_PDAC <- Cibersort_PDAC [,rownames(dataset_survival)]
Cibersort_PDAC_t <- t(Cibersort_PDAC)

data <- data.frame(
  time = as.numeric(dataset_survival[, "time"]),
  status = as.logical(dataset_survival[, "status"]),
  SLIT2 = as.numeric(SLIT2),
  SLIT3 = as.numeric(SLIT3),
  ROBO1 = as.numeric(ROBO1),
  Cibersort_PDAC_t
)
surv_object <- Surv(time = as.numeric(data$time), event = as.numeric(data$status))
cox_model <- coxph(surv_object ~ SLIT3 + Neutrophils, data = data)
summary(cox_model)
coefficients <- coef(cox_model)
data <- data %>%
  mutate(combined_score = coefficients["SLIT3"] * SLIT3 + coefficients["Neutrophils"] * Neutrophils)
data <- data %>%
  mutate(group = cut(combined_score, breaks = quantile(combined_score, probs = seq(0, 1, by = 0.25)), include.lowest = TRUE))
fit_combined <- survfit(surv_object ~ group, data = data)
ggsurvplot(
  fit_combined, 
  data = data, 
  pval = TRUE, 
  title = "Survival Curves",
  risk.table = TRUE, 
  tables.theme = theme_cleantable(),  
  tables.height = 0.3  
)
log_rank_test_combined <- survdiff(surv_object ~ group, data = data)

print(log_rank_test_combined)

neutrophils_median <- median(data$Neutrophils, na.rm = TRUE)
slit3_median <- median(data$SLIT3, na.rm = TRUE)

data <- data %>%
  mutate(
    Neutrophils_Level = ifelse(Neutrophils > neutrophils_median, "Alto", "Basso"),
    SLIT3_Level = ifelse(SLIT3 > slit3_median, "Alto", "Basso")
  )
output_table <- data.frame(
  Sample = rownames(data),
  Time = data$time,
  Combined_Score = data$combined_score,
  Neutrophils = data$Neutrophils,
  SLIT3 = data$SLIT3,
  Neutrophils_Level = data$Neutrophils_Level,
  SLIT3_Level = data$SLIT3_Level
)

#TRAINING DATA E TESTING DATA
set.seed(123)  
sample_indices <- sample(1:nrow(data), size = floor(0.3 * nrow(data)))

training_data <- data[sample_indices, ]
testing_data <- data[-sample_indices, ]

surv_object_train <- Surv(time = as.numeric(training_data$time), event = as.numeric(training_data$status))
cox_model_train <- coxph(surv_object_train ~ SLIT3 + Neutrophils, data = training_data)
summary(cox_model_train)
coefficients_train <- coef(cox_model_train)
testing_data <- testing_data %>%
  mutate(combined_score = coefficients_train["SLIT3"] * SLIT3 + coefficients_train["Neutrophils"] * Neutrophils)

testing_data <- testing_data %>%
  mutate(group = cut(combined_score, breaks = quantile(combined_score, probs = seq(0, 1, by = 0.25)), include.lowest = TRUE))
surv_object_test <- Surv(time = as.numeric(testing_data$time), event = as.numeric(testing_data$status))
fit_combined_test <- survfit(surv_object_test ~ group, data = testing_data)
ggsurvplot(
  fit_combined_test, 
  data = testing_data, 
  pval = TRUE, 
  title = "Survival Curves (Testing Data)",
  risk.table = TRUE, 
  tables.theme = theme_cleantable(),  
  tables.height = 0.3
)

