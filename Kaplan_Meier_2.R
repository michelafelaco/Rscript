library(survival)
library(survminer)
library(readxl)
library(maxstat)

#modello Cox
#importazione dati
dataset <- read.table("C:\\Users\\miche\\Desktop\\CNR lab\\Trascriptomic_R\\Moffit_database\\GSE71729_series_matrix.txt\\GSE71729_series_matrix.txt", header = TRUE, skip = 26, sep = "\t", fill = TRUE)
dataset <- dataset [-nrow(dataset),]
matrix <- as.matrix(dataset[44:nrow(dataset),2:ncol(dataset)])
numeric_matrix_Moffit <- apply(matrix, 2, as.numeric)
rownames(numeric_matrix_Moffit) <- dataset[44:nrow(dataset), 1]
SLIT2 <-  numeric_matrix_Moffit ["SLIT2",]
SLIT3 <- numeric_matrix_Moffit ["SLIT3",]
ROBO1 <- numeric_matrix_Moffit ["ROBO1",]
GEO_ID <- dataset[1,2:ncol(dataset)]
dataset_survival <- as.matrix(read_xlsx("C:\\Users\\miche\\Desktop\\CNR lab\\Trascriptomic_R\\Moffit_database\\SLIT_ROBO\\Santiago Poster\\Tumor_survival.xlsx"))
row.names(dataset_survival) <- as.vector(dataset_survival[,1])
dataset_survival <- dataset_survival [,2:4]
names(SLIT2) <- GEO_ID
SLIT2 <- t(SLIT2 [rownames(dataset_survival)])
names(SLIT3) <- GEO_ID
SLIT3 <- (SLIT3 [rownames(dataset_survival)])
names(ROBO1) <- GEO_ID
ROBO1 <- (ROBO1 [rownames(dataset_survival)])

data <- data.frame(
  time = as.numeric(dataset_survival [,"time"] ),
  status = as.logical(dataset_survival [,"status"] ),
  SLIT2 = as.numeric(SLIT2),
  SLIT3 = as.numeric(SLIT3),
  ROBO1 = as.numeric(ROBO1))
surv_object <- Surv(time = as.numeric(data$time), event = as.numeric(data$status))

cox_model <- coxph(surv_object ~ SLIT2+SLIT3, data = data)
summary(cox_model)

data$risk_score <- predict(cox_model, type = "risk")
median_risk <- median(data$risk_score)
data$risk_group <- ifelse(data$risk_score > median_risk, "High Risk", "Low Risk")
fit <- survfit(surv_object ~ risk_group, data = data)

ggsurvplot(fit, data = data, pval = TRUE, 
           risk.table = TRUE, 
           ggtheme = theme_minimal(), 
           title = "Kaplan-Meier Curve by Risk Group based on Cox Model")

#modello lineare
linear_model <- lm(time ~ ROBO1 + SLIT3, data = data)
data$risk_score <- predict(linear_model, data)
median_risk <- median(data$risk_score, na.rm = TRUE)
data$risk_group <- ifelse(data$risk_score > median_risk, "High Risk", "Low Risk")
surv_object <- Surv(time = as.numeric(data$time), event = as.numeric(data$status))
fit <- survfit(surv_object ~ risk_group, data = data)
ggsurvplot(fit, data = data, pval = TRUE, 
           risk.table = TRUE, 
           ggtheme = theme_minimal(), 
           title = "Kaplan-Meier Curve by Risk Group based on Linear Model")

#ENTRO 500 GG
library(survival)
library(survminer)
library(readxl)
library(maxstat)

# Modello Cox
# Importazione dati
dataset <- read.table("C:\\Users\\miche\\Desktop\\CNR lab\\Trascriptomic_R\\Moffit_database\\GSE71729_series_matrix.txt\\GSE71729_series_matrix.txt", header = TRUE, skip = 26, sep = "\t", fill = TRUE)
dataset <- dataset[-nrow(dataset),]
matrix <- as.matrix(dataset[44:nrow(dataset), 2:ncol(dataset)])
numeric_matrix_Moffit <- apply(matrix, 2, as.numeric)
rownames(numeric_matrix_Moffit) <- dataset[44:nrow(dataset), 1]
SLIT2 <- numeric_matrix_Moffit["SLIT2",]
SLIT3 <- numeric_matrix_Moffit["SLIT3",]
ROBO1 <- numeric_matrix_Moffit["ROBO1",]
GEO_ID <- dataset[1, 2:ncol(dataset)]
dataset_survival <- as.matrix(read_xlsx("C:\\Users\\miche\\Desktop\\CNR lab\\Trascriptomic_R\\Moffit_database\\SLIT_ROBO\\Santiago Poster\\Tumor_survival.xlsx"))
row.names(dataset_survival) <- as.vector(dataset_survival[, 1])
dataset_survival <- dataset_survival[, 2:4]
names(SLIT2) <- GEO_ID
SLIT2 <- SLIT2[rownames(dataset_survival)]
names(SLIT3) <- GEO_ID
SLIT3 <- SLIT3[rownames(dataset_survival)]
names(ROBO1) <- GEO_ID
ROBO1 <- ROBO1[rownames(dataset_survival)]

data <- data.frame(
  time = as.numeric(dataset_survival[, "time"]),
  status = as.logical(dataset_survival[, "status"]),
  SLIT2 = as.numeric(SLIT2),
  SLIT3 = as.numeric(SLIT3),
  ROBO1 = as.numeric(ROBO1)
)

# Censurare artificialmente i dati oltre 500 giorni
data$time <- pmin(data$time, 500)
data$status <- ifelse(data$time == 500, 0, data$status)

# Creare l'oggetto Surv
surv_object <- Surv(time = as.numeric(data$time), event = as.numeric(data$status))

# Fit del modello di regressione di Cox
cox_model <- coxph(surv_object ~ SLIT2 + SLIT3, data = data)
summary(cox_model)

# Calcolare il punteggio di rischio
data$risk_score <- predict(cox_model, type = "risk")
median_risk <- median(data$risk_score)
data$risk_group <- ifelse(data$risk_score > median_risk, "High Risk", "Low Risk")

# Fit del modello di sopravvivenza basato sui gruppi di rischio
fit <- survfit(surv_object ~ risk_group, data = data)

# Plot delle curve di Kaplan-Meier
ggsurvplot(fit, data = data, pval = TRUE, 
           risk.table = TRUE, 
           ggtheme = theme_minimal(), 
           title = "Kaplan-Meier Curve by Risk Group based on Cox Model with 500 days limit")

