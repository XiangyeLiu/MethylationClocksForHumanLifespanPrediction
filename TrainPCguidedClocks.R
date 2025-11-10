#revised from codes provided by "A computational solution for bolstering reliability of epigenetic clocks, implications for clinical trials and longitudinal tracking" (https://doi.org/10.1038/s43587-022-00248-2)

#file requirements
#target methylation sites in a .txt file
#methylation values of target sites for each sample in a .csv file, with target sites in rows and samples in columns
#ages or other phenotypes of samples in a .csv file, with samples in rows and ages in column

library(glmnet)

CpGs <- readLines("shared_methylation_sites.txt")
datMethTrain <- read.table("PC_clock_train_meth.csv",
  header = TRUE, row.names = 1, sep = ",")
datMethTrain <- t(datMethTrain)
datPhenoTrain <- read.table("PC_clock_train_pheno.csv",
  header = TRUE, row.names = 1, sep = ",")

meanimpute <- function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
}
datMethTrain <- apply(datMethTrain, 2, meanimpute)

if (all(colnames(datMethTrain) == CpGs)) {
  message("CpGs are all in order")
} else {
  message(paste("Only", sum(colnames(datMethTrain) == CpGs),
    "CpGs are in order"))
}

if (all(rownames(datMethTrain) == rownames(datPhenoTrain))) {
  message("Samples are all in order")
} else {
  message(paste("Only", sum(rownames(datMethTrain) == rownames(datPhenoTrain)),
    "Samples are in order"))
}

PCA <- prcomp(datMethTrain, scale. = FALSE)
TrainPCData <- PCA$x[, 1:(dim(PCA$x)[2] - 1)]

TrainAge <- datPhenoTrain$Age

cv <- cv.glmnet(TrainPCData, TrainAge, family = "gaussian",
  alpha = 0.5, nfolds = 10)
fit <- glmnet(TrainPCData, TrainAge, family = "gaussian",
  alpha = 0.5, nlambda = 100)
plot(cv)

plot(TrainAge, predict(fit, TrainPCData, s = cv$lambda.min),
  xlab = "Age", ylab = "Predicted Age", main = "Training")
cor(TrainAge, predict(fit, TrainPCData, s = cv$lambda.min))

#optional
plot(TrainAge, predict(fit, TrainPCData, s = cv$lambda.1se),
  xlab = "Age", ylab = "Predicted Age", main = "Training")
cor(TrainAge, predict(fit, TrainPCData, s = cv$lambda.1se))

#CalcPCAge can be changed to other names
CalcPCAge <- vector(mode = "list", length = 0)
temp <- as.matrix(coef(cv, s = cv$lambda.min))
CalcPCAge$model <- temp[temp != 0, ][-1]
CalcPCAge$intercept <- temp[1, 1]
CalcPCAge$center <- PCA$center
CalcPCAge$rotation <- PCA$rotation[, names(CalcPCAge$model)]

save(CalcPCAge, CpGs, file = "CalcPCAge.RData")