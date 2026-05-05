#file requirements
#information of imputed target sites in six .csv files, with rows corresponding to individual sites and columns corresponding to chromosome position, coefficients, and specific imputed values, each clock including one file for the Altai Neandertal and the other for the Denisovan

panmcuseni1 <- read.table("pan_mammalian_clock1_hg19_imputation_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

panmcuseni2 <- read.table("pan_mammalian_clock2_hg19_imputation_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

panmcuseni3 <- read.table("pan_mammalian_clock3_hg19_imputation_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

panmcusedi1 <- read.table("pan_mammalian_clock1_hg19_imputation_Denisovan.csv",
  header = TRUE, sep = ",")

panmcusedi2 <- read.table("pan_mammalian_clock2_hg19_imputation_Denisovan.csv",
  header = TRUE, sep = ",")

panmcusedi3 <- read.table("pan_mammalian_clock3_hg19_imputation_Denisovan.csv",
  header = TRUE, sep = ",")

panmcuse <- function(coefficient, value, type = 1, max = 65.7, gestation = 280 / 365, asm = 13) {
  intercept <- coefficient[1]
  result <- sum(coefficient[-1] * value[-1], intercept)
  if (type == 1) {
    message("Using Clock 1")
    return(exp(result) - 2)
  } else if (type == 2) {
    message("Using Clock 2")
    return(exp(-1 * exp(-1 * result)) * (max + gestation) - gestation)
  } else if (type == 3) {
    message("Using Clock 3")
    m <- 5.0 * (gestation / asm) ^ 0.38
    if (result >= 0) {
      return(m * (asm + gestation) * (result + 1) - gestation)
    } else {
      return(m * (asm + gestation) * exp(result) - gestation)
    }
  } else {
    stop("Wrong Clock, Parameter Should Be 1-3")
  }
}

method <- c("Mean", "Median", rep("KNN(K=1-100)", 100),
  rep("Sample One", 10000), rep("Stochastic Value", 10000),
  "Max Value")

ageni1 <- sapply(panmcuseni1[, -1:-3], function(x) {
  panmcuse(panmcuseni1[, 3], x, 1)
})
ageni2 <- sapply(panmcuseni2[, -1:-3], function(x) {
  panmcuse(panmcuseni2[, 3], x, 2)
})
ageni3 <- sapply(panmcuseni3[, -1:-3], function(x) {
  panmcuse(panmcuseni3[, 3], x, 3)
})

agearchaic1 <- data.frame(Clock = rep(paste("Clock", 1:3, sep = ""),
  each = length(ageni1)), Age = c(ageni1, ageni2, ageni3),
  Type = "Neanderthal", Method = method)

agedi1 <- sapply(panmcusedi1[, -1:-3], function(x) {
  panmcuse(panmcusedi1[, 3], x, 1, max = 69.8)
})
agedi2 <- sapply(panmcusedi2[, -1:-3], function(x) {
  panmcuse(panmcusedi2[, 3], x, 2, max = 69.8)
})
agedi3 <- sapply(panmcusedi3[, -1:-3], function(x) {
  panmcuse(panmcusedi3[, 3], x, 3, max = 69.8)
})

agearchaic2 <- data.frame(Clock = rep(paste("Clock", 1:3, sep = ""),
  each = length(agedi1)), Age = c(agedi1, agedi2, agedi3),
  Type = "Denisovan", Method = method)

datagearchaic <- rbind(agearchaic1, agearchaic2)
datagearchaic$Method <- factor(datagearchaic$Method,
  levels = unique(method))
datagearchaic$Type <- factor(datagearchaic$Type,
  levels = c("Neanderthal", "Denisovan"))
datagearchaic <- datagearchaic[order(datagearchaic$Method), ]

write.table(datagearchaic, "age_pan_mammalian_imputation.csv",
  row.names = FALSE, sep = ",")