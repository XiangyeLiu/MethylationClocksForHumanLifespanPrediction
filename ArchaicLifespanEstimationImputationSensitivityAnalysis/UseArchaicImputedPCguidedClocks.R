#file requirements
#PC-guided clocks in a .RData file or other formats
#information of imputed target sites in six .csv files, with rows corresponding to individual sites and columns corresponding to specific imputed values, each imputation strategy including one file for the Altai Neandertal and the other for the Denisovan

load("CalcPCBigSix.RData")
message("PCClocks Data successfully loaded")
clockname <- c("PCHannum", "PCHorvath1", "PCHorvath2",
  "PCLin", "PCZhang", "PCPhenoAge")

anti.trafo <- function(x, adultage = 20) {
  ifelse(x < 0, (1 + adultage) * exp(x) - 1, (1 + adultage) * x + adultage)
}

calculator <- function (clock, data, trafo = FALSE) {
  if (! all(c("model", "intercept", "center", "rotation") %in% names(clock))) {
    stop("PCClocks Data Missing")
  }
  center <- clock$center
  rotation <- clock$rotation
  modelbeta <- clock$model
  intercept <- clock$intercept
  weightbined <- rotation %*% modelbeta
  if (all(data$Name == row.names(weightbined))) {
    message("CpGs Are All Matched")
    age <- sum(weightbined * (data$Beta_value - center)) + intercept
    if (trafo) {
      age <- anti.trafo(age)
    }
    return(as.numeric(age))
  } else {
    stop(paste("Only", sum(row.names(data) == CpGs),
      "CpGs Are Matched"))
  }
}

agecalstep <- function(data) {
  if (ncol(data) != 5) {
    stop("Data format Does not Meet the Requirements")
  }
  colnames(data)[5] <- "Beta_value"
  age <- numeric()
  age[1] <- calculator(CalcPCHannum, data)
  age[2] <- calculator(CalcPCHorvath1, data, trafo = TRUE)
  age[3] <- calculator(CalcPCHorvath2, data, trafo = TRUE)
  age[4] <- calculator(CalcPCLin, data)
  age[5] <- calculator(CalcPCZhang, data)
  age[6] <- calculator(CalcPCPhenoAge, data)
  return(age)
}

dnamni <- read.table("imputation_pmm_Altai_Neandertal_PC.csv",
  header = TRUE, sep = ",")
dnamdi <- read.table("imputation_pmm_Denisovan_pc.csv",
  header = TRUE, sep = ",")

agearchaic <- agecalstep(dnamni[, -6])
agearchaic <- c(agearchaic, agecalstep(dnamni[, -5]))
agearchaic <- c(agearchaic, agecalstep(dnamdi[, -6]))
agearchaic <- c(agearchaic, agecalstep(dnamdi[, -5]))

datagearchaic1 <- data.frame(Clock = clockname, Imputation = agearchaic,
  Type = rep(c("Neanderthal", "Denisovan"), each = 12),
  Method = rep(c("Mean", "Median"), each = 6))

dnamni <- read.table("imputation_knn_Altai_Neandertal_PC.csv",
  header = TRUE, sep = ",")
dnamdi <- read.table("imputation_knn_Denisovan_PC.csv",
  header = TRUE, sep = ",")

agearchaic <- matrix(nrow = 6, ncol = 200)
for (k in 1:100) {
  agearchaic[, k] <- agecalstep(dnamni[, c(1:4, 4 + k)])
  agearchaic[, 100 + k] <- agecalstep(dnamdi[, c(1:4, 4 + k)])
}

datagearchaic2 <- data.frame(Clock = clockname, Imputation = as.vector(agearchaic),
  Type = rep(c("Neanderthal", "Denisovan"), each = 600), Method = "KNN")

dnamnia <- read.table("imputation_sva_Altai_Neandertal_PC.csv",
  header = TRUE, sep = ",")
dnamdia <- read.table("imputation_sva_Denisovan_PC.csv",
  header = TRUE, sep = ",")

dnamn <- read.table("methylation_wanted_sites_Altai_Neandertal_PC.csv",
  header = TRUE, sep = ",")
dnamd <- read.table("methylation_wanted_sites_Denisovan_PC.csv",
  header = TRUE, sep = ",")

dnamni <- dnamn[, 1:5]
dnamdi <- dnamd[, 1:5]
agearchaic <- matrix(nrow = 6, ncol = 40000)
for (i in 1:20000) {
  dnamni[dnamnia$Number, 5] <- dnamnia[, 2 + i]
  dnamdi[dnamdia$Number, 5] <- dnamdia[, 2 + i]
  agearchaic[, i] <- agecalstep(dnamni)
  agearchaic[, 20000 + i] <- agecalstep(dnamdi)
}

datagearchaic3 <- data.frame(Clock = clockname, Imputation = as.numeric(agearchaic),
  Type = rep(c("Neanderthal", "Denisovan"), each = 120000),
  Method = rep(c("Sample", "Random"), each = 60000))

datagearchaic <- rbind(datagearchaic1, datagearchaic2, datagearchaic3)

write.table(datagearchaic, "imputation_result_long_PC.csv",
  row.names = FALSE, sep = ",")
