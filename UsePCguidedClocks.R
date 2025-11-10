#revised from codes provided by "A computational solution for bolstering reliability of epigenetic clocks, implications for clinical trials and longitudinal tracking" (https://doi.org/10.1038/s43587-022-00248-2)

#file requirements
#PC-guided clocks in a .RData file or other formats
#extremum of target sites among all samples in a .csv file, with rows corresponding to individual sites and columns corresponding to the maximum and minimum
#the order of the target sites should be identical between files

library(tidyverse)

load(file = "CalcPCSix.RData")
message("PC-guided Clocks Data successfully loaded")

#for sample age estimation, use methylation values of the sample instead of extreme methylation values
datMeth <- read.table("extremum_methylation_level.csv", header = TRUE,
  row.names = 1, sep = ",")
cgmax <- datMeth$Max
cgmin <- datMeth$Min

DNAmAge <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(DNAmAge) <- c("Max", "Min")

anti.trafo <- function(x, adult.age = 20) {
  ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age)
}

calculator <- function (clock, data, cgmax, cgmin, trafo = FALSE) {
  if (! all(c("model", "intercept", "center", "rotation") %in% names(clock))) {
    stop("PC-guided Clocks Data missing")
  }
  center <- clock$center
  rotation <- clock$rotation
  modelbeta <- clock$model
  intercept <- clock$intercept
  weightbined <- rotation %*% modelbeta
  if (all(row.names(data) == row.names(weightbined))) {
    message("CpGs are all in order")
    cgoptmax <- ifelse(weightbined > 0, cgmax, cgmin)
    cgoptmin <- ifelse(weightbined < 0, cgmax, cgmin)
    agemax <- sum(weightbined * (cgoptmax - center)) + intercept
    agemin <- sum(weightbined * (cgoptmin - center)) + intercept
    if (trafo) {
      agemax <- anti.trafo(agemax)
      agemin <- anti.trafo(agemin)
    }
    return(list(agemax, agemin))
  } else {
    stop(paste("Only", sum(row.names(data) == CpGs),
      "CpGs are in order"))
  }
}

DNAmAge[1, ] <- calculator(CalcPCHannum, datMeth, cgmax, cgmin)
DNAmAge[2, ] <- calculator(CalcPCHorvath1, datMeth, cgmax, cgmin, trafo = TRUE)
DNAmAge[3, ] <- calculator(CalcPCHorvath2, datMeth, cgmax, cgmin, trafo = TRUE)
DNAmAge[4, ] <- calculator(CalcPCLin, datMeth, cgmax, cgmin)
DNAmAge[5, ] <- calculator(CalcPCZhang, datMeth, cgmax, cgmin)
DNAmAge[6, ] <- calculator(CalcPCPhenoAge, datMeth, cgmax, cgmin)
message("PC-guided Clocks successfully calculated")
row.names(DNAmAge) <- c("PCHannum", "PCHorvath1", "PCHorvath2",
  "PCLin", "PCZhang", "PCPhenoAge")
rm(CalcPCHannum, CalcPCHorvath1, CalcPCHorvath2,
  CalcPCLin, CalcPCZhang, CalcPCPhenoAge)

write.table(DNAmAge, "age_limitation.csv", row.names = TRUE, sep = ",")