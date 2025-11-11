#file requirements
#PC-guided clocks in a .RData file or other formats
#tissue/organ-specific extremum of target sites among samples in two .csv files, with rows corresponding to individual sites and columns corresponding to series, including one file for the maximum and the other for the minimum

library(tidyverse)
library(stringr)

load(file = "CalcPCBigSix.RData")
message("PCClocks Data successfully loaded")

data3 <- read.table("extremum_split_max.csv", header = TRUE,
  row.names = 1, sep = ",")
data4 <- read.table("extremum_split_min.csv", header = TRUE,
  row.names = 1, sep = ",")

agemaxgroup <- data.frame(matrix(nrow = 6, ncol = 0),
  row.names = c("PCHannum", "PCHorvath1", "PCHorvath2",
  "PCLin", "PCZhang", "PCPhenoAge"))
agemingroup <- data.frame(matrix(nrow = 6, ncol = 0),
  row.names = c("PCHannum", "PCHorvath1", "PCHorvath2",
  "PCLin", "PCZhang", "PCPhenoAge"))

calculator <- function (clock, cgmax, cgmin, trafo = FALSE) {
  if (! all(c("model", "intercept", "center", "rotation") %in% names(clock))) {
    stop("PCClocks Data missing")
  }
  center <- clock$center
  rotation <- clock$rotation
  modelbeta <- clock$model
  intercept <- clock$intercept
  weightbined <- rotation %*% modelbeta
  cgoptmax <- ifelse(weightbined > 0, cgmax, cgmin)
  cgoptmin <- ifelse(weightbined < 0, cgmax, cgmin)
  agemax <- sum(weightbined * (cgoptmax - center)) + intercept
  agemin <- sum(weightbined * (cgoptmin - center)) + intercept
  if (trafo) {
    agemax <- anti.trafo(agemax)
    agemin <- anti.trafo(agemin)
  }
  return(list(agemax, agemin))
}

if (all(row.names(data3) == row.names(data4))) {
  message("Data are all matched")
} else {
  message("Data matching problem")
}

for (i in seq_along(data3)) {
  cgmax <- data3[, i]
  cgmin <- data4[, i]
  agemaxgroup[1, i] <- calculator(CalcPCHannum, cgmax, cgmin)[[1]]
  agemaxgroup[2, i] <- calculator(CalcPCHorvath1, cgmax, cgmin, TRUE)[[1]]
  agemaxgroup[3, i] <- calculator(CalcPCHorvath2, cgmax, cgmin, TRUE)[[1]]
  agemaxgroup[4, i] <- calculator(CalcPCLin, cgmax, cgmin)[[1]]
  agemaxgroup[5, i] <- calculator(CalcPCZhang, cgmax, cgmin)[[1]]
  agemaxgroup[6, i] <- calculator(CalcPCPhenoAge, cgmax, cgmin)[[1]]
  agemingroup[1, i] <- calculator(CalcPCHannum, cgmax, cgmin)[[2]]
  agemingroup[2, i] <- calculator(CalcPCHorvath1, cgmax, cgmin, TRUE)[[2]]
  agemingroup[3, i] <- calculator(CalcPCHorvath2, cgmax, cgmin, TRUE)[[2]]
  agemingroup[4, i] <- calculator(CalcPCLin, cgmax, cgmin)[[2]]
  agemingroup[5, i] <- calculator(CalcPCZhang, cgmax, cgmin)[[2]]
  agemingroup[6, i] <- calculator(CalcPCPhenoAge, cgmax, cgmin)[[2]]
}
message("PC Clocks successfully calculated")
colnames(agemaxgroup) <- colnames(data3)
colnames(agemingroup) <- colnames(data4)
rm(CalcPCHannum, CalcPCHorvath1, CalcPCHorvath2,
  CalcPCLin, CalcPCZhang, CalcPCPhenoAge)

write.table(agemaxgroup, "age_split_max.csv", row.names = TRUE, sep = ",")
write.table(agemingroup, "age_split_min.csv", row.names = TRUE, sep = ",")
