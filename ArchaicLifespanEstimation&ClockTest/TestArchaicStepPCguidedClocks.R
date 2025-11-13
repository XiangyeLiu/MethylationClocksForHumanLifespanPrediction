#file requirements
#PC-guided clocks in a .RData file or other formats
#information of target sites in two .csv files, with rows corresponding to individual sites and columns corresponding to chromosome position, specific values, and whether directly matched with the methylation map or not, including one file for the Altai Neandertal and the other for the Denisovan

library(tidyverse)
library(ggplot2)
library(patchwork)

load(file = "CalcPCSix.RData")
message("PCClocks Data successfully loaded")

calculator <- function (clock, data) {
  if (! all(c("model", "intercept", "center", "rotation") %in% names(clock))) {
    stop("PCClocks Data missing")
  }
  center <- clock$center
  rotation <- clock$rotation
  modelbeta <- clock$model
  intercept <- clock$intercept
  weightbined <- rotation %*% modelbeta
  if (all(data$Name == row.names(weightbined))) {
    message("CpGs are all in order")
  }
  mp <- weightbined * (data$Beta_value - center)
  sumg <- numeric(length(mp) + 1)
  sumg[1] <- intercept
  for (i in 2:length(sumg)) {
    sumg[i] <- sumg[i - 1] + mp[i - 1]
  }
  return(sumg)
}

mnp <- read.table("methylation_wanted_sites_Altai_Neandertal_PC.csv",
  header = TRUE, sep = ",")

npg <- data.frame("Number" = 0:nrow(mnp), "Site" = c("(Intercept)", mnp[, 1]))
npg[, 3] <- calculator(CalcPCHannum, mnp)
npg[, 4] <- anti.trafo(calculator(CalcPCHorvath1, mnp))
npg[, 5] <- anti.trafo(calculator(CalcPCHorvath2, mnp))
npg[, 6] <- calculator(CalcPCLin, mnp)
npg[, 7] <- calculator(CalcPCZhang, mnp)
npg[, 8] <- calculator(CalcPCPhenoAge, mnp)

colnames(npg)[3:ncol(npg)] <- c("PCHannum", "PCHorvath1",
  "PCHorvath2", "PCLin", "PCZhang", "PCPhenoAge")
write.table(npg, "cumulative_sum_age_prediction_Neandertal_PC.csv",
  row.names = FALSE, sep = ",")

npgl <- pivot_longer(npg[, -8], "PCHannum":"PCZhang",
  names_to = "Clock", values_to = "Value")
write.table(npgl, "cumulative_sum_age_long_Neandertal_PC.csv",
  row.names = FALSE, sep = ",")

mdp <- read.table("methylation_wanted_sites_Denisovan_PC.csv",
  header = TRUE, sep = ",")

dpg <- data.frame("Number" = 0:nrow(mdp), "Site" = c("(Intercept)", mdp[, 1]))
dpg[, 3] <- calculator(CalcPCHannum, mdp)
dpg[, 4] <- anti.trafo(calculator(CalcPCHorvath1, mdp))
dpg[, 5] <- anti.trafo(calculator(CalcPCHorvath2, mdp))
dpg[, 6] <- calculator(CalcPCLin, mdp)
dpg[, 7] <- calculator(CalcPCZhang, mdp)
dpg[, 8] <- calculator(CalcPCPhenoAge, mdp)

colnames(dpg)[3:ncol(dpg)] <- c("PCHannum", "PCHorvath1",
  "PCHorvath2", "PCLin", "PCZhang", "PCPhenoAge")
write.table(dpg, "cumulative_sum_age_prediction_Denisovan_PC.csv",
  row.names = FALSE, sep = ",")

dpgl <- pivot_longer(dpg[, -8], "PCHannum":"PCZhang",
  names_to = "Clock", values_to = "Value")
write.table(dpgl, "cumulative_sum_age_long_Denisovan_PC.csv",
  row.names = FALSE, sep = ",")