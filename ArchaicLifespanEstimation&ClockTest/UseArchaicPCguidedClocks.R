#file requirements
#PC-guided clocks in a .RData file or other formats
#information of target sites in two .csv files, with rows corresponding to individual sites and columns corresponding to chromosome position, specific values, and whether directly matched with the methylation map or not, including one file for the Altai Neandertal and the other for the Denisovan

load("CalcPCSix.RData")
message("PCClocks Data successfully loaded")
clockname <- c("PCHannum", "PCHorvath1", "PCHorvath2",
  "PCLin", "PCZhang", "PCPhenoAge")

datamn <- read.table("methylation_wanted_sites_Altai_Neandertal_PC.csv",
  header = TRUE, sep = ",")
datamd <- read.table("methylation_wanted_sites_Denisovan_PC.csv",
  header = TRUE, sep = ",")

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

datagearchaic <- data.frame(matrix(nrow = 6, ncol = 2), row.names = clockname)
colnames(datagearchaic) <- c("Neandertal", "Denisovan")

datagearchaic[1, 1] <- calculator(CalcPCHannum, datamn)
datagearchaic[2, 1] <- calculator(CalcPCHorvath1, datamn, trafo = TRUE)
datagearchaic[3, 1] <- calculator(CalcPCHorvath2, datamn, trafo = TRUE)
datagearchaic[4, 1] <- calculator(CalcPCLin, datamn)
datagearchaic[5, 1] <- calculator(CalcPCZhang, datamn)
datagearchaic[6, 1] <- calculator(CalcPCPhenoAge, datamn)

datagearchaic[1, 2] <- calculator(CalcPCHannum, datamd)
datagearchaic[2, 2] <- calculator(CalcPCHorvath1, datamd, trafo = TRUE)
datagearchaic[3, 2] <- calculator(CalcPCHorvath2, datamd, trafo = TRUE)
datagearchaic[4, 2] <- calculator(CalcPCLin, datamd)
datagearchaic[5, 2] <- calculator(CalcPCZhang, datamd)
datagearchaic[6, 2] <- calculator(CalcPCPhenoAge, datamd)

write.table(datagearchaic, "age_archaic_human_PC.csv",
  row.names = TRUE, sep = ",")