#file requirements
#information of target sites from 450K Beadchip to compute unmatched sites in two .RDS file, including individual sites and coefficients, excerpted from "Universal DNA methylation age across mammalian tissues" (https://doi.org/10.1038/s43587-023-00462-6)
#quality control using relation between each unmatched site and its target sites from 450K Beadchip in a .RDS file, obtained from "Universal DNA methylation age across mammalian tissues" (https://doi.org/10.1038/s43587-023-00462-6)
#information of target sites from 450K Beadchip in six .csv files, with rows corresponding to individual sites and columns corresponding to chromosome position, specific values, and whether directly matched with the methylation map or not, each clock including one file for the Altai Neandertal and the other for the Denisovan
#based on Array Converter Algorithm

fitkeepn <- readRDS("fit_EPICtoMM_v2_450k_Altai_Neandertal.RDS")

fitkeepd <- readRDS("fit_EPICtoMM_v2_450k_Denisovan.RDS")

qc <- read.table("MeasuresTEST_EPICtoMM_v2_450k.csv", header = TRUE, sep = ",")
colnames(qc)[1] <- "Name"

confirmvalue <- function(cgwv, test) {
  if (is.na(cgwv)) {
    return(0.5)
  } else if (test < 0.6) {
    return(0.5)
  } else {
    return(cgwv)
  }
}

panmcusen1 <- read.table("pan_mammalian_clock1_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

panmcusen2 <- read.table("pan_mammalian_clock2_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

panmcusen3 <- read.table("pan_mammalian_clock3_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

cgneedvn <- read.table("array_convert_value_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

cgunmatched1 <- cgunmatched2 <- cgunmatched3 <- NULL
cgunmatched4 <- cgunmatched5 <- NULL
for (i in seq_along(fitkeepn)) {
  name <- names(fitkeepn[i])
  cgqc <- qc[which(qc$Name == name), ]
  test <- min(c(cgqc$bicor_test, cgqc$cor_test))
  intercept <- fitkeepn[[i]][1, 2]
  cgneed <- merge(fitkeepn[[i]], cgneedvn[, c(1, 3)],
    by.x = "EPIC_CG", by.y = "Name")
  cgna <- which(is.na(cgneed$Beta_value))

  cgwv1 <- sum(cgneed$coef * cgneed$Beta_value, intercept)
  cgwv1 <- exp(cgwv1) / (exp(cgwv1) + 1)
  cgwv1 <- confirmvalue(cgwv1, test)
  names(cgwv1) <- name
  cgunmatched1 <- c(cgunmatched1, cgwv1)

  cgneed$Beta_value[cgna] <- 0
  cgwv2 <- sum(cgneed$coef * cgneed$Beta_value, intercept)
  cgwv2 <- exp(cgwv2) / (exp(cgwv2) + 1)
  cgwv2 <- confirmvalue(cgwv2, test)
  names(cgwv2) <- name
  cgunmatched2 <- c(cgunmatched2, cgwv2)

  cgneed$Beta_value[cgna] <- 0.5
  cgwv3 <- sum(cgneed$coef * cgneed$Beta_value, intercept)
  cgwv3 <- exp(cgwv3) / (exp(cgwv3) + 1)
  cgwv3 <- confirmvalue(cgwv3, test)
  names(cgwv3) <- name
  cgunmatched3 <- c(cgunmatched3, cgwv3)
  
  cgneed$Beta_value[cgna] <- ifelse(cgneed$coef[cgna] > 0, 1, 0)
  cgwv4 <- sum(cgneed$coef * cgneed$Beta_value, intercept)
  cgwv4 <- exp(cgwv4) / (exp(cgwv4) + 1)
  cgwv4 <- confirmvalue(cgwv4, test)
  names(cgwv4) <- name
  cgunmatched4 <- c(cgunmatched4, cgwv4)

  cgneed$Beta_value[cgna] <- ifelse(cgneed$coef[cgna] > 0, 0, 1)
  cgwv5 <- sum(cgneed$coef * cgneed$Beta_value, intercept)
  cgwv5 <- exp(cgwv5) / (exp(cgwv5) + 1)
  cgwv5 <- confirmvalue(cgwv5, test)
  names(cgwv5) <- name
  cgunmatched5 <- c(cgunmatched5, cgwv5)

  cgneed$Beta_value[cgna] <- NA
}

pumn1 <- which(panmcusen1$Merge == FALSE & panmcusen1$Coefficient_1 > 0)
numn1 <- which(panmcusen1$Merge == FALSE & panmcusen1$Coefficient_1 < 0)
umn1 <- sort(c(pumn1, numn1))
cgpumn1 <- panmcusen1$Name[pumn1]
cgnumn1 <- panmcusen1$Name[numn1]
cgumn1 <- sort(c(cgpumn1, cgnumn1))
panmcusen1 <- within(panmcusen1, {
  Beta_value1 <- Beta_value2 <- Beta_value3 <- Beta_value4 <- Beta_value
  Beta_value1[umn1] <- cgunmatched1[which(names(cgunmatched1) %in% cgumn1)]
  Beta_value2[umn1] <- cgunmatched2[which(names(cgunmatched2) %in% cgumn1)]
  Beta_value3[umn1] <- cgunmatched3[which(names(cgunmatched3) %in% cgumn1)]
  Beta_value4[pumn1] <- cgunmatched4[which(names(cgunmatched4) %in% cgpumn1)]
  Beta_value4[numn1] <- cgunmatched5[which(names(cgunmatched5) %in% cgnumn1)]
})
panmcusen4 <- panmcusen1[, c(1:3, 8:11)]
write.table(panmcusen4, "pan_mammalian_clock1_hg19_convert_Altai_Neandertal.csv",
  row.names = FALSE, sep = ",")

pumn2 <- which(panmcusen2$Merge == FALSE & panmcusen2$Coefficient_2 > 0)
numn2 <- which(panmcusen2$Merge == FALSE & panmcusen2$Coefficient_2 < 0)
umn2 <- sort(c(pumn2, numn2))
cgpumn2 <- panmcusen2$Name[pumn2]
cgnumn2 <- panmcusen2$Name[numn2]
cgumn2 <- sort(c(cgpumn2, cgnumn2))
panmcusen2 <- within(panmcusen2, {
  Beta_value1 <- Beta_value2 <- Beta_value3 <- Beta_value4 <- Beta_value
  Beta_value1[umn2] <- cgunmatched1[which(names(cgunmatched1) %in% cgumn2)]
  Beta_value2[umn2] <- cgunmatched2[which(names(cgunmatched2) %in% cgumn2)]
  Beta_value3[umn2] <- cgunmatched3[which(names(cgunmatched3) %in% cgumn2)]
  Beta_value4[pumn2] <- cgunmatched4[which(names(cgunmatched4) %in% cgpumn2)]
  Beta_value4[numn2] <- cgunmatched5[which(names(cgunmatched5) %in% cgnumn2)]
})
panmcusen5 <- panmcusen2[, c(1:3, 8:11)]
write.table(panmcusen5, "pan_mammalian_clock2_hg19_convert_Altai_Neandertal.csv",
  row.names = FALSE, sep = ",")

pumn3 <- which(panmcusen3$Merge == FALSE & panmcusen3$Coefficient_3 > 0)
numn3 <- which(panmcusen3$Merge == FALSE & panmcusen3$Coefficient_3 < 0)
umn3 <- sort(c(pumn3, numn3))
cgpumn3 <- panmcusen3$Name[pumn3]
cgnumn3 <- panmcusen3$Name[numn3]
cgumn3 <- sort(c(cgpumn3, cgnumn3))
panmcusen3 <- within(panmcusen3, {
  Beta_value1 <- Beta_value2 <- Beta_value3 <- Beta_value4 <- Beta_value
  Beta_value1[umn3] <- cgunmatched1[which(names(cgunmatched1) %in% cgumn3)]
  Beta_value2[umn3] <- cgunmatched2[which(names(cgunmatched2) %in% cgumn3)]
  Beta_value3[umn3] <- cgunmatched3[which(names(cgunmatched3) %in% cgumn3)]
  Beta_value4[pumn3] <- cgunmatched4[which(names(cgunmatched4) %in% cgpumn3)]
  Beta_value4[numn3] <- cgunmatched5[which(names(cgunmatched5) %in% cgnumn3)]
})
panmcusen6 <- panmcusen3[, c(1:3, 8:11)]
write.table(panmcusen6, "pan_mammalian_clock3_hg19_convert_Altai_Neandertal.csv",
  row.names = FALSE, sep = ",")

panmcused1 <- read.table("pan_mammalian_clock1_hg19_Denisovan.csv",
  header = TRUE, sep = ",")

panmcused2 <- read.table("pan_mammalian_clock2_hg19_Denisovan.csv",
  header = TRUE, sep = ",")

panmcused3 <- read.table("pan_mammalian_clock3_hg19_Denisovan.csv",
  header = TRUE, sep = ",")

cgneedvd <- read.table("array_convert_value_Denisovan.csv",
  header = TRUE, sep = ",")

cgunmatched1 <- cgunmatched2 <- cgunmatched3 <- NULL
cgunmatched4 <- cgunmatched5 <- NULL
for (i in seq_along(fitkeepd)) {
  name <- names(fitkeepd[i])
  cgqc <- qc[which(qc$Name == name), ]
  test <- min(c(cgqc$bicor_test, cgqc$cor_test))
  intercept <- fitkeepd[[i]][1, 2]
  cgneed <- merge(fitkeepd[[i]], cgneedvd[, c(1, 3)],
    by.x = "EPIC_CG", by.y = "Name")
  cgna <- which(is.na(cgneed$Beta_value))

  cgwv1 <- sum(cgneed$coef * cgneed$Beta_value, intercept)
  cgwv1 <- exp(cgwv1) / (exp(cgwv1) + 1)
  cgwv1 <- confirmvalue(cgwv1, test)
  names(cgwv1) <- name
  cgunmatched1 <- c(cgunmatched1, cgwv1)

  cgneed$Beta_value[cgna] <- 0
  cgwv2 <- sum(cgneed$coef * cgneed$Beta_value, intercept)
  cgwv2 <- exp(cgwv2) / (exp(cgwv2) + 1)
  cgwv2 <- confirmvalue(cgwv2, test)
  names(cgwv2) <- name
  cgunmatched2 <- c(cgunmatched2, cgwv2)

  cgneed$Beta_value[cgna] <- 0.5
  cgwv3 <- sum(cgneed$coef * cgneed$Beta_value, intercept)
  cgwv3 <- exp(cgwv3) / (exp(cgwv3) + 1)
  cgwv3 <- confirmvalue(cgwv3, test)
  names(cgwv3) <- name
  cgunmatched3 <- c(cgunmatched3, cgwv3)
  
  cgneed$Beta_value[cgna] <- ifelse(cgneed$coef[cgna] > 0, 1, 0)
  cgwv4 <- sum(cgneed$coef * cgneed$Beta_value, intercept)
  cgwv4 <- exp(cgwv4) / (exp(cgwv4) + 1)
  cgwv4 <- confirmvalue(cgwv4, test)
  names(cgwv4) <- name
  cgunmatched4 <- c(cgunmatched4, cgwv4)

  cgneed$Beta_value[cgna] <- ifelse(cgneed$coef[cgna] > 0, 0, 1)
  cgwv5 <- sum(cgneed$coef * cgneed$Beta_value, intercept)
  cgwv5 <- exp(cgwv5) / (exp(cgwv5) + 1)
  cgwv5 <- confirmvalue(cgwv5, test)
  names(cgwv5) <- name
  cgunmatched5 <- c(cgunmatched5, cgwv5)

  cgneed$Beta_value[cgna] <- NA
}

pumd1 <- which(panmcused1$Merge == FALSE & panmcused1$Coefficient_1 > 0)
numd1 <- which(panmcused1$Merge == FALSE & panmcused1$Coefficient_1 < 0)
umd1 <- sort(c(pumd1, numd1))
cgpumd1 <- panmcused1$Name[pumd1]
cgnumd1 <- panmcused1$Name[numd1]
cgumd1 <- sort(c(cgpumd1, cgnumd1))
panmcused1 <- within(panmcused1, {
  Beta_value1 <- Beta_value2 <- Beta_value3 <- Beta_value4 <- Beta_value
  Beta_value1[umd1] <- cgunmatched1[which(names(cgunmatched1) %in% cgumd1)]
  Beta_value2[umd1] <- cgunmatched2[which(names(cgunmatched2) %in% cgumd1)]
  Beta_value3[umd1] <- cgunmatched3[which(names(cgunmatched3) %in% cgumd1)]
  Beta_value4[pumd1] <- cgunmatched4[which(names(cgunmatched4) %in% cgpumd1)]
  Beta_value4[numd1] <- cgunmatched5[which(names(cgunmatched5) %in% cgnumd1)]
})
panmcused4 <- panmcused1[, c(1:3, 8:11)]
write.table(panmcused4, "pan_mammalian_clock1_hg19_convert_Denisovan.csv",
  row.names = FALSE, sep = ",")

pumd2 <- which(panmcused2$Merge == FALSE & panmcused2$Coefficient_2 > 0)
numd2 <- which(panmcused2$Merge == FALSE & panmcused2$Coefficient_2 < 0)
umd2 <- sort(c(pumd2, numd2))
cgpumd2 <- panmcused2$Name[pumd2]
cgnumd2 <- panmcused2$Name[numd2]
cgumd2 <- sort(c(cgpumd2, cgnumd2))
panmcused2 <- within(panmcused2, {
  Beta_value1 <- Beta_value2 <- Beta_value3 <- Beta_value4 <- Beta_value
  Beta_value1[umd2] <- cgunmatched1[which(names(cgunmatched1) %in% cgumd2)]
  Beta_value2[umd2] <- cgunmatched2[which(names(cgunmatched2) %in% cgumd2)]
  Beta_value3[umd2] <- cgunmatched3[which(names(cgunmatched3) %in% cgumd2)]
  Beta_value4[pumd2] <- cgunmatched4[which(names(cgunmatched4) %in% cgpumd2)]
  Beta_value4[numd2] <- cgunmatched5[which(names(cgunmatched5) %in% cgnumd2)]
})
panmcused5 <- panmcused2[, c(1:3, 8:11)]
write.table(panmcused5, "pan_mammalian_clock2_hg19_convert_Denisovan.csv",
  row.names = FALSE, sep = ",")

pumd3 <- which(panmcused3$Merge == FALSE & panmcused3$Coefficient_3 > 0)
numd3 <- which(panmcused3$Merge == FALSE & panmcused3$Coefficient_3 < 0)
umd3 <- sort(c(pumd3, numd3))
cgpumd3 <- panmcused3$Name[pumd3]
cgnumd3 <- panmcused3$Name[numd3]
cgumd3 <- sort(c(cgpumd3, cgnumd3))
panmcused3 <- within(panmcused3, {
  Beta_value1 <- Beta_value2 <- Beta_value3 <- Beta_value4 <- Beta_value
  Beta_value1[umd3] <- cgunmatched1[which(names(cgunmatched1) %in% cgumd3)]
  Beta_value2[umd3] <- cgunmatched2[which(names(cgunmatched2) %in% cgumd3)]
  Beta_value3[umd3] <- cgunmatched3[which(names(cgunmatched3) %in% cgumd3)]
  Beta_value4[pumd3] <- cgunmatched4[which(names(cgunmatched4) %in% cgpumd3)]
  Beta_value4[numd3] <- cgunmatched5[which(names(cgunmatched5) %in% cgnumd3)]
})
panmcused6 <- panmcused3[, c(1:3, 8:11)]
write.table(panmcused6, "pan_mammalian_clock3_hg19_convert_Denisovan.csv",
  row.names = FALSE, sep = ",")
