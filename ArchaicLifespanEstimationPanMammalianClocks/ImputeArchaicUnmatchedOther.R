#file requirements
#information of target sites in six .csv files, with rows corresponding to individual sites and columns corresponding to chromosome position, specific values, and whether directly matched with the methylation map or not, each clock including one file for the Altai Neandertal and the other for the Denisovan

library(tidyverse)

imputation1 <- function(data) {
  unmatched <- which(is.na(data$Beta_value))[-1]
  valueclean <- data$Beta_value[c(-1, -1 * unmatched)]
  data <- within(data[, 1:4], {
    Mean <- Median <- Beta_value
    Mean[unmatched] <- mean(valueclean)
    Median[unmatched] <- median(valueclean)
  })
  data <- data[, -4]
  return(data)
}

imputation2 <- function(data, kmax = 100) {
  unmatched <- which(is.na(data$Beta_value))[-1]
  value <- data$Beta_value
  datanew <- data.frame(matrix(nrow = nrow(data), ncol = kmax))
  for (k in 1:kmax) {
    valuechange <- NULL
    for (i in unmatched) {
      datafilt <- filter(data, chr == data$chr[i] & Merge)
      valuefilt <- datafilt$Beta_value
      range <- abs(datafilt$Position_start - data$Position_start[i])
      valuefilt <- valuefilt[order(range)]
      knew <- ifelse(k <= length(valuefilt), k, length(valuefilt))
      valuechange <- c(valuechange, mean(valuefilt[1:knew]))
    }
    value[unmatched] <- valuechange
    datanew[, k] <- value
    colnames(datanew)[k] <- paste("k", k, sep = "_")
  }
  data <- cbind(data[, 1:3], datanew)
  return(data)
}

imputation3 <- function(data, time = 10000) {
  unmatched <- which(is.na(data$Beta_value))[-1]
  value <- data$Beta_value
  valueclean <- data$Beta_value[c(-1, -1 * unmatched)]
  datanew <- data.frame(matrix(nrow = nrow(data), ncol = 2 * time))
  set.seed(123)
  for (i in 1:time) {
    value[unmatched] <- sample(valueclean, length(unmatched),
      replace = FALSE)
    datanew[, i] <- value
    colnames(datanew)[i] <- paste("sample", i, sep = "_")
    value[unmatched] <- runif(length(unmatched), min(valueclean),
      max(valueclean))
    datanew[, time + i] <- value
    colnames(datanew)[time + i] <- paste("runif", i, sep = "_")
  }
  data <- cbind(data[, 1:3], datanew)
  return(data)
}

imputation4 <- function(data) {
  unmatched <- which(is.na(data$Beta_value))[-1]
  positive <- unmatched[which(data[unmatched, 3] > 0)]
  negative <- unmatched[which(data[unmatched, 3] < 0)]
  data <- within(data[, 1:4], {
    Beta_value[positive] <- 1
    Beta_value[negative] <- 0
  })
  colnames(data)[4] <- "Max"
  return(data)
}

panmcusen1 <- read.table("pan_mammalian_clock1_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")
data1 <- imputation1(panmcusen1)
data2 <- imputation2(panmcusen1)
data3 <- imputation3(panmcusen1)
data4 <- imputation4(panmcusen1)
data <- cbind(data1, data2[, -1:-3], data3[, -1:-3], data4[, -1:-3])
colnames(data)[ncol(data)] <- "Max"
write.table(data, "pan_mammalian_clock1_hg19_imputation_Altai_Neandertal.csv",
  row.names = FALSE, sep = ",")

panmcusen2 <- read.table("pan_mammalian_clock2_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")
data1 <- imputation1(panmcusen2)
data2 <- imputation2(panmcusen2)
data3 <- imputation3(panmcusen2)
data4 <- imputation4(panmcusen2)
data <- cbind(data1, data2[, -1:-3], data3[, -1:-3], data4[, -1:-3])
colnames(data)[ncol(data)] <- "Max"
write.table(data, "pan_mammalian_clock2_hg19_imputation_Altai_Neandertal.csv",
  row.names = FALSE, sep = ",")

panmcusen3 <- read.table("pan_mammalian_clock3_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")
data1 <- imputation1(panmcusen3)
data2 <- imputation2(panmcusen3)
data3 <- imputation3(panmcusen3)
data4 <- imputation4(panmcusen3)
data <- cbind(data1, data2[, -1:-3], data3[, -1:-3], data4[, -1:-3])
colnames(data)[ncol(data)] <- "Max"
write.table(data, "pan_mammalian_clock3_hg19_imputation_Altai_Neandertal.csv",
  row.names = FALSE, sep = ",")

panmcused1 <- read.table("pan_mammalian_clock1_hg19_Denisovan.csv",
  header = TRUE, sep = ",")
data1 <- imputation1(panmcused1)
data2 <- imputation2(panmcused1)
data3 <- imputation3(panmcused1)
data4 <- imputation4(panmcused1)
data <- cbind(data1, data2[, -1:-3], data3[, -1:-3], data4[, -1:-3])
colnames(data)[ncol(data)] <- "Max"
write.table(data, "pan_mammalian_clock1_hg19_imputation_Denisovan.csv",
  row.names = FALSE, sep = ",")

panmcused2 <- read.table("pan_mammalian_clock2_hg19_Denisovan.csv",
  header = TRUE, sep = ",")
data1 <- imputation1(panmcused2)
data2 <- imputation2(panmcused2)
data3 <- imputation3(panmcused2)
data4 <- imputation4(panmcused2)
data <- cbind(data1, data2[, -1:-3], data3[, -1:-3], data4[, -1:-3])
colnames(data)[ncol(data)] <- "Max"
write.table(data, "pan_mammalian_clock2_hg19_imputation_Denisovan.csv",
  row.names = FALSE, sep = ",")

panmcused3 <- read.table("pan_mammalian_clock3_hg19_Denisovan.csv",
  header = TRUE, sep = ",")
data1 <- imputation1(panmcused3)
data2 <- imputation2(panmcused3)
data3 <- imputation3(panmcused3)
data4 <- imputation4(panmcused3)
data <- cbind(data1, data2[, -1:-3], data3[, -1:-3], data4[, -1:-3])
colnames(data)[ncol(data)] <- "Max"
write.table(data, "pan_mammalian_clock3_hg19_imputation_Denisovan.csv",
  row.names = FALSE, sep = ",")