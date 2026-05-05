#file requirements
#information of target sites in four .csv files, with rows corresponding to individual sites and columns corresponding to chromosome position, specific values, and whether directly matched with the methylation map or not, each clock type including one file for the Altai Neandertal and the other for the Denisovan

datamn <- read.table("methylation_wanted_sites_Altai_Neandertal.csv",
  header = TRUE, sep = ",")
datamd <- read.table("methylation_wanted_sites_Denisovan.csv",
  header = TRUE, sep = ",")

datamnpc <- read.table("methylation_wanted_sites_Altai_Neandertal_PC.csv",
  header = TRUE, sep = ",")
datamdpc <- read.table("methylation_wanted_sites_Denisovan_PC.csv",
  header = TRUE, sep = ",")

imputation1 <- function(data) {
  bv <- data$Beta_value
  meanv <- mean(bv)
  medianv <- median(bv)
  unmatched <- which(data$Merge == FALSE)
  imput1 <- imput2 <- bv
  imput1[unmatched] <- meanv
  imput2[unmatched] <- medianv
  datai <- cbind(data[, 1:4], imput1, imput2)
  colnames(datai)[5:6] <- paste(rep("Beta_value", 2), 1:2, sep = "_")
  return(datai)
}

datamni <- imputation1(datamn)
datamdi <- imputation1(datamd)
datamnpci <- imputation1(datamnpc)
datamdpci <- imputation1(datamdpc)

write.table(datamni, "imputation_pmm_Altai_Neandertal.csv",
  row.names = FALSE, sep = ",")
write.table(datamdi, "imputation_pmm_Denisovan.csv",
  row.names = FALSE, sep = ",")
write.table(datamnpci, "imputation_pmm_Altai_Neandertal_PC.csv",
  row.names = FALSE, sep = ",")
write.table(datamdpci, "imputation_pmm_Denisovan_PC.csv",
  row.names = FALSE, sep = ",")

imputation2 <- function(data, kmax) {
  unmatched <- which(data$Merge == FALSE)
  bv <- data$Beta_value
  set <- data.frame(matrix(ncol = kmax, nrow = nrow(data)))
  colnames(set) <- paste("k", 1:kmax, sep = "=")
  for (k in 1:kmax) {
    bvi <- bv
    for (i in unmatched) {
      chr <- data$chr[i]
      poss <- data$Position_start[i]
      candidate <- which(data$chr == chr & data$Merge)
      candiposs <- data$Position_start[candidate]
      dis <- abs(candiposs - poss)
      if (k < length(dis)) {
        near <- order(dis)[1:k]
      } else {
        near <- order(dis)
      }
      nn <- candidate[near]
      nndis <- dis[near]
      nnbv <- data$Beta_value[nn]
      nnweight <- 1 / nndis
      nnsweight <- nnweight / sum(nnweight)
      vi <- sum(nnsweight * nnbv)
      bvi[i] <- vi
    }
    set[, k] <- bvi
  }
  datai <- cbind(data[, 1:4], set)
  return(datai)
}

datamni <- imputation2(datamn, 100)
datamdi <- imputation2(datamd, 100)
datamnpci <- imputation2(datamnpc, 100)
datamdpci <- imputation2(datamdpc, 100)

write.table(datamni, "imputation_knn_Altai_Neandertal.csv",
  row.names = FALSE, sep = ",")
write.table(datamdi, "imputation_knn_Denisovan.csv",
  row.names = FALSE, sep = ",")
write.table(datamnpci, "imputation_knn_Altai_Neandertal_PC.csv",
  row.names = FALSE, sep = ",")
write.table(datamdpci, "imputation_knn_Denisovan_PC.csv",
  row.names = FALSE, sep = ",")

imputation3 <- function(data, time) {
  set.seed(123)
  unmatched <- which(data$Merge == FALSE)
  bv <- data$Beta_value[-unmatched]
  set1 <- set2 <- data.frame(matrix(nrow = length(unmatched) , ncol = time))
  for (i in 1:time) {
    set1[, i] <- sample(bv, length(unmatched), replace = FALSE)
    set2[, i] <- runif(unmatched, min = min(bv), max = max(bv))
  }
  datai <- cbind(data$Name[unmatched], unmatched, set1, set2)
  colnames(datai) <- c("Name", "Number", paste(rep(c("sample", "runif"), each = time), 1:time, sep = "_"))
  return(datai)
}

datamni <- imputation3(datamn, 10000)
datamdi <- imputation3(datamd, 10000)
datamnpci <- imputation(datamnpc, 10000)
datamdpci <- imputation(datamdpc, 10000)

write.table(datamni, "imputation_sva_Altai_Neandertal.csv",
  row.names = FALSE, sep = ",")
write.table(datamdi, "imputation_sva_Denisovan.csv",
  row.names = FALSE, sep = ",")
write.table(datamnpci, "imputation_sva_Altai_Neandertal_PC.csv",
  row.names = FALSE, sep = ",")
write.table(datamdpci, "imputation_sva_Denisovan_PC.csv",
  row.names = FALSE, sep = ",")
