#file requirements
#target sites in a .txt file
#extremum of each series  in a .csv file

filegroup <- list.files(path = ".", "*_twisted.csv")
cggroup <- readLines("shared_methylation_sites.txt")
data4 <- data.frame("ID_REF" = cggroup, "Max" = 0, "Min" = 1)

for (file in filegroup) {
  data3 <- read.table(file, header = FALSE, row.names = 1, sep = ",")
  remove <- c(which(data3[, 1] > 1), which(data3[, 2] < 0))
  if (length(remove)) {
    data3 <- data3[-remove, ]
  }
  if (all(row.names(data3) %in% cggroup)) {
    message(paste("All Settled", length(row.names(data3)), "Matched"))
    row <- which(data4[, 1] %in% row.names(data3))
    orimax <- data4[row, 2]
    orimin <- data4[row, 3]
    data4[row, 2] <- ifelse(orimax > data3[, 1], orimax, data3[, 1])
    data4[row, 3] <- ifelse(orimin < data3[, 2], orimin, data3[, 2])
  } else {
    stop("Wrong Happened")
  }
}

write.table(data4, "extremum_methylation_level.csv", row.names = FALSE, sep = ",")