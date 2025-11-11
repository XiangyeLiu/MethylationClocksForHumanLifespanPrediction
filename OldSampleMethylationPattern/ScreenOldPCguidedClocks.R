#file requirements

#target methylation sites in a .txt file
#sample series in multiple .csv files, with rows corresponding to individual sites and columns corresponding to samples (annotated with ages)

library(stringr)

cggroup <- readLines("shared_methylation_sites.txt")
data4 <- data.frame(matrix(nrow = 0, ncol = 3))

filegroup <- list.files(path = ".", "*.csv")

for (i in seq_along(filegroup)) {
  if (grepl(".csv", filegroup[i])) {
    data3 <- read.table(filegroup[i], header = TRUE, sep = ",")
    filename <- str_replace_all(filegroup[i], ".csv", "")
  } else if (grepl(".txt", filegroup[i])) {
    data3 <- read.table(filegroup[i], header = TRUE, sep = "\t")
    filename <- str_replace_all(filegroup[i], ".txt", "")
  }
  if (sum(grepl("yearsold", colnames(data3)))) {
    colnames(data3) <- str_replace_all(colnames(data3), "yearsold", "_")
  }
  remove <- which(grepl("_NA_", colnames(data3)))
  if (length(remove)) {
    data3 <- data3[, -remove]
  }
  colselect <- which(grepl("_\\d+\\.?\\d*_", colnames(data3)))
  samplegroup <- colnames(data3)[colselect]
  if (length(samplegroup)) {
    agegroup <- NULL
    for (sample in samplegroup) {
      age <- regmatches(sample, gregexpr("_\\d+\\.?\\d*_", sample))[[1]]
      age <- as.numeric(str_replace_all(age, "_", ""))
      agegroup <- c(agegroup, age)
    }
    datain <- data.frame("Age" = agegroup,
      "Sample" = samplegroup, "Name" = filename)
    data4 <- rbind(data4, datain)
    message(paste("Added", nrow(datain), "Samples", filegroup[i]))
  } else {
    message(paste("Added 0 Sample", filegroup[i]))
  }
}

write.table(data4, "age_display.csv", header = TRUE,
  row.names = FALSE, sep = ",")

data5 <- data4[order(data4[, 1], decreasing = TRUE), ]

data6 <- data5[which(data5[, 1] > 90), ]
data7 <- data.frame(matrix(ncol = 0, nrow = length(cggroup)))
for (filename in unique(data6[, 3])) {
  agegroup <- data6[which(data6[, 3] == filename), 1]
  colgroup <- data6[which(data6[, 3] == filename), 2]
  datain <- read.table(paste(filename, ".csv", sep = ""), header = TRUE,
    row.names = 1, sep = ",")
  if (sum(grepl("yearsold", colnames(datain)))) {
    colnames(datain) <- str_replace_all(colnames(datain), "yearsold", "_")
  }
  datain <- datain[which(row.names(datain) %in% cggroup), ]
  namerow <- row.names(datain)
  datain <- datain[, which(colnames(datain) %in% colgroup)]
  if (is.vector(datain)) {
    datain <- as.data.frame(datain, row.names = namerow)
    colnames(datain)[1] <- colgroup
  }
  datain <- data.frame(lapply(datain, as.numeric),
    row.names = row.names(datain))

  cgadd <- cggroup[! cggroup %in% row.names(datain)]
  if (length(cgadd)) {
    add <- data.frame(matrix(rep(colMeans(datain), length(cgadd)),
      nrow = length(cgadd), byrow = TRUE), row.names = cgadd)
    colnames(add) <- colnames(datain)
    datain <- rbind(datain, add)
  }
  datain <- datain[order(row.names(datain)), ]
  data7 <- cbind(data7, datain)
}

write.table(data7, "age_old.csv", row.names = TRUE, sep = ",")