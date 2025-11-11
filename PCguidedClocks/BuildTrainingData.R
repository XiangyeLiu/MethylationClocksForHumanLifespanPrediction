#file requirements

#target methylation sites in a .txt file
#training series in multiple .csv files, with rows corresponding to individual sites and columns corresponding to samples (annotated with ages)

library(stringr)

cggroup <- readLines("shared_methylation_sites.txt")
data4 <- data.frame(matrix(ncol = 0, nrow = length(cggroup)))

filegroup <- list.files(path = ".", "*.csv")
for (file in filegroup) {
  data3 <- read.table(file, header = TRUE, row.names = 1, sep = ",")
  if (sum(grepl("yearsold", colnames(data3)))) {
    colnames(data3) <- str_replace_all(colnames(data3), "yearsold", "_")
  }
  remove <- which(grepl("_NA_", colnames(data3)))
  if (length(remove)) {
    data3 <- data3[, -remove]
  }
  data3 <- data3[which(row.names(data3) %in% cggroup), ]
  cgadd <- cggroup[! cggroup %in% row.names(data3)]
  if (length(cgadd)) {
    add <- data.frame(matrix(rep(colMeans(data3), length(cgadd)),
      nrow = length(cgadd), byrow = TRUE), row.names = cgadd)
    colnames(add) <- colnames(data3)
    data3 <- rbind(data3, add)
  }
  data3 <- data3[order(row.names(data3)), ]
  data4 <- cbind(data4, data3)

  samplegroup <- colnames(data3)
  agegroup <- NULL
  for (sample in samplegroup) {
    age <- regmatches(sample, gregexpr("_\\d+\\.?\\d*_", sample))[[1]]
    age <- as.numeric(str_replace_all(age, "_", ""))
    agegroup <- c(agegroup, age)
  }
  data5 <- as.data.frame(agegroup, row.names = samplegroup)
  colnames(data5) <- "Age"
  if (!file.exists("PC_clock_train_pheno.csv")) {
    write.table(data5, "PC_clock_train_pheno.csv", row.names = TRUE, sep = ",")
  } else {
    write.table(data5, "PC_clock_train_pheno.csv", row.names = TRUE, sep = ",",
      col.names = FALSE, append = TRUE)
  }
}

write.table(data4, "PC_clock_train_meth.csv", row.names = TRUE, sep = ",")
