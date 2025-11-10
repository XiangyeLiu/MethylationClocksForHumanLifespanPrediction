#file requirements
#target sites in a .txt file
#all methylation values of samples in each series in a .csv file

library(stringr)

filegroup <- list.files(path = ".", "*.csv")
cggroup <- readLines("shared_methylation_sites.txt")
data4 <- data5 <- data.frame("ID_REF" = cggroup)
newcol <- 1

for (file in filegroup) {
  data3 <- read.table(file, header = TRUE, row.names = 1, sep = ",")
  data3 <- data3[which(row.names(data3) %in% cggroup), ]
  data3 <- data3[order(row.names(data3)), ]
  data3 <- data.frame(lapply(data3, as.numeric), row.names = row.names(data3))
  newmax <- apply(data3, 1, function(x) {
    x <- x[x >= 0 & x <= 1 & !is.na(x)]
    max(x)
  })
  newmin <- apply(data3, 1, function(x) {
    x <- x[x >= 0 & x <= 1 & !is.na(x)]
    min(x)
  })

  cgrow <- which(cggroup %in% row.names(data3))
  message(paste(file, length(cgrow), "Sites Matched"))
  if (length(cgrow) >= 0.9 * length(cggroup)) {
    newcol <- newcol + 1
    data4[cgrow, newcol] <- newmax
    data5[cgrow, newcol] <- newmin
    if (length(cgrow) != length(cggroup)) {
      data4[-cgrow, newcol] <- (max(colMeans(data3, na.rm = TRUE)) + mean(newmax)) / 2
      data5[-cgrow, newcol] <- (min(colMeans(data3, na.rm = TRUE)) + mean(newmin)) / 2
    }
    colnames(data4)[newcol] <- str_replace_all(file, ".csv", "_max")
    colnames(data5)[newcol] <- str_replace_all(file, ".csv", "_min")
  }

  data6 <- data.frame("ID_REF" = row.names(data3),
    "Max" = newmax, "Min" = newmin)
  scfile <- str_replace_all(file, ".csv", "_twisted.csv")
  write.table(data6, scfile, row.names = FALSE, col.names = FALSE, sep = ",")
}

write.table(data4, "extremum_methylation_level_max.csv", row.names = FALSE, sep = ",")
write.table(data5, "extremum_methylation_level_min.csv", row.names = FALSE, sep = ",")