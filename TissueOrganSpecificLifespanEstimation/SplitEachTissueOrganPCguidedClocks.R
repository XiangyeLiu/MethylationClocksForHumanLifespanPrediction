#file requirements
#series-specific extremum of target sites among samples in two .csv files, with rows corresponding to individual sites and columns corresponding to series, including one file for the maximum and the other for the minimum


cgmaxgroup <- read.table("extremum_methylation_level_max.csv", header = TRUE,
  row.names = 1, sep = ",")
cgmingroup <- read.table("extremum_methylation_level_min.csv", header = TRUE,
  row.names = 1, sep = ",")

tissue <- as.vector(sapply(colnames(cgmaxgroup), function(x) {
  strsplit(x, "_")[[1]][3]
}))
tissue1 <- unique(tissue)

for(i in tissue1) {
  cols <- which(grepl(paste("_", i, "_", sep = ""), tissue))
  cgmax <- cgmaxgroup[, cols]
  cgmin <- cgmingroup[, cols]
  if (length(cols) == 1) {
    cgmax <- data.frame(cgmax, row.names = row.names(cgmaxgroup))
    colnames(cgmax) <- colnames(cgmaxgroup)[cols]
    cgmin <- data.frame(cgmin, row.names = row.names(cgmingroup))
    colnames(cgmin) <- colnames(cgmingroup)[cols]
  }
  res <- lapply(1:78464, function(x) {
    maxv <- max(unlist(cgmax[x, ], use.names = FALSE))
    minv <- min(unlist(cgmin[x, ], use.names = FALSE))
    return(c(maxv, minv))
  })
  res <- matrix(unlist(res, use.names = FALSE), ncol = 2, nrow = 78464, byrow = TRUE)
  resmax <- res[, 1]
  resmin <- res[, 2]
  data3 <- cbind(data3, resmax)
  data4 <- cbind(data4, resmin)
}

tissuename <- c("blood and blood cell", "blood", "leukocyte",
  "breast", "brain different parts", "chorionic villi", "liver",
  "buffy coat", "skin", "lung", "pancreas", "thyroid", "colon",
  "saliva", "semen", "vagina", "adipose tissue", "brain",
  "buccal and oral", "fibroblast", "kidney", "fallopian tube",
  "ovary", "neuron", "buccal")

for (i in 1:nrow(data3)) {
  data3[i, 1] <- max(data3[i, c(1, 2, 3, 8)])
  data3[i, 5] <- max(data3[i, c(5, 18)])
  data3[i, 19] <- max(data3[i, c(19, 25)])
}
colnames(data3) <- tissuename
data3 <- data3[, -c(2, 3, 8, 18, 25)]
data3 <- data3[, order(colnames(data3))]

for (i in 1:nrow(data4)) {
  data4[i, 1] <- min(data4[i, c(1, 2, 3, 8)])
  data4[i, 5] <- min(data4[i, c(5, 18)])
  data4[i, 19] <- min(data4[i, c(19, 25)])
}
colnames(data4) <- tissuename
data4 <- data4[, -c(2, 3, 8, 18, 25)]
data4 <- data4[, order(colnames(data4))]

write.table(data3, "extremum_split_max.csv", row.names = TRUE, sep = ",")
write.table(data4, "extremum_split_min.csv", row.names = TRUE, sep = ",")
