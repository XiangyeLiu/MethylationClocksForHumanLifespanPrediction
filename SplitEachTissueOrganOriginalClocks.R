#file requirements
#target methylation sites in a .txt file
#sample series in multiple .csv files, with rows corresponding to individual sites and columns corresponding to samples

cggroup <- readLines("wanted_sites.txt")
filegroup <- list.files(path = ".", ".csv")
tissue <- unique(as.vector(sapply(filegroup, function(x) {
  strsplit(x, "_")[[1]][3]
})))
datamax <- data.frame(matrix(nrow = length(cggroup), ncol = length(tissue) + 1))
datamax[, 1] <- cggroup
colnames(datamax) <- c("sites", tissue)
datamin <- datamax

for (i in seq_along(tissue)) {
  tissueone <- paste("_", tissue[i], "_", sep = "")
  filepart <- filegroup[which(grepl(tissueone, filegroup))]
  allmax <- data.frame("site" = cggroup, "max" = rep(-1, length(cggroup)))
  allmin <- data.frame("site" = cggroup, "min" = rep(2, length(cggroup)))
  for (file in filepart) {
    data <- read.table(file, header = FALSE, sep = ",")
    data <- na.omit(data)
    newmax <- data[, 2]
    newmin <- data[, 3]
    row <- which(cggroup %in% data[, 1])
    oldmax <- allmax[row, 2]
    oldmin <- allmin[row, 2]
    allmax[row, 2] <- ifelse(newmax >= oldmax & newmax <= 1, newmax, oldmax)
    allmin[row, 2] <- ifelse(newmin <= oldmin & newmin >= 0, newmin, oldmin)
  }
  datamax[, i + 1] <- allmax[, 2]
  datamin[, i + 1] <- allmin[, 2]
}

for (i in seq_along(cggroup)) {
  datamax[i, 2] <- max(datamax[i, c(2, 5, 7, 14)])
  datamax[i, 3] <- max(datamax[i, c(3, 19, 20)])
  datamax[i, 6] <- max(datamax[i, c(6, 24)])
  datamax[i, 8] <- max(datamax[i, c(8, 31)])
  datamax[i, 9] <- max(datamax[i, c(9, 21)])

  datamin[i, 2] <- min(datamin[i, c(2, 5, 7, 14)])
  datamin[i, 3] <- min(datamin[i, c(3, 19, 20)])
  datamin[i, 6] <- min(datamin[i, c(6, 24)])
  datamin[i, 8] <- min(datamin[i, c(8, 31)])
  datamin[i, 9] <- min(datamin[i, c(9, 21)])
}

tissuename <- c("blood and blood cells", "brain different parts",
  "mesenchymal stromal cell", "blood", "skin or dermal", "leukocyte", 
  "buccal or oral", "heart", "kidney", "liver", "prostate",
  "leukocyte lymphoblast", "buffy coat", "saliva", "gastric", "uterine cervix",
  "breast", "cerebellum", "frontal cortex", "left ventricular myocardium",
  "chorionic villi", "pancreas", "skin", "lung", "thyroid", "colon",
  "semen", "vagina", "adipose tissue", "oral cavity", "fibroblast",
  "fallopian tube", "ovary", "airway", "neuron")
colnames(datamax) <- colnames(datamin) <- c("site", tissuename)
datamax <- datamax[, -c(5, 7, 14, 19, 20, 21, 24, 31)]
datamin <- datamin[, -c(5, 7, 14, 19, 20, 21, 24, 31)]

datamax[, -1] <- apply(datamax[, -1], 2, function(x) {
  x[which(x == -1)] <- NA
  return(x)
})
datamin[, -1] <- apply(datamin[, -1], 2, function(x) {
  x[which(x == 2)] <- NA
  return(x)
})

write.table(datamax, "extremum_split_max.csv", row.names = FALSE, sep = ",")
write.table(datamin, "extremum_split_min.csv", row.names = FALSE, sep = ",")