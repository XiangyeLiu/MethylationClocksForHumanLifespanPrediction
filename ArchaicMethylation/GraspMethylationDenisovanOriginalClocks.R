#file requirements
#methylation map of the Denisovan in a .txt file, obtained from "Reconstructing the DNA methylation maps of the Neandertal and the Denisovan" (https://doi.org/10.1126/science.1250368)
#target sites chromosome position information in a .csv file

datmet3 <- read.table("Recon_Meth_Denisovan.txt",
  header = FALSE, sep = "\t")
datmet3[, 1] <- paste("chr", datmet3[, 1], sep = "")

datasite <- read.table("methylation_sites_relevant_information.csv",
  header = TRUE, sep = ",")
datasit3 <- datasite[, 1:3]

chrgroup <- paste("chr", c(1:22, "X"), sep = "")
result <- data.frame(matrix(ncol = 6, nrow = 0))
for (chr in chrgroup) {
  datam <- datmet3[which(datmet3[, 1] == chr), ]
  datas <- datasit3[which(datasit3[, 2] == chr), ]
  datac <- merge(datas, datam[, -1], by.x = "pos", by.y = "V2")
  datac <- datac[, c(2, 3, 1, 4, 5)]
  datac <- cbind(datac, TRUE)
  colnames(datac) <- c("Name", "chr", "Position_start",
    "Position_end", "Beta_value", "Merge")
  message(paste(nrow(datac), "in", nrow(datas), "collected", chr))

  row <- which(! datas[, 1] %in% datac[, 1])
  if (length(row)) {
    cgsome <- datas[row, 1]
    posgroup <- datas[row, 3]
    dataa <- data.frame(cbind(cgsome, chr, posgroup, posgroup + 1, 0, FALSE))
    colnames(dataa) <- c("Name", "chr", "Position_start",
      "Position_end", "Beta_value", "Merge")
    datac <- rbind(datac, dataa)
  }

  result <- rbind(result, datac)
}

result[, 5] <- as.numeric(result[, 5]) / 100
result[, 2] <- factor(result[, 2], levels = chrgroup)
result <- result[order(result[, 1]), ]
row.names(result) <- 1:nrow(result)

write.table(result, "methylation_wanted_sites_Denisovan.csv",
  row.names = FALSE, sep = ",")