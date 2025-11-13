#file requirements
#methylation map of the Altai Neandertal in a .txt file, obtained from "Reconstructing the DNA methylation maps of the Neandertal and the Denisovan" (https://doi.org/10.1126/science.1250368)
#target sites chromosome position information in a .csv file

datmet3 <- read.table("Recon_Meth_Altai_Neandertal.txt",
  header = FALSE, sep = "\t")
datmet3[, 1] <- paste("chr", datmet3[, 1], sep = "")

datasitp <- read.table("methylation_sites_relevant_position.csv",
  header = TRUE, sep = ",")

chrgroup <- paste("chr", 1:22, sep = "")
result <- data.frame("Name" = character(0), "chr" = character(0),
  "Position_start" = numeric(0), "Position_end" = numeric(0),
  "Beta_value" = numeric(0), "Merge" = logical(0))
for (chr in chrgroup) {
  datam <- datmet3[which(datmet3[, 1] == chr), ]
  datas <- datasitp[which(datasitp[, 2] == chr), ]
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
    dataa <- data.frame("Name" = cgsome, "chr" = chr,
      "Position_start" = posgroup, "Position_end" = posgroup + 1,
      "Beta_value" = 0, "Merge" = FALSE)
    datac <- rbind(datac, dataa)
  }

  result <- rbind(result, datac)
}

result[, 5] <- result[, 5] / 100
result[, 2] <- factor(result[, 2], levels = chrgroup)
result <- result[order(result[, 1]), ]
row.names(result) <- 1:nrow(result)

write.table(result, "methylation_wanted_sites_Altai_Neandertal_PC.csv",
  row.names = FALSE, sep = ",")