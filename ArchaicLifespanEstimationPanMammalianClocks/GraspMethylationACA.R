#file requirements
#methylation map of the Altai Neandertal in a .txt file, obtained from "Reconstructing the DNA methylation maps of the Neandertal and the Denisovan" (https://doi.org/10.1126/science.1250368)
#target sites from 450K Beadchip to compute unmatched sites and chromosome position information in two .csv files, including one file for the Altai Neandertal and the other for the Denisovan
#based on Array Converter Algorithm

library(tidyverse)

metmapn <- read.table("Recon_Meth_Altai_Neandertal.txt",
  header = FALSE, sep = "\t")
metmapn[, 1] <- paste("chr", metmapn[, 1], sep = "")
colnames(metmapn) <- c("chr", "Position_start",
  "Position_end", "Beta_value_percentage")

metmapd <- read.table("Recon_Meth_Denisovan.txt",
  header = FALSE, sep = "\t")
metmapd[, 1] <- paste("chr", metmapd[, 1], sep = "")
colnames(metmapd) <- c("chr", "Position_start",
  "Position_end", "Beta_value_percentage")

cgneedmn <- read.table("array_convert_annotation_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

cgneedmd <- read.table("array_convert_annotation_Denisovan.csv",
  header = TRUE, sep = ",")

mcgn <- NULL
chrgroup <- sort(unique(cgneedmn$chr))
for (chrsg in chrgroup) {
  datam <- filter(metmapn, chr == chrsg)
  datas <- filter(cgneedmn, chr == chrsg)
  datac <- merge(datam, datas,
    by.x = "Position_start", by.y = "Position")
  datac <- datac[, c(5, 2, 4, 1, 3)]
  datac <- tibble(datac, TRUE)
  datac$Beta_value_percentage <- datac$Beta_value_percentage / 100
  colnames(datac)[c(2, 3, 6)] <- c("chr", "Beta_value", "Merge")

  failed <- which(! datas$Name %in% datac$Name)
  if (length(failed)) {
    name <- datas$Name[failed]
    position <- datas$Position[failed]
    dataa <- tibble(name, chrsg, NA, position, position + 1, FALSE)
    colnames(dataa) <- colnames(datac)
    datac <- bind_rows(datac, dataa)
  }
  
  mcgn <- bind_rows(mcgn, datac)
}
mcgn$chr <- factor(mcgn$chr, levels = paste("chr", c(1:22, "X"), sep = ""))
mcgn <- mcgn[order(mcgn$Name), ]

mcgd <- NULL
chrgroup <- sort(unique(cgneedmd$chr))
for (chrsg in chrgroup) {
  datam <- filter(metmapd, chr == chrsg)
  datas <- filter(cgneedmd, chr == chrsg)
  datac <- merge(datam, datas,
    by.x = "Position_start", by.y = "Position")
  datac <- datac[, c(5, 2, 4, 1, 3)]
  datac <- tibble(datac, TRUE)
  datac$Beta_value_percentage <- datac$Beta_value_percentage / 100
  colnames(datac)[c(2, 3, 6)] <- c("chr", "Beta_value", "Merge")

  failed <- which(! datas$Name %in% datac$Name)
  if (length(failed)) {
    name <- datas$Name[failed]
    position <- datas$Position[failed]
    dataa <- tibble(name, chrsg, NA, position, position + 1, FALSE)
    colnames(dataa) <- colnames(datac)
    datac <- bind_rows(datac, dataa)
  }

  mcgd <- bind_rows(mcgd, datac)
}
mcgd$chr <- factor(mcgd$chr, levels = paste("chr", c(1:22, "X"), sep = ""))
mcgd <- mcgd[order(mcgd$Name), ]

write.table(mcgn, "array_convert_value_Altai_Neandertal.csv",
  row.names = FALSE, sep = ",")

write.table(mcgd, "array_convert_value_Denisovan.csv",
  row.names = FALSE, sep = ",")
