#file requirements
#methylation map of the Denisovan in a .txt file, obtained from "Reconstructing the DNA methylation maps of the Neandertal and the Denisovan" (https://doi.org/10.1126/science.1250368)
#pan-mammalian clocks in three .csv files, one clock for each, with rows corresponding to individual sites and columns corresponding to transformed chromosome position, coefficients and relevant gene information

library(tidyverse)

metmapd <- read.table("Recon_Meth_Denisovan.txt",
  header = FALSE, sep = "\t")
metmapd[, 1] <- paste("chr", metmapd[, 1], sep = "")
colnames(metmapd) <- c("chr", "Position_start",
  "Position_end", "Beta_value_percentage")

panmcdata1 <- read.table("pan_mammalian_clock1_hg19.csv",
  header = TRUE, sep = ",")
intercept <- panmcdata1[1, c(1, 3, 2)]

result1 <- NULL
chrgroup <- unique(panmcdata1[-1, 3])
for (chrsg in chrgroup) {
  datam <- filter(metmapd, chr == chrsg)
  datas <- filter(panmcdata1[, 1:4], chr == chrsg)
  datac <- merge(datam, datas,
    by.x = "Position_start", by.y = "Position_hg19")
  datac <- datac[, c(5, 2, 6, 4, 1, 3)]
  datac <- tibble(datac, TRUE)
  colnames(datac) <- c("Name", "chr", "Coefficient_1", "Beta_value",
    "Position_start", "Position_end",  "Merge")
  message(paste("In ", chrsg, ": ",
    nrow(datac), "/", nrow(datas), " Are Collected", sep = ""))

  failed <- which(! datas$Name %in% datac$Name)
  if (length(failed)) {
    name <- datas[failed, 1]
    coef <- datas[failed, 2]
    position <- datas[failed, 4]
    dataa <- tibble(name, chrsg, coef, NA,
      position, position + 1, FALSE)
    colnames(dataa) <- colnames(datac)
    datac <- bind_rows(datac, dataa)
  }
  
  result1 <- bind_rows(result1, datac)
}

result1$Beta_value <- result1$Beta_value / 100
result1$chr <- factor(result1$chr,
  levels = paste("chr", c(1:22, "X"), sep = ""))
result1 <- result1[order(result1$Name), ]
result1 <- bind_rows(intercept, result1)

write.table(result1, "pan_mammalian_clock1_hg19_Denisovan.csv",
  row.names = FALSE, sep = ",")

panmcdata2 <- read.table("pan_mammalian_clock2_hg19.csv",
  header = TRUE, sep = ",")
intercept <- panmcdata2[1, c(1, 3, 2)]

result2 <- NULL
chrgroup <- unique(panmcdata2[-1, 3])
for (chrsg in chrgroup) {
  datam <- filter(metmapd, chr == chrsg)
  datas <- filter(panmcdata2[, 1:4], chr == chrsg)
  datac <- merge(datam, datas,
    by.x = "Position_start", by.y = "Position_hg19")
  datac <- datac[, c(5, 2, 6, 4, 1, 3)]
  datac <- tibble(datac, TRUE)
  colnames(datac) <- c("Name", "chr", "Coefficient_2", "Beta_value",
    "Position_start", "Position_end",  "Merge")
  message(paste("In ", chrsg, ": ",
    nrow(datac), "/", nrow(datas), " Are Collected", sep = ""))

  failed <- which(! datas$Name %in% datac$Name)
  if (length(failed)) {
    name <- datas[failed, 1]
    coef <- datas[failed, 2]
    position <- datas[failed, 4]
    dataa <- tibble(name, chrsg, coef, NA,
      position, position + 1, FALSE)
    colnames(dataa) <- colnames(datac)
    datac <- bind_rows(datac, dataa)
  }
  
  result2 <- bind_rows(result2, datac)
}

result2$Beta_value <- result2$Beta_value / 100
result2$chr <- factor(result2$chr,
  levels = paste("chr", c(1:22, "X"), sep = ""))
result2 <- result2[order(result2$Name), ]
result2 <- bind_rows(intercept, result2)

write.table(result2, "pan_mammalian_clock2_hg19_Denisovan.csv",
  row.names = FALSE, sep = ",")

panmcdata3 <- read.table("pan_mammalian_clock3_hg19.csv",
  header = TRUE, sep = ",")
intercept <- panmcdata3[1, c(1, 3, 2)]

result3 <- NULL
chrgroup <- unique(panmcdata3[-1, 3])
for (chrsg in chrgroup) {
  datam <- filter(metmapd, chr == chrsg)
  datas <- filter(panmcdata3[, 1:4], chr == chrsg)
  datac <- merge(datam, datas,
    by.x = "Position_start", by.y = "Position_hg19")
  datac <- datac[, c(5, 2, 6, 4, 1, 3)]
  datac <- tibble(datac, TRUE)
  colnames(datac) <- c("Name", "chr", "Coefficient_3", "Beta_value",
    "Position_start", "Position_end",  "Merge")
  message(paste("In ", chrsg, ": ",
    nrow(datac), "/", nrow(datas), " Are Collected", sep = ""))

  failed <- which(! datas$Name %in% datac$Name)
  if (length(failed)) {
    name <- datas[failed, 1]
    coef <- datas[failed, 2]
    position <- datas[failed, 4]
    dataa <- tibble(name, chrsg, coef, NA,
      position, position + 1, FALSE)
    colnames(dataa) <- colnames(datac)
    datac <- bind_rows(datac, dataa)
  }
  
  result3 <- bind_rows(result3, datac)
}

result3$Beta_value <- result3$Beta_value / 100
result3$chr <- factor(result3$chr,
  levels = paste("chr", c(1:22, "X"), sep = ""))
result3 <- result3[order(result3$Name), ]
result3 <- bind_rows(intercept, result3)

write.table(result3, "pan_mammalian_clock3_hg19_Denisovan.csv",
  row.names = FALSE, sep = ",")