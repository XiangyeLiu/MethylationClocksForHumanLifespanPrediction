#file requirements
#pan-mammalian clocks in three .csv files, one clock for each, with rows corresponding to individual sites and columns corresponding to chromosome position, coefficients and relevant gene information
#corresponding relationship between GRCh38 (hg38) and GRCh37 (hg19) in a .over.chain file or other formats

library(tidyverse)
library(rtracklayer)

panmcdata <- read.table("pan_mammalian_clock1.csv",
  header = TRUE, sep = ",")

intercept <- panmcdata[1, ]
colnames(intercept)[4] <- "Position_hg19"
panmcdata <- panmcdata[-1, ]
grhg38 <- GRanges(seqnames = panmcdata$chr,
  ranges = IRanges(start = panmcdata$Position_hg38,
    end = panmcdata$Position_hg38), strand = "*")
mcols(grhg38) <- panmcdata[, c(-3, -4)]

chain <- import.chain("hg38ToHg19.over.chain")
grhg19 <- liftOver(grhg38, chain)

success <- which(! sapply(grhg19, is.null) & lengths(grhg19) > 0)
if(length(success) == length(grhg38)) {
  message("All Sites Transformed")
} else if (length(grhg38) > length(success)) {
  failed <- panmcdata[which(! 1:length(grhg38) %in% success), ]
  failed[, -1:-3] <- NA
  colnames(failed)[4] <- "Position_hg19"
  msg <- paste(length(success), "Sites Transformed,",
    length(grhg38) - length(success), "Sites Left")
  message(msg)
} else {
  message("Wrong Happened")
}
grhg19 <- unlist(grhg19[success])

result <- data.frame(Name = mcols(grhg19)$Name,
  Coefficient_1 = mcols(grhg19)$Coefficient_1,
  chr = as.character(seqnames(grhg19)),
  Position_hg19 = start(grhg19),
  Gene = mcols(grhg19)$Gene,
  Gene_hg19 = mcols(grhg19)$Gene_hg19,
  ENTREZID = mcols(grhg19)$ENTREZID,
  Annotation = mcols(grhg19)$Annotation,
  Gene_region_ID = mcols(grhg19)$Gene_region_ID)

newpanmcdata <- bind_rows(intercept, result)

write.table(newpanmcdata, "pan_mammalian_clock1_hg19.csv",
  row.names = FALSE, sep = ",")

panmcdata <- read.table("pan_mammalian_clock2.csv",
  header = TRUE, sep = ",")

intercept <- panmcdata[1, ]
colnames(intercept)[4] <- "Position_hg19"
panmcdata <- panmcdata[-1, ]
grhg38 <- GRanges(seqnames = panmcdata$chr,
  ranges = IRanges(start = panmcdata$Position_hg38,
    end = panmcdata$Position_hg38), strand = "*")
mcols(grhg38) <- panmcdata[, c(-3, -4)]

chain <- import.chain("hg38ToHg19.over.chain")
grhg19 <- liftOver(grhg38, chain)

success <- which(! sapply(grhg19, is.null) & lengths(grhg19) > 0)
if(length(success) == length(grhg38)) {
  message("All Sites Transformed")
} else if (length(grhg38) > length(success)) {
  failed <- panmcdata[which(! 1:length(grhg38) %in% success), ]
  failed[, -1:-3] <- NA
  colnames(failed)[4] <- "Position_hg19"
  msg <- paste(length(success), "Sites Transformed,",
    length(grhg38) - length(success), "Sites Left")
  message(msg)
} else {
  message("Wrong Happened")
}
grhg19 <- unlist(grhg19[success])

result <- data.frame(Name = mcols(grhg19)$Name,
  Coefficient_2 = mcols(grhg19)$Coefficient_2,
  chr = as.character(seqnames(grhg19)),
  Position_hg19 = start(grhg19),
  Gene = mcols(grhg19)$Gene,
  Gene_hg19 = mcols(grhg19)$Gene_hg19,
  ENTREZID = mcols(grhg19)$ENTREZID,
  Annotation = mcols(grhg19)$Annotation,
  Gene_region_ID = mcols(grhg19)$Gene_region_ID)

newpanmcdata <- bind_rows(intercept, result)

write.table(newpanmcdata, "pan_mammalian_clock2_hg19.csv",
  row.names = FALSE, sep = ",")

panmcdata <- read.table("pan_mammalian_clock3.csv",
  header = TRUE, sep = ",")

intercept <- panmcdata[1, ]
colnames(intercept)[4] <- "Position_hg19"
panmcdata <- panmcdata[-1, ]
grhg38 <- GRanges(seqnames = panmcdata$chr,
  ranges = IRanges(start = panmcdata$Position_hg38,
    end = panmcdata$Position_hg38), strand = "*")
mcols(grhg38) <- panmcdata[, c(-3, -4)]

chain <- import.chain("hg38ToHg19.over.chain")
grhg19 <- liftOver(grhg38, chain)

success <- which(! sapply(grhg19, is.null) & lengths(grhg19) > 0)
if(length(success) == length(grhg38)) {
  message("All Sites Transformed")
} else if (length(grhg38) > length(success)) {
  failed <- panmcdata[which(! 1:length(grhg38) %in% success), ]
  failed[, -1:-3] <- NA
  colnames(failed)[4] <- "Position_hg19"
  msg <- paste(length(success), "Sites Transformed,",
    length(grhg38) - length(success), "Sites Left")
  message(msg)
} else {
  message("Wrong Happened")
}
grhg19 <- unlist(grhg19[success])

result <- data.frame(Name = mcols(grhg19)$Name,
  Coefficient_3 = mcols(grhg19)$Coefficient_3,
  chr = as.character(seqnames(grhg19)),
  Position_hg19 = start(grhg19),
  Gene = mcols(grhg19)$Gene,
  Gene_hg19 = mcols(grhg19)$Gene_hg19,
  ENTREZID = mcols(grhg19)$ENTREZID,
  Annotation = mcols(grhg19)$Annotation,
  Gene_region_ID = mcols(grhg19)$Gene_region_ID)

newpanmcdata <- bind_rows(intercept, result)

write.table(newpanmcdata, "pan_mammalian_clock3_hg19.csv",
  row.names = FALSE, sep = ",")