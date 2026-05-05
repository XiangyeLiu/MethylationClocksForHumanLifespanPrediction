#file requirements
#information of target sites from 450K Beadchip to compute unmatched sites in a .RDS file, including individual sites and coefficients, obtained from "Universal DNA methylation age across mammalian tissues" (https://doi.org/10.1038/s43587-023-00462-6)
#information of target sites in six .csv files, with rows corresponding to individual sites and columns corresponding to chromosome position, coefficients, specific values, and whether directly matched with the methylation map or not, each clock including one file for the Altai Neandertal and the other for the Denisovan
#based on Array Converter Algorithm

library(tidyverse)

fit <- readRDS("fit_EPICtoMM_v2_450k.RDS")

panmcusen1 <- read.table("pan_mammalian_clock1_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

panmcusen2 <- read.table("pan_mammalian_clock2_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

panmcusen3 <- read.table("pan_mammalian_clock3_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

cgunmatched <- sort(unique(c(filter(panmcusen1, ! Merge)$Name,
  filter(panmcusen2, ! Merge)$Name, filter(panmcusen3, ! Merge)$Name)))

fitneedn <- fit[which(names(fit) %in% cgunmatched)]
saveRDS(fitneedn, "fit_EPICtoMM_v2_450k_filtered_Altai_Neandertal.RDS")

cgneedn <- lapply(seq_along(fitneedn), function(x) {
  cgneed <- fitneedn[[x]]$EPIC_CG
  return(cgneed[-1])
})
cgneedn <- sort(unique(unlist(cgneedn, use.names = FALSE)))
writeLines(cgneedn, "array_convert_need_Altai_Neandertal.txt")

panmcused1 <- read.table("pan_mammalian_clock1_hg19_Denisovan.csv",
  header = TRUE, sep = ",")

panmcused2 <- read.table("pan_mammalian_clock2_hg19_Denisovan.csv",
  header = TRUE, sep = ",")

panmcused3 <- read.table("pan_mammalian_clock3_hg19_Denisovan.csv",
  header = TRUE, sep = ",")

cgunmatched <- sort(unique(c(filter(panmcused1, ! Merge)$Name,
  filter(panmcused2, ! Merge)$Name, filter(panmcused3, ! Merge)$Name)))

fitneedd <- fit[which(names(fit) %in% cgunmatched)]
saveRDS(fitneedd, "fit_EPICtoMM_v2_450k_filtered_Denisovan.RDS")

cgneedd <- lapply(seq_along(fitneedd), function(x) {
  cgneed <- fitneedd[[x]]$EPIC_CG
  return(cgneed[-1])
})

cgneedd <- sort(unique(unlist(cgneedd, use.names = FALSE)))
writeLines(cgneedd, "array_convert_need_Denisovan.txt")
