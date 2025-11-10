#file requirements
#cluster and outlier information of target methylation sites in two .RDS files
#methylation values of target sites for each sample (age >= 90) in a .csv file, with rows corresponding to individual sites (the last row corresponding to ages) and columns corresponding to samples

library(tidyverse)

dataold <- read.table("age_old_enhanced_V2.csv", header = TRUE, sep = ",")
dataold <- dataold[-nrow(dataold), ]

cglist <- readRDS("detailed_methylation_level_enhanced_PC.RDS")
boundarylist <- readRDS("boundary_methylation_level_enhanced_PC.RDS")

cggroup <- dataold[, 1]

cgrelation <- data.frame("site" = cggroup, "reach_outlier" = NA,
  "reach_extremum" = NA, "all_outlier" = NA)

cgrelation[, 2] <- as.vector(sapply(seq_along(cggroup), function(x) {
  cgv <- unlist(dataold[x, -1], use.names = FALSE)
  boundary <- boundarylist[[x]]
  len <- length(boundary)
  if (min(cgv) <= boundary[1]) {
    reach <- "left"
  } else if (max(cgv) >= boundary[len]) {
    reach <- "right"
  } else if (min(cgv) <= boundary[1] && max(cgv) >= boundray[len]) {
    reach <- "both"
  } else {
    reach  <- "no"
  }
  return(reach)
}))

cgrelation[, 3] <- as.vector(sapply(seq_along(cggroup), function(x) {
  cgv <- unlist(dataold[x, -1], use.names = FALSE)
  cgsome <- cglist[[x]]
  if (min(cgv) <= min(cgsome)) {
    reach <- "min"
  } else if (max(cgv) >= max(cgsome)) {
    reach <- "max"
  } else if (min(cgv) <= min(cgsome) && max(cgv) >= max(cgsome)) {
    reach <- "both"
  } else {
    reach  <- "no"
  }
  return(reach)
}))

cgrelation[, 4] <- as.vector(sapply(seq_along(cggroup), function(x) {
  cgv <- unlist(dataold[x, -1], use.names = FALSE)
  ro <- cgrelation[x, 2]
  boundary <- boundarylist[[x]]
  len <- length(boundary)
  if (ro == "left" && max(cgv) <= boundary[1]) {
    reach <- "leftall"
  } else if (ro == "right" && min(cgv) >= boundary[len]) {
    reach <- "rightall"
  } else {
    reach  <- "no"
  }
  return(reach)
}))

bb <- data.frame(value = c("No Outlier / Extremum", "Outlier / Extremum 1 Only",
  "Outlier / Extremum 2 Only", "Both Outlier / Extremum"))
freq1 <- as.data.frame(table(cgrelation$reach_outlier))[, 2]
freq1 <- c(freq1[c(2, 1, 3)], 0)
bb$freq1 <- freq1
freq2 <- as.data.frame(table(cgrelation$reach_extremum))[, 2]
freq2 <- c(rev(freq2), 0)
bb$freq2 <- freq2
freq3 <- as.data.frame(table(cgrelation$all_outlier))[, 2]
freq3 <- c(freq3[1], 0, freq3[2], 0)
bb$freq3 <- freq3
bbl <- pivot_longer(bb, freq1:freq3, names_to = "freq", values_to = "count")
bbl$value <- factor(bbl$value, levels = bb$value)