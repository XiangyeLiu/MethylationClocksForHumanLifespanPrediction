#file requirements
#original clocks in a .csv file, with rows corresponding to individual sites and columns corresponding to the coefficients and which original clock sites are used in
##pan-mammalian clocks in three .csv files, one clock for each, with rows corresponding to individual sites and columns corresponding to transformed chromosome position, coefficients and relevant gene information

library(ggplot2)
library(ggvenn)
library(patchwork)

orimc <- read.table("methylation_sites_relevant_information.csv",
  header = TRUE, sep = ",")
oripos <- orimc[, 1:3]

panmc1 <- read.table("pan_mammalian_clock1_hg19.csv",
  header = TRUE, sep = ",")
panpos1 <- panmc1[-1, c(1, 3, 4)]
panpos1$Name[227] <- "cg18473521"

panmc2 <- read.table("pan_mammalian_clock2_hg19.csv",
  header = TRUE, sep = ",")
panpos2 <- panmc2[-1, c(1, 3, 4)]

panmc3 <- read.table("pan_mammalian_clock3_hg19.csv",
  header = TRUE, sep = ",")
panpos3 <- panmc3[-1, c(1, 3, 4)]

oripanvenn1 <- list("Original Clock" = oripos$Name, "Clock1" = panpos1$Name)
oripanvenn2 <- list("Original Clock" = oripos$Name, "Clock2" = panpos2$Name)
oripanvenn3 <- list("Original Clock" = oripos$Name, "Clock3" = panpos3$Name)
oripanvenn4 <- list("Original Clock" = oripos$Name, "Clock1" = panpos1$Name,
  "Clock2" = panpos2$Name, "Clock3" = panpos3$Name)