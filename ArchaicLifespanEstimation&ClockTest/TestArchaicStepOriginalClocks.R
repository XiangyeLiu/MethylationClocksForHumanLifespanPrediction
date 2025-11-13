#file requirements
#original clocks in a .csv file, with rows corresponding to individual sites and columns corresponding to the coefficients and which original clock sites are used in
#information of target sites in two .csv files, with rows corresponding to individual sites and columns corresponding to chromosome position, specific values, and whether directly matched with the methylation map or not, including one file for the Altai Neandertal and the other for the Denisovan

library(tidyverse)
library(ggplot2)
library(patchwork)

cgdata <- read.table("wanted_sites_arrangement.csv", header = TRUE, sep = ",")
cgpartgroup <- split(cgdata, cgdata$Number)
names(cgpartgroup) <- c("Hannum-71 Clock (2013)",
  "Horvath-353 Clock (2013)", "Horvath-391 Clock (2018)",
  "Weidner-102 Clock (2014)", "Weidner-3 Clock (2014)",
  "Lin-99 Clock (2016)", "Lin-3 Clock (2016)",
  "Zhang-514 Clock (2019)", "McEwen-94 Clock (2019)",
  "DNAm PhenoAge")

anti.trafo <- function(x, adultage = 20) {
  ifelse(x < 0, (1 + adultage) * exp(x) - 1, (1 + adultage) * x + adultage)
}

mno<- read.table("methylation_wanted_sites_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

nog <- data.frame("Number" = 0:nrow(mno), "Site" = c("(Intercept)", mno[, 1]))
for (i in seq_along(cgpartgroup)) {
  cgpart <- cgpartgroup[[i]]
  intercept <- cgpart[1, 2]
  mergedcm <- merge(cgpart, mno, by.x = "Site", by.y = "Name")
  mp1 <- mergedcm[, 2] * mergedcm[, 4]
  names(mp1) <- mergedcm$Site
  cgleft <- mno[(! mno$Name %in% mergedcm$Site), 1]
  mp2 <- numeric(length(cgleft))
  names(mp2) <- cgleft
  mp <- c(mp1, mp2)
  mp <- mp[order(names(mp))]
  if (length(mp) != nrow(mno)) {
    stop("Arrange Problem")
  }
  sumg <- numeric(length(mp) + 1)
  sumg[1] <- intercept
  for (j in 2:length(sumg)) {
    sumg[j] <- sumg[j - 1] + mp[j - 1]
  }
  nog[, i + 2] <- sumg
}

colnames(nog)[3:ncol(nog)] <- names(cgpartgroup)

nog[, 4] <- anti.trafo(nog[, 4])
nog[, 5] <- anti.trafo(nog[, 5])
nog[, 11] <- anti.trafo(nog[, 11])
write.table(nog, "cumulative_sum_age_prediction_Neandertal.csv",
  row.names = FALSE, sep = ",")

nogl <- pivot_longer(nog[, -ncol(nog)], "Hannum-71 Clock (2013)":"McEwen-94 Clock (2019)",
  names_to = "Clock", values_to = "Value")
write.table(nogl, "cumulative_sum_age_long_Neandertal.csv",
  row.names = FALSE, sep = ",")

mdo<- read.table("methylation_wanted_sites_Denisovan.csv",
  header = TRUE, sep = ",")

dog <- data.frame("Number" = 0:nrow(mdo), "Site" = c("(Intercept)", mdo[, 1]))
for (i in seq_along(cgpartgroup)) {
  cgpart <- cgpartgroup[[i]]
  intercept <- cgpart[1, 2]
  mergedcm <- merge(cgpart, mdo, by.x = "Site", by.y = "Name")
  mp1 <- mergedcm[, 2] * mergedcm[, 4]
  names(mp1) <- mergedcm$Site
  cgleft <- mdo[(! mdo$Name %in% mergedcm$Site), 1]
  mp2 <- numeric(length(cgleft))
  names(mp2) <- cgleft
  mp <- c(mp1, mp2)
  mp <- mp[order(names(mp))]
  if (length(mp) != nrow(mdo)) {
    stop("Arrange Problem")
  }
  sumg <- numeric(length(mp) + 1)
  sumg[1] <- intercept
  for (j in 2:length(sumg)) {
    sumg[j] <- sumg[j - 1] + mp[j - 1]
  }
  dog[, i + 2] <- sumg
}

colnames(dog)[3:ncol(dog)] <- names(cgpartgroup)

dog[, 4] <- anti.trafo(dog[, 4])
dog[, 5] <- anti.trafo(dog[, 5])
dog[, 11] <- anti.trafo(dog[, 11])
write.table(dog, "cumulative_sum_age_prediction_Denisovan.csv",
  row.names = FALSE, sep = ",")

dogl <- pivot_longer(dog[, -ncol(dog)], "Hannum-71 Clock (2013)":"McEwen-94 Clock (2019)",
  names_to = "Clock", values_to = "Value")
write.table(dogl, "cumulative_sum_age_long_Denisovan.csv",
  row.names = FALSE, sep = ",")