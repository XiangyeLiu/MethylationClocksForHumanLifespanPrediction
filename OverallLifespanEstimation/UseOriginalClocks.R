#file requirements
#Original clocks in a .csv file, with rows corresponding to individual sites and columns corresponding to the coefficients and which original clock sites are used in
#extremum of target sites among all samples in a .csv file, with rows corresponding to individual sites and columns corresponding to the maximum and minimum

library(ggplot2)

mldata <- read.csv("extremum_methylation_level.csv", header = FALSE)
sitedata <- read.csv("wanted_sites_arrangement.csv", header = TRUE)
colnames(mldata) <- c("Site", "Max", "Min")

number <- unique(sitedata$Number)
barrier <- 20
limitation <- data.frame()

for (i in number) {
  partdata <- sitedata[which(sitedata$Number == i), ]
  intercept <- partdata[1, 2]
  cgdata <- merge(partdata[, -3], mldata, by = "Site")
  pnsite <- which(cgdata$Coefficient > 0)
  ngsite <- which(cgdata$Coefficient < 0)
  max <- sum(cgdata[pnsite, 2] * cgdata[pnsite, 3],
    cgdata[ngsite, 2] * cgdata[ngsite, 4], intercept)
  min <- sum(cgdata[pnsite, 2] * cgdata[pnsite, 4],
    cgdata[ngsite, 2] * cgdata[ngsite, 3], intercept)
  if (i %in% c(2, 3, 7)) {
    if (max >= 0) {
      max <- max * (barrier + 1) + barrier
    } else {
      max <- exp(max) * (barrier + 1) - 1
    }
    if (min >= 0) {
      min <- min * (barrier + 1) + barrier
    } else {
      min <- exp(min) * (barrier + 1) - 1
    }
  }
  limitation <- rbind(limitation, c(max, min))
}
name <- c("Hannum-71 Clock (2013)",
  "Horvath-353 Clock (2013)", "Horvath-391 Clock (2018)",
  "Weidner-102 Clock (2014)", "Weidner-3 Clock (2014)",
  "Lin-99 Clock (2016)", "Lin-3 Clock (2016)",
  "Zhang-514 Clock (2019)", "McEwen-94 Clock (2019)",
  "DNAm PhenoAge")
limitation <- cbind(name, limitation)
colnames(limitation) <- c("Name", "Max", "Min")
limitationplot <- data.frame(Name = rep(name, 2), Extremum = c(limitation[, 2], limitation[, 3]))
limitationplot[, 1] <- factor(limitationplot[, 1], levels = rev(name))