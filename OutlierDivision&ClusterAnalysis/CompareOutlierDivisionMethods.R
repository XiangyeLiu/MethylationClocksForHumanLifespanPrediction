#file requirements
#original clocks in a .csv file, with rows corresponding to individual sites and columns corresponding to the coefficients and which original clock sites are used in
#all values of target sites in a .RData file

library(moments)

cgfile <- read.table("wanted_sites_arrangement.csv", header = TRUE, sep = ",")
cgpartgroup <- split(cgfile, cgfile$Number)
names(cgpartgroup) <- c("Hannum-71 Clock (2013)",
  "Horvath-353 Clock (2013)", "Horvath-391 Clock (2018)",
  "Weidner-102 Clock (2014)", "Weidner-3 Clock (2014)",
  "Lin-99 Clock (2016)", "Lin-3 Clock (2016)",
  "Zhang-514 Clock (2019)", "McEwen-94 Clock (2019)",
  "DNAm PhenoAge")

load(file = "detailed_methylation_level_enhanced.RData")

outlierrange <- function(cglist) {
  data3 <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(data3) <- c("ID_REF", "Method", "Outlier_down", "Outlier_up")

  for (i in seq_along(cglist)) {
    cg <- cglist[[i]]

    q1 <- cg[round(length(cg) / 4, 0)]
    q3 <- cg[round(length(cg) / 4 * 3, 0)]
    iqr <- q3 - q1
    outd1 <- q1 - 1.5 * iqr
    outu1 <- q3 + 1.5 * iqr
    outd1 <- ifelse(outd1 >= 0, outd1, 0)
    outu1 <- ifelse(outu1 <= 1, outu1, 1)

    mediancg <- median(cg)
    absgroup <- abs(cg - mediancg)
    mad <- median(absgroup)
    if (skewness(cg) < -0.5) {
      outd2 <- mediancg - 3 * mad
      outu2 <- mediancg + 5 * mad
    } else if (skewness(cg) > 0.5) {
      outd2 <- mediancg - 5 * mad
      outu2 <- mediancg + 3 * mad
    } else {
      outd2 <- mediancg - 3 * mad
      outu2 <- mediancg + 3 * mad
    }
    outd2 <- ifelse(outd2 >= 0, outd2, 0)
    outu2 <- ifelse(outu2 <= 1, outu2, 1)

    outd3 <- cg[round(length(cg) * 0.05, 0)]
    outu3 <- cg[round(length(cg) * 0.95, 0)]

    data4 <- data.frame("ID_REF" = rep(names(cglist[i]), 3),
      "Method" = c("IQR", "MAD", "PCT"),
      "Outlier_down" = c(outd1, outd2, outd3),
      "Outlier_up" = c(outu1, outu2, outu3))
    data3 <- rbind(data3, data4)
  }

  return(data3)
}

outlier <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(outlier) <- c("ID_REF", "Method",
  "Outlier_down", "Outlier_up", "Clock")
for (j in seq_along(cgpartgroup)) {
  cgpart <- cgpartgroup[[j]]
  cgpartlist <- cglist[which(names(cglist) %in% cgpart[, 1])]
  data3 <- outlierrange(cgpartlist)
  data3$Clock <- names(cgpartgroup[j])
  outlier <- rbind(outlier, data3)
}

write.table(outlier, "outlier_methylation_density.csv",
  row.names = FALSE, sep = ",")