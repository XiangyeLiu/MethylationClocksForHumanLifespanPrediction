#file requirements
#original clocks in a .csv file, with rows corresponding to individual sites and columns corresponding to the coefficients and which original clock sites are used in
#tissue/organ-specific extremum of target sites among samples in two .csv files, with rows corresponding to individual sites and columns corresponding samples, including one file for the maximum and the other for the minimum

library(stringr)

cgdata <- read.table("wanted_sites_arrangement.csv", header = TRUE, sep = ",")
cgpartgroup <- split(cgdata, cgdata$Number)
names(cgpartgroup) <- c("Hannum-71 Clock (2013)",
  "Horvath-353 Clock (2013)", "Horvath-391 Clock (2018)",
  "Weidner-102 Clock (2014)", "Weidner-3 Clock (2014)",
  "Lin-99 Clock (2016)", "Lin-3 Clock (2016)",
  "Zhang-514 Clock (2019)", "McEwen-94 Clock (2019)",
  "DNAm PhenoAge")

datamax <- read.table("extremum_split_max.csv", header = TRUE, sep = ",")
datamin <- read.table("extremum_split_min.csv", header = TRUE, sep = ",")
tissuename <- str_replace_all(colnames(datamax)[-1], "\\.", " ")

anti.trafo <- function(x, adultage = 20) {
  ifelse(x < 0, (1 + adultage) * exp(x) - 1, (1 + adultage) * x + adultage)
}

agemax <- agemin <- data.frame(matrix(nrow = length(cgpartgroup),
  ncol = length(tissuename) + 1))
colnames(agemax) <- colnames(agemin) <- c("Clock", tissuename)
agemax[, 1] <- agemin[, 1] <- names(cgpartgroup)

for (i in seq_along(tissuename) + 1) {
  dmax <- datamax[, c(1, i)]
  dmin <- datamin[, c(1, i)]
  dall <- cbind(dmax, dmin[, -1])
  result <- lapply(seq_along(cgpartgroup), function(x) {
    cgpart <- cgpartgroup[[x]]
    intercept <- cgpart[1, 2]
    coef <- cgpart[-1, 2]
    cgsome <- cgpart[-1, 1]
    cgsome <- cgsome[order(cgsome)]

    dall <- dall[which(dall[, 1] %in% cgsome), ]
    rowmax <- which(is.na(dall[, 2]))
    rowmin <- which(is.na(dall[, 3]))
    if (length(rowmax)) {
      dall[rowmax, 2] <- mean(dall[, 2], na.rm = TRUE)
    }
    if (length(rowmin)) {
      dall[rowmin, 3] <- mean(dall[, 3], na.rm = TRUE)
    }
    dall[, 4] <- ifelse(coef >= 0, dall[, 2], dall[, 3])
    dall[, 5] <- ifelse(coef <= 0, dall[, 2], dall[, 3])

    maxpre <- sum(coef * dall[, 4], intercept)
    minpre <- sum(coef * dall[, 5], intercept)
    if(x %in% c(2, 3, 9)) {
      maxpre <- anti.trafo(maxpre)
      minpre <- anti.trafo(minpre)
    }
    return(c(maxpre, minpre))
  })
  result <- unlist(result, use.names = FALSE)
  maxresult <- result[seq(1, length(result), 2)]
  minresult <- result[seq(2, length(result), 2)]
  agemax[, i] <- maxresult
  agemin[, i] <- minresult
}

write.table(agemax, "age_split_max.csv", row.names = FALSE, sep = ",")
write.table(agemin, "age_split_min.csv", row.names = FALSE, sep = ",")
