#file requirements
#original clocks in a .csv file, with rows corresponding to individual sites and columns corresponding to the coefficients and which original clock sites are used in
#extremum of target sites among all samples in a .csv file, with rows corresponding to individual sites and columns corresponding to the maximum and minimum
#cluster and outlier information of all combined target methylation sites in two .RData files

cgdata <- read.table("wanted_sites_arrangement.csv", header = TRUE, sep = ",")
cgpartgroup <- split(cgdata, cgdata$Number)
names(cgpartgroup) <- c("Hannum-71 Clock (2013)",
  "Horvath-353 Clock (2013)", "Horvath-391 Clock (2018)",
  "Weidner-102 Clock (2014)", "Weidner-3 Clock (2014)",
  "Lin-99 Clock (2016)", "Lin-3 Clock (2016)",
  "Zhang-514 Clock (2019)", "McEwen-94 Clock (2019)",
  "DNAm PhenoAge")

mldata <- read.table("extremum_methylation_level.csv", header = FALSE,
  row.names = 1, sep = ",")

load("detailed_methylation_level_enhanced.RData")
load("whole_methylation_level_enhanced.RData")

cgsome <- sort(unlist(cglist, use.names = FALSE))

anti.trafo <- function(x, adultage = 20) {
  ifelse(x < 0, (1 + adultage) * exp(x) - 1, (1 + adultage) * x + adultage)
}

wholecluster <- function(sample, wholelist) {
  boundary <- sort(wholelist[[1]])
  cluster <- matrix(nrow = length(sample), ncol = 2,
    dimnames = list(seq_along(sample), c("ID_REF", "Cluster")))

  for (i in seq_along(sample)) {
    interval <- cut(sample[i], breaks = c(0, boundary, 1),
      include.lowest = TRUE, right = TRUE)
    table <- table(interval)
    name <- names(table)
    value <- as.numeric(table)
    name <- name[which(value != 0)]
    cluster[i, ] <- c(names(sample[i]), name)
  }

  return(cluster)
}

samplewhole <- function(cgsome, mldata, wholelist, times = 100) {
  if (times < 100) {
    stop("Repeat Times Might Be Too Few")
  }

  sampledata <- matrix(nrow = nrow(mldata), ncol = times,
    dimnames = list(row.names(mldata), paste("sample", 1:times, sep = "")))

  set.seed(123)

  boundary <- sort(wholelist[[1]])
  percentage <- wholelist[[2]]
  percentagecount <- round(percentage * times, 0)

  pos1 <- percentagecount[1]
  pos2 <- percentagecount[1] + percentagecount[5]
  pos3 <- sum(percentagecount[c(1, 2, 5)])
  pos4 <- sum(percentagecount[c(1, 2, 3, 5)])
  pos5 <- sum(percentagecount)

  for (i in seq_along(row.names(mldata))) {
    wholesample <- numeric(times)
    wholesample[1:pos1] <- rep(mldata[i, 2], percentagecount[1])
    wholesample[(pos1 + 1):pos2] <- rep(mldata[i, 1],
      percentagecount[5])
    interval <- cgsome[cgsome > boundary[1] & cgsome < boundary[2]]
    wholesample[(pos2 + 1):pos3] <- interval[sample.int(length(interval), percentagecount[2])]
    interval <- cgsome[cgsome > boundary[2] & cgsome < boundary[3]]
    wholesample[(pos3 + 1):pos4] <- interval[sample.int(length(interval), percentagecount[3])]
    interval <- cgsome[cgsome > boundary[3] & cgsome < boundary[4]]
    wholesample[(pos4 + 1):pos5] <- interval[sample.int(length(interval), percentagecount[4])]

    if (length(wholesample) > times) {
      remove <- sample((pos2 + 1):pos5,
        length(wholesample) - times)
      wholesample <- wholesample[-remove]
    }
    wholesample <- sample(wholesample)

    sampledata[i, ] <- wholesample
  }

  return(sampledata)
}

#optional, at least 10000
times <- 10000
sampleagegroup <-matrix(nrow = length(cgpartgroup), ncol = times,
  dimnames = list(names(cgpartgroup), paste("sample", 1:times, sep = "")))
clustergroup <- matrix(nrow = length(cgpartgroup), ncol = 2,
  dimnames = list(names(cgpartgroup), c("Outlier1", "Outlier2")))
result <- samplewhole(cgsome, mldata, wholelist, times)
for (i in seq_along(cgpartgroup)) {
  cgpart <- cgpartgroup[[i]]
  intercept <- cgpart[1, 2]
  coef <- cgpart[-1, 2]
  cgv <- sort(cgpart[-1, 1])
  resultpart <- result[which(row.names(result) %in% cgv), ]
  if (! all(row.names(resultpart) == cgv)) {
    stop("Arrangement Problem")
  }
  sampleage <- coef %*% resultpart + intercept
  sampleage <- as.vector(sampleage)
  if (i %in% c(2, 3, 9)) {
    sampleage <- anti.trafo(sampleage)
  }

  position <- which.max(sampleage)
  positionsample <- resultpart[, position]
  names(positionsample) <- row.names(resultpart)
  cluster <- wholecluster(positionsample, wholelist)

  sampleagegroup[i, ] <- sampleage
  clustergroup[i, ] <- c(sum(cluster[, 2] == "Outlier_down") / nrow(cluster),
    sum(cluster[, 2] == "Outlier_up") / nrow(cluster))
  message(paste(i, "Finished"))
}