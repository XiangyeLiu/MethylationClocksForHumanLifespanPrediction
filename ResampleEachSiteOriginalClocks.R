#file requirements
#original clocks in a .csv file, with rows corresponding to individual sites and columns corresponding to the coefficients and which original clock sites are used in
#extremum of target sites among all samples in a .csv file, with rows corresponding to individual sites and columns corresponding to the maximum and minimum
#cluster and outlier information of target methylation sites in three .RData files

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

load("detailed_methylation_level.RData")
load("boundary_methylation_level.RData")
load("distribution_old_methylation_level.RData")

anti.trafo <- function(x, adultage = 20) {
  ifelse(x < 0, (1 + adultage) * exp(x) - 1, (1 + adultage) * x + adultage)
}

getcluster <- function(sample, boundarylist) {
  cggroup <- row.names(sample)
  boundarylistp <- boundarylist[which(names(boundarylist) %in% cggroup)]
  outd <- which(sapply(boundarylistp, function(x) {x[1] == 0}))
  outu <- which(sapply(boundarylistp, function(x) {x[length(x)] == 1}))

  cluster <- data.frame(matrix(nrow = length(cggroup), ncol = 2))
  colnames(cluster) <- c("ID_REF", "Cluster")

  for (i in seq_along(cggroup)) {
    cgvalue <- sample[i, 1]
    breaks <- unique(c(0, boundarylistp[[i]], 1))
    interval <- cut(cgvalue, breaks = breaks, include.lowest = TRUE,
      right = TRUE)
    table <- table(interval)
    name <- names(table)
    if (! i %in% outd) {
      name[1] <- "Outlier_down"
    }
    if (! i %in% outu) {
      name[length(name)] <- "Outlier_up"
    }
    value <- as.vector(table)
    name <- name[which(value != 0)]
    cluster[i, ] <- c(cggroup[i], name)
  }

  return(cluster)
}

samplelate <- function(cgpart, mldata, distlist, cglist, times = 10000) {
  if (times < 10000) {
    stop("Repeat Times Might Be Too Few")
  }
  intercept <- cgpart[1, 2]
  coef <- cgpart[-1, 2]
  cggroup <- cgpart[-1, 1]
  cggroup <- cggroup[order(cggroup)]
  distlistp <- distlist[which(names(distlist) %in% cggroup)]
  cglistp <- cglist[which(names(cglist) %in% cggroup)]
  if (! all(names(distlistp) == cggroup) || ! all(names(cglistp) == cggroup)) {
    stop("Arrangement Problem")
  }

  sampledata <- data.frame(matrix(nrow = length(cggroup), ncol = times),
    row.names = cggroup)
  colnames(sampledata) <- paste("sample", 1:times, sep = "")

  for (i in seq_along(cggroup)) {
    set.seed(123)
    cgsingle <- cggroup[i]
    name <- names(distlistp[[i]])
    percentage <- round(distlistp[[i]], 4)
    cgv <- cglistp[[i]]

    minsample <- midsample <- maxsample <- NULL

    if (name[1] == "Outlier_down") {
      min <- mldata[which(row.names(mldata) == cgsingle), 2]
      minsample <- rep(min, percentage[1] * times)
      name <- name[-1]
      percentage <- percentage[-1]
    }

    if (name[length(name)] == "Outlier_up") {
      max <- mldata[which(row.names(mldata) == cgsingle), 1]
      maxsample <- rep(max, percentage[length(percentage)] * times)
      name <- name[-length(name)]
      percentage <- percentage[-length(percentage)]
    }

    boundary <- unique(as.numeric(sapply(name, function(x) {
      unlist(regmatches(x, gregexpr("\\d+\\.?\\d*", x)))
    })))
    if (length(boundary) != length(percentage) + 1) {
      stop(paste(cgsingle, "Boundary Percentage Matching Problem"))
    }
    for (j in seq_along(percentage)) {
      interval <- cgv[cgv >= boundary[j] & cgv < boundary[j + 1]]
      midsample <- c(midsample,
        sample(interval, percentage[j] * times, replace = TRUE))
    }

    wholesample <- c(minsample, midsample, maxsample)
    if (times > length(wholesample)) {
      wholesample <- c(wholesample,
        sample(midsample, times - length(wholesample)))
    } else if (length(wholesample) > times) {
      remove <- sample(which(wholesample %in% midsample),
        length(wholesample) - times)
      wholesample <- wholesample[-remove]
    }
    wholesample <- sample(wholesample)

    sampledata[i, ] <- wholesample
  }

  sampleage <- coef %*% as.matrix(sampledata) + intercept
  sampleage <- as.vector(sampleage)

  return(list(sampleage, sampledata))
}

#optional, at least 10000
times <- 10000
sampleagegroup <- data.frame(matrix(nrow = length(cgpartgroup), ncol = times),
  row.names = names(cgpartgroup))
colnames(sampleagegroup) <- paste("sample", 1:times, sep = "")
clustergroup <- data.frame(matrix(nrow = length(cgpartgroup), ncol = 2),
  row.names = names(cgpartgroup))
colnames(clustergroup) <- c("Outlier_down", "Outlier_up")
for (i in seq_along(cgpartgroup)) {
  cgpart <- cgpartgroup[[i]]
  result <- samplelate(cgpart, mldata, distlist, cglist, times)
  sampleage <- result[[1]]
  if (i %in% c(2, 3, 9)) {
    sampleage <- anti.trafo(sampleage)
  }

  position <- which.max(sampleage)
  positionsample <- data.frame(result[[2]][, position],
    row.names = row.names(result[[2]]))
  cluster <- getcluster(positionsample, boundarylist)

  sampleagegroup[i, ] <- sampleage
  clustergroup[i, ] <- c(sum(cluster[, 2] == "Outlier_down") / nrow(cluster),
    sum(cluster[, 2] == "Outlier_up") / nrow(cluster))
  message(paste(i, "Finished"))
}