#file requirements
#PC-guided clocks in a .RData file or other formats
#extremum of target sites among all samples in a .csv file, with rows corresponding to individual sites and columns corresponding to the maximum and minimum
#cluster and outlier information of target methylation sites in three .RDS files
#the order of the target sites should be identical between files

library(data.table)

mldata <- fread("extremum_methylation_level.csv", header = TRUE,
  sep = ",", data.table = TRUE)
setkey(mldata, "ID_REF")

load("CalcPCSix.RData")
message("PCClocks Data successfully loaded")
clockname <- c("PCHannum", "PCHorvath1", "PCHorvath2",
  "PCLin", "PCZhang", "PCPhenoAge")

cglist <- readRDS("detailed_methylation_level_enhanced_PC.RDS")
boundarylist <- readRDS("boundary_methylation_level_enhanced_PC.RDS")
distlist <- readRDS("distribution_old_methylation_level_enhanced_PC.RDS")

outd <- which(sapply(boundarylist, function(x) {x[1] == 0}))
outu <- which(sapply(boundarylist, function(x) {x[length(x)] == 1}))

getcluster <- function(CpGs, sample, boundarylist) {
  cluster <- matrix(nrow = length(CpGs), ncol = 2,
    dimnames = list(NULL, c("ID_REF", "Cluster")))

  if(! all(names(sample) == CpGs)) {
    stop("Sample Arrange Problem")
  }

  namelist <- lapply(seq_along(CpGs), function (x) {
    cgvalue <- sample[x]
    breaks <- unique(c(0, boundarylist[[x]], 1))
    interval <- findInterval(cgvalue, breaks, rightmost.closed = TRUE)
    if (interval == 1 && (! pos %in% outd)) {
      name <- "Outlier_down"
    } else if (interval == length(breaks) - 1 && (! pos %in% outu)) {
      name <- "Outlier_up"
    } else {
      lower <- round(breaks[interval], 4)
      upper <- round(breaks[interval + 1], 4)
      name <- paste("(", lower, " - ", upper, "]", sep = "")
    }
    return(name)
  })

  cluster[, 1] <- CpGs
  cluster[, 2] <- unlist(namelist, use.names = FALSE)
  return(cluster)
}

samplelate <- function(mldata, distlist, boundarylist, cglist, times = 1000) {
  if (times < 100) {
    stop("Repeat Times Might Be Too Few")
  }
  if (! all(names(distlist) == CpGs) || ! all(names(cglist) == CpGs)) {
    stop("Arrangement Problem")
  }

  sampledata <- matrix(nrow = length(CpGs), ncol = times,
    dimnames = list(CpGs, paste("sample", 1:times, sep = "")))
  
  set.seed(123)
  for (i in seq_along(CpGs)) {
    cgsingle <- CpGs[i]
    name <- names(distlist[[i]])
    percentage <- distlist[[i]]
    count <- round(percentage * times, 0)
    cgv <- cglist[[i]]
    boundary <- boundarylist[[i]]
    wholesample <- numeric(sum(count))
    currentpos <- 1L

    if (name[1] == "Outlier_down") {
      minv <- mldata[.(cgsingle), "Min", with = FALSE][[1]]
      wholesample[currentpos:(currentpos + count[1] - 1)] <- minv
      currentpos <- currentpos + count[1]
      name <- name[-1]
      count <- count[-1]
    }

    if (length(name) > 0 && name[length(name)] == "Outlier_up") {
      maxv <- mldata[.(cgsingle), "Max", with = FALSE][[1]]
      wholesample[currentpos:(currentpos + count[length(count)] - 1)] <- maxv
      currentpos <- currentpos + count[length(count)]
      name <- name[-length(name)]
      count <- count[-length(count)]
    }

    if (length(boundary) != length(count) + 1) {
      stop(paste(cgsingle, "Boundary Count Matching Problem"))
    }
    for (j in seq_along(count)) {
      lower <- boundary[j]
      upper <- boundary[j + 1]
      startpos <- ifelse(lower == 0, 1, findInterval(lower, cgv) + 1)
      endpos <- ifelse(upper == 1, length(cgv), findInterval(upper, cgv))
      interval <- cgv[startpos:endpos]
      if(count[j] > 0) {
        if (length(interval) == 1) {
          midv <- rep(interval, count[j])
        } else if (length(interval) < count[j]) {
          midv <- interval[sample.int(length(interval), count[j], replace = TRUE)]
        } else {
          midv <- interval[sample.int(length(interval), count[j])]
        }
        wholesample[currentpos:(currentpos + count[j] - 1)] <- midv
        currentpos <- currentpos + count[j]
      }
    }

    if (currentpos <= times) {
      wlower <- boundary[1]
      wupper <- boundary[length(boundary)]
      wmidv <- wholesample[wholesample > wlower & wholesample < wupper]
      wholesample[currentpos:times] <- sample(wmidv, times - currentpos + 1)
    } else if (currentpos > times + 1) {
      remove <- (times + 1):(currentpos - 1)
      wholesample <- wholesample[-remove]
    }
    wholesample <- sample(wholesample)

    sampledata[i, ] <- wholesample
  }

  return(sampledata)
}

calculator <- function (clock, data, trafo = FALSE) {
  if (! all(c("model", "intercept", "center", "rotation") %in% names(clock))) {
    stop("PCClocks Data Missing")
  }
  center <- clock$center
  rotation <- clock$rotation
  modelbeta <- clock$model
  intercept <- clock$intercept
  weightbined <- rotation %*% modelbeta
  if (all(row.names(data) == row.names(weightbined))) {
    message("CpGs Are All Matched")
    age <- apply(data, 2, function(x) {
      sum(weightbined * (x - center)) + intercept
    })
    if (trafo) {
      age <- anti.trafo(age)
    }
    return(as.numeric(age))
  } else {
    stop(paste("Only", sum(row.names(data) == CpGs),
      "CpGs Are Matched"))
  }
}

#optional, at least 1000
times <- 1000
sampleagegroup <- matrix(nrow = length(clockname), ncol = times,
  dimnames = list(clockname, paste("sample", 1:times, sep = "")))
clustergroup <- matrix(nrow = length(clockname), ncol = 2,
  dimnames = list(clockname, c("Outlier 1", "Outlier 2")))

result <- samplelate(mldata, distlist, boundarylist, cglist, times)
sampleagegroup[1, ] <- calculator(CalcPCHannum, result)
sampleagegroup[2, ] <- calculator(CalcPCHorvath1, result, trafo = TRUE)
sampleagegroup[3, ] <- calculator(CalcPCHorvath2, result, trafo = TRUE)
sampleagegroup[4, ] <- calculator(CalcPCLin, result)
sampleagegroup[5, ] <- calculator(CalcPCZhang, result)
sampleagegroup[6, ] <- calculator(CalcPCPhenoAge, result)

positiongroup <- max.col(sampleagegroup, ties.method = "first")
for (i in seq_along(positiongroup)) {
  position <- positiongroup[i]
  positionsample <- result[, position]
  names(positionsample) <- row.names(result)
  cluster <- getcluster(CpGs, positionsample, boundarylist)

  clustergroup[i, ] <- c(sum(cluster[, 2] == "Outlier_down") / nrow(cluster),
    sum(cluster[, 2] == "Outlier_up") / nrow(cluster))
  message(paste(i, "Finished"))

}
