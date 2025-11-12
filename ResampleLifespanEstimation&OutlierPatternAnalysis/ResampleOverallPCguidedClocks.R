#file requirements
#PC-guided clocks in a .RData file or other formats
#extremum of target sites among all samples in a .csv file, with rows corresponding to individual sites and columns corresponding to the maximum and minimum
#cluster and outlier information of all combined target methylation sites in twwo .RDS files
#the order of the target sites should be identical between files

library(data.table)

mldata <- fread("extremum_methylation_level.csv", header = TRUE,
  sep = ",", data.table = TRUE)
setkey(mldata, "ID_REF")

load("CalcPCBigSix.RData")
message("PCClocks Data successfully loaded")
clockname <- c("PCHannum", "PCHorvath1", "PCHorvath2",
  "PCLin", "PCZhang", "PCPhenoAge")

cglist <- readRDS("detailed_methylation_level_enhanced_PC.RDS")
wholelist <- readRDS("whole_old_methylation_level_enhanced_PC.RDS")

cgsome <- sort(unlist(cglist, use.names = FALSE))

boundary <- wholelist[[1]]
lower <- boundary[1]
upper <- boundary[2]
startpos <- findInterval(lower, cgsome) + 1
endpos <- findInterval(upper, cgsome)
interval2 <- cgsome[startpos:endpos]
lower <- boundary[2]
upper <- boundary[3]
startpos <- findInterval(lower, cgsome) + 1
endpos <- findInterval(upper, cgsome)
interval3 <- cgsome[startpos:endpos]
lower <- boundary[3]
upper <- boundary[4]
startpos <- findInterval(lower, cgsome) + 1
endpos <- findInterval(upper, cgsome)
interval4 <- cgsome[startpos:endpos]
  
getcluster <- function(CpGs, sample, wholelist) {
  cluster <- matrix(nrow = length(CpGs), ncol = 2,
    dimnames = list(NULL, c("ID_REF", "Cluster")))
  breaks <- unique(c(0, wholelist[[1]], 1))
  
  if(! all(names(sample) == CpGs)) {
    stop("Sample Arrange Problem")
  }

  namelist <- lapply(seq_along(CpGs), function (x) {
    cgvalue <- sample[x]
    interval <- findInterval(cgvalue, breaks, rightmost.closed = TRUE)
    if (interval == 1) {
      name <- "Outlier1"
    } else if (interval == length(breaks) - 1) {
      name <- "Outlier2"
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

samplelate <- function(mldata, wholelist, cgsome, times = 1000) {
  if (times < 1000) {
    stop("Repeat Times Might Be Too Few")
  }

  sampledata <- matrix(nrow = length(CpGs), ncol = times,
    dimnames = list(CpGs, paste("sample", 1:times, sep = "")))
  
  count <- round(wholelist[[2]] * times, 0)

  set.seed(123)
  for (i in seq_along(CpGs)) {
    cgsingle <- CpGs[i]
    wholesample <- numeric(sum(count))
    currentpos <- 1L

    minv <- mldata[.(cgsingle), "Min", with = FALSE][[1]]
    wholesample[currentpos:(currentpos + count[1] - 1)] <- minv
    currentpos <- currentpos + count[1]

    maxv <- mldata[.(cgsingle), "Max", with = FALSE][[1]]
    wholesample[currentpos:(currentpos + count[5] - 1)] <- maxv
    currentpos <- currentpos + count[5]

    midv <- interval2[sample.int(length(interval2), count[2])]
    wholesample[currentpos:(currentpos + count[2] - 1)] <- midv
    currentpos <- currentpos + count[2]
    midv <- interval3[sample.int(length(interval3), count[3])]
    wholesample[currentpos:(currentpos + count[3] - 1)] <- midv
    currentpos <- currentpos + count[3]
    midv <- interval4[sample.int(length(interval4), count[4])]
    wholesample[currentpos:(currentpos + count[4] - 1)] <- midv
    currentpos <- currentpos + count[4]

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
  dimnames = list(clockname, c("Outlier1", "Outlier2")))

result <- samplelate(mldata, wholelist, cgsome, times)
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
  cluster <- getcluster(CpGs, positionsample, wholelist)

  clustergroup[i, ] <- c(sum(cluster[, 2] == "Outlier1") / nrow(cluster),
    sum(cluster[, 2] == "Outlier2") / nrow(cluster))
  message(paste(i, "Finished"))
}