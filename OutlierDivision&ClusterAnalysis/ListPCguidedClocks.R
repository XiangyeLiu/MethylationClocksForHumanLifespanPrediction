#file requirements
#target sites in a .txt file
#all methylation values of samples in each series in a .csv file

options(scipen = 100)
filename <- list.files(path = ".", "*.csv")
cggroup <- readLines("shared_methylation_sites.txt")

cglist <-list()
for (i in seq_along(cggroup)) {
  cglist[[i]] <- 2
}
names(cglist) <- cggroup

cggrab <- function(data, cggroup, cglist) {
  row <- which(row.names(data) %in% cggroup)
  if (! length(row)) {
    stop("No Sites Matched or Row Names Should Be Sites")
  }
  data <- data[row, ]
  for (i in seq_along(row.names(data))) {
    cgname <- row.names(data)[i]
    cgadd <- as.numeric(unlist(data[i, ], use.names = FALSE))
    if (cglist[[cgname]][1] == 2) {
      cglist[[cgname]] <- c(cglist[[cgname]][-1], cgadd[cgadd >= 0 & cgadd <= 1 & ! is.na(cgadd)])
    } else {
      cglist[[cgname]] <- c(cglist[[cgname]], cgadd[cgadd >= 0 & cgadd <= 1 & ! is.na(cgadd)])
    }
  }
  return(cglist)
}

for (file in filename) {
  data3 <- read.table(file, header = TRUE, row.names = 1, sep = ",")
  cglist <- cggrab(data3, cggroup, cglist)
  message(paste(file, "Finished"))
}

saveRDS(cglist, file = "detailed_methylation_level_enhanced_PC.RDS", compress = "xz")


kgroup <- sapply(cglist, function(x) {
  cgsome <- x
  dcgsome <- density(cgsome)
  set.seed(123)
  wss <- sapply(1:10, function(y) {
    kmeans(dcgsome$y, centers = y, nstart = 25)$tot.withinss
  })
  secdifference <- diff(diff(wss))
  if (all(secdifference <= 1000)) {
    k <- which.max(secdifference) + 1
  } else {
    candidate <- which(secdifference > 1000)
    k  <- candidate[length(candidate)] + 1
  }
  return(k)
})

boundarylist <- lapply(seq_along(kgroup), function(x) {
  cgsome <- cglist[[x]]
  dcgsome <- density(cgsome)
  cluster <- kmeans(dcgsome$y, centers = kgroup[x])$cluster
  boundary <- dcgsome$x[which(diff(cluster) != 0)]
  num <- length(boundary) / 2 + 1
  name <- paste("boundary", c(seq(2, num, 1), rev(seq(2, num, 1))), sep = "")

  q1 <- cgsome[round(length(cgsome) / 4, 0)]
  q3 <- cgsome[round(length(cgsome) / 4 * 3, 0)]
  iqr <- q3 - q1
  outd <- q1 - 1.5 * iqr
  outu <- q3 + 1.5 * iqr
  outd <- ifelse(outd > 0, outd, 0)
  outu <- ifelse(outu < 1, outu, 1)

  remove1 <- which(boundary < 0)
  remove2 <- which(boundary > 1)
  remove <- c(remove1, remove2)
  if (length(remove)) {
    boundary <- boundary[-remove]
    name <- name[-remove]
  }
  boundary <- sort(c(outd, boundary, outu))
  name <- c("Outlier_down", name, "Outlier_up")
  names(boundary) <- name
  return(boundary)
})
names(boundarylist) <- names(cglist)

saveRDS(boundarylist, file = "boundary_methylation_level_enhanced_PC.RDS")


data3 <- read.table("age_old_enhanced_V2.csv",
  header = TRUE, row.names = 1, sep = ",")
agegroup <- unlist(data3[nrow(data3), ], use.names = FALSE)
data3 <- data3[-nrow(data3), ]

if (! all(row.names(data3) == names(boundarylist))) {
  stop("Arrangement Problem")
}

outd <- which(sapply(boundarylist, function(x) {x[1] == 0}))
outu <- which(sapply(boundarylist, function(x) {x[length(x)] == 1}))

distlist <- lapply(seq_along(row.names(data3)), function(x) {
  cgvalue <- na.omit(unlist(data3[x, ], use.names = FALSE))
  breaks <- unique(c(0, boundarylist[[x]], 1))
  interval <- cut(cgvalue, breaks = breaks,
    include.lowest = TRUE, right = FALSE)
  level <- levels(interval)
  pattern <- "(\\d+\\.?\\d*)e\\-0(\\d+)"
  change <- which(grepl(pattern, level))
  if (length(change)) {
    for (i in change) {
      number <- regmatches(level[i], gregexpr(pattern, level[i]))[[1]]
      number <- as.character(as.numeric(number))
      levels(interval)[i] <- str_replace_all(level[i], pattern, number)
    }
  }
  table <- prop.table(table(interval))
  name <- names(table)
  if (! x %in% outd) {
    name[1] <- "Outlier_down"
  }
  if (! x %in% outu) {
    name[length(name)] <- "Outlier_up"
  }
  value <- as.vector(table)
  names(value) <- name
  return(value)
})
names(distlist) <- names(boundarylist)

saveRDS(distlist, file = "distribution_old_methylation_level_enhanced_PC.RDS")


wholelist <- list()
cgsome <- unlist(cglist, use.names = FALSE)
cgsome <- cgsome[order(cgsome)]
dcgsome <- density(cgsome)
cluster <- kmeans(dcgsome$y, centers = 2)$cluster
boundary <- dcgsome$x[which(diff(cluster) != 0)]
name <- paste("boundary", rep(2, 2), sep = "")

outd <- cgsome[round(length(cgsome) * 0.05, 0)]
outu <- cgsome[round(length(cgsome) * 0.95, 0)]

boundary <- c(outd, boundary, outu)
name <- c("Outlier_down", name, "Outlier_up")
names(boundary) <- name

wholelist[[1]] <- boundary

cgvalue <- as.vector(na.omit(unlist(data3, use.names = FALSE)))
cgvalue <- cgvalue[order(cgvalue)]
breaks <- c(0, boundary, 1)
interval <- cut(cgvalue, breaks = breaks, include.lowest = TRUE, right = FALSE)
table <- prop.table(table(interval))
tablename <- names(table)
tablename[1] <- "Outlier_down"
tablename[length(tablename)] <- "Outlier_up"
value <- as.vector(table)
names(value) <- tablename

wholelist[[2]] <- value

saveRDS(wholelist, file = "whole_old_methylation_level_enhanced_PC.RDS")
