#file requirements
#original clocks in a .csv file, with rows corresponding to individual sites and columns corresponding to the coefficients and which original clock sites are used in
#information of imputed target sites in six .csv files, with rows corresponding to individual sites and columns corresponding to specific imputed values, each imputation strategy including one file for the Altai Neandertal and the other for the Denisovan

cgdata <- read.table("wanted_sites_arrangement.csv", header = TRUE, sep = ",")
cgpartgroup <- split(cgdata, cgdata$Number)
names(cgpartgroup) <- c("Hannum-71 Clock",
  "Horvath-353 Clock", "Horvath-391 Clock",
  "Weidner-102 Clock", "Weidner-3 Clock",
  "Lin-99 Clock", "Lin-3 Clock",
  "Zhang-514 Clock", "McEwen-94 Clock",
  "DNAm PhenoAge")

anti.trafo <- function(x, adultage = 20) {
  ifelse(x < 0, (1 + adultage) * exp(x) - 1, (1 + adultage) * x + adultage)
}

dnamni <- read.table("imputation_pmm_Altai_Neandertal.csv",
  header = TRUE, sep = ",")
dnamdi <- read.table("imputation_pmm_Denisovan.csv",
  header = TRUE, sep = ",")

agearchaic1 <- lapply(seq_along(cgpartgroup), function(x) {
  cgpart <- cgpartgroup[[x]]
  cgsome <- cgpart[-1, 1]
  coef <- cgpart[-1, 2]
  intercept <- cgpart[1, 2]
  cgvaluen1 <- dnamni[which(dnamni[, 1] %in% cgsome), 5]
  cgvaluen2 <- dnamni[which(dnamni[, 1] %in% cgsome), 6]
  cgvalued1 <- dnamdi[which(dnamdi[, 1] %in% cgsome), 5]
  cgvalued2 <- dnamdi[which(dnamdi[, 1] %in% cgsome), 6]
  agen1 <- sum(coef * cgvaluen1, intercept)
  agen2 <- sum(coef * cgvaluen2, intercept)
  aged1 <- sum(coef * cgvalued1, intercept)
  aged2 <- sum(coef * cgvalued2, intercept)
  if (x %in% c(2, 3, 9)) {
    agen1 <- anti.trafo(agen1)
    agen2 <- anti.trafo(agen2)
    aged1 <- anti.trafo(aged1)
    aged2 <- anti.trafo(aged2)
  }
  return(c(agen1, agen2, aged1, aged2))
})

datagearchaic1 <- data.frame(matrix(unlist(agearchaic1, use.names = FALSE),
  nrow = 40, ncol = 1, byrow = TRUE))
colnames(datagearchaic1) <- "Imputation"
datagearchaic1$Clock <- rep(names(cgpartgroup), each = 4)
datagearchaic1$Type <- rep(rep(c("Neanderthal", "Denisovan"), each = 2), 10)
datagearchaic1$Method <- rep(c("PopMean", "PopMedian"), 20)
datagearchaic1 <- datagearchaic1[, c(2, 1, 3, 4)]

dnamni <- read.table("imputation_knn_Altai_Neandertal.csv",
  header = TRUE, sep = ",")
dnamdi <- read.table("imputation_knn_Denisovan.csv",
  header = TRUE, sep = ",")

agearchaic2 <- lapply(seq_along(cgpartgroup), function(x) {
  cgpart <- cgpartgroup[[x]]
  cgsome <- cgpart[-1, 1]
  coef <- cgpart[-1, 2]
  intercept <- cgpart[1, 2]
  agengroup <- apply(dnamni[, -1:-4], 2, function(y) {
    cgvaluen <- y[which(dnamni[, 1] %in% cgsome)]
    agen <- sum(coef * cgvaluen, intercept)
    if (x %in% c(2, 3, 9)) {
      agen <- anti.trafo(agen)
    }
    return(agen)
  })
  agedgroup <- apply(dnamdi[, -1:-4], 2, function(y) {
    cgvalued <- y[which(dnamdi[, 1] %in% cgsome)]
    aged <- sum(coef * cgvalued, intercept)
    if (x %in% c(2, 3, 9)) {
      aged <- anti.trafo(aged)
    }
    return(aged)
  })
  return(c(agengroup, agedgroup))
})
datagearchaic2 <- data.frame(Clock = rep(names(cgpartgroup), each = 200),
  Imputation = unlist(agearchaic2, use.names = FALSE),
  Type = rep(rep(c("Neanderthal", "Denisovan"), each = 100), 10),
  Method = rep("KNN", 2000))

dnamnia <- read.table("imputation_sva_Altai_Neandertal.csv",
  header = TRUE, sep = ",")
dnamdia <- read.table("imputation_sva_Denisovan.csv",
  header = TRUE, sep = ",")

agearchaic3 <- lapply(seq_along(cgpartgroup), function(x) {
  cgpart <- cgpartgroup[[x]]
  cgsome <- cgpart[-1, 1]
  coef <- cgpart[-1, 2]
  intercept <- cgpart[1, 2]
  cgrdyn <- dnamni[, 5]
  agengroup <- apply(dnamnia[, -1:-2], 2, function(y) {
    cgrdyn[dnamnia[, 2]] <- y
    cgvaluen <- cgrdyn[which(dnamni[, 1] %in% cgsome)]
    agen <- sum(coef * cgvaluen, intercept)
    if (x %in% c(2, 3, 9)) {
      agen <- anti.trafo(agen)
    }
    return(agen)
  })
  cgrdyd <- dnamdi[, 5]
  agedgroup <- apply(dnamdia[, -1:-2], 2, function(y) {
    cgrdyd[dnamdia[, 2]] <- y
    cgvalued <- cgrdyd[which(dnamdi[, 1] %in% cgsome)]
    aged <- sum(coef * cgvalued, intercept)
    if (x %in% c(2, 3, 9)) {
      aged <- anti.trafo(aged)
    }
    return(aged)
  })
  return(c(agengroup, agedgroup))
})
datagearchaic3 <- data.frame(Clock = rep(names(cgpartgroup), each = 40000),
  Imputation = unlist(agearchaic3, use.names = FALSE),
  Type = rep(rep(c("Neanderthal", "Denisovan"), each = 20000), 10),
  Method = rep(rep(c("Sample", "Random"), each = 10000), 20))

datagearchaic <- rbind(datagearchaic1, datagearchaic2, datagearchaic3)

write.table(datagearchaic, "imputation_result_long.csv",
  row.names = FALSE, sep = ",")
