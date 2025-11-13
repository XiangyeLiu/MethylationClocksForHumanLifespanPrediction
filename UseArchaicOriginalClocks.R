#file requirements
#original clocks in a .csv file, with rows corresponding to individual sites and columns corresponding to the coefficients and which original clock sites are used in
#information of target sites in two .csv files, with rows corresponding to individual sites and columns corresponding to chromosome position, specific values, and whether directly matched with the methylation map or not, including one file for the Altai Neandertal and the other for the Denisovan

cgdata <- read.table("wanted_sites_arrangement.csv", header = TRUE, sep = ",")
cgpartgroup <- split(cgdata, cgdata$Number)
names(cgpartgroup) <- c("Hannum-71 Clock (2013)",
  "Horvath-353 Clock (2013)", "Horvath-391 Clock (2018)",
  "Weidner-102 Clock (2014)", "Weidner-3 Clock (2014)",
  "Lin-99 Clock (2016)", "Lin-3 Clock (2016)",
  "Zhang-514 Clock (2019)", "McEwen-94 Clock (2019)",
  "DNAm PhenoAge")

datamn <- read.table("methylation_wanted_sites_Altai_Neandertal.csv",
  header = TRUE, sep = ",")
datamd <- read.table("methylation_wanted_sites_Denisovan.csv",
  header = TRUE, sep = ",")

anti.trafo <- function(x, adultage = 20) {
  ifelse(x < 0, (1 + adultage) * exp(x) - 1, (1 + adultage) * x + adultage)
}

agearchaic <- lapply(seq_along(cgpartgroup), function(x) {
  cgpart <- cgpartgroup[[x]]
  cgsome <- cgpart[-1, 1]
  coef <- cgpart[-1, 2]
  intercept <- cgpart[1, 2]
  cgvaluen <- datamn[which(datamn[, 1] %in% cgsome), 5]
  cgvalued <- datamd[which(datamd[, 1] %in% cgsome), 5]
  agen <- sum(coef * cgvaluen, intercept)
  aged <- sum(coef * cgvalued, intercept)
  if (x %in% c(2, 3, 9)) {
    agen <- anti.trafo(agen)
    aged <- anti.trafo(aged)
  }
  return(c(agen, aged))
})

datagearchaic <- data.frame(matrix(unlist(agearchaic, use.names = FALSE),
  nrow = 10, ncol = 2, byrow = TRUE), row.names = names(cgpartgroup))
colnames(datagearchaic) <- c("Neandertal", "Denisovan")

write.table(datagearchaic, "age_archaic_human.csv",
  row.names = TRUE, sep = ",")