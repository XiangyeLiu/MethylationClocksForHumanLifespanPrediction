#file requirements
#information of imputed target sites in six .csv files, with rows corresponding to individual sites and columns corresponding to chromosome position, coefficients, and specific imputed values, each clock including one file for the Altai Neandertal and the other for the Denisovan

panmcusenc1 <- read.table("pan_mammalian_clock1_hg19_convert_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

panmcusenc2 <- read.table("pan_mammalian_clock2_hg19_convert_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

panmcusenc3 <- read.table("pan_mammalian_clock3_hg19_convert_Altai_Neandertal.csv",
  header = TRUE, sep = ",")

panmcusedc1 <- read.table("pan_mammalian_clock1_hg19_convert_Denisovan.csv",
  header = TRUE, sep = ",")

panmcusedc2 <- read.table("pan_mammalian_clock2_hg19_convert_Denisovan.csv",
  header = TRUE, sep = ",")

panmcusedc3 <- read.table("pan_mammalian_clock3_hg19_convert_Denisovan.csv",
  header = TRUE, sep = ",")

panmcuse <- function(coefficient, value, type = 1, max = 65.7, gestation = 280 / 365, asm = 13) {
  intercept <- coefficient[1]
  result <- sum(coefficient[-1] * value[-1], intercept)
  if (type == 1) {
    message("Using Clock 1")
    return(exp(result) - 2)
  } else if (type == 2) {
    message("Using Clock 2")
    return(exp(-1 * exp(-1 * result)) * (max + gestation) - gestation)
  } else if (type == 3) {
    message("Using Clock 3")
    m <- 5.0 * (gestation / asm) ^ 0.38
    if (result >= 0) {
      return(m * (asm + gestation) * (result + 1) - gestation)
    } else {
      return(m * (asm + gestation) * exp(result) - gestation)
    }
  } else {
    stop("Wrong Clock, Parameter Should Be 1-3")
  }
}

agearchaic1 <- data.frame(Clock = rep(paste("Clock", 1:3, sep = ""), each = 4),
  Age = NA, Type = "Neanderthal", Convert = 1:4)

agearchaic1$Age[1:4] <- sapply(panmcusenc1[, -1:-3], function(x) {
  panmcuse(panmcusenc1$Coefficient_1, x, 1)
})
agearchaic1$Age[5:8] <- sapply(panmcusenc2[, -1:-3], function(x) {
  panmcuse(panmcusenc2$Coefficient_2, x, 2)
})
agearchaic1$Age[9:12] <- sapply(panmcusenc3[, -1:-3], function(x) {
  panmcuse(panmcusenc3$Coefficient_3, x, 3)
})

agearchaic2 <- data.frame(Clock = rep(paste("Clock", 1:3, sep = ""), each = 4),
  Age = NA, Type = "Denisovan", Convert = 1:4)

agearchaic2$Age[1:4] <- sapply(panmcusedc1[, -1:-3], function(x) {
  panmcuse(panmcusedc1$Coefficient_1, x, 1, 69.8)
})
agearchaic2$Age[5:8] <- sapply(panmcusedc2[, -1:-3], function(x) {
  panmcuse(panmcusedc2$Coefficient_2, x, 2, 69.8)
})
agearchaic2$Age[9:12] <- sapply(panmcusedc3[, -1:-3], function(x) {
  panmcuse(panmcusedc3$Coefficient_3, x, 3, 69.8)
})

datagearchaic <- rbind(agearchaic1, agearchaic2)
datagearchaic$Type <- factor(datagearchaic$Type,
  levels = c("Neanderthal", "Denisovan"))

write.table(datagearchaic, "age_pan_mammalian_convert.csv",
  row.names = FALSE, sep = ",")