#file requirements
#information of target sites in six .csv files, with rows corresponding to individual sites and columns corresponding to chromosome position, coefficients, specific values, and whether directly matched with the methylation map or not, each clock including one file for the Altai Neandertal and the other for the Denisovan

panmcage <- data.frame(Clock = rep(paste("Clock", 1:3, sep = ""), each = 2),
  Age = NA, Type = c("Neanderthal", "Denisovan"), Method = "No Imputations")

panmcuse <- function(data, type = 1, max = 65.7, gestation = 280 / 365, asm = 13) {
  if (type == 1) {
    message("Using Clock 1")
    intercept <- data$Coefficient_1[1]
    result <- sum(data$Coefficient_1[-1] * data$Beta_value[-1], intercept)
    return(exp(result) - 2)
  } else if (type == 2) {
    message("Using Clock 2")
    intercept <- data$Coefficient_2[1]
    result <- sum(data$Coefficient_2[-1] * data$Beta_value[-1], intercept)
    return(exp(-1 * exp(-1 * result)) * (max + gestation) - gestation)
  } else if (type == 3) {
    message("Using Clock 3")
    intercept <- data$Coefficient_3[1]
    result <- sum(data$Coefficient_3[-1] * data$Beta_value[-1], intercept)
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

panmcusen1 <- read.table("pan_mammalian_clock1_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")
unmatched <- which(is.na(panmcusen1$Beta_value))[-1]
panmcusen1$Beta_value[unmatched] <- 0
panmcage$Age[1] <- panmcuse(panmcusen1, type = 1)

panmcusen2 <- read.table("pan_mammalian_clock2_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")
unmatched <- which(is.na(panmcusen2$Beta_value))[-1]
panmcusen2$Beta_value[unmatched] <- 0
panmcage$Age[3] <- panmcuse(panmcusen2, type = 2)

panmcusen3 <- read.table("pan_mammalian_clock3_hg19_Altai_Neandertal.csv",
  header = TRUE, sep = ",")
unmatched <- which(is.na(panmcusen3$Beta_value))[-1]
panmcusen3$Beta_value[unmatched] <- 0
panmcage$Age[5] <- panmcuse(panmcusen3, type = 3)

panmcused1 <- read.table("pan_mammalian_clock1_hg19_Denisovan.csv",
  header = TRUE, sep = ",")
unmatched <- which(is.na(panmcused1$Beta_value))[-1]
panmcused1$Beta_value[unmatched] <- 0
panmcage$Age[2] <- panmcuse(panmcused1, type = 1)

panmcused2 <- read.table("pan_mammalian_clock2_hg19_Denisovan.csv",
  header = TRUE, sep = ",")
unmatched <- which(is.na(panmcused2$Beta_value))[-1]
panmcused2$Beta_value[unmatched] <- 0
panmcage$Age[4] <- panmcuse(panmcused2, type = 2, max = 69.8)

panmcused3 <- read.table("pan_mammalian_clock3_hg19_Denisovan.csv",
  header = TRUE, sep = ",")
unmatched <- which(is.na(panmcused3$Beta_value))[-1]
panmcused3$Beta_value[unmatched] <- 0
panmcage$Age[6] <- panmcuse(panmcused3, type = 3)

write.table(panmcage, "age_pan_mammalian_direct.csv",
  row.names = FALSE, sep = ",")