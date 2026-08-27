#Convert published thresholds in umol/kg into kPa
library(respR)

temp <- 12
pcrit <- c(4, 44)
in_units <- "umol/kg"

convert_DO(x = pcrit, t = pcrit, S = 32, from = in_units, to = "kPa")


temp <- 10
pcrit <- c(52, 58)
in_units <- "umol/kg"

convert_DO(x = pcrit, t = pcrit, S = 32, from = in_units, to = "kPa")

temp <- 8
pcrit <- 60
in_units <- "umol/kg"
convert_DO(x = pcrit, t = pcrit, S = 32, from = in_units, to = "kPa")