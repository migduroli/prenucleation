#install.packages("ConfigParser")

library(ConfigParser)

config <- ConfigParser$new()
config$read("config.ini")

d <-  config$getfloat("hard_sphere", NA, "thermo")
n_points <- config$getfloat("points", NA, "path")
r_min <- config$getfloat("r_min", NA, "path")

kT <- config$getfloat("kT", NA, "bulk")
b_beta <- 1 / kT
b_supsat <- config$getfloat("supsat", NA, "bulk")