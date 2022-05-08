library(odin)
library(odin.dust)
library(dust)
library(mcstate)
library(socialmixr)
library(coda)
library(scales)
library(qs)
library(data.table)
library(readxl)
library(ggplot2)
library(MASS)
library(abind)

source("R/utils.R")
source("R/vaccination.R")
source("R/covid_multi_strain.R")
source("R/fit_covid_multi_strain.R")
source("R/pmcmc.R")

# Process FP data
source("R/process_FP_data.R")

# Fit covid_multi_strain to FP data
# u <- 1:9 # all parameters
# u <- c(1:4,6:9) # beta parameters, seed date, strain seed date, IHR scaling, IFR scaling
u <- c(1:4,6:8) # beta parameters, seed date, strain seed date, IHR scaling
# u <- 1:4 # only beta parameters
n_iters <- 1e3 #1e3 #2e4
run <- 37
deterministic <- T # flag for whether to use "deterministic particle filter" or not
thinning <- 10
fit_covid_multi_strain(u,n_iters,run,deterministic,thinning)

