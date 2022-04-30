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
# source("R/cpp11.R")
# source("R/dust.R")
# source("R/seirhdagevaxmultistrainserotimedepbeta.R")
source("R/covid_multi_strain.R")
source("R/fit_covid_multi_strain.R")
source("R/pmcmc.R")

# Process FP data
source("R/process_FP_data.R")

# Fit covid_multi_strain to FP data
u <- c(1:4,6:8) #c(1,6:7) #c(1:4,6:9) #1:9 # all parameters 1:5 # only update beta parameters
n_iters <- 2e4
run <- 36
fit_covid_multi_strain(u,n_iters,run)