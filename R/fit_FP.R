library(odin)
library(odin.dust)
library(dust)
library(mcstate)
library(socialmixr)
library(coda)
library(scales)
library(qs)
library(data.table)
library(MASS)
library(abind)

source("R/utils.R")
source("R/vaccination.R")
# source("R/seirhdagevaxmultistrainserotimedepbeta.R")
source("R/covid_multi_strain.R")
source("R/fit_covid_multi_strain.R")
source("R/pmcmc.R")

# Process FP data
source("R/process_FP_data.R")

# Fit covid_multi_strain to FP data
n_iters <- 100
run <- 3
fit_covid_multi_strain(n_iters,run)