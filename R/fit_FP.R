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
library(lubridate)
library(ISOweek)
library(ggplot2)
library(MASS)
library(abind)

source("R/utils.R")
source("R/date.R")
source("R/vaccination.R")
source("R/covid_multi_strain.R")
source("R/simulate.R")
source("R/fit_covid_multi_strain.R")
source("R/pmcmc.R")
source("R/plot_fit.R")

# Process FP data
source("R/process_FP_data.R")

# Fit covid_multi_strain to FP data
# u <- 1:9 # all parameters
# u <- c(1:4,6:9) # beta parameters, seed date, strain seed date, IHR scaling, IFR scaling
u <- c(1:6,8:10,13:15) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, observation parameters for case data
# u <- 1:4 # only beta parameters
n_iters <- 1e4 #5e4 #1e3 #2e4 #
run <- 55
deterministic <- T # flag for whether to use "deterministic particle filter" or not
Rt <- T # flag for whether to return variables needed for calculating Rt in "state"
thinning <- 10
fit_covid_multi_strain(u,n_iters,run,deterministic,Rt,thinning)

# Plot fit
plot_fit(paste0("output/MCMCoutput",run,".RData"),run)
