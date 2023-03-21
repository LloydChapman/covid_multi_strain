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
library(GGally)
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
source("R/process_fit.R")

# Process FP data
source("R/process_FP_data.R")

# Fit covid_multi_strain to FP data
# u <- 1:9 # all parameters
# u <- 1:4 # only beta parameters
# u <- c(1:4,6:9) # beta parameters, seed date, strain seed date, IHR scaling, IFR scaling
# u <- c(1:6,8:10,13) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, observation parameters for case data
# u <- c(1:5,7:9,12) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, observation parameters for case data
u <- c(1:5,7:9,12:13) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, observation parameters for case data
# u <- c(1:10,12:14,17) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, observation parameters for case data
# u <- c(1:7,9:11,14) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, observation parameters for case data
n_iters <- 5e4 #1e3 #2e4 #1e4 #
run <- 76
deterministic <- T # flag for whether to use "deterministic particle filter" or not
Rt <- T #F # flag for whether to return variables needed for calculating Rt in "state"
thinning <- 10
fit_covid_multi_strain(u,n_iters,run,deterministic,Rt,thinning)

# Set burn-in
burnin <- 1500

# Set whether to plot moving average of data
moving_avg <- F

# Set output to use
output <- paste0("output/MCMCoutput",run,".RData")

# Plot fit
plot_fit(output,run,burnin,moving_avg)

# Process fit
pars_qntls <- calculate_parameter_quantiles(output,burnin = burnin)
write.csv(pars_qntls,paste0("output/parameter_quantiles",run,".csv"), row.names = F)
