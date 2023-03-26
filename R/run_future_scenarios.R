library(odin)
library(odin.dust)
library(dust)
library(socialmixr)
library(qs)
library(data.table)
library(readxl)
library(ggplot2)


source("R/utils.R")
source("R/date.R")
source("R/vaccination.R")
source("R/covid_multi_strain.R")
source("R/process_fit.R")
source("R/simulate.R")

## Process FP data
source("R/process_FP_data.R")

# Load MCMC output
# TODO: Save only samples of fitted parameters, probabilities, final model state, 
# trajectories, and predict and info objects in MCMC output (see process_fit()) 
# not full workspace inside fit_covid_multi_strain(), and just load those here
load("output/MCMCoutput42.RData")

# Create dust model generator
covid_multi_strain <- odin_dust("inst/odin/covid_multi_strain.R")

# Date of end of fit
date <- end_date #"2021-11-21" # TODO: check if this should be following day
# Number of parameter samples to use for simulations
n_smpls <- 10
# End date for simulations
end_date <- "2022-01-01"
# Number of threads to use for simulations
n_threads <- 4

# # Run MCMC
# samples <- mcmc(...)

# Process MCMC output
dat <- process_fit(res,init_pars,idx,n_steps,dt,transform,index,filter,date,
                   beta_date,mean_days_between_doses,schedule,data)
saveRDS(dat,"output/fit.rds")

# Load MCMC output and process for forward simulations
dat <- load_fit("output/fit.rds")

# Get onward object for forward simulations
onward <- dat$onward

# Prepare input for simulations
onward1 <- simulate_prepare(onward,n_smpls)

# Create list of simulation parameters
# N.B. Do this just as a list here for now and look at setting up scenario grids
# as in spimalot later (see 
# https://github.com/mrc-ide/spimalot/blob/main/R/simulate.R#L79-L137
# and
# https://github.com/mrc-ide/sarscov2-roadmap-england/blob/main/src/vaccine_simulation/script.R#L249-L360)
# TODO: Check start date for booster schedule - is it too early by 1 day?
vaccine_booster_daily_doses <- rep(5000/7,as.integer(as.Date(end_date) - as.Date(date)) + 1)
names(vaccine_booster_daily_doses) <- as.character(seq(as.Date(date), as.Date(end_date), by = 1))
args <- list(end_date = as.Date(end_date),
             seed = NULL,
             n_threads = n_threads,
             output_keep = names(idx$state), #c("hosps","deaths","sero_pos_1"),
             output_time_series = TRUE,
             output_vaccination = FALSE, # N.B. Not implemented yet, so set to FALSE
             vaccine_uptake = rep(0.9,8),
             vaccine_eligibility = rep(1,8),
             vaccine_delay_multiplier = 1,
             vaccine_daily_dose = 0,
             vaccine_lag_groups = NULL,
             vaccine_lag_days = NULL,
             vaccine_booster_daily_doses = vaccine_booster_daily_doses,
             vaccine_booster_eligibility = rep(1,8),
             strain_transmission = NULL,
             strain_seed_date = NULL,
             strain_seed_size = NULL,
             strain_seed_pattern = NULL,
             strain_cross_immunity = NULL,
             waning_rate = NULL)

# Simulate future scenarios
ret <- simulate_future_scenario(args,onward1)
