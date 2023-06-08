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
library(cowplot)
library(MASS)
library(abind)

source("R/utils.R")
source("R/date.R")
source("R/vaccination.R")
source("R/parameters.R")
source("R/pars.R")
source("R/covid_multi_strain.R")
source("R/simulate.R")
source("R/fit.R")
source("R/pmcmc.R")
source("R/plot_fit.R")
source("R/fit_process.R")
source("parameters/transform.R")

# Load data
data_raw <- fread("data/data_cases_hosps_deaths_serology.csv")
# data_raw <- data_raw[day < covid_multi_strain_date(as.Date("2021-06-01")),]
# cols <- setdiff(names(data_raw),"day")
# data_raw <- data_raw[!(day <= 213 | day >= covid_multi_strain_date(as.Date("2020-09-15"))),(cols):= NA]
vax <- fread("data/data_vaccination.csv", colClasses = c(number = "numeric"))
pop <- fread("data/population.csv")

# Set assumption for booster waning rate
assumptions <- "central" #-log(67.7/82.8)/(105-25) # (Stowe Nat Comm 2022 Table S11)
# assumptions <- "optimistic" # -log(0.923)/140 (Barnard Nat Com 2022 Table S4)

## Load parameters
# Output pars is a list containing:
# info - the loaded info.csv
# prior - the loaded prior.csv
# proposal - the loaded proposal.csv
# transform - the functions in transform.R
# raw - exact output of first 3 csvs, without any treatment
# base - the base.rds object, which contains the fixed parameters
# mcmc - initialisation object built from the above to pass to the mcmc
pars <- fit_pars_load("parameters",assumptions)

# Plot vaccination coverage by age
# vaccination_coverage_plot(pars$base$vaccine_schedule,pars$base$age_groups,vax,pop)
dir.create("output")
# ggsave("output/vax_cov_by_dose.pdf",width = 9,height = 3)

# Fit covid_multi_strain to FP data
u <- c(1:6,8:10,11:13,15,16,17:21) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, reporting rate for confirmed cases
# u <- c(1:3,8,10,13,16:19)
# u <- c(2:3,8,10,13,16:19)
n_iters <- 5e4 #4e4 #2e4 #1e4 #3e4 #
# Change run number for different assumption on booster waning rate
# run <- 77
# run <- 78
run <- 101
deterministic <- T # flag for whether to use deterministic model or not
Rt <- T #F # flag for whether to return variables needed for calculating Rt in "state" object
initial_date <- pars$info$min[pars$info$name == "start_date"] - 1
    
# Construct particle filter
filter <- covid_multi_strain_particle_filter(data_raw,pars,deterministic,Rt,initial_date)

# Run fitting
thinning <- 10
samples <- fit_run(pars,filter,u,n_iters,deterministic,Rt,thinning)
saveRDS(samples,paste0("output/MCMCsamples",run,".RDS"))

## Post processing
# Set burn-in
burnin <- 3000 #1500 #500 #

# Process MCMC output
dat <- fit_process(samples,pars,data_raw,filter,burnin,simulate_object = T)

# Save output
saveRDS(dat,paste0("output/MCMCoutput",run,".RDS"))


# Set whether to plot moving average of data
moving_avg <- F

# Set number of posterior samples for age-decomposition plots
n_smpls <- 1000 #500 #

# Plot fit
plot_fit(dat,run,pop,u,moving_avg,n_smpls)

# Process fit
pars_qntls <- calculate_parameter_quantiles(dat)
write.csv(pars_qntls,paste0("output/parameter_quantiles",run,".csv"), row.names = F)
