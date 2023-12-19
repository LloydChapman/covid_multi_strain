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
library(parallel)
library(latex2exp)
library(stringr)
library(ggh4x)

source("R/utils.R")
source("R/date.R")
source("R/plot_timeline.R")
source("R/vaccination.R")
source("R/parameters.R")
source("R/pars.R")
source("parameters/transform.R")
source("R/covid_multi_strain.R")
source("R/run_fitting.R")
source("R/fit.R")
source("R/pmcmc.R")
source("R/chains.R")
source("R/fit_process.R")
source("R/plot_fit.R")
source("R/run_simulations.R")
source("R/simulate.R")
source("R/plot_simulations.R")
source("R/plot_sensitivity_analysis.R")

# Set MCMC run number
runs <- 130:132
# Set counterfactual simulations run number
sim_runs <- 28:30

# # Process input data
# source("R/process_FP_data.R")
# 
# # Create parameters
# assumptions <- "central"
# # assumptions <- "alt_contact_matrix"
# source("R/create_parameters.R")

# Parameter assumptions
assumption <- c("pessimistic","central","optimistic")

# Load data
data_raw <- fread("data/data_cases_hosps_deaths_serology.csv")
vax <- fread("data/data_vaccination.csv", colClasses = c(number = "numeric"))
pop <- fread("data/population.csv")
hosps_by_age <- fread("data/hosps_by_age.csv")
deaths_by_age <- fread("data/deaths_by_age.csv")

# Set parameters to fit
# 1-5: transmission rate parameters (beta1-beta5)
# 6: relative transmissibility of Delta vs wildtype (rel_strain_transmission)
# 7: start date of outbreak (start_date)
# 8: Delta seeding date (strain_seed_date)
# 9: maximum probability of hospitalisation given symptomatic infection across all age groups (p_H_max)
# 10: maximum probability of death given hospitalisation across all age groups on and before 2021-06-11 (p_D)
# 11: maximum probability of death given hospitalisation across all age groups on 2021-08-15 (p_D_2)
# 12: maximum probability of death given hospitalisation across all age groups on and after 2021-11-01 (p_D_3)
# 13: relative transmissibility of Omicron BA.1/BA.2 vs Delta (rel_strain_transmission1)
# 14: Omicron BA.1 seeding date (strain_seed_date1)
# 15: relative risk of severe disease for Delta vs wildtype (strain_rel_p_hosp_if_sympt)
# 16: Symptomatic case reporting rate (phi_cases)
# 17: Overdispersion parameter for negative binomial observation process for cases (alpha_cases)
# 18: Overdispersion parameter for negative binomial observation process for hospitalisations (alpha_hosp)
# 19: Overdispersion parameter for negative binomial observation process for hospital deaths (alpha_death)
u <- c(1:5,7:9,10:12,14,15:19)
# Set number of MCMC iterations
n_iters <- 5e4
# Set number of chains to run
n_chains <- 4
# Set flag for whether to use deterministic model or not
deterministic <- T
# Set flag for whether to return variables needed for calculating Rt in "state" object
Rt <- T
# Set flag for whether to use "fixed" initial parameter values for MCMC chains 
# or not (i.e. values slightly perturbed from initial values in parameters/info.csv)
fixed_initial <- T
# Set factor by which to thin MCMC samples
thinning <- 10
# Set burn-in for post processing
burnin <- 4000
# Set number of posterior samples for age-decomposition plots
n_smpls <- 1000
# Set number of samples for counterfactual simulations
n_sims <- 500

# Create output directory
dir.create("output")
    
for (j in seq_along(runs)){
    # Run fitting
    run_fitting(runs[j],assumption[j],data_raw,deterministic,Rt,pop,u,n_iters,n_chains,fixed_initial,thinning,burnin,n_smpls)
    
    # Run counterfactual simulations
    run_simulations(runs[j],sim_runs[j],assumption[j],deterministic,Rt,n_sims)
}

# Plot timeline of epidemic in French Polynesia
plot_timeline(data_raw,pop,hosps_by_age,deaths_by_age)
ggsave("output/timeline_plot.pdf",width = 14,height = 6)

# Plot vaccination coverage by age
base <- readRDS("parameters/base.rds")
vaccination_coverage_plot(base$vaccine_schedule,base$age_groups,vax,pop)
ggsave("output/vax_cov_by_dose.pdf",width = 7,height = 2.33)

# Plot sensitivity analysis results
plot_sensitivity_analysis(sim_runs,assumption)
