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
source("R/fit_FP.R")
source("R/fit.R")
source("R/pmcmc.R")
source("R/chains.R")
source("R/fit_process.R")
source("R/plot_fit.R")
source("R/run_simulations.R")
source("R/simulate.R")
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

# Plot timeline of epidemic in French Polynesia
plot_timeline()
ggsave("output/timeline_plot.pdf",width = 14,height = 6)

# Parameter assumptions
assumption <- c("pessimistic","central","optimistic")

u <- c(1:5,7:9,10:12,14,15:19) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, reporting rate for confirmed cases
n_iters <- 5e4 #3e4 #1e4 #3e4 #100 #
n_chains <- 4 #2 #1 #
    
for (j in seq_along(runs)){
    # Run fitting
    run_fitting(runs[j],assumption[j],u,n_iters,n_chains)
    
    # Run counterfactual simulations
    run_simulations(runs[j],sim_runs[j],assumption[j])
}

# Plot sensitivity analysis results
plot_sensitivity_analysis(sim_runs,assumption)
