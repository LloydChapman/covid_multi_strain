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
source("R/covid_multi_strain.R")
source("R/simulate.R")
source("R/fit_covid_multi_strain.R")
source("R/pmcmc.R")
source("R/plot_fit.R")
source("R/process_fit.R")

# # Process FP data
# source("R/process_FP_data.R")

# Load data
data_raw <- fread("data/data_cases_hosps_deaths_serology.csv")
vax <- fread("data/data_vaccination.csv", colClasses = c(number = "numeric"))
pop <- fread("data/population.csv")

# Set end date for data
end_date <- as.Date("2022-05-06") # last death date in data files

# Set age groups
age_groups <- c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+")

# Create vaccination schedule
# Set delays for immune response to different vaccine doses
delay_dose1 <- 28
delay_dose2 <- 14

# Get vaccination age groups
age_groups_vax <- vax[,unique(age_group)]

# Matrix of uptake rates (age group x dose)
uptake <- matrix(1,nrow = length(age_groups),ncol = vax[,length(unique(dose))])

# Make vaccine schedule
schedule <- vaccination_data(vax,delay_dose1,delay_dose2,pop,age_groups_vax,
                             age_groups,end_date,uptake)

# Plot vaccination coverage by age
vaccination_coverage_plot(schedule,age_groups,vax,pop)
ggsave("output/vax_cov_by_dose.pdf",width = 9,height = 3)

# Fit covid_multi_strain to FP data
# u <- 1:9 # all parameters
# u <- 1:4 # only beta parameters
# u <- c(1:4,6:9) # beta parameters, seed date, strain seed date, IHR scaling, IFR scaling
# u <- c(1:6,8:10,13) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, observation parameters for case data
# u <- c(1:5,7:9,12) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, observation parameters for case data
u <- c(1:5,7:9,12:13) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, observation parameters for case data
# u <- c(1:10,12:14,17) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, observation parameters for case data
# u <- c(1:7,9:11,14) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, observation parameters for case data
n_iters <- 5e4 #1e3 #2e4 #1e4 #200 #
run <- 77 #75 #76 #
deterministic <- T # flag for whether to use "deterministic particle filter" or not
Rt <- T #F # flag for whether to return variables needed for calculating Rt in "state"
thinning <- 10
fit_covid_multi_strain(data_raw,schedule,pop,age_groups,u,n_iters,run,deterministic,Rt,thinning)

# Post processing
# Set output to use
output <- paste0("output/MCMCoutput",run,".RData")

# Set burn-in
burnin <- 1500 #1 #

# Set whether to plot moving average of data
moving_avg <- F

# Set number of posterior samples for age-decomposition plots
n_smpls <- 1000 #10 #

# Set seed for immune status plot
seed <- 1L

# Plot fit
plot_fit(output,run,pop,burnin,moving_avg,n_smpls)

# Process fit
pars_qntls <- calculate_parameter_quantiles(output,burnin = burnin)
write.csv(pars_qntls,paste0("output/parameter_quantiles",run,".csv"), row.names = F)
