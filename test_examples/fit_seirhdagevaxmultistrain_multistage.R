library(odin)
library(odin.dust)
library(dust)
library(mcstate)
library(socialmixr)
library(coda)
library(scales)
library(qs)
library(data.table)

source("R/utils.R")
source("R/vaccination.R")
source("R/seirhdagevaxmultistrain.R")



#### Set up model and parameters ####


# Age-structured SEIR model with deaths and vaccination
gen_seirhd_age_vax_multistrain <- odin_dust("test_examples/seirhdagevaxmultistrain.R")

# Get age-structured contact matrix
data(polymod, package = "socialmixr")

age.limits <- seq(0,70,10)
contact <- socialmixr::contact_matrix(survey = polymod,countries = "United Kingdom",age.limits = age.limits,symmetric = T)

## Transform the matrix to the (symmetrical) transmission matrix
## rather than the contact matrix. This transmission matrix is
## weighted by the population in each age band.
# population <- contact$demography$population
population <- round(contact$demography$population/1000)
transmission <- contact$matrix/rep(population, each = ncol(contact$matrix))
# for (i in 1:nrow(transmission)){
#   transmission[i,i] <- transmission[2,2]  
# }
transmission

# Get age-dependent symptomatic fraction, IHR and IFR
sympt_frac <- qread("~/covidm/newcovid3/fitting_data/2-linelist_both_fit_fIa0.5-rbzvih.qs")
p_C <- unname(colMeans(sympt_frac[,13:20]))
prob_death_given_hosp <- fread("data/prob_death_given_hosp_salje.csv")
p_D <- prob_death_given_hosp[,median_perc_mean/100]
IHR <- fread("../covid_remaining_burden/data/ihr_salje.csv")
ihr <- IHR[,median_perc_mean/100]
IFR_salje <- fread("data/ifr_salje.csv")
ifr_salje <- IFR_salje[,median_perc_mean/100]
# # O'Driscoll IFR
# IFR <- fread("data/ifr_odriscoll.csv")
# get_min_age = function(x) as.numeric(sub("-.*","",sub("\\+|<","-",x)))
# IFR[,age_low:=get_min_age(IFR$Age_group)]
# IFR[,age_group:=cut(age_low,c(age.limits,Inf),labels=sub(",","-",gsub("\\[|)","",contact$demography$age.group)),right=F)]
# agg_IFR <- IFR[,lapply(.SD,mean),.SDcols=setdiff(names(IFR),c("Age_group","age_low","age_group")),by=.(age_group)]
# ifr <- agg_IFR[,Median_perc_mean/100]
# For the time being just assume that IFR including non-hospital deaths is 10% higher
ifr <- 1.1*ifr_salje
p_G <- (ifr - p_D*ihr)/((1-p_D)*ihr + ifr)
p_H <- ihr/(p_C*(1-p_G))

# IFR <- readRDS("~/UCSF/COVIDVaccineModelling/Data/IFR_by_age_ODriscoll.RDS")
# p_death <- IFR$median_perc[19:27]/100
# p_death[length(p_death)-1] <- (p_death[length(p_death)-1]+p_death[length(p_death)])/2
# p_death <- p_death[1:(length(p_death)-1)]

# Generate model
dt <- 0.25
n_age <- length(age.limits)
n_vax <- 3

# Transmission and natural history parameters
gamma_E <- 0.5
gamma_A <- 0.2
gamma_H <- 0.1
gamma_G <- 1/3
theta_A <- 0.5

# relative transmissibilities of 1st and 2nd strains
strain_transmission <- c(1,2)

# Seeding parameters
# 1st strain
# start_date <- 1L
initial_seed_size <- 10
initial_seed_pattern <- 1

# 2nd strain
# strain_seed_date <- 70
strain_seed_size <- 100
strain_seed_pattern <- 1

# Vaccination parameters
vaccine_progression_rate <- c(0,0,1/(26*7))

daily_doses <- rep(0,365) #rep(1000,365) #rep(5,365)
mean_days_between_doses <- 12*7
# index_dose <- c(1L,2L)
# 
# # Build vaccination progression rate matrix
# p$vaccine_progression_rate_base <- build_vaccine_progression_rate(
#     vaccine_progression_rate, n_age, n_vax, index_dose)

# Make example vaccine schedule
pop_mat <- matrix(rep(population,1),nrow = length(population))
schedule <- vaccine_schedule_future(30, daily_doses, mean_days_between_doses, pop_mat)

vaccine_index_dose2 <- 2L
vaccine_catchup_fraction <- 1
n_doses <- 2L

# Relative probabilities of symptoms, hospitalisation and death for different strains
strain_rel_p_sympt <- 1
strain_rel_p_hosp_if_sympt <- 1
strain_rel_p_death <- 1

# Parameters for impact of vaccination of susceptibility and infectiousness
rel_susceptibility <- c(1,0.8,0.5) # relative susceptibility to infection in each vaccine stratum
rel_p_sympt <- c(1,0.6,0.3) # relative risk of symptoms in each vaccine stratum
rel_p_hosp_if_sympt <- c(1,0.95,0.95) # relative risk of hospitalisation given infection in each vaccine stratum
rel_p_death <- c(1,0.9,0.9) # relative risk of death in each vaccine stratum
rel_infectivity <- c(1,0.5,0.5) # relative infectiousness of infected individuals in each vaccine stratum

# Waning parameters
waning_rate <- 1/365

# Cross immunity parameters
cross_immunity <- 0.5 # 1 #

# Construct parameters object
p <- parameters(dt,
                n_age,
                n_vax,
                transmission,
                beta = 0.08,
                gamma_E,
                gamma_P = 0.4,
                gamma_A,
                gamma_C = 0.4,
                gamma_H,
                gamma_G,
                theta_A,
                p_C,
                p_H,
                p_G,
                p_D,
                population,
                start_date = 1L,
                initial_seed_size,
                initial_seed_pattern,
                strain_transmission,
                strain_seed_date = 70L,
                strain_seed_size,
                strain_seed_pattern,
                strain_rel_p_sympt,
                strain_rel_p_hosp_if_sympt,
                strain_rel_p_death,
                rel_susceptibility,
                rel_p_sympt,
                rel_p_hosp_if_sympt,
                rel_p_death,
                rel_infectivity,
                vaccine_progression_rate,
                schedule,
                vaccine_index_dose2,
                vaccine_catchup_fraction = vaccine_catchup_fraction,
                n_doses = n_doses,
                waning_rate = waning_rate,
                cross_immunity = cross_immunity)

# Create instance of model
seirhd_age_vax_multistrain <- gen_seirhd_age_vax_multistrain$new(p,step = 0,n_particles = 1,n_threads = 1,seed = 1)

info <- seirhd_age_vax_multistrain$info()
initial_state <- initial(info, NULL, p)

seirhd_age_vax_multistrain$update_state(state = initial_state)

# Run epidemic forward
n_steps <- 800



#### Create data to be fitted to ####


# Create an array to contain outputs after looping the model.
# Array contains XX rows = Total S, E, I, R (4), and
# in each age compartment (XX) as well as the cumulative incidence (XX)
x <- array(NA, dim = c(info$len, 1, n_steps+1))

# For loop to run the model iteratively
x[ , ,1] <- seirhd_age_vax_multistrain$state()
for (t in seq_len(n_steps)) {
    # print(t)
    x[ , ,t+1] <- seirhd_age_vax_multistrain$run(t)
}
time <- x[1,1,-1]

# Drop time row
x <- x[-1, , ,drop=FALSE]

n_strains <- 4

# Plot trajectories
plot_trajectories(time,x,n_age,n_strains,n_vax)


# Rotate strains so that strain 1 becomes the sum of strains 1 & 2 and strain 2
# is empty
state1 <- seirhd_age_vax_multistrain$state()
state2 <- rotate_strains(state1,info)

# Parameters for strain 3
# relative transmissibility of strain 3 is 4 times that of strain 1
strain_transmission1 <- c(1,4)
# lower cross-immunity from previous infection with strains 1 or 2, but can't be
# infected with strains 1 or 2 after being infected with strain 3
cross_immunity1 <- c(0.2,1) 
p1 <- parameters(dt,
                 n_age,
                 n_vax,
                 transmission,
                 beta = 0.08,
                 gamma_E,
                 gamma_P = 0.4,
                 gamma_A,
                 gamma_C = 0.4,
                 gamma_H,
                 gamma_G,
                 theta_A,
                 p_C,
                 p_H,
                 p_G,
                 p_D,
                 population,
                 start_date = n_steps*dt,
                 initial_seed_size = 0,
                 initial_seed_pattern,
                 strain_transmission1,
                 strain_seed_date = 210L,
                 strain_seed_size,
                 strain_seed_pattern,
                 strain_rel_p_sympt,
                 strain_rel_p_hosp_if_sympt,
                 strain_rel_p_death,
                 rel_susceptibility,
                 rel_p_sympt,
                 rel_p_hosp_if_sympt,
                 rel_p_death,
                 rel_infectivity,
                 vaccine_progression_rate,
                 schedule,
                 vaccine_index_dose2,
                 vaccine_catchup_fraction = vaccine_catchup_fraction,
                 n_doses = n_doses,
                 waning_rate = waning_rate,
                 cross_immunity = cross_immunity1)

# Update model parameters and state
# seirhd_age_vax_multistrain1 <- gen_seirhd_age_vax_multistrain$new(p1,step = 0,n_particles = 1,n_threads = 1,seed = 1)
# 
# # Initialise with final state from first epoch
# seirhd_age_vax_multistrain1$update_state(state = state2)
seirhd_age_vax_multistrain$update_state(pars = p1,state = state2)

n_steps1 <- 1200 # number of steps to run model up to

# Create data to be fitted to
# Create an array to contain outputs after looping the model.
# Array contains XX rows = Total S, E, I, R (4), and
# in each age compartment (XX) as well as the cumulative incidence (XX)
x1 <- array(NA, dim = c(info$len, 1, n_steps1-n_steps+1))

# For loop to run the model iteratively
x1[ , ,1] <- seirhd_age_vax_multistrain$state()
for (t in seq_len(n_steps1-n_steps)) {
    # print(t)
    x1[ , ,t+1] <- seirhd_age_vax_multistrain$run(n_steps+t)
}
time1 <- x1[1,1,-1]

# Drop time row
x1 <- x1[-1, , ,drop=FALSE]

# Join with first epoch
time2 <- c(time,time1)
x2 <- array_bind(x,x1[,,-1,drop=FALSE])

# Plot trajectories
plot_trajectories(time2,x2,n_age,n_strains,n_vax)

# Extract true history of model states
true_history <- x2[ , ,seq(0,n_steps1+1,by=1/dt)+1,drop=F]


index(seirhd_age_vax_multistrain$info())

# Add noise to simulated data
hosps <- true_history[3:7, ,-1]
deaths <- true_history[8:12, ,-1]

par(mfrow = c(1,1))
days <- seq(1,n_steps1*dt)
matplot(days,t(hosps),type="l",xlab="Day",ylab="Hospitalisations")
matplot(days,t(deaths),type="l",xlab="Day",ylab="Deaths")

add_noise <- function(x,f){
    noise <- apply(x,2,function(y) round(rnorm(length(y),0,f*y)))
    x <- x + noise
    x <- pmax(x,0)
}

set.seed(0)
hosps <- add_noise(hosps,0.2)
rownames(hosps) <- paste0("hosps_",age.limits[c(1,5:length(age.limits))],"_",c(as.character(age.limits[5:length(age.limits)]-1),"plus"))
deaths <- add_noise(deaths,0.3)
rownames(deaths) <- paste0("deaths_",age.limits[c(1,5:length(age.limits))],"_",c(as.character(age.limits[5:length(age.limits)]-1),"plus"))

matplot(days,t(hosps),type="l",xlab="Day",ylab="Hospitalisations")
matplot(days,t(deaths),type="l",xlab="Day",ylab="Deaths")

# Create "observed" data
data_raw <- data.frame(t(rbind(hosps,deaths)))
data_raw$day <- days
# Add empty columns for total cases and deaths
data_raw$hosps <- NA
data_raw$deaths <- NA
# Convert to required format
data <- particle_filter_data(data_raw,"day",1/dt)

# plot(data_raw$day,data$hosps_70_plus,type="l")



#### Fit to multiple age-stratified data streams ####


# Create multistage parameters object
epochs <- list(
    multistage_epoch(n_steps*dt, pars = p1, transform_state = transform_state))
pars <- multistage_parameters(p,epochs)

# Create particle filter object
n_particles <- 200
filter <- particle_filter$new(data, gen_seirhd_age_vax_multistrain, n_particles, 
                              compare, index, initial)
filter$run(
    pars = pars,
    save_history = TRUE)

# Check variable names in particle filter history
dimnames(filter$history())

# Plot filtered trajectories
plot_particle_filter(filter$history(),true_history,data_raw$day)
plot_hosps_and_deaths_age(filter$history(),data,data_raw$day,n_age,n_vax,n_strains)

# Infer parameters by pMCMC
beta <- pmcmc_parameter("beta",0.08,min = 0,
                        prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
# R0 <- pmcmc_parameter("R0",2,min = 0)
gamma <- pmcmc_parameter("gamma",0.4,min = 0,
                         prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
start_date <- pmcmc_parameter("start_date",1L,min = 1L,max = 10L,discrete = TRUE)
strain_seed_date <- pmcmc_parameter("strain_seed_date",65L,min = 60L,max = 80L,discrete = TRUE)
strain_seed_date1 <- pmcmc_parameter("strain_seed_date1",200L,min = 200L,max = 250L,discrete = TRUE)

# proposal <- matrix(c(0.01^2,0,0,0.01^2),nrow = 2,ncol = 2,byrow = TRUE)
# proposal <- matrix(c(0.01^2,0,0,0,0.01^2,0,0,0,2),nrow = 3,ncol = 3,byrow = TRUE)
# proposal <- matrix(c(0.01^2,0,0,0,
#                      0,0.01^2,0,0,
#                      0,0,2,0,
#                      0,0,0,2),nrow = 4,ncol = 4,byrow = TRUE)
proposal <- matrix(c(0.01^2,0,0,0,0,
                     0,0.01^2,0,0,0,
                     0,0,2,0,0,
                     0,0,0,2,0,
                     0,0,0,0,2),nrow = 5,ncol = 5,byrow = TRUE)
# transform <- make_transform(dt,
#                             n_age,
#                             n_vax,
#                             transmission,
#                             gamma_E,
#                             gamma_A,
#                             gamma_H,
#                             gamma_G,
#                             theta_A,
#                             p_C,
#                             p_H,
#                             p_G,
#                             p_D,
#                             population,
#                             # start_date,
#                             initial_seed_size,
#                             initial_seed_pattern,
#                             strain_transmission,
#                             # strain_seed_date,
#                             strain_seed_size,
#                             strain_seed_pattern,
#                             strain_rel_p_sympt,
#                             strain_rel_p_hosp_if_sympt,
#                             strain_rel_p_death,
#                             rel_susceptibility,
#                             rel_p_sympt,
#                             rel_p_hosp_if_sympt,
#                             rel_p_death,
#                             rel_infectivity,
#                             vaccine_progression_rate,
#                             schedule,
#                             vaccine_index_dose2,
#                             vaccine_catchup_fraction,
#                             n_doses,
#                             waning_rate,
#                             cross_immunity)
transform <- make_transform_multistage(dt,
                                       n_age,
                                       n_vax,
                                       transmission,
                                       gamma_E,
                                       gamma_A,
                                       gamma_H,
                                       gamma_G,
                                       theta_A,
                                       p_C,
                                       p_H,
                                       p_G,
                                       p_D,
                                       population,
                                       # start_date,
                                       initial_seed_size,
                                       initial_seed_pattern,
                                       strain_transmission,
                                       # strain_seed_date,
                                       strain_seed_size,
                                       strain_seed_pattern,
                                       strain_rel_p_sympt,
                                       strain_rel_p_hosp_if_sympt,
                                       strain_rel_p_death,
                                       rel_susceptibility,
                                       rel_p_sympt,
                                       rel_p_hosp_if_sympt,
                                       rel_p_death,
                                       rel_infectivity,
                                       vaccine_progression_rate,
                                       schedule,
                                       vaccine_index_dose2,
                                       vaccine_catchup_fraction,
                                       n_doses,
                                       waning_rate,
                                       cross_immunity,
                                       start_date1 = n_steps * dt,
                                       strain_transmission1,
                                       cross_immunity1)
mcmc_pars <- pmcmc_parameters$new(list(beta = beta,gamma = gamma,
                                       start_date = start_date,
                                       strain_seed_date = strain_seed_date,
                                       strain_seed_date1 = strain_seed_date1),
                                  proposal,transform = transform)
# mcmc_pars <- pmcmc_parameters$new(list(R0 = R0,gamma = gamma),
#                                   proposal,transform = transform)

# Run MCMC
n_samples <- 10 #1000
control <- pmcmc_control(
    n_samples,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE)
pmcmc_run <- pmcmc(mcmc_pars,filter,control = control)


dimnames(pmcmc_run$trajectories$state)

plot_particle_filter(pmcmc_run$trajectories$state,true_history,data_raw$day)
plot_hosps_and_deaths_age(pmcmc_run$trajectories$state,data,data_raw$day,n_age,n_vax,n_strains)

# # Plot MCMC output
# mcmc_out <- as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars))
# summary(mcmc_out)
# plot(mcmc_out)
# 
# # Calculate effective sample size
# effectiveSize(mcmc_out)
# 1 - rejectionRate(mcmc_out)
# 
# # Pairwise correlation plot
# par(mfrow = c(1,1))
# pairs(pmcmc_run$pars[])
# 
# # Tune MCMC
# beta <- pmcmc_parameter("beta",mean(pmcmc_run$pars[(round(n_samples/5):n_samples)+1,1]),min = 0,
#                         prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
# gamma <- pmcmc_parameter("gamma",mean(pmcmc_run$pars[(round(n_samples/5):n_samples)+1,2]),min = 0,
#                          prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
# # beta <- pmcmc_parameter("beta",0.2,min = 0,
# #                         prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
# # gamma <- pmcmc_parameter("gamma",0.5,min = 0,
# #                          prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
# proposal1 <- cov(pmcmc_run$pars)
# mcmc_pars1 <- pmcmc_parameters$new(list(beta = beta,gamma = gamma),
#                                    proposal1,transform = transform)
# 
# pmcmc_run1 <- pmcmc(mcmc_pars1,filter,control = control)
# 
# plot_particle_filter(pmcmc_run$trajectories$state[,101:1001,],true_history,data_raw$day)
# plot_hosps_and_deaths_age(pmcmc_run1$trajectories$state[,101:1001,],data,data_raw$day,n_age,n_vax)
# 
# # Plot MCMC output
# mcmc_out1 <- as.mcmc(cbind(pmcmc_run1$probabilities, pmcmc_run1$pars))
# summary(mcmc_out1)
# plot(mcmc_out1)
# 
# # Calculate effective sample size
# effectiveSize(mcmc_out1)
# 1 - rejectionRate(mcmc_out1)
# 
# # Pairwise correlation plot
# par(mfrow = c(1,1))
# pairs(pmcmc_run1$pars[])
