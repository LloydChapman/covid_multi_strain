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

source("R/utils.R")
source("R/vaccination.R")
source("R/seirhdagevaxmultistrainserotimedepbeta.R")
source("R/pmcmc.R")

#### Set up model and parameters ####


# Age-structured SEIR model with deaths and vaccination
gen_seirhd_age_vax_multistrain_sero_time_dep_beta <- odin_dust("test_examples/seirhdagevaxmultistrainserotimedepbeta.R")

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
n_age <- length(age_groups)
n_vax <- 5

# Transmission and natural history parameters
intvtn_dates <- as.Date(c("2020-08-09","2020-08-27","2020-10-24","2021-06-01","2021-08-12"))
beta_date <- as.integer(intvtn_dates - min(intvtn_dates))
beta_value_sim <- 7/8*c(0.04,0.03,0.02,0.04,0.02)
beta_type <- "piecewise-constant"
gamma_E <- 0.5
gamma_A <- 0.2
gamma_H <- 0.1
gamma_G <- 1/3
gamma_pre_1 <- 1/13
gamma_P_1 <- 1/200
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
strain_seed_size <- 10
strain_seed_pattern <- 1

# Vaccination parameters
vaccine_progression_rate <- c(0,0,1/(26*7),0,1/(26*7))


vaccine_index_dose2 <- 2L
vaccine_index_booster <- 4L
vaccine_catchup_fraction <- 1
n_doses <- 3L #2L

# Relative probabilities of symptoms, hospitalisation and death for different strains
strain_rel_p_sympt <- 1
strain_rel_p_hosp_if_sympt <- 1
strain_rel_p_death <- 1

# Parameters for impact of vaccination of susceptibility and infectiousness
rel_susceptibility <- c(1,0.8,0.5,0.8,0.5) # relative susceptibility to infection in each vaccine stratum
rel_p_sympt <- c(1,0.6,0.3,0.6,0.3) # relative risk of symptoms in each vaccine stratum
rel_p_hosp_if_sympt <- c(1,0.95,0.95,0.95,0.95) # relative risk of hospitalisation given infection in each vaccine stratum
rel_p_death <- c(1,0.9,0.9,0.9,0.9) # relative risk of death in each vaccine stratum
rel_infectivity <- c(1,0.5,0.5,0.5,0.5) # relative infectiousness of infected individuals in each vaccine stratum

# Waning parameters
waning_rate <- 1/365

# Cross immunity parameters
cross_immunity <- 0.95 # 1 #

# Sensitivity and specificity of serological tests
sero_sensitivity_1 <- 0.9
sero_specificity_1 <- 0.99

# Construct parameters object for 1st epoch
p <- parameters(dt,
                n_age,
                n_vax,
                transmission,
                beta_date,
                beta_value = beta_value_sim,
                beta_type,
                gamma_E,
                gamma_P = 0.4,
                gamma_A,
                gamma_C = 0.4,
                gamma_H,
                gamma_G,
                gamma_pre_1,
                gamma_P_1,
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
                strain_seed_date = 315L,
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
                vaccine_index_booster,
                vaccine_catchup_fraction,
                n_doses,
                waning_rate,
                cross_immunity,
                sero_sensitivity_1,
                sero_specificity_1)

# Plot p.w. linear beta to check it looks right
beta_t <- seq(0, beta_date[length(beta_date)], by = dt)
par(mfrow = c(1,1))
plot(beta_t, p$beta_step, type="o", cex = 0.25)
points(beta_date, beta_value_sim, pch = 19, col = "red")

# Number of steps in 1st epoch
n_steps <- max(data_raw$day)/dt

# # Parameters for 2nd epoch (in which new strain is introduced)
# # relative transmissibility of strain 3 is 4 times that of strain 1
# strain_transmission1 <- c(1,4)
# # lower cross-immunity from previous infection with strains 1 or 2, but can't be
# # infected with strains 1 or 2 after being infected with strain 3
# cross_immunity1 <- c(0.2,1) 
# p1 <- parameters(dt,
#                  n_age,
#                  n_vax,
#                  transmission,
#                  beta_date,
#                  beta_value = beta_value_sim,
#                  beta_type,
#                  gamma_E,
#                  gamma_P = 0.4,
#                  gamma_A,
#                  gamma_C = 0.4,
#                  gamma_H,
#                  gamma_G,
#                  gamma_pre_1,
#                  gamma_P_1,
#                  theta_A,
#                  p_C,
#                  p_H,
#                  p_G,
#                  p_D,
#                  population,
#                  start_date = n_steps*dt,
#                  initial_seed_size = 0,
#                  initial_seed_pattern,
#                  strain_transmission1,
#                  strain_seed_date = 210L,
#                  strain_seed_size,
#                  strain_seed_pattern,
#                  strain_rel_p_sympt,
#                  strain_rel_p_hosp_if_sympt,
#                  strain_rel_p_death,
#                  rel_susceptibility,
#                  rel_p_sympt,
#                  rel_p_hosp_if_sympt,
#                  rel_p_death,
#                  rel_infectivity,
#                  vaccine_progression_rate,
#                  schedule,
#                  vaccine_index_dose2,
#                  vaccine_index_booster,
#                  vaccine_catchup_fraction,
#                  n_doses,
#                  waning_rate,
#                  cross_immunity,
#                  sero_sensitivity_1,
#                  sero_specificity_1)
# 
# n_steps1 <- 1200 # number of steps to run model up to



#### Try simulating some data to see if the parameters give the right general pattern ####


out <- simulate_data(gen_seirhd_age_vax_multistrain_sero_time_dep_beta, p, n_steps)

# Drop time row
x <- out$x
time <- x[1,1,-1]
x <- x[-1, , ,drop=FALSE]

n_strains <- 4

# Plot trajectories
plot_trajectories(time,x,n_age,n_strains,n_vax)

# Extract true history of model states
true_history <- x[ , ,seq(0,n_steps+1,by=1/dt)+1,drop=F]

# Add noise to simulated data
info <- out$info
idx <- index(info)
hosps <- true_history[idx$run[grep("hosps_",names(idx$run))]-1, ,-1]
deaths <- true_history[idx$run[grep("deaths_",names(idx$run))]-1, ,-1]
sero_pos <- true_history[idx$run[grep("sero_pos_1_",names(idx$run))]-1, ,-1]
cases <- matrix(true_history[idx$run[match("cases",names(idx$run))]-1, ,-1],nrow = 1)
cases_non_variant <- matrix(true_history[idx$run[grep("cases_non_variant",names(idx$run))]-1, ,-1],nrow = 1)

par(mfrow = c(1,1))
days <- seq(1,n_steps*dt)
matplot(days,t(hosps),type="l",xlab="Day",ylab="Hospitalisations")
matplot(days,t(deaths),type="l",xlab="Day",ylab="Deaths")
matplot(days,t(sero_pos),type="l",xlab="Day",ylab="Seropositive")
matplot(days,t(cases),type="l",xlab="Day",ylab="Cases")
matplot(days,t(cases_non_variant),type="l",xlab="Day",ylab="Non-variant cases")



#### Fit to multiple age-stratified data streams ####


# Convert raw data to required format for particle filter
data <- particle_filter_data(data_raw,"day",1/dt)

# Create particle filter object
n_particles <- 200
filter <- particle_filter$new(data, gen_seirhd_age_vax_multistrain_sero_time_dep_beta, n_particles, 
                              compare, index, initial)
filter$run(
    pars = p,
    save_history = TRUE)

# Plot filtered trajectories
# plot_particle_filter(filter$history(),true_history,data_raw$day,idx)
plot_hosps_and_deaths_age(filter$history(),data,data_raw$day,n_age,n_vax,n_strains)
plot_sero(filter$history(),data,data_raw$day)
plot_cases(filter$history(),data,data_raw$day)

# Infer parameters by pMCMC
beta_value_list <- list(name = c("beta1","beta2","beta3","beta4","beta5"), initial = c(0.04,0.03,0.02,0.04,0.02),
                        min = rep(0,5), max = rep(Inf,5), discrete = rep(F,5),
                        prior = replicate(5,function(x) dgamma(x, shape = 1, scale = 1, log = TRUE)))
beta_value <- Map(pmcmc_parameter,beta_value_list$name,beta_value_list$initial,
                  beta_value_list$min,beta_value_list$max,
                  beta_value_list$discrete,beta_value_list$prior)
# R0 <- pmcmc_parameter("R0",2,min = 0)
gamma <- pmcmc_parameter("gamma",0.4,min = 0,
                         prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
start_date <- pmcmc_parameter("start_date",1L,min = 1L,max = 10L,discrete = TRUE)
strain_seed_date <- pmcmc_parameter("strain_seed_date",300L,min = 290L,max = 330L,discrete = TRUE)
# strain_seed_date1 <- pmcmc_parameter("strain_seed_date1",200L,min = 200L,max = 250L,discrete = TRUE)

# proposal <- matrix(c(0.01^2,0,0,0.01^2),nrow = 2,ncol = 2,byrow = TRUE)
# proposal <- matrix(c(0.01^2,0,0,0,0.01^2,0,0,0,2),nrow = 3,ncol = 3,byrow = TRUE)
# proposal <- matrix(c(0.01^2,0,0,0,
#                      0,0.01^2,0,0,
#                      0,0,2,0,
#                      0,0,0,2),nrow = 4,ncol = 4,byrow = TRUE)
# proposal <- matrix(c(0.01^2,0,0,0,0,
#                      0,0.01^2,0,0,0,
#                      0,0,2,0,0,
#                      0,0,0,2,0,
#                      0,0,0,0,2),nrow = 5,ncol = 5,byrow = TRUE)
proposal <- diag(c(rep(0.01^2,6),rep(2^2,2)))
transform <- make_transform(dt,
                            n_age,
                            n_vax,
                            transmission,
                            beta_date,
                            beta_type,
                            gamma_E,
                            gamma_A,
                            gamma_H,
                            gamma_G,
                            gamma_pre_1,
                            gamma_P_1,
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
                            vaccine_index_booster,
                            vaccine_catchup_fraction,
                            n_doses,
                            waning_rate,
                            cross_immunity)
pars_mcmc <- c(beta_value,
               list(gamma = gamma,
                    start_date = start_date,
                    strain_seed_date = strain_seed_date))
# transform <- make_transform_multistage(dt,
#                                        n_age,
#                                        n_vax,
#                                        transmission,
#                                        beta_date,
#                                        beta_type,
#                                        gamma_E,
#                                        gamma_A,
#                                        gamma_H,
#                                        gamma_G,
#                                        gamma_pre_1,
#                                        gamma_P_1,
#                                        theta_A,
#                                        p_C,
#                                        p_H,
#                                        p_G,
#                                        p_D,
#                                        population,
#                                        # start_date,
#                                        initial_seed_size,
#                                        initial_seed_pattern,
#                                        strain_transmission,
#                                        # strain_seed_date,
#                                        strain_seed_size,
#                                        strain_seed_pattern,
#                                        strain_rel_p_sympt,
#                                        strain_rel_p_hosp_if_sympt,
#                                        strain_rel_p_death,
#                                        rel_susceptibility,
#                                        rel_p_sympt,
#                                        rel_p_hosp_if_sympt,
#                                        rel_p_death,
#                                        rel_infectivity,
#                                        vaccine_progression_rate,
#                                        schedule,
#                                        vaccine_index_dose2,
#                                        vaccine_index_booster,
#                                        vaccine_catchup_fraction,
#                                        n_doses,
#                                        waning_rate,
#                                        cross_immunity,
#                                        start_date1 = n_steps * dt,
#                                        strain_transmission1,
#                                        cross_immunity1,
#                                        sero_sensitivity_1,
#                                        sero_specificity_1)
# pars_mcmc <- c(beta_value,
#                list(gamma = gamma,
#                     start_date = start_date,
#                     strain_seed_date = strain_seed_date,
#                     strain_seed_date1 = strain_seed_date1))
mcmc_pars <- pmcmc_parameters$new(pars_mcmc,proposal,transform = transform)
# mcmc_pars <- pmcmc_parameters$new(list(R0 = R0,gamma = gamma),
#                                   proposal,transform = transform)

# Run MCMC
# Using mcstate
n_samples <- 100 #1000
control <- pmcmc_control(
    n_samples,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE)
pmcmc_run <- pmcmc(mcmc_pars,filter,control = control)

# Check variable names in saved state trajectories
dimnames(pmcmc_run$trajectories$state)

plot_particle_filter(pmcmc_run$trajectories$state,true_history,data_raw$day,idx)
plot_hosps_and_deaths_age(pmcmc_run$trajectories$state,data,data_raw$day,n_age,n_vax,n_strains)
plot_sero(pmcmc_run$trajectories$state,data,data_raw$day)
plot_cases(pmcmc_run$trajectories$state,data,data_raw$day)

# Plot MCMC output
mcmc_out <- as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars))
summary(mcmc_out)
plot(mcmc_out)

# Using custom accelerated adaptive MCMC algorithm
init_pars <- c(beta_value_list$initial,gamma$initial,start_date$initial,strain_seed_date$initial)
priors <- c(beta_value_list$prior,gamma$prior,replicate(2,function(x) 0))
n_iters <- 100
scaling_factor_start <- 1
pars_min <- c(beta_value_list$min,gamma$min,start_date$min,strain_seed_date$min)
pars_max <- c(beta_value_list$max,gamma$max,start_date$max,strain_seed_date$max)
iter0 <- 10
discrete <- c(rep(F,6),rep(T,2))
names(init_pars) <- names(priors) <- names(pars_min) <- names(pars_max) <- names(discrete) <- c(beta_value_list$name,"gamma","start_date","strain_seed_date")

tstart <- Sys.time()
res <- mcmc(transform,filter,init_pars,priors,n_particles,n_iters,scaling_factor_start,proposal,pars_min,pars_max,iter0,discrete)
tend <- Sys.time()
print(tend - tstart)
# Time difference of 26.56045 mins

# Trace plots
for (i in seq_along(init_pars)){
    par(mfrow = c(1,1))
    plot(res$pars[,i],type="l")
}

# Scaling factor
par(mfrow = c(1,1))
plot(res$scaling_factor,type="l")


dimnames(res$states) <- dimnames(pmcmc_run$trajectories$state)

# Plot fitted hospitalisations and deaths against data
plot_hosps_and_deaths_age(res$states,data,data_raw$day,n_age,n_vax,n_strains)
