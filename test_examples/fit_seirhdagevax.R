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
source("R/seirhdagevax.R")

#### Age-structured SEIR model with deaths and vaccination
gen_seirhd_age_vax <- odin_dust("test_examples/seirhdagevax.R")

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
# S_ini <- population
# E_ini <- c(0,0,0,0,0,0,0,0)
# I_ini <- c(1,1,2,2,2,1,1,0)

# # Transmission and natural history parameters
gamma_E <- 0.5
gamma_A <- 0.2
gamma_H <- 0.1
gamma_G <- 1/3
theta_A <- 0.5
# p <- list(dt = dt,n_age = n_age,n_vax = n_vax,m = transmission,beta = 0.04,
#           gamma_E = 0.5,gamma_P = 0.4,gamma_A = 0.2,gamma_C = 0.4,gamma_H = 0.1, 
#           gamma_G = 1/3,theta_A = 0.5,p_C = p_C,p_H = p_H,p_G = p_G,p_D = p_D)
# 
# # Total population
# p$N_tot <- population
# 
# Seeding parameters
start_date <- 1L
initial_seed_size <- 10
initial_seed_pattern <- 1

# start_step <- start_date/dt
# seed_step_start <- floor(start_step)
# seed_value <- seed_over_steps(start_step, initial_seed_pattern) * initial_seed_size
# 
# p$seed_step_start <- seed_step_start
# p$seed_value <- seed_value

# Vaccination parameters
vaccine_progression_rate <- c(0,0,1/(26*7))

daily_doses <- rep(1000,365) #rep(5,365)
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

# Parameters for impact of vaccination of susceptibility and infectiousness
rel_susceptibility <- c(1,0.8,0.5) # relative susceptibility to infection in each vaccine stratum
rel_p_sympt <- c(1,0.6,0.3) # relative risk of symptoms in each vaccine stratum
rel_p_hosp_if_sympt <- c(1,0.95,0.95) # relative risk of hospitalisation given infection in each vaccine stratum
rel_p_death <- c(1,0.9,0.9) # relative risk of death in each vaccine stratum
rel_infectivity <- c(1,0.5,0.5) # relative infectiousness of infected individuals in each vaccine stratum

# # Construct vaccination parameters
# vaccination <- vaccination_parameters(p$N_tot,
#                                       p$dt,
#                                       rel_susceptibility,
#                                       rel_p_sympt,
#                                       rel_p_hosp_if_sympt,
#                                       rel_p_death,
#                                       rel_infectivity,
#                                       vaccine_progression_rate,
#                                       schedule,
#                                       vaccine_index_dose2 = 2L)
# 
# # Concatenate parameters into one list
# p <- c(p,vaccination)

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
                start_date,
                initial_seed_size,
                initial_seed_pattern,
                rel_susceptibility,
                rel_p_sympt,
                rel_p_hosp_if_sympt,
                rel_p_death,
                rel_infectivity,
                vaccine_progression_rate,
                schedule,
                vaccine_index_dose2,
                vaccine_catchup_fraction = vaccine_catchup_fraction,
                n_doses = n_doses)

# Create instance of model
seirhd_age_vax <- gen_seirhd_age_vax$new(p,step = 0,n_particles = 1,n_threads = 1,seed = 1)

info <- seirhd_age_vax$info()
initial_state <- initial(info, NULL, p)

seirhd_age_vax$update_state(state = initial_state)

# Run epidemic forward
n_steps <- 400

# Create data to be fitted to
# Create an array to contain outputs after looping the model.
# Array contains XX rows = Total S, E, I, R (4), and
# in each age compartment (XX) as well as the cumulative incidence (XX)
x <- array(NA, dim = c(seirhd_age_vax$info()$len, 1, n_steps+1))

# For loop to run the model iteratively
x[ , ,1] <- seirhd_age_vax$state()
for (t in seq_len(n_steps)) {
    x[ , ,t+1] <- seirhd_age_vax$run(t)
}
time <- x[1,1,-1]

# Drop time row
x <- x[-1, , ,drop=FALSE]

# Plot trajectories
for (j in 1:n_vax) {
    par(mfrow = c(2,4), oma=c(2,3,0,0))
    for (i in 1:n_age) {
        par(mar = c(3, 4, 2, 0.5))
        cols <- c(S = "#8c8cd9", E = "#ffff00", I_P = "#cc0044", I_A = "green", I_C = "blue", R = "#999966", D = "#000000")
        matplot(time, x[12 + i + (j-1)*n_age, ,-1], type = "l", # Offset to access numbers in age compartment
                xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
                col = cols[["S"]], lty = 1, ylim=range(x[-1:-4,,]))
        matlines(time, x[12 + i + (j-1)*n_age + n_age*n_vax, ,-1], col = cols[["E"]], lty = 1)
        matlines(time, x[12 + i + (j-1)*n_age + 2*n_age*n_vax, ,-1], col = cols[["I_P"]], lty = 1)
        matlines(time, x[12 + i + (j-1)*n_age + 3*n_age*n_vax, ,-1], col = cols[["I_A"]], lty = 1)
        matlines(time, x[12 + i + (j-1)*n_age + 4*n_age*n_vax, ,-1], col = cols[["I_C"]], lty = 1)
        matlines(time, x[12 + i + (j-1)*n_age + 5*n_age*n_vax, ,-1], col = cols[["R"]], lty = 1)
        matlines(time, x[12 + i + (j-1)*n_age + 8*n_age*n_vax, ,-1], col = cols[["D"]], lty = 1)
        legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
        axis(2, las = 2)
    }
    mtext("Number of individuals", side=2, line=1, outer=T, las=0)
    mtext("Time", side = 1, line = 0, outer =T)    
}

# Extract true history of model states
true_history <- x[ , ,seq(0,n_steps+1,by=1/dt)+1,drop=F]



#### Fit to multiple age-stratified data streams
index(seirhd_age_vax$info())

# Add noise to simulated data
hosps <- true_history[3:7, ,-1]
deaths <- true_history[8:12, ,-1]

par(mfrow = c(1,1))
days <- seq(1,n_steps*dt)
matplot(days,t(hosps),type="l",xlab="Day",ylab="Hospitalisations")

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

# Create particle filter object
n_particles <- 200
filter <- particle_filter$new(data, gen_seirhd_age_vax, n_particles, 
                              compare, index, initial)
filter$run(
    save_history = TRUE,
    pars = p)

# Check variable names in particle filter history
dimnames(filter$history())

# Plot filtered trajectories
plot_particle_filter(filter$history(),true_history,data_raw$day)
plot_hosps_and_deaths_age(filter$history(),data,data_raw$day,n_age,n_vax) #[,seq(0,1e4,by=10),] [,seq(0,5e3,by=5),]

# Infer parameters by pMCMC
beta <- pmcmc_parameter("beta",0.1,min = 0,
                        prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
# R0 <- pmcmc_parameter("R0",2,min = 0)
gamma <- pmcmc_parameter("gamma",0.5,min = 0,
                         prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))

proposal <- matrix(c(0.01^2,0,0,0.01^2),nrow = 2,ncol = 2,byrow = TRUE)
transform <- make_transform(dt,
                            n_age,
                            n_vax,
                            transmission,
                            gamma_E,
                            gamma_P,
                            gamma_A,
                            gamma_C,
                            gamma_H,
                            gamma_G,
                            theta_A,
                            p_C,
                            p_H,
                            p_G,
                            p_D,
                            population,
                            start_date,
                            initial_seed_size,
                            initial_seed_pattern,
                            rel_susceptibility,
                            rel_p_sympt,
                            rel_p_hosp_if_sympt,
                            rel_p_death,
                            rel_infectivity,
                            vaccine_progression_rate,
                            schedule,
                            vaccine_index_dose2,
                            vaccine_catchup_fraction,
                            n_doses)
mcmc_pars <- pmcmc_parameters$new(list(beta = beta,gamma = gamma),
                                  proposal,transform = transform)
# mcmc_pars <- pmcmc_parameters$new(list(R0 = R0,gamma = gamma),
#                                   proposal,transform = transform)

# Run MCMC
n_samples <- 1000
control <- pmcmc_control(
    n_samples,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE)
pmcmc_run <- pmcmc(mcmc_pars,filter,control = control)


dimnames(pmcmc_run$trajectories$state)

plot_particle_filter(pmcmc_run$trajectories$state[,401:1001,],true_history,data_raw$day)
plot_hosps_and_deaths_age(pmcmc_run$trajectories$state[,401:1001,],data,data_raw$day,n_age,n_vax)

# Plot MCMC output
mcmc_out <- as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars))
summary(mcmc_out)
plot(mcmc_out)

# Calculate effective sample size
effectiveSize(mcmc_out)
1 - rejectionRate(mcmc_out)

# Pairwise correlation plot
par(mfrow = c(1,1))
pairs(pmcmc_run$pars[])

# Tune MCMC
beta <- pmcmc_parameter("beta",mean(pmcmc_run$pars[(round(n_samples/5):n_samples)+1,1]),min = 0,
                        prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
gamma <- pmcmc_parameter("gamma",mean(pmcmc_run$pars[(round(n_samples/5):n_samples)+1,2]),min = 0,
                         prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
# beta <- pmcmc_parameter("beta",0.2,min = 0,
#                         prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
# gamma <- pmcmc_parameter("gamma",0.5,min = 0,
#                          prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
proposal1 <- cov(pmcmc_run$pars)
mcmc_pars1 <- pmcmc_parameters$new(list(beta = beta,gamma = gamma),
                                   proposal1,transform = transform)

pmcmc_run1 <- pmcmc(mcmc_pars1,filter,control = control)

plot_particle_filter(pmcmc_run$trajectories$state[,101:1001,],true_history,data_raw$day)
plot_hosps_and_deaths_age(pmcmc_run1$trajectories$state[,101:1001,],data,data_raw$day,n_age,n_vax)

# Plot MCMC output
mcmc_out1 <- as.mcmc(cbind(pmcmc_run1$probabilities, pmcmc_run1$pars))
summary(mcmc_out1)
plot(mcmc_out1)

# Calculate effective sample size
effectiveSize(mcmc_out1)
1 - rejectionRate(mcmc_out1)

# Pairwise correlation plot
par(mfrow = c(1,1))
pairs(pmcmc_run1$pars[])
