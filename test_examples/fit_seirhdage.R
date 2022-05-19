library(odin)
library(odin.dust)
library(dust)
library(mcstate)
library(socialmixr)
library(coda)
library(scales)
library(qs)
library(data.table)
library(extraDistr)

#### Age-structured SEIR model with deaths
gen_seirhd_age <- odin.dust::odin_dust("test_examples/seirhdage.R")

# Get age-structured contact matrix
data(polymod, package = "socialmixr")

age.limits <- seq(0,70,10)
contact <- socialmixr::contact_matrix(survey = polymod,countries = "United Kingdom",age.limits = age.limits,symmetric = T)

## Transform the matrix to the (symmetrical) transmission matrix
## rather than the contact matrix. This transmission matrix is
## weighted by the population in each age band.
# population <- contact$demography$population
population <- round(contact$demography$population/10000)
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

# # Get max values
# p_D_max0 <- p_D[length(p_D)]
# p_H_max0 <- p_H[length(p_H)]
# # Normalise
# p_D <- p_D/max(p_D)
# p_H <- p_H/max(p_H)

# IFR <- readRDS("~/UCSF/COVIDVaccineModelling/Data/IFR_by_age_ODriscoll.RDS")
# p_death <- IFR$median_perc[19:27]/100
# p_death[length(p_death)-1] <- (p_death[length(p_death)-1]+p_death[length(p_death)])/2
# p_death <- p_death[1:(length(p_death)-1)]

# Generate model
n_age <- length(age.limits)
dt <- 0.25
S_ini <- population
E_ini <- c(0,0,0,0,0,0,0,0)
I_ini <- c(1,1,2,2,2,1,1,0)
seirhd_age <- gen_seirhd_age$new(
    list(dt = dt,n_age = n_age,S_ini = S_ini,E_ini = E_ini,I_ini = I_ini,
         m = transmission,beta = 0.04,sigma = 0.5,gamma_P = 0.4,gamma_A = 0.2,
         gamma_C = 0.4,gamma_H = 0.1, gamma_G = 1/3, 
         p_C = p_C,p_H = p_H,p_G = p_G,p_D = p_D),
         # p_C = p_C,p_H = p_H_max*p_H,p_G = p_G,p_D = p_D_max*p_D),
    step = 0,n_particles = 1,n_threads = 1,seed = 1)

seirhd_age$info()

# Run epidemic forward
n_steps <- 400

# Create data to be fitted to
# Create an array to contain outputs after looping the model.
# Array contains XX rows = Total S, E, I, R (4), and
# in each age compartment (XX) as well as the cumulative incidence (XX)
x <- array(NA, dim = c(seirhd_age$info()$len, 1, n_steps+1))

# For loop to run the model iteratively
x[ , ,1] <- seirhd_age$state()
for (t in seq_len(n_steps)) {
    x[ , ,t+1] <- seirhd_age$run(t)
}
time <- x[1,1,-1]

# Drop time row
x <- x[-1, , ,drop=FALSE]

# Plot trajectories
par(mfrow = c(2,4), oma=c(2,3,0,0))
for (i in 1:n_age) {
    par(mar = c(3, 4, 2, 0.5))
    cols <- c(S = "#8c8cd9", E = "#ffff00", I_P = "#cc0044", I_A = "green", I_C = "blue", R = "#999966", D = "#000000")
    matplot(time, x[i + 12, ,-1], type = "l", # Offset to access numbers in age compartment
            xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
            col = cols[["S"]], lty = 1, ylim=range(x[-1:-4,,]))
    matlines(time, x[i + 12 + n_age, ,-1], col = cols[["E"]], lty = 1)
    matlines(time, x[i + 12 + 2*n_age, ,-1], col = cols[["I_P"]], lty = 1)
    matlines(time, x[i + 12 + 3*n_age, ,-1], col = cols[["I_A"]], lty = 1)
    matlines(time, x[i + 12 + 4*n_age, ,-1], col = cols[["I_C"]], lty = 1)
    matlines(time, x[i + 12 + 5*n_age, ,-1], col = cols[["R"]], lty = 1)
    matlines(time, x[i + 12 + 8*n_age, ,-1], col = cols[["D"]], lty = 1)
    legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
    axis(2, las = 2)
}
mtext("Number of individuals", side=2,line=1, outer=T)
mtext("Time", side = 1, line = 0, outer =T)

# Extract true history of model states
true_history <- x[ , ,seq(0,n_steps+1,by=1/dt)+1,drop=F]



#### Fit to multiple age-stratified data streams
index <- function(info){
    index <- info$index
    
    index_run <- c(hosps = index[["H_inc"]],
                   hosps_0_39 = index[["H_inc_0_39"]],
                   hosps_40_49 = index[["H_inc_40_49"]],
                   hosps_50_59 = index[["H_inc_50_59"]],
                   hosps_60_69 = index[["H_inc_60_69"]],
                   hosps_70_plus = index[["H_inc_70_plus"]],
                   # hosps_70_79 = index[["H_inc_70_79"]],
                   # hosps_80_plus = index[["H_inc_80_plus"]],
                   deaths = index[["D_inc"]],
                   deaths_0_39 = index[["D_inc_0_39"]],
                   deaths_40_49 = index[["D_inc_40_49"]],
                   deaths_50_59 = index[["D_inc_50_59"]],
                   deaths_60_69 = index[["D_inc_60_69"]],
                   deaths_70_plus = index[["D_inc_70_plus"]]
                   # deaths_70_79 = index[["D_inc_70_79"]],
                   # deaths_80_plus = index[["D_inc_80_plus"]]
    )
    
    list(run = index_run,
         state = c(S = index[["S"]],
                   E = index[["E"]],
                   I_P = index[["I_P"]],
                   I_A = index[["I_A"]],
                   I_C = index[["I_C"]],
                   R = index[["R"]],
                   H = index[["H"]],
                   G = index[["G"]],
                   D = index[["D"]],
                   index_run))
}

index(seirhd_age$info())

# log-likelihood of Poisson count
ll_pois <- function(obs, model, exp_noise) {
    if (is.na(obs)) {
        # Creates vector of zeros in ll with same length, if no data
        return(numeric(length(model)))
    } 
    lambda <- model + rexp(n = length(model), rate = exp_noise)
    dpois(x = obs, lambda = lambda, log = TRUE)
}

# Define negative binomial log-likelihood
ll_nbinom <- function(data, model, kappa, exp_noise){
    if(is.na(data)) {
        return(numeric(length(model)))
    }
    mu <- model + rexp(n = length(model), rate = exp_noise)
    dnbinom(data, kappa, mu = mu, log = TRUE)
}

ll_dirmnom <- function(data, model, size, exp_noise){
    if(any(is.na(data))){
        return(numeric(length(model)))
    }
    model <- model + rexp(n = length(model), rate = exp_noise)
    alpha <- size * model/rowSums(model)
    ddirmnom(data, rowSums(data), alpha, log = TRUE)
}

# Define comparison function for age-stratified data
compare <- function(state, observed, pars = NULL){
    if (is.null(pars$kappa_hosp)){
        kappa_hosp <- 2
    } else {
        kappa_hosp <- pars$kappa_hosp
    }
    if (is.null(pars$kappa_death)){
        kappa_death <- 2
    } else {
        kappa_death <- pars$kappa_death
    }
    if (is.null(pars$size_hosp)){
        size_hosp <- 10
    } else {
        size_hosp <- pars$size_hosp
    }
    if (is.null(pars$size_death)){
        size_death <- 10
    } else {
        size_death <- pars$size_death
    }
    if (is.null(pars$exp_noise)){
        exp_noise <- 1e6
    } else {
        exp_noise <- pars$exp_noise
    }
    
    model_hosps <- state["hosps", ] # for locations without age-stratified data
    model_hosps_0_39 <- state["hosps_0_39", ]
    model_hosps_40_49 <- state["hosps_40_49", ]
    model_hosps_50_59 <- state["hosps_50_59", ]
    model_hosps_60_69 <- state["hosps_60_69", ]
    model_hosps_70_plus <- state["hosps_70_plus", ]
    # model_hosps_70_79 <- state["hosps_70_79", ]
    # model_hosps_80_plus <- state["hosps_80_plus", ]
    model_deaths <- state["deaths", ] # for locations without age-stratified data
    model_deaths_0_39 <- state["deaths_0_39", ]
    model_deaths_40_49 <- state["deaths_40_49", ]
    model_deaths_50_59 <- state["deaths_50_59", ]
    model_deaths_60_69 <- state["deaths_60_69", ]
    model_deaths_70_plus <- state["deaths_70_plus", ]
    # model_deaths_70_79 <- state["deaths_70_79", ]
    # model_deaths_80_plus <- state["deaths_80_plus", ]
    
    # Log-likelihoods for hospitalisations
    # ll_hosps <- ll_nbinom(observed$hosps,model_hosps,5*kappa_hosp,exp_noise)
    # ll_hosps_0_39 <- ll_nbinom(observed$hosps_0_39,model_hosps_0_39,kappa_hosp,exp_noise)
    # ll_hosps_40_49 <- ll_nbinom(observed$hosps_40_49,model_hosps_40_49,kappa_hosp,exp_noise)
    # ll_hosps_50_59 <- ll_nbinom(observed$hosps_50_59,model_hosps_50_59,kappa_hosp,exp_noise)
    # ll_hosps_60_69 <- ll_nbinom(observed$hosps_60_69,model_hosps_60_69,kappa_hosp,exp_noise)
    # ll_hosps_70_plus <- ll_nbinom(observed$hosps_70_plus,model_hosps_70_plus,kappa_hosp,exp_noise)
    # # ll_hosps_70_79 <- ll_nbinom(observed$hosps_70_79,model_hosps_70_79,kappa_hosp,exp_noise)
    # # ll_hosps_80_plus <- ll_nbinom(observed$hosps_80_plus,model_hosps_80_plus,kappa_hosp,exp_noise)

    hosps_by_age <- matrix(c(observed$hosps_0_39,observed$hosps_40_49,observed$hosps_50_59,observed$hosps_60_69,observed$hosps_70_plus),nrow = 1)
    model_hosps_by_age <- cbind(model_hosps_0_39,model_hosps_40_49,model_hosps_50_59,model_hosps_60_69,model_hosps_70_plus)
    ll_hosps <- ll_dirmnom(hosps_by_age,model_hosps_by_age,size_hosp,exp_noise)
    
    # Log-likelihoods for deaths
    # ll_deaths <- ll_nbinom(observed$deaths,model_deaths,5*kappa_death,exp_noise)
    # ll_deaths_0_39 <- ll_nbinom(observed$deaths_0_39,model_deaths_0_39,kappa_death,exp_noise)
    # ll_deaths_40_49 <- ll_nbinom(observed$deaths_40_49,model_deaths_40_49,kappa_death,exp_noise)
    # ll_deaths_50_59 <- ll_nbinom(observed$deaths_50_59,model_deaths_50_59,kappa_death,exp_noise)
    # ll_deaths_60_69 <- ll_nbinom(observed$deaths_60_69,model_deaths_60_69,kappa_death,exp_noise)
    # ll_deaths_70_plus <- ll_nbinom(observed$deaths_70_plus,model_deaths_70_plus,kappa_death,exp_noise)
    # # ll_deaths_70_79 <- ll_nbinom(observed$deaths_70_79,model_deaths_70_79,kappa_death,exp_noise)
    # # ll_deaths_80_plus <- ll_nbinom(observed$deaths_80_plus,model_deaths_80_plus,kappa_death,exp_noise)
    
    deaths_by_age <- matrix(c(observed$deaths_0_39,observed$deaths_40_49,observed$deaths_50_59,observed$deaths_60_69,observed$deaths_70_plus),nrow = 1)
    model_deaths_by_age <- cbind(model_deaths_0_39,model_deaths_40_49,model_deaths_50_59,model_deaths_60_69,model_deaths_70_plus)
    ll_deaths <- ll_dirmnom(deaths_by_age,model_deaths_by_age,size_death,exp_noise)

    # ll_hosps <- ll_pois(observed$hosps,model_hosps,exp_noise)
    # ll_hosps_0_39 <- ll_pois(observed$hosps_0_39,model_hosps_0_39,exp_noise)
    # ll_hosps_40_49 <- ll_pois(observed$hosps_40_49,model_hosps_40_49,exp_noise)
    # ll_hosps_50_59 <- ll_pois(observed$hosps_50_59,model_hosps_50_59,exp_noise)
    # ll_hosps_60_69 <- ll_pois(observed$hosps_60_69,model_hosps_60_69,exp_noise)
    # ll_hosps_70_plus <- ll_pois(observed$hosps_70_plus,model_hosps_70_plus,exp_noise)
    # # ll_hosps_70_79 <- ll_pois(observed$hosps_70_79,model_hosps_70_79,exp_noise)
    # # ll_hosps_80_plus <- ll_pois(observed$hosps_80_plus,model_hosps_80_plus,exp_noise)
    # 
    # # Log-likelihoods for deaths
    # ll_deaths <- ll_pois(observed$deaths,model_deaths,exp_noise)
    # ll_deaths_0_39 <- ll_pois(observed$deaths_0_39,model_deaths_0_39,exp_noise)
    # ll_deaths_40_49 <- ll_pois(observed$deaths_40_49,model_deaths_40_49,exp_noise)
    # ll_deaths_50_59 <- ll_pois(observed$deaths_50_59,model_deaths_50_59,exp_noise)
    # ll_deaths_60_69 <- ll_pois(observed$deaths_60_69,model_deaths_60_69,exp_noise)
    # ll_deaths_70_plus <- ll_pois(observed$deaths_70_plus,model_deaths_70_plus,exp_noise)
    # # ll_deaths_70_79 <- ll_pois(observed$deaths_70_79,model_deaths_70_79,exp_noise)
    # # ll_deaths_80_plus <- ll_pois(observed$deaths_80_plus,model_deaths_80_plus,exp_noise)
    
    # Calculate total log-likelihood
    # ll_hosps + ll_hosps_0_39 + ll_hosps_40_49 + ll_hosps_50_59 + ll_hosps_60_69 + ll_hosps_70_plus + #ll_hosps_70_79 + ll_hosps_80_plus +
    #     ll_deaths + ll_deaths_0_39 + ll_deaths_40_49 + ll_deaths_50_59 + ll_deaths_60_69 + ll_deaths_70_plus #+ ll_deaths_70_79 + ll_deaths_80_plus
    # ll_hosps_0_39 + ll_hosps_40_49 + ll_hosps_50_59 + ll_hosps_60_69 + ll_hosps_70_plus + #ll_hosps_70_79 + ll_hosps_80_plus +
    #     ll_deaths_0_39 + ll_deaths_40_49 + ll_deaths_50_59 + ll_deaths_60_69 + ll_deaths_70_plus #+ ll_deaths_70_79 + ll_deaths_80_plus
    # ll_hosps_0_39 + ll_deaths_0_39
    # ll_hosps_70_plus + ll_deaths_70_plus
    ll_hosps + ll_deaths
}

# Add noise to simulated data
hosps <- true_history[3:7, ,-1]
deaths <- true_history[8:12, ,-1]

par(mfrow = c(1,1))
days <- seq(1,n_steps*dt)
matplot(days,t(hosps),type="l",xlab="Day",ylab="Hospitalisations")
matplot(days,t(deaths),type="l",xlab="Day",ylab="Deaths")

add_noise <- function(x,f){
    noise <- apply(x,2,function(y) round(rnorm(length(y),0,f*y)))
    x <- x + noise
    x <- pmax(x,0)
}

add_nbinom_noise <- function(x,size){
    x <- apply(x,2,function(y) rnbinom(length(y),size = size,mu = y))
}

set.seed(0)
# hosps <- add_noise(hosps,0.2)
hosps <- add_nbinom_noise(hosps,100)
rownames(hosps) <- paste0("hosps_",age.limits[c(1,5:length(age.limits))],"_",c(as.character(age.limits[5:length(age.limits)]-1),"plus"))
# deaths <- add_noise(deaths,0.3)
deaths <- add_nbinom_noise(deaths,10)
rownames(deaths) <- paste0("deaths_",age.limits[c(1,5:length(age.limits))],"_",c(as.character(age.limits[5:length(age.limits)]-1),"plus"))

matplot(days,t(hosps),type="l",xlab="Day",ylab="Hospitalisations")
matplot(days,t(deaths),type="l",xlab="Day",ylab="Deaths")

# Create "observed" data
data_raw <- data.frame(t(rbind(hosps,deaths)))
data_raw$day <- days
# Add empty columns for total cases and deaths
# data_raw$hosps <- NA
# data_raw$deaths <- NA
data_raw$hosps <- data$hosps_0_39 + data$hosps_40_49 + data$hosps_50_59 + data$hosps_60_69 + data$hosps_70_plus
data_raw$deaths <- data$deaths_0_39 + data$deaths_40_49 + data$deaths_50_59 + data$deaths_60_69 + data$deaths_70_plus
# Convert to required format
data <- particle_filter_data(data_raw,"day",1/dt)

# Function for plotting fitted trajectories
plot_particle_filter <- function(history, true_history, times, obs_end = NULL) {
    if (is.null(obs_end)) {
        obs_end = max(times)
    }
    
    # par(mfrow = c(2,4), oma=c(2,3,0,0))
    par(mfrow = c(1,1), oma=c(2,3,0,0))
    for (i in 1:n_age){
        par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
        cols <- c(S = "#8c8cd9", E = "#ffff00", I_P = "#cc0044", I_A = "green", I_C = "blue", R = "#999966", D = "#000000")
        # matplot(times, t(history[i, ,-1]), type = "l", # Offset to access numbers in age compartment
        #         xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
        #         col = alpha(cols[["S"]],0.1), lty = 1, ylim=range(history))
        matplot(times, t(history[i + n_age, ,-1]), type = "l", # Offset to access numbers in age compartment
                xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
                col = alpha(cols[["E"]],0.1), lty = 1, ylim=range(true_history[i + 12 + seq(n_age,4*n_age,by = n_age), ,-1]))
        # matlines(times, t(history[i + n_age, ,-1]), col = alpha(cols[["E"]],0.1), lty = 1)
        matlines(times, t(history[i + 2*n_age, ,-1]), col = alpha(cols[["I_P"]],0.1), lty = 1)
        matlines(times, t(history[i + 3*n_age, ,-1]), col = alpha(cols[["I_A"]],0.1), lty = 1)
        matlines(times, t(history[i + 4*n_age, ,-1]), col = alpha(cols[["I_C"]],0.1), lty = 1)
        # matlines(times, t(history[i + 5*n_age, ,-1]), col = alpha(cols[["R"]],0.1), lty = 1)
        # matlines(times, t(history[i + 8*n_age, ,-1]), col = alpha(cols[["D"]],0.1), lty = 1)
        # matpoints(times[1:obs_end], t(true_history[seq(i+12,i+12+8*n_age,by = n_age), ,-1]), pch = 19,
        #           col = cols)
        matpoints(times[1:obs_end], t(true_history[i + 12 + seq(n_age,4*n_age,by = n_age), ,-1]), pch = 19,
                  col = cols[2:5])
        # legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
        legend("right", lwd = 1, col = cols[2:5], legend = names(cols[2:5]), bty = "n")
        axis(2, las = 2)
    }
}

plot_hosps_and_deaths_age <- function(incidence_modelled, incidence_observed, times){
    par(mfrow = c(2,4), oma=c(2,3,0,0))
    for (i in 1:5){
        par(mar = c(3, 4, 2, 0.5))
        matplot(times, t(incidence_modelled[73+i, ,-1]),
                type="l",col = alpha("black",0.1),xlab = "Day",ylab = "Hospitalisations",
                main = paste0("Age ", rownames(incidence_modelled)[73+i]))
        points(times, incidence_observed[[4+i]],pch=19,col="red")
        axis(2, las = 2)
    }
    par(mfrow = c(2,4), oma=c(2,3,0,0))
    for (i in 1:5){
        par(mar = c(3, 4, 2, 0.5))
        matplot(times, t(incidence_modelled[79+i, ,-1]),
                type="l",col = alpha("black",0.1),xlab = "Day",ylab = "Deaths",
                main = paste0("Age ", rownames(incidence_modelled)[79+i]))
        points(times, incidence_observed[[9+i]],pch=19,col="red")
        axis(2, las = 2)
    }
}

# Create particle filter object
n_particles <- 200
filter <- particle_filter$new(data, gen_seirhd_age, n_particles, 
                              compare, index, seed = 1)
filter$run(
    save_history = TRUE,
    pars = list(dt = dt,n_age = n_age,S_ini = S_ini,E_ini = E_ini,I_ini = I_ini,
                m = transmission,beta = 0.04,sigma = 0.5,gamma_P = 0.4,gamma_A = 0.2,
                gamma_C = 0.4,gamma_H = 0.1, gamma_G = 1/3, 
                p_C = p_C,p_H = p_H,p_G = p_G,p_D = p_D,kappa_hosp = 100,kappa_death = 10))

# Check variable names in particle filter history
dimnames(filter$history())

# Plot filtered trajectories
plot_particle_filter(filter$history(),true_history,data_raw$day)
plot_hosps_and_deaths_age(filter$history(),data,data_raw$day)

# Infer parameters by pMCMC
beta <- pmcmc_parameter("beta",0.04,min = 0,
                        prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
gamma <- pmcmc_parameter("gamma",0.5,min = 0,
                         prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
# p_H_max <- pmcmc_parameter("p_H_max",p_H_max0,min = 0,max = 1,
#                            prior = function(x) dbeta(x,1,1,log = TRUE))
# p_D_max <- pmcmc_parameter("p_D_max",p_D_max0,min = 0,max = 1,
#                            prior = function(x) dbeta(x,1,1,log = TRUE))
# alpha_H <- pmcmc_parameter("alpha_H",0.01,min = 0,max = 1,
#                            prior = function(x) dbeta(x,1,1,log = TRUE))
# alpha_D <- pmcmc_parameter("alpha_D",0.1,min = 0,max = 1,
#                            prior = function(x) dbeta(x,1,1,log = TRUE))

parameter_transform <- function(pars,dt,n_age,S_ini,E_ini,I_ini,m,sigma,gamma_A,gamma_H,gamma_G,p_C,p_H,p_G,p_D){
    beta <- pars[["beta"]]
    gamma <- pars[["gamma"]]
    # p_H_max <- pars[["p_H_max"]]
    # p_D_max <- pars[["p_D_max"]]
    # alpha_H <- pars[["alpha_H"]]
    # alpha_D <- pars[["alpha_D"]]
    list(beta = beta,dt = dt,n_age = n_age,S_ini = S_ini,E_ini = E_ini,
         I_ini = I_ini,m = m,sigma = sigma,gamma_P = gamma,gamma_A = gamma_A,gamma_C = gamma,
         gamma_H = gamma_H,gamma_G = gamma_G,p_C = p_C,p_H = p_H,p_G = p_G,p_D = p_D)
    # list(beta = beta,dt = dt,n_age = n_age,S_ini = S_ini,E_ini = E_ini,
    #      I_ini = I_ini,m = m,sigma = sigma,gamma_P = gamma,gamma_A = gamma_A,gamma_C = gamma,
    #      gamma_H = gamma_H,gamma_G = gamma_G,p_C = p_C,p_H = p_H_max*p_H,p_G = p_G,p_D = p_D_max*p_D)#,
    #      # kappa_hosp = 1/alpha_H,kappa_death = 1/alpha_D)
}

transform <- function(pars){
    parameter_transform(pars,dt = dt,n_age = n_age,S_ini = S_ini,E_ini = E_ini,
                        I_ini = I_ini,m = transmission,sigma = 0.5,
                        gamma_A = 0.2,gamma_H = 0.1,gamma_G = 1/3,
                        p_C = p_C,p_H = p_H,p_G = p_G,p_D = p_D)
}

proposal <- diag(c(1e-6,1e-5),2)
# proposal <- matrix(c(0.01^2,0,0,0.01^2),nrow = 2,ncol = 2,byrow = TRUE)
# proposal <- diag(c(1e-6,1e-5,1e-5,1e-5,1e-7,1e-5),6)
# proposal <- diag(c(1e-6,1e-5,1e-5,1e-5),4)
mcmc_pars <- pmcmc_parameters$new(list(beta = beta,gamma = gamma),
                                       # p_H_max = p_H_max, p_D_max = p_D_max),#,
                                       # alpha_H = alpha_H, alpha_D = alpha_D),
                                  proposal,transform = transform)

# Run MCMC
control <- pmcmc_control(
    1000,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE)
pmcmc_run <- pmcmc(mcmc_pars,filter,control = control)

plot_particle_filter(pmcmc_run$trajectories$state[,301:1000,],true_history,data_raw$day)
plot_hosps_and_deaths_age(pmcmc_run$trajectories$state[,301:1000,],data,data_raw$day)

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
beta <- pmcmc_parameter("beta",pmcmc_run$pars[nrow(pmcmc_run$pars),1],min = 0,
                        prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
gamma <- pmcmc_parameter("gamma",pmcmc_run$pars[nrow(pmcmc_run$pars),2],min = 0,
                        prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
# alpha_H <- pmcmc_parameter("alpha_H",pmcmc_run$pars[nrow(pmcmc_run$pars),3],min = 0,max = 1,
#                            prior = function(x) dbeta(x,1,1,log = TRUE))
# alpha_D <- pmcmc_parameter("alpha_D",pmcmc_run$pars[nrow(pmcmc_run$pars),4],min = 0,max = 1,
#                            prior = function(x) dbeta(x,1,1,log = TRUE))
proposal1 <- cov(pmcmc_run$pars)
mcmc_pars1 <- pmcmc_parameters$new(list(beta = beta,gamma = gamma),
                                        # alpha_H = alpha_H,alpha_D = alpha_D),
                                   proposal1,transform = transform)

pmcmc_run1 <- pmcmc(mcmc_pars1,filter,control = control)

plot_particle_filter(pmcmc_run1$trajectories$state[,101:1000,],true_history,data_raw$day)
plot_hosps_and_deaths_age(pmcmc_run1$trajectories$state[,101:1000,],data,data_raw$day)

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

# Save workspace
save.image("output/MCMC_output_seirhd_age_dirmnom_lklhd.RData")
