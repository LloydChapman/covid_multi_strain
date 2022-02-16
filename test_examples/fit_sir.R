library(odin)
library(odin.dust)
library(dust)
library(mcstate)
library(socialmixr)
library(coda)
library(scales)

## SIR example from "sir_models" vignette in mcstate package
gen_sir <- odin_dust("test_examples/sir.R")

# Simulate some data
dt <- 0.25
S_ini <- 1000
I_ini <- 10
sir <- gen_sir$new(
    list(dt = dt,S_ini = S_ini,I_ini = I_ini,
         beta = 0.2,gamma = 0.1),
    step = 0,n_particles = 1,n_threads = 1,seed = 1)

# Run epidemic forward
n_steps <- 400

# Create data to be fitted to
# Create an array to contain outputs after looping the model.
# Array contains XX rows = Total S, E, I, R (4), and
# in each age compartment (XX) as well as the cumulative incidence (XX)
x <- array(NA, dim = c(sir$info()$len, 1, n_steps+1))

# For loop to run the model iteratively
x[ , ,1] <- sir$state()
for (t in seq_len(n_steps)) {
    x[ , ,t+1] <- sir$run(t)
}
time <- x[1,1,-1]

# Drop time row
x <- x[-1, , ,drop=FALSE]

# Extract true history of model states
true_history <- x[ , ,seq(0,n_steps+1,by=1/dt)+1,drop=F]

# Get case incidence
incidence <- true_history[4, ,-1]
plot(incidence,type="l")

# Add noise
set.seed(0)
noise <- round(rnorm(n = length(true_history[2,1,-1]),mean = 0,
                     sd = 0.2 * sqrt(true_history[2,1,-1])))
incidence <- incidence + noise
incidence <- sapply(incidence,max,0)

# Convert to required format
data_raw <- data.frame(day = seq(1,n_steps*dt),cases = incidence)
data <- particle_filter_data(data_raw,"day",1/dt)

plot(data_raw$day,data_raw$cases,type = "l",xlab = "Day",ylab = "Cases")

# Define comparison function
case_compare <- function(state, observed, pars = NULL) {
    exp_noise <- 1e6
    
    incidence_modelled <- state[5, , drop = TRUE]
    incidence_observed <- observed$cases
    lambda <- incidence_modelled +
        rexp(n = length(incidence_modelled), rate = exp_noise)
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
}

# Create particle filter
filter <- particle_filter$new(data = data,
                              model = gen_sir,
                              n_particles = 100,
                              compare = case_compare,
                              seed = 1)

filter$run(save_history = TRUE, pars = list(dt = dt))

plot_particle_filter <- function(history, true_history, times, obs_end = NULL) {
    if (is.null(obs_end)) {
        obs_end = max(times)
    }
    
    par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
    cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966")
    matplot(times, t(history[2, , -1]), type = "l",
            xlab = "Time", ylab = "Number of individuals",
            col = alpha(cols[["S"]],0.1), lty = 1, ylim = range(history))
    matlines(times, t(history[3, , -1]), col = alpha(cols[["I"]],0.1), lty = 1)
    matlines(times, t(history[4, , -1]), col = alpha(cols[["R"]],0.1), lty = 1)
    matpoints(times[1:obs_end], t(true_history[1:3, , -1]), pch = 19,
              col = cols)
    legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
}

plot_particle_filter(filter$history(), true_history, data_raw$day)

# Set up MCMC
# beta <- pmcmc_parameter("beta", 0.2, min = 0)
R0 <- pmcmc_parameter("R0", 2, min = 0)
gamma <- pmcmc_parameter("gamma", 0.1, min = 0,
                         prior = function(p)
                             dgamma(p, shape = 1, scale = 0.2, log = TRUE))

parameter_transform <- function(pars,dt,S_ini,I_ini){
    # beta <- pars[["beta"]]
    beta <- pars[["R0"]] * pars[["gamma"]]
    gamma <- pars[["gamma"]]
    list(beta = beta,gamma = gamma,dt = dt,S_ini = S_ini,I_ini = I_ini)
}

transform <- function(pars){
    parameter_transform(pars,dt = dt,S_ini = S_ini,I_ini = I_ini)
}

proposal_matrix <- matrix(c(0.01^2, 0, 0, 0.01^2), nrow = 2, ncol = 2, byrow = TRUE)
# mcmc_pars <- pmcmc_parameters$new(list(beta = beta, gamma = gamma), proposal_matrix, transform = transform)
mcmc_pars <- pmcmc_parameters$new(list(R0 = R0, gamma = gamma), proposal_matrix, transform = transform)

# Run MCMC
control <- pmcmc_control(
    500,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE)
pmcmc_run <- pmcmc(mcmc_pars, filter, control = control)

# Plot fitted trajectories
plot_particle_filter(pmcmc_run$trajectories$state, true_history, data_raw$day)

# processed_chains <- pmcmc_thin(pmcmc_run, burnin = 200, thin = 2)
# parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
# parameter_mean_hpd

mcmc1 <- as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars))
summary(mcmc1)
plot(mcmc1)

# Calculate effective sample size and acceptance rate
effectiveSize(mcmc1)
1 - rejectionRate(mcmc1)

## Tune MCMC
proposal_matrix <- cov(pmcmc_run$pars)
# mcmc_pars <- pmcmc_parameters$new(list(beta = beta, gamma = gamma),
#                                   proposal_matrix,transform = transform)
mcmc_pars <- pmcmc_parameters$new(list(R0 = R0, gamma = gamma),
                                  proposal_matrix,transform = transform)

# Run tuned MCMC
pmcmc_tuned_run <- pmcmc(mcmc_pars, filter, control = control)

# Plot fitted trajectories
plot_particle_filter(pmcmc_tuned_run$trajectories$state, true_history, data_raw$day)

mcmc2 <- as.mcmc(cbind(pmcmc_tuned_run$probabilities, pmcmc_tuned_run$pars))
summary(mcmc2)
plot(mcmc2)

# Calculate effective sample size and acceptance rate
effectiveSize(mcmc2)
1 - rejectionRate(mcmc2)

