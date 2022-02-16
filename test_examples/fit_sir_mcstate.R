library(odin)
library(odin.dust)
library(dust)
library(mcstate)
library(socialmixr)
library(coda)
library(scales)

## SIR example from "sir_models" vignette in mcstate package
gen_sir <- dust_example("sir")

incidence <- read.csv(system.file("sir_incidence.csv", package = "mcstate"))
dt <- 0.25
sir_data <- mcstate::particle_filter_data(data = incidence,
                                          time = "day",
                                          rate = 1 / dt)

case_compare <- function(state, observed, pars = NULL) {
    exp_noise <- 1e6
    
    incidence_modelled <- state[5, , drop = TRUE]
    incidence_observed <- observed$cases
    lambda <- incidence_modelled +
        rexp(n = length(incidence_modelled), rate = exp_noise)
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
}

filter <- particle_filter$new(data = sir_data,
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
    matplot(times, t(history[1, , -1]), type = "l",
            xlab = "Time", ylab = "Number of individuals",
            col = alpha(cols[["S"]],0.1), lty = 1, ylim = range(history))
    matlines(times, t(history[2, , -1]), col = alpha(cols[["I"]],0.1), lty = 1)
    matlines(times, t(history[3, , -1]), col = alpha(cols[["R"]],0.1), lty = 1)
    matpoints(times[1:obs_end], t(true_history[1:3, , -1]), pch = 19,
              col = cols)
    legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
}

beta <- pmcmc_parameter("beta", 0.2, min = 0)
gamma <- pmcmc_parameter("gamma", 0.1, min = 0,
                         prior = function(p)
                             dgamma(p, shape = 1, scale = 0.2, log = TRUE))

proposal_matrix <- matrix(c(0.01^2, 0, 0, 0.01^2), nrow = 2, ncol = 2, byrow = TRUE)
mcmc_pars <- pmcmc_parameters$new(list(beta = beta, gamma = gamma), proposal_matrix)

control <- pmcmc_control(
    500,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE)
pmcmc_run <- pmcmc(mcmc_pars, filter, control = control)

true_history <- readRDS("test_examples/sir_true_history.rds")
plot_particle_filter(pmcmc_run$trajectories$state, true_history, incidence$day)

processed_chains <- pmcmc_thin(pmcmc_run, burnin = 200, thin = 2)
parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
parameter_mean_hpd

mcmc1 <- as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars))
summary(mcmc1)
plot(mcmc1)
