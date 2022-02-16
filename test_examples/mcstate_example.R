library(dust)
library(mcstate)

# A basic SIR model included in the dust package
gen <- dust::dust_example("sir")

# Some data that we will fit to, using 1 particle:
sir <- gen$new(pars = list(), step = 0, n_particles = 1)
dt <- 1 / 4
day <- seq(1, 100)
incidence <- rep(NA, length(day))
true_history <- array(NA_real_, c(5, 1, 101))
true_history[, 1, 1] <- sir$state()
for (i in day) {
    state_start <- sir$state()
    sir$run(i / dt)
    state_end <- sir$state()
    true_history[, 1, i + 1] <- state_end
    # Reduction in S
    incidence[i] <- state_start[1, 1] - state_end[1, 1]
}

# Convert this into our required format:
data_raw <- data.frame(day = day, incidence = incidence)
data <- particle_filter_data(data_raw, "day", 4)

# A comparison function
compare <- function(state, observed, pars = NULL) {
    if (is.null(pars$exp_noise)) {
        exp_noise <- 1e6
    } else {
        exp_noise <- pars$exp_noise
    }
    incidence_modelled <- state[5,]
    incidence_observed <- observed$incidence
    lambda <- incidence_modelled +
        rexp(length(incidence_modelled), exp_noise)
    dpois(incidence_observed, lambda, log = TRUE)
}

# Construct the particle_filter object with 100 particles
p <- particle_filter$new(data, gen, 100, compare)
p$run(save_history = TRUE)

# Our simulated trajectories, with the "real" data superimposed
history <- p$history()
matplot(data_raw$day, t(history[1, , -1]), type = "l",
        xlab = "Time", ylab = "State",
        col = "#ff000022", lty = 1, ylim = range(history))
matlines(data_raw$day, t(history[2, , -1]), col = "#ffff0022", lty = 1)
matlines(data_raw$day, t(history[3, , -1]), col = "#0000ff22", lty = 1)
matpoints(data_raw$day, t(true_history[1:3, , -1]), pch = 19,
          col = c("red", "yellow", "blue"))