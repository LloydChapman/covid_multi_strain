library(odin)
library(odin.dust)
library(dust)
library(mcstate)

gen_seir <- odin.dust::odin_dust("test_examples/seir.R")

seir_model <- gen_seir$new(list(),step = 0,n_particles = 1,n_threads = 1,seed = 2)

# seir_model$state()

# Run epidemic forward
n_particles <- 10
n_steps <- 200

x <- array(NA, dim = c(seir_model$info()$len, n_particles, n_steps))
for (t in seq_len(n_steps)) {
  x[ , , t] <- seir_model$run(t)
}
time <- x[1, 1, ]
x <- x[-1, , ]

# Plot trajectories
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
cols <- c(S = "#8c8cd9", E = "#ffff00", I = "#cc0044", R = "#999966")
matplot(time, t(x[1, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = range(x))
matlines(time, t(x[2, , ]), col = cols[["E"]], lty = 1)
matlines(time, t(x[3, , ]), col = cols[["I"]], lty = 1)
matlines(time, t(x[4, , ]), col = cols[["R"]], lty = 1)
legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
