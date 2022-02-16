library(odin)
library(odin.dust)
library(dust)
library(mcstate)
library(socialmixr)

gen_seird_age <- odin.dust::odin_dust("test_examples/seirdage.R")

# Get age-structured contact matrix
data(polymod, package = "socialmixr")

age.limits <- seq(0,70,10)
contact <- socialmixr::contact_matrix(survey = polymod,countries = "United Kingdom",age.limits = age.limits,symmetric = T)

## Transform the matrix to the (symmetrical) transmission matrix
## rather than the contact matrix. This transmission matrix is
## weighted by the population in each age band.
transmission <- contact$matrix/rep(contact$demography$population, each = ncol(contact$matrix))
# for (i in 1:nrow(transmission)){
#   transmission[i,i] <- transmission[2,2]  
# }
transmission

# Get death probability
IFR <- readRDS("~/UCSF/COVIDVaccineModelling/Data/IFR_by_age_ODriscoll.RDS")
p_death <- IFR$median_perc[19:27]/100
p_death[length(p_death)-1] <- (p_death[length(p_death)-1]+p_death[length(p_death)])/2
p_death <- p_death[1:(length(p_death)-1)]

# Generate model
N_age <- length(age.limits)
n_particles <- 100
dt <- 0.25
seird_age_model <- gen_seird_age$new(
  list(dt = dt,S_ini = contact$demography$population,
       E_ini = c(0,0,0,0,0,0,0,0),I_ini = c(0,10,0,0,0,0,0,0),
       beta = 0.05,sigma = 0.5,m = transmission,p_death = p_death,N_age = N_age),
  step = 0,n_particles = n_particles,
  n_threads = 1,seed = 1)

# Run epidemic forward
n_steps <- 1000

# Create an array to contain outputs after looping the model.
# Array contains XX rows = Total S, E, I, R (4), and
# in each age compartment (XX) as well as the cumulative incidence (XX)
x <- array(NA, dim = c(seird_age_model$info()$len, n_particles, n_steps))

# For loop to run the model iteratively
for (t in seq_len(n_steps)) {
  x[ , , t] <- seird_age_model$run(t)
}
time <- x[1, 1, ]
x <- x[-1, , ]

# Plot trajectories
par(mfrow = c(2,4), oma=c(2,3,0,0))
for (i in 1:N_age) {
  par(mar = c(3, 4, 2, 0.5))
  cols <- c(S = "#8c8cd9", E = "#ffff00", I = "#cc0044", R = "#999966", D = "#000000")
  matplot(time, t(x[i + 5,, ]), type = "l", # Offset to access numbers in age compartment
          xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
          col = cols[["S"]], lty = 1, ylim=range(x[-1:-4,,]))
  matlines(time, t(x[i + 5 + N_age, , ]), col = cols[["E"]], lty = 1)
  matlines(time, t(x[i + 5 + 2*N_age, , ]), col = cols[["I"]], lty = 1)
  matlines(time, t(x[i + 5 + 3*N_age, , ]), col = cols[["R"]], lty = 1)
  matlines(time, t(x[i + 5 + 4*N_age, , ]), col = cols[["D"]], lty = 1)
  legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
  axis(2, las =2)
}
mtext("Number of individuals", side=2,line=1, outer=T)
mtext("Time", side = 1, line = 0, outer =T)
