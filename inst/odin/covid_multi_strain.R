## Compartments E to R stratified by i, j and k with
## i for the age group
## j for the strain
## j = 1: infection with strain 1
## j = 2: infection with strain 2
## j = 3: infection with strain 1 followed by strain 2
## j = 4: infection with strain 2 followed by strain 1
## k for the vaccine group

## Number of age and vaccination groups
n_age <- user()
n_strains <- user()
n_vax <- user()

## Definition of the time-step and output as "time" 
dt <- user()
initial(time) <- 0
update(time) <- (step + 1) * dt
steps_per_day <- 1/dt

## Core equations for transitions between compartments:
update(S[, ]) <- S[i,j] + sum(n_RS[i, ,j]) - sum(n_SE[i, ,j]) - n_SV[i,j] + 
    (if (j == 1) n_SV[i,n_vax] else n_SV[i,j-1]) - 
    (if (vacc_skip_to[j] > 0) n_SV_skip[i,j] else 0) + 
    (if (vacc_skip_from[j] > 0) n_SV_skip[i,vacc_skip_from[j]] else 0)
update(E[, , ]) <- E[i,j,k] + n_SE[i,j,k] + (if (j > 2) n_RE[i,j-2,k] else 0) - 
    n_EI[i,j,k] - n_EV[i,j,k] + 
    (if (k == 1) n_EV[i,j,n_vax] else n_EV[i,j,k-1]) - 
    (if (vacc_skip_to[k] > 0) n_EV_skip[i,j,k] else 0) +
    (if (vacc_skip_from[k] > 0) n_EV_skip[i,j,vacc_skip_from[k]] else 0)
update(I_P[, , ]) <- new_I_P[i,j,k]
update(I_A[, , ]) <- new_I_A[i,j,k]
update(I_C[, , ]) <- new_I_C[i,j,k]
update(R[, , ]) <- R[i,j,k] + n_I_AR[i,j,k] + n_I_CR[i,j,k] + n_HR[i,j,k] - 
    n_Rx[i,j,k] - n_RV[i,j,k] + 
    (if (k == 1) n_RV[i,j,n_vax] else n_RV[i,j,k-1]) - 
    (if (vacc_skip_to[k] > 0) n_RV_skip[i,j,k] else 0) + 
    (if (vacc_skip_from[k] > 0) n_RV_skip[i,j,vacc_skip_from[k]] else 0)
update(H[, , ]) <- H[i,j,k] + n_I_CH[i,j,k] - n_HR[i,j,k] - n_HD[i,j,k]
update(G[, , ]) <- G[i,j,k] + n_I_CG[i,j,k] - n_GD[i,j,k]
update(D[, , ]) <- D[i,j,k] + n_HD[i,j,k] + n_GD[i,j,k]
update(T_pre_1[, , ]) <- T_pre_1[i,j,k] + n_EI[i,j,k] - n_T_pre_1x[i,j,k]
update(T_P_1[, , ]) <- T_P_1[i,j,k] + n_T_pre_1T_P_1[i,j,k] - n_T_P_1T_N_1[i,j,k]
update(T_N_1[, , ]) <- T_N_1[i,j,k] + n_T_pre_1x[i,j,k] - n_T_pre_1T_P_1[i,j,k] + n_T_P_1T_N_1[i,j,k]

update(N_tot[]) <- sum(S[i,]) + sum(E[i, , ]) + sum(I_P[i, , ]) + 
    sum(I_A[i, , ]) + sum(I_C[i, , ]) + sum(R[i, , ]) + sum(H[i, , ]) + 
    sum(G[i, , ]) + sum(D[i, , ])

## Individual probabilities of transition:
# infection progression:
p_SE[, ] <- 1 - exp(-sum(lambda_sus[i, ,j]) * dt) # S to E, sum susceptibility-weighted force of infection over strains
p_EI <- 1 - exp(-gamma_E * dt) # E to I
p_I_PI_C <- 1 - exp(-gamma_P * dt) # I_P to I_C
p_I_Cx <- 1 - exp(-gamma_C * dt) # I_C to x
p_I_AR <- 1 - exp(-gamma_A * dt) # I_A to R
p_Hx <- 1 - exp(-gamma_H * dt) # H to x
p_GD <- 1 - exp(-gamma_G * dt) # G to D
p_T_pre_1x <- 1 - exp(-gamma_pre_1 * dt) # T_pre_1 to x
p_T_P_1T_N_1 <- 1 - exp(-gamma_P_1 * dt) # T_P_1 to T_N_1

# vaccination:
p_SV[, ] <- p_V[i,j]
p_EV[, , ] <- p_V[i,k]
p_I_AV[, , ] <- p_V[i,k]
p_I_PV[, , ] <- p_V[i,k]
p_RV[, , ] <- p_V[i,k]

p_SV_skip[, ] <- p_V_skip[i,j]
p_EV_skip[, , ] <- p_V_skip[i,k]
p_I_PV_skip[, , ] <- p_V_skip[i,k]
p_I_AV_skip[, , ] <- p_V_skip[i,k]
p_RV_skip[, , ] <- p_V_skip[i,k]

p_C[] <- user()
p_H[] <- user()
p_G[] <- user()
p_D[] <- user()

## Compute the force of infection
# Multiply numbers of infected individuals for each strain j in each vaccine 
# stratum k by their relative infectiousness
I_rel_inf[, , ] <- rel_infectivity[i,j,k] * strain_transmission[j] * 
    (theta_A*I_A[i,j,k] + I_P[i,j,k] + I_C[i,j,k])
# Calculate contact rate for individual in age group i with infectious 
# individuals in age group j infected with strain k
s_ij[, , ] <- m[i, j] * sum(I_rel_inf[j,k, ])
# Calculate force of infection on individual in age group i from individuals 
# infected with strain j
# P(Strain = 1) := P(Strain = Only 1) + P(Strain = 2->1), same for Strain = 2
lambda[, ] <- beta * (if (n_real_strains == 1) sum(s_ij[i, ,1]) else
    (if (j == 1) sum(s_ij[i, ,1]) + sum(s_ij[i, ,4]) else
        sum(s_ij[i, ,2]) + sum(s_ij[i, ,3])))
# Multiply by relative susceptibility to get susceptibility-weighted force of 
# infection on individual in age group i in vaccine stratum k from individuals 
# infected with strain j
lambda_sus[, , ] <- rel_susceptibility[i,j,k] * lambda[i,j]

## Draws from binomial distributions for numbers changing between 
## compartments:

# Flow out of S:
# infection
# Compute the new infections with multiple strains using nested binomials
# No one can move from S to E3 or E4
n_SE_tot[, ] <- rbinom(S[i,j], p_SE[i,j])
n_SE[, , ] <- if (j == 1 || n_real_strains == 1)
    rbinom(n_SE_tot[i,k], rel_foi_strain[i,j,k]) else
        (if (j == 2) n_SE_tot[i,k] - n_SE[i,1,k] else 0)

# rel_foi_strain is probability of an infection in group i, vaccination class k
# being of strain j
# min() here is to avoid rel_foi_strain going above 1 due to rounding errors
rel_foi_strain[, , ] <- (if (sum(lambda_sus[i, ,k]) == 0)
    (if (j == 1) 1 else 0) else # don't quite understand why rel_foi_strain should be 1 when FOI=0 if j=1? Is it because rel_foi_strain should always sum to 1 (for each i and k)?
        min(lambda_sus[i,j,k]/sum(lambda_sus[i, ,k]),as.numeric(1)))

# vaccination
n_SV[, ] <- rbinom(S[i,j] - sum(n_SE[i, ,j]), p_SV[i,j])
n_SV_skip[, ] <- rbinom(S[i,j] - sum(n_SE[i, ,j]) - n_SV[i,j], p_SV_skip[i,j])

# Flow out of E:
# progression
n_EI[, , ] <- rbinom(E[i,j,k], p_EI)
n_EI_P[, , ] <- rbinom(n_EI[i,j,k], rel_p_sympt[i,j,k]*
                           strain_rel_p_sympt[j]*p_C[i])
n_EI_A[, , ] <- n_EI[i,j,k] - n_EI_P[i,j,k]
# vaccination
n_EV[, , ] <- rbinom(E[i,j,k] - n_EI[i,j,k], p_EV[i,j,k])
n_EV_skip[, , ] <- rbinom(E[i,j,k] - n_EI[i,j,k] - n_EV[i,j,k], p_EV_skip[i,j,k])

# Flow out of I_P:
# progression
n_I_PI_C[, , ] <- rbinom(I_P[i,j,k], p_I_PI_C)
# vaccination
n_I_PV[, , ] <- rbinom(I_P[i,j,k] - n_I_PI_C[i,j,k], p_I_PV[i,j,k])
n_I_PV_skip[, , ] <- rbinom(I_P[i,j,k] - n_I_PI_C[i,j,k], p_I_PV_skip[i,j,k])

# Flow out of I_A:
# progression
n_I_AR[, , ] <- rbinom(I_A[i,j,k], p_I_AR)
# vaccination
n_I_AV[, , ] <- rbinom(I_A[i,j,k] - n_I_AR[i,j,k], p_I_AV[i,j,k])
n_I_AV_skip[, , ] <- rbinom(I_A[i,j,k] - n_I_AR[i,j,k], p_I_AV_skip[i,j,k])

# Flow out of I_C:
# progression
n_I_Cx[, , ] <- rbinom(I_C[i,j,k], p_I_Cx)
n_I_CR[, , ] <- rbinom(n_I_Cx[i,j,k], 1 - rel_p_hosp_if_sympt[i,j,k]*
                           strain_rel_p_hosp_if_sympt[j]*p_H[i]) #
n_I_CHG[, , ] <- n_I_Cx[i,j,k] - n_I_CR[i,j,k]
n_I_CH[, , ] <- rbinom(n_I_CHG[i,j,k], 1 - rel_p_death[i,j,k]*
                           strain_rel_p_death[j]*p_G[i])
n_I_CG[, , ] <- n_I_CHG[i,j,k] - n_I_CH[i,j,k]

# Flow out of H:
n_Hx[, , ] <- rbinom(H[i,j,k], p_Hx)
n_HD[, , ] <- rbinom(n_Hx[i,j,k], rel_p_death[i,j,k]*
                         strain_rel_p_death[j]*p_D[i])
n_HR[, , ] <- n_Hx[i,j,k] - n_HD[i,j,k]

# Flow out of G:
n_GD[, , ] <- rbinom(G[i,j,k], p_GD)

# Flow out of R:
# progression
# rate of progressing from R w/o superinfection is just waning_rate
# if n_strains is 1 then can only go to S w.r. waning_rate
# if in R3 or R4 then can only go to S w.r. waning_rate
# otherwise:
#  R1 and R2 can progress to S (w.r. waning_rate[i]) or
#  R1 can progress to E3 (w.r. strain 2 (3 - 1))
#  R2 can progress to E4 (w.r. strain 1 (3 - 2))
#
# Note that (if n_real_strains == 2)
# cross_immunity[1] is the cross immunity of strain 1 against strain 2
# cross_immunity[2] is the cross immunity of strain 2 against strain 1
rate_Rx[, , ] <- waning_rate +
    if (n_strains == 1 || j > 2) 0 else
        lambda_sus[i,3-j,k] * (1 - cross_immunity[j])

p_Rx[, , ] <- 1 - exp(-rate_Rx[i,j,k] * dt)

# number leaving R for S or E
n_Rx[, , ] <- rbinom(R[i,j,k], p_Rx[i,j,k])

# number going from R to S
# In one-strain model, all progress to S only
# In multi-strain model, R3 and R4 progress to S only
# If not modelling super infection then all progress to S only
#
# In multi-strain model, number of R1 and R2 to S is binomial w.p. waning over
# waning plus prob strain
p_RS[, , ] <- if (n_strains == 1 || j > 2) 1 else
    (if (waning_rate == 0) 0 else
        (waning_rate/(waning_rate + lambda_sus[i,3-j,k] * (1 - cross_immunity[j]))))
n_RS[, , ] <- rbinom(n_Rx[i,j,k], p_RS[i,j,k])

# n_RE[i, j, k] is the number going from
# R[age i, strain j, vacc class k] to
# E[age i, strain j + 2, vacc class k]
# Movement possible to E3 and E4 are from R1 and R2 respectively
n_RE[, ,] <- n_Rx[i,j,k] - n_RS[i,j,k]

# vaccination
n_RV[, , ] <- rbinom(R[i,j,k] - n_Rx[i,j,k], p_RV[i,j,k])
# n_RV[, , ] <- rbinom(R[i,j,k], p_RV[i,j,k])
n_RV_skip[, , ] <- rbinom(R[i,j,k] - n_Rx[i,j,k] - n_RV[i,j,k], p_RV_skip[i,j,k])

# Number vaccinated by age and vaccine stage

# Note that if an individual makes a vaccine skip move, then in the following
# they are counted as making as having moved through each of the intermediate
# classes they skip over, even though they do not spend any time in them.
# Note this should not affect model dynamics, the assumption is purely for
# bookkeeping purposes
n_S_vaccinated[, ] <- n_SV[i,j] + 
    (if (vacc_skipped[j] > 0) n_SV_skip[i,vacc_skipped[j]] else 0)
n_E_vaccinated[, ] <- sum(n_EV[i, ,j]) + 
    (if (vacc_skipped[j] > 0) sum(n_EV[i, ,vacc_skipped[j]]) else 0)
n_I_P_vaccinated[, ] <- sum(n_I_PV[i, ,j]) + 
    (if (vacc_skipped[j] > 0) sum(n_I_PV_skip[i, ,vacc_skipped[j]]) else 0)
n_I_A_vaccinated[, ] <- sum(n_I_AV[i, ,j]) + 
    (if (vacc_skipped[j] > 0) sum(n_I_AV_skip[i, ,vacc_skipped[j]]) else 0)
n_R_vaccinated[, ] <- sum(n_RV[i, ,j]) + 
    (if (vacc_skipped[j] > 0) sum(n_RV_skip[i, ,vacc_skipped[j]]) else 0)

n_vaccinated[, ] <- 
    n_S_vaccinated[i,j] + 
    n_E_vaccinated[i,j] + 
    n_I_P_vaccinated[i,j] +
    n_I_A_vaccinated[i,j] + 
    n_R_vaccinated[i,j]

# Parallel flow for serological testing
n_T_pre_1x[, , ] <- rbinom(T_pre_1[i,j,k], p_T_pre_1x)
n_T_pre_1T_P_1[, , ] <- rbinom(n_T_pre_1x[i,j,k], p_P_1)
n_T_P_1T_N_1[, , ] <- rbinom(T_P_1[i,j,k], p_T_P_1T_N_1)

# number of new I_A
new_I_A[, , ] <- I_A[i,j,k] + n_EI_A[i,j,k] - n_I_AR[i,j,k] - n_I_AV[i,j,k] + 
    (if (k == 1) n_I_AV[i,j,n_vax] else n_I_AV[i,j,k-1]) - 
    (if (vacc_skip_to[k] > 0) n_I_AV_skip[i,j,k] else 0) + 
    (if (vacc_skip_from[k] > 0) n_I_AV_skip[i,j,vacc_skip_from[k]] else 0)

# number of new I_P
new_I_P[, , ] <- I_P[i,j,k] + n_EI_P[i,j,k] - n_I_PI_C[i,j,k] - n_I_PV[i,j,k] +
    (if (k == 1) n_I_PV[i,j,n_vax] else n_I_PV[i,j,k-1]) - 
    (if (vacc_skip_to[k] > 0) n_I_PV_skip[i,j,k] else 0) + 
    (if (vacc_skip_from[k] > 0) n_I_PV_skip[i,j,vacc_skip_from[k]] else 0)

# number of new I_C
new_I_C[, , ] <- I_C[i,j,k] + n_I_PI_C[i,j,k] - n_I_Cx[i,j,k]

## Initial conditions 
initial(S[, ]) <- 0
initial(E[, , ]) <- 0
initial(I_P[, , ]) <- 0
initial(I_A[, , ]) <- 0
initial(I_C[, , ]) <- 0
initial(R[, , ]) <- 0
initial(H[, , ]) <- 0
initial(G[, , ]) <- 0
initial(D[, , ]) <- 0
initial(T_pre_1[, , ]) <- 0
initial(T_P_1[, , ]) <- 0
initial(T_N_1[, , ]) <- 0

initial(N_tot[]) <- 0

## Seeding
# First strain
seed_step_end <- seed_step_start + length(seed_value)
seed_rate <- if (step >= seed_step_start && step < seed_step_end)
    seed_value[as.integer(step - seed_step_start + 1)] else 0
seed <- rpois(seed_rate)

n_SE[seed_age,1,1] <- n_SE[seed_age,1,1] + min(S[seed_age,1] - n_SE_tot[seed_age,1], seed)

# Introduction of new strain
strain_seed_step_end <- strain_seed_step_start + length(strain_seed_value)
strain_seed_rate <- if (step >= strain_seed_step_start && step < strain_seed_step_end)
    strain_seed_value[as.integer(step - strain_seed_step_start + 1)] else 0
strain_seed <- rpois(strain_seed_rate)

n_SE[seed_age,2:n_strains,1] <- 
    if (j < 3) min(n_SE[i,j,k] + S[i,k] - sum(n_SE[i, ,k]), # can't move more individuals than are in this S category
                   n_SE[i,j,k] + strain_seed) else 0 # ensure no individuals move from S to strain E3 or E4

## Store hospitalisation and death incidence
initial(H_inc) <- 0
initial(D_inc) <- 0
update(H_inc) <- if (step %% steps_per_day == 0) sum(n_I_CH) else H_inc + sum(n_I_CH)
update(D_inc) <- if (step %% steps_per_day == 0) sum(n_HD) else D_inc + sum(n_HD)

## Store hospitalisation and death incidence by age
initial(H_inc_0_39) <- 0
update(H_inc_0_39) <- if (step %% steps_per_day == 0) sum(n_I_CH[1:4, , ]) else H_inc_0_39 + sum(n_I_CH[1:4, , ])
initial(H_inc_40_49) <- 0
update(H_inc_40_49) <- if (step %% steps_per_day == 0) sum(n_I_CH[5, , ]) else H_inc_40_49 + sum(n_I_CH[5, , ])
initial(H_inc_50_59) <- 0
update(H_inc_50_59) <- if (step %% steps_per_day == 0) sum(n_I_CH[6, , ]) else H_inc_50_59 + sum(n_I_CH[6, , ])
initial(H_inc_60_69) <- 0
update(H_inc_60_69) <- if (step %% steps_per_day == 0) sum(n_I_CH[7, , ]) else H_inc_60_69 + sum(n_I_CH[7, , ])
initial(H_inc_70_plus) <- 0
update(H_inc_70_plus) <- if (step %% steps_per_day == 0) sum(n_I_CH[8, , ]) else H_inc_70_plus + sum(n_I_CH[8, , ])

initial(D_inc_0_39) <- 0
update(D_inc_0_39) <- if (step %% steps_per_day == 0) sum(n_HD[1:4, , ]) else D_inc_0_39 + sum(n_HD[1:4, , ])
initial(D_inc_40_49) <- 0
update(D_inc_40_49) <- if (step %% steps_per_day == 0) sum(n_HD[5, , ]) else D_inc_40_49 + sum(n_HD[5, , ])
initial(D_inc_50_59) <- 0
update(D_inc_50_59) <- if (step %% steps_per_day == 0) sum(n_HD[6, , ]) else D_inc_50_59 + sum(n_HD[6, , ])
initial(D_inc_60_69) <- 0
update(D_inc_60_69) <- if (step %% steps_per_day == 0) sum(n_HD[7, , ]) else D_inc_60_69 + sum(n_HD[7, , ])
initial(D_inc_70_plus) <- 0
update(D_inc_70_plus) <- if (step %% steps_per_day == 0) sum(n_HD[8, , ]) else D_inc_70_plus + sum(n_HD[8, , ])

# Store seropositive number
initial(sero_pos_1) <- 0
update(sero_pos_1) <- sum(T_P_1) + sum(n_T_pre_1T_P_1) - sum(n_T_P_1T_N_1)

# Store seropositive numbers by age
initial(sero_pos_1_20_29) <- 0
update(sero_pos_1_20_29) <- sum(T_P_1[3, , ]) + sum(n_T_pre_1T_P_1[3, , ]) - sum(n_T_P_1T_N_1[3, , ])
initial(sero_pos_1_30_39) <- 0
update(sero_pos_1_30_39) <- sum(T_P_1[4, , ]) + sum(n_T_pre_1T_P_1[4, , ]) - sum(n_T_P_1T_N_1[4, , ])
initial(sero_pos_1_40_49) <- 0
update(sero_pos_1_40_49) <- sum(T_P_1[5, , ]) + sum(n_T_pre_1T_P_1[5, , ]) - sum(n_T_P_1T_N_1[5, , ])
initial(sero_pos_1_50_59) <- 0
update(sero_pos_1_50_59) <- sum(T_P_1[6, , ]) + sum(n_T_pre_1T_P_1[6, , ]) - sum(n_T_P_1T_N_1[6, , ])
initial(sero_pos_1_60_69) <- 0
update(sero_pos_1_60_69) <- sum(T_P_1[7, , ]) + sum(n_T_pre_1T_P_1[7, , ]) - sum(n_T_P_1T_N_1[7, , ])
initial(sero_pos_1_70_plus) <- 0
update(sero_pos_1_70_plus) <- sum(T_P_1[8, , ]) + sum(n_T_pre_1T_P_1[8, , ]) - sum(n_T_P_1T_N_1[8, , ])

# Case incidence 
new_cases <- sum(n_EI_P)
initial(cases_inc) <- 0
update(cases_inc) <- (
    if (step %% steps_per_day == 0) new_cases
    else cases_inc + new_cases)

new_cases_non_variant <- sum(n_EI_P[,1, ]) +
    (if (n_real_strains == 2) sum(n_EI_P[,4, ]) else 0)
initial(cases_non_variant_inc) <- 0
update(cases_non_variant_inc) <- (
    if (step %% steps_per_day == 0) new_cases_non_variant
    else cases_non_variant_inc + new_cases_non_variant)

new_cases_0_9 <- sum(n_EI_P[1, , ])
initial(cases_inc_0_9) <- 0
update(cases_inc_0_9) <- ( 
    if (step %% steps_per_day == 0) new_cases_0_9
    else cases_inc_0_9 + new_cases_0_9)

new_cases_10_19 <- sum(n_EI_P[2, , ])
initial(cases_inc_10_19) <- 0
update(cases_inc_10_19) <- ( 
    if (step %% steps_per_day == 0) new_cases_10_19
    else cases_inc_10_19 + new_cases_10_19)

new_cases_20_29 <- sum(n_EI_P[3, , ])
initial(cases_inc_20_29) <- 0
update(cases_inc_20_29) <- ( 
    if (step %% steps_per_day == 0) new_cases_20_29
    else cases_inc_20_29 + new_cases_20_29)

new_cases_30_39 <- sum(n_EI_P[4, , ])
initial(cases_inc_30_39) <- 0
update(cases_inc_30_39) <- ( 
    if (step %% steps_per_day == 0) new_cases_30_39
    else cases_inc_30_39 + new_cases_30_39)

new_cases_40_49 <- sum(n_EI_P[5, , ])
initial(cases_inc_40_49) <- 0
update(cases_inc_40_49) <- ( 
    if (step %% steps_per_day == 0) new_cases_40_49
    else cases_inc_40_49 + new_cases_40_49)

new_cases_50_59 <- sum(n_EI_P[6, , ])
initial(cases_inc_50_59) <- 0
update(cases_inc_50_59) <- ( 
    if (step %% steps_per_day == 0) new_cases_50_59
    else cases_inc_50_59 + new_cases_50_59)

new_cases_60_69 <- sum(n_EI_P[7, , ])
initial(cases_inc_60_69) <- 0
update(cases_inc_60_69) <- ( 
    if (step %% steps_per_day == 0) new_cases_60_69
    else cases_inc_60_69 + new_cases_60_69)

new_cases_70_plus <- sum(n_EI_P[8, , ])
initial(cases_inc_70_plus) <- 0
update(cases_inc_70_plus) <- ( 
    if (step %% steps_per_day == 0) new_cases_70_plus
    else cases_inc_70_plus + new_cases_70_plus)

## Model parameters (default in parenthesis) 
# Transmission and progression
beta_step[] <- user()
beta <- if (as.integer(step) >= length(beta_step))
    beta_step[length(beta_step)] else beta_step[step + 1]
gamma_E <- user(0.5)
gamma_P <- user(0.4)
gamma_C <- user(0.4)
gamma_A <- user(0.2)
gamma_H <- user(0.1)
gamma_G <- user(1/3)
gamma_pre_1 <- user(1/13)
gamma_P_1 <- user(1/200)
p_P_1 <- user(0.85)
waning_rate <- user()
theta_A <- user(0.5)
m[, ] <- user() # age-structured contact matrix
rel_infectivity[, , ] <- user() # average vaccine effectiveness against infection by age group and vaccine stratum
strain_transmission[] <- user() # relative transmissibility of different strains
rel_susceptibility[, , ] <- user() # infectiousness by age group and vaccine stratum
cross_immunity[] <- user() # cross-immunity between 1st and 2nd strain

# Seeding
seed_age <- user(4) # 30-39 years band
seed_step_start <- user() # step at which to start seeding
seed_value[] <- user() # number of infections seeded into seed_age age group over length(seed_value) steps
strain_seed_step_start <- user() # step at which to start seeding of 2nd strain
strain_seed_value[] <- user() # number of infections seeded into seed_age age group over length(strain_seed_value) steps
n_real_strains <- if (n_strains == 4) 2 else 1 # number of real strains for history-based multi-strain model of infection

# Vaccination
rel_p_sympt[, , ] <- user()
rel_p_hosp_if_sympt[, , ] <- user()
rel_p_death[, , ] <- user()
strain_rel_p_sympt[] <- user()
strain_rel_p_hosp_if_sympt[] <- user()
strain_rel_p_death[] <- user()

n_doses <- user()
index_dose[] <- user(integer = TRUE)
index_dose_inverse[] <- user(integer = TRUE)
vaccine_dose_step[, , ] <- user() # n_age x n_doses x n_time

# Number of vaccination candidates
vaccine_n_candidates[, ] <-
    S[i,index_dose[j]] +
    sum(E[i, ,index_dose[j]]) +
    sum(I_A[i, ,index_dose[j]]) +
    sum(I_P[i, ,index_dose[j]]) +
    sum(R[i, ,index_dose[j]])

vacc_skip_n_candidates[, ] <-
    (if (vacc_skip_dose[j] > 0)
        S[i,vacc_skip_dose[j]] +
         sum(E[i, ,vacc_skip_dose[j]]) +
         sum(I_A[i, ,vacc_skip_dose[j]]) +
         sum(I_P[i, ,vacc_skip_dose[j]]) +
         sum(R[i, ,vacc_skip_dose[j]])
     else 0)

# Work out the vaccination probability via doses, driven by the
# schedule
vaccine_probability_doses[, ] <- min(
    if (vaccine_n_candidates[i, j] > 0)
        vaccine_attempted_doses[i, j] / vaccine_n_candidates[i, j] else 0,
    as.numeric(1))

# Work out the total attempted doses
total_attempted_doses[, ] <- vaccine_missed_doses[i, j] + (
    if (as.integer(step) >= dim(vaccine_dose_step, 3)) 0
    else vaccine_dose_step[i, j, step + 1])

# Now we work out the split of the total attempted doses, firstly for the
# next vaccine class moves, then for the vaccine skip moves.
#
# Note attempted doses for next vaccine class moves competing with a vaccine
# skip are weighted here, with remaining doses made available to the vaccine
# skip move below.
#
# Note that since we require vacc_skip_dose_weight <= 1, it is not possible for
# there to be an excess of doses for vaccine skip moves while having not
# enough doses for the next vaccine class candidates.
vaccine_attempted_doses[, ] <-
    (if (vaccine_n_candidates[i, j] == 0) 0
     else
         min(vaccine_n_candidates[i, j] /
                 (vaccine_n_candidates[i, j] +
                      vacc_skip_dose_weight[j] * vacc_skip_n_candidates[i, j])
             * total_attempted_doses[i, j],
             vaccine_n_candidates[i, j]))

vacc_skip_attempted_doses[, ] <-
    (if (vacc_skip_dose_weight[j] > 0)
        (if (vacc_skip_dose[j] > 0)
            total_attempted_doses[i, j] -
             vaccine_attempted_doses[i, j]
         else 0)
     else 0)

# Calculate catch-up vaccine doses, i.e. doses not delivered according to
# schedule, e.g. because too many individuals are in the I or H states, that
# are postponed till later
initial(vaccine_missed_doses[, ]) <- 0
update(vaccine_missed_doses[, ]) <-
    vaccine_catchup_fraction *
    max(total_attempted_doses[i, j] - n_vaccinated[i, index_dose[j]],
        as.numeric(0))
vaccine_catchup_fraction <- user(1)

# Calculate probability of vaccination with a certain dose based on progression
# at a constant rate (vaccine_progression_rate_base) or the supplied
# time-varying probabilities (vaccine_probability_doses)
p_V[, ] <- 
    (if (index_dose_inverse[j] > 0)
        vaccine_probability_doses[i, index_dose_inverse[j]] 
     else
        1 - exp(-vaccine_progression_rate_base[i, j] * dt))
vaccine_progression_rate_base[, ] <- user()

p_V_skip[, ] <- 
    (if (vacc_skip_dose_inverse[j] > 0)
        (if (vacc_skip_n_candidates[i, vacc_skip_dose_inverse[j]] > 0)
            min(vacc_skip_attempted_doses[i, vacc_skip_dose_inverse[j]]/
                    vacc_skip_n_candidates[i, vacc_skip_dose_inverse[j]], 
                as.numeric(1))
         else 0)
     else
         1 - exp(-vacc_skip_progression_rate_base[j] * dt))

# Vaccine skip inputs
# 1. vacc_skip_to[j] is the vaccine stratum that the vaccine skip move
#    from stratum j goes to (0 represents no vaccine skip move from j)
# 2. vacc_skip_from[j] is the vaccine stratum that the vaccine skip move
#    to stratum j comes from (0 represents no vaccine skip move to j)
# 3. vacc_skip_progression_rate_base[j] is the progression rate used for
#    the vaccine skip from stratum j (unless the vaccine skip move is
#    controlled by doses)
# 4. vacc_skip_dose[j] is the vaccine stratum that the vaccine skip move comes
#    from that is controlled by dose j (0 represents no vaccine skip
#    move controlled by dose j)
# 5. vacc_skip_dose_inverse[j] is the dose that the vaccine skip move from
#    stratum j is controlled by (0 represents that no dose controls the vaccine
#    skip move)
# 6. vacc_skip_dose_weight[j] represents how much vaccine skip candidates are
#    weighted for distribution of dose j relative to standard vaccine
#    candidates for that dose
# 7. vacc_skipped is used for bookkeeping - vacc_skipped[j] is the vaccine
#    stratum a vaccine skip move goes from that either starts at j or skips
#    over j (if there is no such move then the vacc_skipped[j] is 0)
vacc_skip_to[] <- user(integer = TRUE)
vacc_skip_from[] <- user(integer = TRUE)
vacc_skip_progression_rate_base[] <- user()
vacc_skip_dose[] <- user(integer = TRUE)
vacc_skip_dose_inverse[] <- user(integer = TRUE)
vacc_skip_dose_weight[] <- user()
vacc_skipped[] <- user(integer = TRUE)

# Calculation of I_weighted: infectious-weighted number of infectious individuals
new_I_weighted[, , ] <- 
    theta_A * new_I_A[i,j,k] +
    new_I_P[i,j,k] +
    new_I_C[i,j,k]
sum_new_I_weighted <- sum(new_I_weighted)
initial(I_weighted[, , ]) <- 0
# If there are zero infectious individuals default to putting weight
# in group 4/strain 1/vaccine stratum 1 to avoid NAs in Rt
update(I_weighted[, , ]) <- 
    (if (sum_new_I_weighted == 0)
        (if (i == seed_age && j == 1 && k == 1) 1 else 0)
     else new_I_weighted[i,j,k])

# prob_strain is proportion of total I_weighted in each strain
# If there are zero infectious individuals, default to full weight on strain 1
# to avoid NAs in Rt
prob_strain_1 <- if (n_real_strains == 1 || sum_new_I_weighted == 0) 1 else
    (sum(new_I_weighted[, 1, ]) + sum(new_I_weighted[, 4, ])) /
    sum_new_I_weighted
initial(prob_strain[1:n_real_strains]) <- 0
initial(prob_strain[1]) <- 1
update(prob_strain[]) <- if (i == 1) prob_strain_1 else 1 - prob_strain_1

# Calculate effective susceptibles to each strain
# Weight each person in S/R by their relative susceptibility
# Note that for those in R we further account for cross immunity
# to strains. Those recovered from strain 3 - j will be (partially)
# susceptible to strain j
eff_S[, , ] <- if (n_real_strains == 1)
    S[i, k] * rel_susceptibility[i, j, k] else
        (S[i, k] + (1 - cross_immunity[3 - j]) * R[i, 3 - j, k]) *
    rel_susceptibility[i, j, k]

initial(effective_susceptible[, ]) <- 0
update(effective_susceptible[, ]) <- sum(eff_S[i, j, ])

# Output new infections
initial(new_infections[, , ]) <- 0
update(new_infections[, , ]) <- n_SE[i, j, k]

# Output new reinfections
initial(new_reinfections[, , ]) <- 0
update(new_reinfections[, , ]) <- n_RE[i, j, k]

## Array dimensions
dim(S) <- c(n_age,n_vax)
dim(E) <- c(n_age,n_strains,n_vax)
dim(I_P) <- c(n_age,n_strains,n_vax)
dim(I_A) <- c(n_age,n_strains,n_vax)
dim(I_C) <- c(n_age,n_strains,n_vax)
dim(R) <- c(n_age,n_strains,n_vax)
dim(H) <- c(n_age,n_strains,n_vax)
dim(G) <- c(n_age,n_strains,n_vax)
dim(D) <- c(n_age,n_strains,n_vax)
dim(T_pre_1) <- c(n_age,n_strains,n_vax)
dim(T_P_1) <- c(n_age,n_strains,n_vax)
dim(T_N_1) <- c(n_age,n_strains,n_vax)
dim(N_tot) <- n_age
dim(n_SE_tot) <- c(n_age,n_vax)
dim(n_SE) <- c(n_age,n_strains,n_vax)
dim(n_EI) <- c(n_age,n_strains,n_vax)
dim(n_EI_P) <- c(n_age,n_strains,n_vax)
dim(n_EI_A) <- c(n_age,n_strains,n_vax)
dim(n_I_PI_C) <- c(n_age,n_strains,n_vax)
dim(n_I_AR) <- c(n_age,n_strains,n_vax)
dim(n_I_Cx) <- c(n_age,n_strains,n_vax)
dim(n_I_CR) <- c(n_age,n_strains,n_vax)
dim(n_I_CHG) <- c(n_age,n_strains,n_vax)
dim(n_I_CH) <- c(n_age,n_strains,n_vax)
dim(n_I_CG) <- c(n_age,n_strains,n_vax)
dim(n_Hx) <- c(n_age,n_strains,n_vax)
dim(n_HD) <- c(n_age,n_strains,n_vax)
dim(n_HR) <- c(n_age,n_strains,n_vax)
dim(n_GD) <- c(n_age,n_strains,n_vax)
dim(n_Rx) <- c(n_age,n_strains,n_vax)
dim(n_RS) <- c(n_age,n_strains,n_vax)
dim(n_RE) <- c(n_age,n_strains,n_vax)
dim(n_SV) <- c(n_age,n_vax)
dim(n_EV) <- c(n_age,n_strains,n_vax)
dim(n_I_PV) <- c(n_age,n_strains,n_vax)
dim(n_I_AV) <- c(n_age,n_strains,n_vax)
dim(n_RV) <- c(n_age,n_strains,n_vax)
dim(n_SV_skip) <- c(n_age,n_vax)
dim(n_EV_skip) <- c(n_age,n_strains,n_vax)
dim(n_I_PV_skip) <- c(n_age,n_strains,n_vax)
dim(n_I_AV_skip) <- c(n_age,n_strains,n_vax)
dim(n_RV_skip) <- c(n_age,n_strains,n_vax)
dim(n_T_pre_1x) <- c(n_age,n_strains,n_vax)
dim(n_T_pre_1T_P_1) <- c(n_age,n_strains,n_vax)
dim(n_T_P_1T_N_1) <- c(n_age,n_strains,n_vax)
dim(n_S_vaccinated) <- c(n_age,n_vax)
dim(n_E_vaccinated) <- c(n_age,n_vax)
dim(n_I_P_vaccinated) <- c(n_age,n_vax)
dim(n_I_A_vaccinated) <- c(n_age,n_vax)
dim(n_R_vaccinated) <- c(n_age,n_vax)
dim(n_vaccinated) <- c(n_age,n_vax)
dim(new_I_A) <- c(n_age,n_strains,n_vax)
dim(new_I_P) <- c(n_age,n_strains,n_vax)
dim(new_I_C) <- c(n_age,n_strains,n_vax)
dim(p_SE) <- c(n_age,n_vax)
dim(lambda) <- c(n_age,n_real_strains)
dim(lambda_sus) <- c(n_age,n_real_strains,n_vax)
dim(beta_step) <- user()
dim(m) <- c(n_age, n_age)
dim(s_ij) <- c(n_age,n_age,n_strains)
dim(rel_infectivity) <- c(n_age,n_strains,n_vax)
dim(strain_transmission) <- n_strains
dim(I_rel_inf) <- c(n_age,n_strains,n_vax)
dim(rel_susceptibility) <- c(n_age,n_strains,n_vax)
dim(rel_foi_strain) <- c(n_age,n_real_strains,n_vax)
dim(cross_immunity) <- n_real_strains
dim(p_C) <- n_age
dim(p_H) <- n_age
dim(p_G) <- n_age
dim(p_D) <- n_age
dim(rel_p_sympt) <- c(n_age,n_strains,n_vax)
dim(rel_p_hosp_if_sympt) <- c(n_age,n_strains,n_vax)
dim(rel_p_death) <- c(n_age,n_strains,n_vax)
dim(strain_rel_p_sympt) <- n_strains
dim(strain_rel_p_hosp_if_sympt) <- n_strains
dim(strain_rel_p_death) <- n_strains
dim(rate_Rx) <- c(n_age,n_strains,n_vax)
dim(p_Rx) <- c(n_age,n_strains,n_vax)
dim(p_RS) <- c(n_age,n_strains,n_vax)
dim(p_V) <- c(n_age,n_vax)
dim(p_V_skip) <- c(n_age,n_vax)
dim(vaccine_progression_rate_base) <- c(n_age,n_vax)
dim(p_SV) <- c(n_age,n_vax)
dim(p_EV) <- c(n_age,n_strains,n_vax)
dim(p_I_PV) <- c(n_age,n_strains,n_vax)
dim(p_I_AV) <- c(n_age,n_strains,n_vax)
dim(p_RV) <- c(n_age,n_strains,n_vax)
dim(p_SV_skip) <- c(n_age,n_vax)
dim(p_EV_skip) <- c(n_age,n_strains,n_vax)
dim(p_I_PV_skip) <- c(n_age,n_strains,n_vax)
dim(p_I_AV_skip) <- c(n_age,n_strains,n_vax)
dim(p_RV_skip) <- c(n_age,n_strains,n_vax)
dim(seed_value) <- user()
dim(strain_seed_value) <- user()
dim(index_dose) <- n_doses
dim(index_dose_inverse) <- n_vax
dim(vaccine_dose_step) <- user()
dim(vaccine_n_candidates) <- c(n_age,n_doses)
dim(vacc_skip_n_candidates) <- c(n_age,n_doses)
dim(vaccine_attempted_doses) <- c(n_age,n_doses)
dim(vacc_skip_attempted_doses) <- c(n_age,n_doses)
dim(vaccine_probability_doses) <- c(n_age,n_doses)
dim(total_attempted_doses) <- c(n_age,n_doses)
dim(vaccine_missed_doses) <- c(n_age,n_doses)
dim(vacc_skip_to) <- n_vax
dim(vacc_skip_from) <- n_vax
dim(vacc_skip_progression_rate_base) <- n_vax
dim(vacc_skip_dose) <- n_doses
dim(vacc_skip_dose_inverse) <- n_vax
dim(vacc_skip_dose_weight) <- n_doses
dim(vacc_skipped) <- n_vax
dim(new_I_weighted) <- c(n_age,n_strains,n_vax)
dim(I_weighted) <- c(n_age,n_strains,n_vax)
dim(prob_strain) <- n_real_strains
dim(eff_S) <- c(n_age,n_real_strains,n_vax)
dim(effective_susceptible) <- c(n_age,n_real_strains)
dim(new_infections) <- c(n_age,n_real_strains,n_vax)
dim(new_reinfections) <- c(n_age,n_real_strains,n_vax)