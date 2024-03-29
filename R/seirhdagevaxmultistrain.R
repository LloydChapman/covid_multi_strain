initial <- function(info, n_particles, pars) {
    index <- info$index
    state <- numeric(info$len)
    
    index_S <- index[["S"]]
    index_S_no_vacc <- index_S[seq_len(length(pars$N_tot))]
    index_N_tot <- index[["N_tot"]]
    
    # S0 is the population totals
    initial_S <- pars$N_tot
    
    state[index_S_no_vacc] <- initial_S
    state[index_N_tot] <- pars$N_tot
    
    state
}


parameters <- function(dt,
                       n_age,
                       n_vax,
                       m,
                       beta,
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
                       start_date = 1L,
                       initial_seed_size = 10,
                       initial_seed_pattern = 1,
                       strain_transmission = 1,
                       strain_seed_date = NULL,
                       strain_seed_size = NULL,
                       strain_seed_pattern = NULL,
                       strain_rel_p_sympt = 1,
                       strain_rel_p_hosp_if_sympt = 1,
                       strain_rel_p_death = 1,
                       rel_susceptibility = 1,
                       rel_p_sympt = 1,
                       rel_p_hosp_if_sympt = 1,
                       rel_p_death = 1,
                       rel_infectivity = 1,
                       vaccine_progression_rate = NULL,
                       vaccine_schedule = NULL,
                       vaccine_index_dose2 = NULL,
                       vaccine_index_booster = NULL,
                       vaccine_catchup_fraction = 1,
                       n_doses = 2L,
                       waning_rate = 0,
                       cross_immunity = 1) {
    
    n_real_strains <- length(strain_transmission)
    
    if (n_real_strains > 2) {
        stop("Only 1 or 2 strains valid ('strain_transmission' too long)'.")
    }
    
    # Transmission and progression parameters
    p <- list(dt = dt,
              n_age = n_age,
              n_vax = n_vax,
              m = m,
              beta = beta,
              gamma_E = gamma_E,
              gamma_P = gamma_P,
              gamma_A = gamma_A,
              gamma_C = gamma_C,
              gamma_H = gamma_H, 
              gamma_G = gamma_G,
              theta_A = theta_A,
              p_C = p_C,
              p_H = p_H,
              p_G = p_G,
              p_D = p_D)
    
    # Waning
    p$waning_rate <- waning_rate
    
    # Total population
    p$N_tot <- population
    
    # Seeding parameters
    start_step <- start_date/dt
    p$seed_step_start <- floor(start_step)
    p$seed_value <- seed_over_steps(start_step, initial_seed_pattern) * initial_seed_size
    
    # Cross-immunity
    p$cross_immunity <- recycle(cross_immunity,n_real_strains)
    
    # Number of strains and relative transmissibility
    strain <- strain_parameters(strain_transmission, strain_seed_date, 
                                strain_seed_size, strain_seed_pattern,p$dt)
    
    # Make example vaccine schedule
    pop_mat <- matrix(rep(population,1),nrow = length(population))
    schedule <- vaccine_schedule_future(0, daily_doses, mean_days_between_doses, pop_mat)
    
    # Construct vaccination parameters
    vaccination <- vaccination_parameters(p$N_tot,
                                          p$dt,
                                          rel_susceptibility,
                                          rel_p_sympt,
                                          rel_p_hosp_if_sympt,
                                          rel_p_death,
                                          rel_infectivity,
                                          vaccine_progression_rate,
                                          vaccine_schedule,
                                          vaccine_index_dose2,
                                          vaccine_index_booster,
                                          strain$n_strains, # now up to 2 strains
                                          vaccine_catchup_fraction,
                                          n_doses)
    
    # Relative probabilities for severe outcomes for strains
    p$strain_rel_p_sympt <- process_strain_rel_p(strain_rel_p_sympt,
                                                 strain$n_strains,
                                                 n_real_strains)
    
    p$strain_rel_p_hosp_if_sympt <- process_strain_rel_p(strain_rel_p_hosp_if_sympt,
                                                  strain$n_strains,
                                                  n_real_strains)
    
    p$strain_rel_p_death <- process_strain_rel_p(strain_rel_p_death,
                                                 strain$n_strains,
                                                 n_real_strains)
    
    # Concatenate parameters into one list
    p <- c(p,strain,vaccination)
    
}


process_strain_rel_p <- function(p, n_strains, n_real_strains) {
    if (length(p) < n_strains) {
        p <- recycle(p, n_real_strains)
        if (n_real_strains > 1) {
            p <- mirror_strain(p)
        }
    }
    p
}


vaccination_parameters <- function(N_tot,
                                   dt,
                                   rel_susceptibility = 1,
                                   rel_p_sympt = 1,
                                   rel_p_hosp_if_sympt = 1,
                                   rel_p_death = 1,
                                   rel_infectivity = 1,
                                   vaccine_progression_rate = NULL,
                                   vaccine_schedule = NULL,
                                   vaccine_index_dose2 = NULL,
                                   vaccine_index_booster = NULL,
                                   n_strains = 1,
                                   vaccine_catchup_fraction = 1,
                                   n_doses = 2L) {
    n_groups <- length(N_tot)
    calc_n_vacc_classes <- function(x) {
        if (is.array(x)) nlayer(x) else length(x)
    }
    
    # assert_proportion(rel_susceptibility)
    
    rel_params <- list(rel_susceptibility = rel_susceptibility,
                       rel_p_sympt = rel_p_sympt,
                       rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
                       rel_p_death = rel_p_death,
                       rel_infectivity = rel_infectivity)
    
    n <- vnapply(rel_params, calc_n_vacc_classes)
    
    if (any(n > 1) && length(unique(n[n > 1])) != 1) {
        msg1 <- paste(names(rel_params), collapse = ", ")
        msg2 <- "should have the same dimension"
        stop(paste(msg1, msg2))
    }
    n_vacc_classes <- max(n)
    
    ret <- Map(function(value, name)
        build_rel_param(value, n_groups, n_strains, n_vacc_classes, name),
        rel_params, names(rel_params))
    
    if (is.null(vaccine_schedule)) {
        if (!is.null(vaccine_index_dose2) && vaccine_index_dose2 != 1L) {
            stop("'vaccine_index_dose2' set without schedule")
        }
        ret$vaccine_dose_step <- array(0, c(n_groups, n_doses, 1))
        ret$index_dose <- rep(1L, n_doses)
    } else {
        # assert_is(vaccine_schedule, "vaccine_schedule")
        vaccine_index_dose2 <- vaccine_index_dose2 %||% 1L
        if (vaccine_index_dose2 > n_vacc_classes) {
            stop(sprintf(
                "Invalid value for 'vaccine_index_dose2', must be in [1, %d]",
                n_vacc_classes))
        }
        stopifnot(vaccine_schedule$n_doses == n_doses)
        
        if (is.null(vaccine_index_booster)) {
            if (n_doses != 2L) {
                stop("'n_doses' must be 2 as boosters not used")
            }
        } else {
            if (n_doses != 3L) {
                stop("'n_doses' must be 3 as boosters are used")
            }
            if (vaccine_index_booster > n_vacc_classes) {
                stop(sprintf(
                    "Invalid value for 'vaccine_index_booster', must be in [1, %d]",
                    n_vacc_classes))
            }
        }
        
        n_days <- dim(vaccine_schedule$doses)[[3]]
        i <- rep(seq_len(n_days), each = 1 / dt)
        len <- vaccine_schedule$date / dt
        ret$index_dose <- c(1L, vaccine_index_dose2, vaccine_index_booster)
        
        ret$vaccine_dose_step <- mcstate::array_bind(
            array(0, c(n_groups, n_doses, len)),
            (vaccine_schedule$doses * dt)[, , i])
    }
    
    ret$index_dose_inverse <- create_index_dose_inverse(n_vacc_classes,
                                                        ret$index_dose)
    
    ret$n_vacc_classes <- n_vacc_classes
    ret$vaccine_progression_rate_base <- build_vaccine_progression_rate(
        vaccine_progression_rate, n_groups, n_vacc_classes, ret$index_dose)
    
    
    # assert_scalar(vaccine_catchup_fraction)
    # assert_proportion(vaccine_catchup_fraction)
    ret$vaccine_catchup_fraction <- vaccine_catchup_fraction
    
    ret$n_doses <- n_doses
    
    ret
}


create_index_dose_inverse <- function(n_vacc_classes, index_dose) {
    index_dose_inverse <- integer(n_vacc_classes)
    index_dose_inverse[index_dose] <- seq_along(index_dose)
    index_dose_inverse
}


strain_parameters <- function(strain_transmission, strain_seed_date,
                              strain_seed_size, strain_seed_pattern,
                              dt) {
    if (length(strain_transmission) == 0) {
        stop("At least one value required for 'strain_transmission'")
    }
    if (length(strain_transmission) > 2) {
        stop(paste(
            "Only 1 or 2 strains valid ('strain_transmission' too long)'.",
            "See 'n_S_progress' in the odin code to fix this"))
    }
    
    # assert_non_negative(strain_transmission)
    
    if (is.null(strain_seed_date)) {
        if (!is.null(strain_seed_size)) {
            stop(paste("As 'strain_seed_date' is NULL, expected 'strain_seed_size'",
                       "to be NULL"))
        }
        if (!is.null(strain_seed_pattern)) {
            stop(paste("As 'strain_seed_date' is NULL, expected",
                       "'strain_seed_pattern' to be NULL"))
        }
        strain_seed_step_start <- 0
        strain_seed_value <- 0
    } else {
        if (length(strain_transmission) == 1L) {
            stop("Can't use 'strain_seed_date' if only using one strain")
        }
        if (length(strain_seed_date) != 1L) {
            stop("'strain_seed_date' must be a single date")
        }
        if (length(strain_seed_size) != 1L) {
            stop("'strain_seed_size' must be a single value")
        }
        # assert_covid_multi_strain_date(strain_seed_date)
        # assert_non_negative(strain_seed_size)
        # assert_positive(strain_seed_pattern)
        
        strain_seed_step <- strain_seed_date / dt
        
        strain_seed_value <- strain_seed_size *
            seed_over_steps(strain_seed_step, strain_seed_pattern)
        
        strain_seed_step_start <- floor(strain_seed_step)
    }
    
    if (length(strain_transmission) == 2) {
        strain_transmission <- mirror_strain(strain_transmission)
    }
    
    list(n_strains = length(strain_transmission),
         strain_transmission = strain_transmission,
         strain_seed_step_start = strain_seed_step_start,
         strain_seed_value = strain_seed_value)
}

seed_over_steps <- function(start_step, weights) {
    ## The weights vector must be over steps, not dates
    weights <- weights / sum(weights)
    p <- start_step %% 1
    if (p == 0) {
        ret <- weights
    } else{
        ret <- p * c(0, weights) + (1 - p) * c(weights, 0)
    }
    ret
}


plot_trajectories <- function(time,x,n_age,n_strains,n_vax){
    # Plot trajectories
    for (k in 1) { #:n_vax){
        for (j in 1:n_strains) {
            par(mfrow = c(2,4), oma=c(2,3,0,0))
            for (i in 1:n_age) {
                par(mar = c(3, 4, 2, 0.5))
                cols <- c(S = "#8c8cd9", E = "#ffff00", I_P = "#cc0044", I_A = "green", I_C = "blue", R = "#999966", D = "#000000")
                matplot(time, x[20 + i + (k-1)*n_age, ,-1], type = "l", # Offset to access numbers in age compartment
                        xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
                        col = cols[["S"]], lty = 1, ylim=range(x[-1:-4,,]))
                matlines(time, x[60 + i + (j-1)*n_age + (k-1)*n_age*n_strains, ,-1], col = cols[["E"]], lty = 1)
                matlines(time, x[60 + i + (j-1)*n_age + (k-1)*n_age*n_strains + n_age*n_strains*n_vax, ,-1], col = cols[["I_P"]], lty = 1)
                matlines(time, x[60 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 2*n_age*n_strains*n_vax, ,-1], col = cols[["I_A"]], lty = 1)
                matlines(time, x[60 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 3*n_age*n_strains*n_vax, ,-1], col = cols[["I_C"]], lty = 1)
                matlines(time, x[60 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 4*n_age*n_strains*n_vax, ,-1], col = cols[["R"]], lty = 1)
                matlines(time, x[60 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 7*n_age*n_strains*n_vax, ,-1], col = cols[["D"]], lty = 1)
                legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
                axis(2, las = 2)
            }
            mtext("Number of individuals", side=2, line=1, outer=T, las=0)
            mtext("Time", side = 1, line = 0, outer =T)    
        }    
    }
}


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


# log-likelihood of Poisson count
ll_pois <- function(data, model, exp_noise) {
    if (is.na(data)) {
        # Creates vector of zeros in ll with same length, if no data
        return(numeric(length(model)))
    } 
    lambda <- model + rexp(n = length(model), rate = exp_noise)
    dpois(x = data, lambda = lambda, log = TRUE)
}


# Define negative binomial log-likelihood
ll_nbinom <- function(data, model, kappa, exp_noise){
    if(is.na(data)) {
        return(numeric(length(model)))
    }
    mu <- model + rexp(length(model),rate = exp_noise)
    dnbinom(data, kappa, mu = mu, log = TRUE)
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
    
    # Log-likelihoods for deaths
    ll_hosps <- ll_nbinom(observed$hosps,model_hosps,kappa_hosp,exp_noise)
    ll_hosps_0_39 <- ll_nbinom(observed$hosps_0_39,model_hosps_0_39,kappa_hosp,exp_noise)
    ll_hosps_40_49 <- ll_nbinom(observed$hosps_40_49,model_hosps_40_49,kappa_hosp,exp_noise)
    ll_hosps_50_59 <- ll_nbinom(observed$hosps_50_59,model_hosps_50_59,kappa_hosp,exp_noise)
    ll_hosps_60_69 <- ll_nbinom(observed$hosps_60_69,model_hosps_60_69,kappa_hosp,exp_noise)
    ll_hosps_70_plus <- ll_nbinom(observed$hosps_70_plus,model_hosps_70_plus,kappa_hosp,exp_noise)
    # ll_hosps_70_79 <- ll_nbinom(observed$hosps_70_79,model_hosps_70_79,kappa_hosp,exp_noise)
    # ll_hosps_80_plus <- ll_nbinom(observed$hosps_80_plus,model_hosps_80_plus,kappa_hosp,exp_noise)
    
    # Log-likelihoods for deaths
    ll_deaths <- ll_nbinom(observed$deaths,model_deaths,kappa_death,exp_noise)
    ll_deaths_0_39 <- ll_nbinom(observed$deaths_0_39,model_deaths_0_39,kappa_death,exp_noise)
    ll_deaths_40_49 <- ll_nbinom(observed$deaths_40_49,model_deaths_40_49,kappa_death,exp_noise)
    ll_deaths_50_59 <- ll_nbinom(observed$deaths_50_59,model_deaths_50_59,kappa_death,exp_noise)
    ll_deaths_60_69 <- ll_nbinom(observed$deaths_60_69,model_deaths_60_69,kappa_death,exp_noise)
    ll_deaths_70_plus <- ll_nbinom(observed$deaths_70_plus,model_deaths_70_plus,kappa_death,exp_noise)
    # ll_deaths_70_79 <- ll_nbinom(observed$deaths_70_79,model_deaths_70_79,kappa_death,exp_noise)
    # ll_deaths_80_plus <- ll_nbinom(observed$deaths_80_plus,model_deaths_80_plus,kappa_death,exp_noise)
    
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
    # ll_deaths + ll_deaths_0_39 + ll_deaths_40_49 + ll_deaths_50_59 + ll_deaths_60_69 + ll_deaths_70_plus #+ ll_deaths_70_79 + ll_deaths_80_plus
    # ll_hosps_70_plus + ll_deaths_70_plus
    ll_deaths_70_plus
}


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
        matplot(times, t(history[i + n_age*n_vax + n_age, ,-1]), type = "l", # Offset to access numbers in age compartment
                xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
                col = alpha(cols[["E"]],0.1), lty = 1, ylim=range(true_history[60 + i, ,-1]))
        # # matlines(times, t(history[i + n_age*n_vax, ,-1]), col = alpha(cols[["E"]],0.01), lty = 1)
        # matlines(times, t(history[i + 2*n_age*n_vax, ,-1]), col = alpha(cols[["I_P"]],0.01), lty = 1)
        # matlines(times, t(history[i + 3*n_age*n_vax, ,-1]), col = alpha(cols[["I_A"]],0.01), lty = 1)
        # matlines(times, t(history[i + 4*n_age*n_vax, ,-1]), col = alpha(cols[["I_C"]],0.01), lty = 1)
        # # matlines(times, t(history[i + 5*n_age*n_vax, ,-1]), col = alpha(cols[["R"]],0.01), lty = 1)
        # # matlines(times, t(history[i + 8*n_age*n_vax, ,-1]), col = alpha(cols[["D"]],0.01), lty = 1)
        # # matpoints(times[1:obs_end], t(true_history[seq(i+12,i+12+8*n_age*n_vax,by = n_age*n_vax), ,-1]), pch = 19,
        # #           col = cols)
        # matpoints(times[1:obs_end], t(true_history[i + 12 + seq(n_age*n_vax,4*n_age*n_vax,by = n_age*n_vax), ,-1]), pch = 19,
        #           col = cols[2:5])
        points(times[1:obs_end], true_history[60 + n_age + i, ,-1], pch = 19,
               col = cols[2])
        # legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
        legend("right", lwd = 1, col = cols[2:5], legend = names(cols[2:5]), bty = "n")
        axis(2, las = 2)
    }
}


make_transform <- function(dt,
                           n_age,
                           n_vax,
                           m,
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
                           cross_immunity){
    
    function(pars){
        beta <- pars[["beta"]]
        # ngm <- t(t(transmission)*(p_C*(1/pars[["gamma"]]+1/pars[["gamma"]])+(1-p_C)*theta_A/gamma_A)*population)
        # beta <- pars[["R0"]]/eigen(ngm)$value[1]
        gamma <- pars[["gamma"]]
        start_date <- pars[["start_date"]]
        strain_seed_date <- pars[["strain_seed_date"]]
        
        p <- parameters(dt,
                        n_age,
                        n_vax,
                        m,
                        beta = beta,
                        gamma_E,
                        gamma_P = gamma,
                        gamma_A,
                        gamma_C = gamma,
                        gamma_H,
                        gamma_G,
                        theta_A,
                        p_C,
                        p_H,
                        p_G,
                        p_D,
                        population,
                        start_date = start_date,
                        initial_seed_size,
                        initial_seed_pattern,
                        strain_transmission,
                        strain_seed_date = strain_seed_date,
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
        p
    }
    
}

make_transform_multistage <- function(dt,
                                      n_age,
                                      n_vax,
                                      m,
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
                                      start_date1,
                                      strain_transmission1,
                                      cross_immunity1){
    
    function(pars){
        beta <- pars[["beta"]]
        # ngm <- t(t(transmission)*(p_C*(1/pars[["gamma"]]+1/pars[["gamma"]])+(1-p_C)*theta_A/gamma_A)*population)
        # beta <- pars[["R0"]]/eigen(ngm)$value[1]
        gamma <- pars[["gamma"]]
        start_date <- pars[["start_date"]]
        strain_seed_date <- pars[["strain_seed_date"]]
        strain_seed_date1 <- pars[["strain_seed_date1"]]
        
        # Parameters for 1st epoch
        p <- parameters(dt,
                        n_age,
                        n_vax,
                        m,
                        beta = beta,
                        gamma_E,
                        gamma_P = gamma,
                        gamma_A,
                        gamma_C = gamma,
                        gamma_H,
                        gamma_G,
                        theta_A,
                        p_C,
                        p_H,
                        p_G,
                        p_D,
                        population,
                        start_date = start_date,
                        initial_seed_size,
                        initial_seed_pattern,
                        strain_transmission,
                        strain_seed_date = strain_seed_date,
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
        
        # Parameters for 2nd epoch
        p1 <- parameters(dt,
                         n_age,
                         n_vax,
                         transmission,
                         beta = beta,
                         gamma_E,
                         gamma_P = gamma,
                         gamma_A,
                         gamma_C = gamma,
                         gamma_H,
                         gamma_G,
                         theta_A,
                         p_C,
                         p_H,
                         p_G,
                         p_D,
                         population,
                         start_date = start_date1,
                         initial_seed_size = 0,
                         initial_seed_pattern,
                         strain_transmission = strain_transmission1,
                         strain_seed_date = strain_seed_date1,
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
        
        epochs <- list(
            multistage_epoch(start_date1, pars = p1, transform_state = transform_state)
        )
        p_multistage <- multistage_parameters(p, epochs)
        p_multistage
    }
    
}


plot_hosps_and_deaths_age <- function(incidence_modelled, incidence_observed, times, n_age, n_vax, n_strains){
    par(mfrow = c(2,4), oma=c(2,3,0,0))
    idx1 <- n_age*n_vax + 8*n_age*n_strains*n_vax + 1
    for (i in 1:5){
        par(mar = c(3, 4, 2, 0.5))
        matplot(times, t(incidence_modelled[idx1+i, ,-1]),
                type="l",col = alpha("black",0.1),xlab = "Day",ylab = "Hospitalisations",
                main = paste0("Age ", rownames(incidence_modelled)[idx1+i]))
        points(times, incidence_observed[[4+i]],pch=19,col="red")
        axis(2, las = 2)
    }
    par(mfrow = c(2,4), oma=c(2,3,0,0))
    for (i in 1:5){
        par(mar = c(3, 4, 2, 0.5))
        matplot(times, t(incidence_modelled[idx1+6+i, ,-1]),
                type="l",col = alpha("black",0.1),xlab = "Day",ylab = "Deaths",
                main = paste0("Age ", rownames(incidence_modelled)[idx1+6+i]))
        points(times, incidence_observed[[9+i]],pch=19,col="red")
        axis(2, las = 2)
    }
}

##' Rotate strains, so that strain 1 becomes the sum of strains 1 and
##' 2 and strain 2 is empty. Use this to allow sequential replacement
##' of strains.
##'
##' @title Rotate strains
##'
##' @param state Model state
##'
##' @param info Model info
##'
##' @param ... Additional arguments, ignored. This exists so that this
##'   function can be used as an argument to
##'   [mcstate::multistage_epoch].  Practically your two model
##'   informations in this case would be equivalent.
##'
##' @export
rotate_strains <- function(state, info, ...) {
    if (is.null(dim(state))) {
        stop("Expected a matrix or array for 'state'")
    }
    if (nrow(state) != info$len) {
        stop(sprintf("Expected a matrix with %d rows for 'state'",
                     info$len))
    }
    
    ## Push all state down into a common rank to avoid lots of boring
    ## dimension arithmetic (e.g., we might have this structured by
    ## region, or a single particle's state etc, but what matters is the
    ## structure by particle only).
    dim_orig <- dim(state)
    if (length(dim(state)) > 2) {
        state <- matrix(state, nrow(state))
    }
    
    ## Currently we only support one sort of move: Anyone who has been
    ## infected in the past is moved to now be indexed by strain 1 (no
    ## matter whether they had multiple infections before), and strain 2
    ## will be completely empty.
    strain_from_idx <- c(2, 3, 4)
    strain_to_idx <- 1
    
    n_particle <- ncol(state)
    for (i in seq_along(rotate_strain_compartments)) {
        name <- rotate_strain_compartments[[i]]
        dim <- info$dim[[name]]
        state_i <- array(state[info$index[[name]], ], c(dim, n_particle))
        
        if (length(dim) == 4) {
            for (j in strain_from_idx) {
                tomove <- state_i[, j, , , , drop = FALSE]
                state_i[, strain_to_idx, , , ] <-
                    state_i[, strain_to_idx, , , , drop = FALSE] + tomove
                state_i[, j, , , ] <- 0
            }
        } else if (length(dim) == 3) {
            for (j in strain_from_idx) {
                tomove <- state_i[, j, , , drop = FALSE]
                state_i[, strain_to_idx, , ] <-
                    state_i[, strain_to_idx, , , drop = FALSE] + tomove
                state_i[, j, , ] <- 0
            }
        } else if (length(dim) == 1) {
            ## this loop range ensures that the move can still happen when
            ## the object has dimension n_real_strains not n_strains
            ## e.g. for prob_strain
            stopifnot(dim %in% c(2, 4))
            for (j in strain_from_idx[strain_from_idx <= dim]) {
                tomove <- state_i[j, , drop = FALSE]
                state_i[strain_to_idx, ] <-
                    state_i[strain_to_idx, , drop = FALSE] + tomove
                state_i[j, ] <- 0
            }
        } else {
            ## This is unreachable unless the model changes to include
            ## something that has a rank-2 variable that needs transforming.
            stop(sprintf("Unexpected dimensions (%d) in rotate_strain", # nocov
                         length(dim)))                                  # nocov
        }
        state[info$index[[name]], ] <- state_i
    }
    
    dim(state) <- dim_orig
    state
}


rotate_strain_compartments <- c(
    ## those with dimension c(n_groups, n_strains, n_vacc_classes):
    "E", "I_A", "I_P", "I_C", "R", "G", "H", "D")

transform_state <- function(state, info_old, info_new){
    rotate_strains(state,info_old)
}
