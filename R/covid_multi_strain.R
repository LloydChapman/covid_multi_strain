transmission_matrix <- function(country, pop, age_groups){
    # Get contact matrix of country
    contact <- contact_matrix(country)
    
    # Convert contact matrix for country to one with population age structure pop
    contact1 <- contact_matrix_density_correction(contact,country,pop)
    
    # Aggregate contacts over age groups
    contact2 <- contact_matrix_aggregate_age_groups(contact1,age_groups,pop)
    
    # Convert the contact matrix to the "transmission matrix" (the contact matrix weighted by the population in each age group)
    contact_matrix <- matrix(contact2[,contacts],nrow = length(age_groups),ncol = length(age_groups), byrow = T) # N.B. individuals in rows, contacts in columns
    population <- pop[,.(sum(total)),by = .(age_group)][,V1]
    transmission <- contact_matrix/rep(population, each = ncol(contact_matrix))
    
    transmission
}


contact_matrix <- function(country){
    # Load contact matrices
    contact_matrices <- fread("data/synthetic_contacts_2020.csv")
    setnames(contact_matrices,"age_cotactee","age_contactee")
    
    # Get contact matrix for desired country
    contact <- contact_matrices[iso3c == country & setting == "overall" & location_contact == "all"]
    age_cols <- names(contact)[grep("age",names(contact))]
    contact[,(paste0("min_",age_cols)):=lapply(.SD,function(x) as.numeric(sub("\\+","",sub(" to.*","",x)))),.SDcols = age_cols] 
    # Convert age group columns to factors to maintain ordering 
    contact[,`:=`(age_contactor = factor(age_contactor,levels = unique(age_contactor)),
                  age_contactee = factor(age_contactee, levels = unique(age_contactee)))]
    
    contact
}


contact_matrix_density_correction <- function(contact,cntry,pop){
    min_ages_contact <- contact[,unique(min_age_contactee)]
    age_groups_contact <- contact[,unique(age_contactee)]
    
    # Get population age structure of country with contact matrix
    pop_all <- qread("data/unwpp_data.qs")
    
    pop_contact <- pop_all[iso3 == cntry & year == 2020,]
    pop_contact[,age_group_contact := cut(age,c(min_ages_contact,Inf),labels = age_groups_contact,right = F)]
    pop_contact <- pop_contact[,.(population = sum(total)),by = .(age_group_contact)]
    pop_contact <- pop_contact[,pop_prop := population/sum(population)]
    
    # Get population age structure of country without contact matrix (by contact age groups)
    pop1 <- copy(pop)
    pop1[,age_group_contact := cut(age,c(min_ages_contact,Inf),labels = age_groups_contact,right = F)]
    pop1 <- pop1[,.(population = sum(total)),by = .(age_group_contact)]
    pop1 <- pop1[,pop_prop := population/sum(population)]
    
    # Merge populations of country with contact matrix and country without with contact data
    contact1 <- merge(contact,pop_contact,by.x = "age_contactee",by.y = "age_group_contact")
    contact1 <- merge(contact1,pop1,by.x = "age_contactee",by.y = "age_group_contact")
    
    # Correct mean number of contacts by relative population densities
    contact1[,mean_number_of_contacts := mean_number_of_contacts*pop_prop.y/pop_prop.x]
    
    contact1
}


contact_matrix_aggregate_age_groups <- function(contact,age_groups,pop){
    # Sum mean numbers of contacts over contact age groups being aggregated
    min_ages <- get_min_age(age_groups)
    contact[,age_contactee := cut(min_age_contactee,c(min_ages,Inf),labels = age_groups,right = F)]
    contact1 <- contact[,.(contacts = sum(mean_number_of_contacts)),by = .(age_contactor,age_contactee,min_age_contactor)]
    
    min_ages_contact <- contact[,unique(min_age_contactor)]
    age_groups_contact <- contact[,unique(age_contactor)]
    
    # Get population age structure of country without contact matrix (by contactor age groups)
    pop1 <- copy(pop)
    pop1[,age_group_contact := cut(age,c(min_ages_contact,Inf),labels = age_groups_contact,right = F)]
    pop1 <- pop1[,.(population = sum(total)),by = .(age_group_contact)]
    pop1 <- pop1[,pop_prop := population/sum(population)]
    
    contact1 <- merge(contact1,pop1,by.x = "age_contactor",by.y = "age_group_contact")
    contact1[,age_contactor := cut(min_age_contactor,c(min_ages,Inf),labels = age_groups,right = F)]
    
    # Take population-weighted average of mean number of contacts over contactor age groups being aggregated
    contact2 <- contact1[,.(contacts = sum(contacts * population)/sum(population)), by = .(age_contactor,age_contactee)]
    
    # # Plot
    # ggplot(contact2,aes(x = age_contactee,y = age_contactor,fill = contacts)) + geom_tile()
    
    contact2
}


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
                       beta_date,
                       beta_value,
                       beta_type,
                       gamma_E,
                       gamma_P,
                       gamma_A,
                       gamma_C,
                       gamma_H,
                       gamma_G,
                       gamma_pre_1,
                       gamma_P_1,
                       theta_A,
                       severity,
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
                       strain_rel_gamma_E = 1,
                       strain_rel_gamma_P = 1,
                       strain_rel_gamma_C = 1,
                       strain_rel_gamma_A = 1,
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
                       vacc_skip_progression_rate = NULL,
                       vacc_skip_to = NULL,
                       vacc_skip_weight = NULL,
                       waning_rate = 0,
                       cross_immunity = 1,
                       observation,
                       sero_sensitivity_1 = 0.9,
                       sero_specificity_1 = 0.99,
                       test_sensitivity = 0.99,
                       test_specificity = 0.99) {
    
    n_real_strains <- length(strain_transmission)
    
    if (n_real_strains > 2) {
        stop("Only 1 or 2 strains valid ('strain_transmission' too long)'.")
    }
    
    # Construct time-varying transmission coefficient
    if (beta_type == "piecewise-linear") {
        beta_step <- parameters_piecewise_linear(beta_date, 
                                                 beta_value %||% 0.1, dt)
    } else if (beta_type == "piecewise-constant") {
        beta_step <- parameters_piecewise_constant(beta_date, 
                                                   beta_value %||% 0.1, dt)
    } else {
        stop("'beta_type' must be 'piecewise-linear' or 'piecewise-constant'")
    }
    
    # Transmission and progression parameters
    p <- list(dt = dt,
              n_age = n_age,
              n_vax = n_vax,
              m = m,
              beta_step = beta_step,
              gamma_E = gamma_E,
              gamma_P = gamma_P,
              gamma_A = gamma_A,
              gamma_C = gamma_C,
              gamma_H = gamma_H, 
              gamma_G = gamma_G,
              gamma_pre_1 = gamma_pre_1,
              gamma_P_1 = gamma_P_1,
              theta_A = theta_A)
    
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
    strain <- parameters_strain(strain_transmission, strain_seed_date, 
                                strain_seed_size, strain_seed_pattern,p$dt)
    
    # Construct vaccination parameters
    vaccination <- parameters_vaccination(p$N_tot,
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
    
    # Vaccine stratum skip parameters
    vacc_skip_to <- vacc_skip_to %||% integer(vaccination$n_vacc_classes)
    vacc_skip_progression_rate <-
        vacc_skip_progression_rate %||% rep(0, vaccination$n_vacc_classes)
    vacc_skip_weight <- vacc_skip_weight %||% rep(0, vaccination$n_vacc_classes)
    vacc_skip <- parameters_vacc_skip(vacc_skip_to,
                                      vacc_skip_progression_rate,
                                      vacc_skip_weight,
                                      vaccination$n_vacc_classes,
                                      n_doses,
                                      vaccination$index_dose_inverse)
    
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
    
    # Strain-dependent progression rates
    strain_rel_gammas <- list(E = strain_rel_gamma_E,
                              P = strain_rel_gamma_P,
                              C = strain_rel_gamma_C,
                              A = strain_rel_gamma_A)
    for (name in names(strain_rel_gammas)) {
        rel_gamma <- strain_rel_gammas[[name]]
        rel_gamma_name <- paste0("rel_gamma_", name)
        if (is.null(rel_gamma)) {
            p[[rel_gamma_name]] <- rep(1, strain$n_strains)
        } else {
            # rel_gamma <- recycle(assert_non_negative(rel_gamma),
            #                      n_real_strains)
            rel_gamma <- recycle(rel_gamma,n_real_strains)
            if (length(rel_gamma) == 2) {
                p[[rel_gamma_name]] <- mirror_strain(rel_gamma)
            } else {
                p[[rel_gamma_name]] <- rep(1, strain$n_strains)
            }
        }
    }
    
    # # Observation parameters
    # p$phi_cases <- phi_cases
    # p$kappa_cases <- kappa_cases
    # p$kappa_hosp <- kappa_hosp
    # p$kappa_death <- kappa_death
    
    # Sensitivity and specificity of serological tests
    p$sero_sensitivity_1 <- sero_sensitivity_1
    p$sero_specificity_1 <- sero_specificity_1
    
    # Sensitivity and specificity of PCR and rapid Ag tests
    p$test_sensitivity <- test_sensitivity
    p$test_specificity <- test_specificity
    
    # Concatenate parameters into one list
    p <- c(p,severity,strain,vaccination,vacc_skip,observation)
    
    p
}


parameters_piecewise_linear <- function(date, value, dt) {
    ## for one strain coerce numeric to matrix, will coerce back later
    if (!inherits(value, "matrix")) {
        value <- matrix(value, ncol = 1)
    }
    if (is.null(date)) {
        if (nrow(value) != 1L) {
            stop("As 'date' is NULL, expected single value")
        }
        ## coerce back
        if (ncol(value) == 1) {
            value <- as.numeric(value)
        }
        return(value)
    }
    if (length(date) != nrow(value)) {
        stop("'date' and 'value' must have the same length")
    }
    if (length(date) < 2) {
        stop("Need at least two dates and values for a varying piecewise linear")
    }
    # assert_covid_multi_strain_date(date)
    # assert_increasing(date)
    
    if (date[[1]] != 0) {
        date <- c(0, date)
        value <- rbind(value[1, ], value)
    }
    
    value <- apply(value, 2, function(x) {
        stats::approx(date, x, seq(0, date[[length(date)]], by = dt))$y
    })
    
    ## coerce back
    if (ncol(value) == 1) {
        value <- as.numeric(value)
    }
    
    value
}


parameters_piecewise_constant <- function(date, value, dt) {
    if (is.null(date)) {
        if (length(value) != 1L) {
            stop("As 'date' is NULL, expected single value")
        }
        return(value)
    }
    if (length(date) != length(value)) {
        stop("'date' and 'value' must have the same length")
    }
    # assert_covid_multi_strain_date(date)
    # assert_increasing(date)
    if (!is.null(date)) {
        if (date[1L] != 0) {
            stop("As 'date' is not NULL, first date should be 0")
        }
    }
    
    stats::approx(date, value, method = "constant",
                  xout = seq(0, date[[length(date)]], by = dt))$y
    
}


##' Expand `value_step` based on a series of `step`s.  Use this to
##' convert between the values passed to
##' [parameters_piecewise_linear()] and the actual values
##' for a given set of steps.
##'
##' @title Expand beta steps
##'
##' @param step A vector of steps
##'
##' @param value_step A vector of values
##'
##' @return A numeric vector the same length as `step`
##'
##' @export
parameters_expand_step <- function(step, value_step) {
    value_step[pmin(step, length(value_step) - 1L) + 1L]
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


parameters_severity <- function(dt,
                                severity,
                                p_D) {
    
    time_varying_severity <- list(D = p_D)
    
    get_p_step <- function(x, name) {
        
        p_name <- paste0("p_", name)
        p <- x[[p_name]]
        time_vary <- time_varying_severity[[name]]
        
        if (is.null(time_vary)) {
            p_value <- NULL
            p_date <- NULL
        } else {
            p_value <- time_vary$value
            if ("date" %in% names(time_vary)) {
                p_date <- time_vary$date
            } else {
                p_date <- NULL
            }
        }
        
        if (!is.null(p_value)) {
            # assert_proportion(p_value, p_name)
            if (length(p_value) == 1L) {
                if (length(p_date) != 0) {
                    stop(sprintf(
                        "As '%s' has a single 'value', expected NULL or missing 'date'",
                        p_name))
                }
            } else if (length(p_date) != length(p_value)) {
                stop(sprintf("'date' and 'value' for '%s' must have the same length",
                             p_name))
            }
        }
        
        if (all(p == 0)) {
            psi <- p
        } else {
            psi <- p / max(p)
        }
        
        if (is.null(p_value)) {
            p_step <- max(p)
        } else {
            # p_step <- parameters_piecewise_constant(p_date, p_value, dt)
            p_step <- parameters_piecewise_linear(p_date, p_value, dt)
        }
        
        p_step <- outer(p_step, psi)
        
        x[[paste0(p_name, "_step")]] <- p_step
        x[[paste0(p_name)]] <- NULL
        x[[paste0("n_", p_name, "_steps")]] <- dim(p_step)[1]
        
        x
    }
    
    for (name in names(time_varying_severity)) {
        severity <- get_p_step(severity, name)
    }
    
    severity
}


parameters_vaccination <- function(N_tot,
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


parameters_vacc_skip <- function(vacc_skip_to,
                                 vacc_skip_progression_rate,
                                 vacc_skip_weight,
                                 n_vacc_classes,
                                 n_doses,
                                 index_dose_inverse) {
    
    vacc_classes <- seq_len(n_vacc_classes)
    if (length(vacc_skip_to) != n_vacc_classes) {
        stop(sprintf("There are %s vaccine classes so 'vacc_skip_to' must be of
                 length %s", n_vacc_classes, n_vacc_classes))
    }
    if (length(vacc_skip_weight) != n_vacc_classes) {
        stop(sprintf("There are %s vaccine classes so 'vacc_skip_weight' must be of
                 length %s", n_vacc_classes, n_vacc_classes))
    }
    if (length(vacc_skip_progression_rate) != n_vacc_classes) {
        stop(sprintf("There are %s vaccine classes so 'vacc_skip_progression_rate'
                 must be of length %s", n_vacc_classes, n_vacc_classes))
    }
    if (!all(vacc_skip_to %in% vacc_classes | vacc_skip_to == 0)) {
        stop(sprintf("There are %s vaccine classes so the values in 'vacc_skip_to'
                 must be 0 or one of: %s",
                     n_vacc_classes,
                     paste(vacc_classes, collapse = ", ")))
    }
    if (any(vacc_skip_to == 0 & vacc_skip_weight != 0)) {
        stop("Require 0 values in 'vacc_skip_weight' for vaccine strata that have 0
         values in 'vacc_skip_to'")
    }
    if (any(vacc_skip_to == 0 & vacc_skip_progression_rate != 0)) {
        stop("Require 0 values in 'vacc_skip_progression_rate' for vaccine strata
         that have 0 values in 'vacc_skip_to'")
    }
    if (any(vacc_skip_to != 0 & vacc_skip_to <= vacc_classes + 1)) {
        stop("Require vacc_skip_to[j] = 0 or vacc_skip_to[j] > j + 1")
    }
    vacc_skip_moves <- which(vacc_skip_to > 0)
    vacc_skip_to_moves <- vacc_skip_to[vacc_skip_moves]
    if (length(unique(vacc_skip_to_moves)) != length(vacc_skip_to_moves)) {
        stop("Cannot have more than one vaccine skip move to the same stratum")
    }
    
    # assert_proportion(vacc_skip_weight)
    
    
    vacc_skip_from <- integer(n_vacc_classes)
    vacc_skip_from[vacc_skip_to[vacc_skip_moves]] <- vacc_skip_moves
    vacc_skip_dose_inverse <- integer(n_vacc_classes)
    vacc_skip_dose_inverse[vacc_skip_moves] <-
        index_dose_inverse[vacc_skip_to[vacc_skip_moves] - 1]
    vacc_skip_dose <- integer(n_doses)
    vacc_skip_dose[vacc_skip_dose_inverse[vacc_skip_moves]] <- vacc_skip_moves
    vacc_skipped <- integer(n_vacc_classes)
    if (length(vacc_skip_moves) > 0) {
        for (i in vacc_classes) {
            skipped_from <- which(vacc_classes <= i & vacc_skip_to > i)
            if (length(skipped_from) > 1) {
                stop("Cannot have overlapping vaccine skip moves")
            } else if (length(skipped_from) == 1) {
                vacc_skipped[i] <- skipped_from
            }
        }
    }
    vacc_skip_dose_weight <- rep(0, n_doses)
    vacc_skip_dose_weight[vacc_skip_dose_inverse[vacc_skip_moves]] <-
        vacc_skip_weight[vacc_skip_moves]
    
    list(vacc_skip_to = vacc_skip_to,
         vacc_skip_from = vacc_skip_from,
         vacc_skip_weight = vacc_skip_weight,
         vacc_skip_dose_weight = vacc_skip_dose_weight,
         vacc_skip_dose = vacc_skip_dose,
         vacc_skip_dose_inverse = vacc_skip_dose_inverse,
         vacc_skipped = vacc_skipped,
         vacc_skip_progression_rate_base = vacc_skip_progression_rate)
}


parameters_strain <- function(strain_transmission, strain_seed_date,
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


index <- function(info, min_ages = seq(0,70,by = 10), Rt = TRUE){
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
                   deaths_70_plus = index[["D_inc_70_plus"]],
                   # deaths_70_79 = index[["D_inc_70_79"]],
                   # deaths_80_plus = index[["D_inc_80_plus"]]
                   cases = index[["cases_inc"]],
                   # cases_week = index[["cases_inc_week"]],
                   cases_non_variant = index[["cases_non_variant_inc"]],
                   cases_0_9 = index[["cases_inc_0_9"]],
                   cases_10_19 = index[["cases_inc_10_19"]],
                   cases_20_29 = index[["cases_inc_20_29"]],
                   cases_30_39 = index[["cases_inc_30_39"]],
                   cases_40_49 = index[["cases_inc_40_49"]],
                   cases_50_59 = index[["cases_inc_50_59"]],
                   cases_60_69 = index[["cases_inc_60_69"]],
                   cases_70_plus = index[["cases_inc_70_plus"]],
                   sero_pos_1 = index[["sero_pos_1"]],
                   sero_pos_1_20_29 = index[["sero_pos_1_20_29"]],
                   sero_pos_1_30_39 = index[["sero_pos_1_30_39"]],
                   sero_pos_1_40_49 = index[["sero_pos_1_40_49"]],
                   sero_pos_1_50_59 = index[["sero_pos_1_50_59"]],
                   sero_pos_1_60_69 = index[["sero_pos_1_60_69"]],
                   sero_pos_1_70_plus = index[["sero_pos_1_70_plus"]]
    )
    
    suffix <- paste0("_",min_ages)
    
    n_vacc_classes <- info$dim$S[[2]]
    n_strains <- info$dim$prob_strain
    if (n_strains == 2){
        n_tot_strains <- 4
    } else {
        n_tot_strains <- 1
    }
    
    # S states (age x vacc class)
    index_S <- calculate_index(index, "S", list(n_vacc_classes), suffix)
    
    # Strain weights for Rt calculation - relative probability of infection with each strain
    index_prob_strain <- calculate_index(index, "prob_strain", list(n_strains))
    
    # Effective susceptibles (age x strain)
    index_effective_susceptible <- 
        calculate_index(index, "effective_susceptible", list(n_strains),
                        suffix, "effective_susceptible")
    
    # # New infections (age x strain x vacc class)
    # index_new_infections <- 
    #     calculate_index(index, "new_infections", 
    #                     list(S = n_strains, V = n_vacc_classes), suffix)
    # 
    # # New reinfections (age x strain x vacc class)
    # index_new_reinfections <- 
    #     calculate_index(index, "new_reinfections", 
    #                     list(S = n_strains, V = n_vacc_classes), suffix)
    
    # R states (age x (total) strain x vacc class)
    index_R <- calculate_index(index, "R", list(S = n_tot_strains, 
                                                V = n_vacc_classes), suffix)
    # index_N_tot <- index[["N_tot"]]
    
    # index_state <- c(index_run, index_effective_susceptible, index_N_tot)
    index_state <- c(index_run, index_effective_susceptible)
    # index_state <- c(index_run, index_effective_susceptible,
    #                  index_new_infections, index_new_reinfections)
    
    if (Rt){
        index_state <- c(index_state, index_S, index_R, index_prob_strain)
    }
    
    # list(run = index_run,
    #      state = c(S = index[["S"]],
    #                E = index[["E"]],
    #                I_P = index[["I_P"]],
    #                I_A = index[["I_A"]],
    #                I_C = index[["I_C"]],
    #                R = index[["R"]],
    #                H = index[["H"]],
    #                G = index[["G"]],
    #                D = index[["D"]],
    #                T_pre_1 = index[["T_pre_1"]],
    #                T_P_1 = index[["T_P_1"]],
    #                T_N_1 = index[["T_N_1"]],
    #                index_run))
    list(run = index_run,
         state = index_state)
}


# log-likelihood of binomial count
ll_binom <- function(data_x, data_size, model_prob) {
    if (is.na(data_x) || is.na(data_size)) {
        return(numeric(length(model_prob)))
    }
    dbinom(data_x, data_size, model_prob, log = TRUE)
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


# log-likelihood of binomial count where success probability is beta-distributed
ll_betabinom <- function(data_x, data_size, model_prob, rho) {
    if (is.na(data_x) || is.na(data_size)) {
        return(numeric(length(model_prob)))
    }
    dbetabinom(data_x, data_size, model_prob, rho, log = TRUE)
}


dbetabinom <- function(x, size, prob, rho, log = FALSE) {
    
    a <- prob * (1 / rho - 1)
    b <- (1 - prob) * (1 / rho - 1)
    
    out <- lchoose(size, x) + lbeta(x + a, size - x + b) - lbeta(a, b)
    
    if (!log) {
        out <- exp(out)
    }
    
    out
}


ll_dirmnom <- function(data, model, size, exp_noise){
    if(any(is.na(data))){
        return(numeric(length(model)))
    }
    model <- model + rexp(n = length(model), rate = exp_noise)
    alpha <- size * model/rowSums(model)
    ddirmnom(data, rowSums(data), alpha, log = TRUE)
}


test_prob_pos <- function(pos, neg, sensitivity, specificity, exp_noise) {
    
    ## We add some exponential noise to the number of positives and negatives
    ## to help ensure prob_pos is not 0 or 1. If e.g. prob_pos were 0 and there
    ## were individuals who tested positive, this would result in a weight of 0
    ## for a particle. If all particles have weights of 0, the particle filter
    ## breaks. The exponential noise produces small non-zero weights in these
    ## circumstances to prevent the particle filter from breaking.
    
    pos <- pos + rexp(length(pos), exp_noise)
    neg <- neg + rexp(length(neg), exp_noise)
    
    prob_pos <- (sensitivity * pos + (1 - specificity) * neg) / (pos + neg)
    prob_pos
}


# Define comparison function for age-stratified data
compare <- function(state, observed, pars){
    if (is.null(pars$phi_cases)){
        phi_cases <- 0.5
    } else {
        phi_cases <- pars$phi_cases
    }
    if (is.null(pars$kappa_cases)){
        kappa_cases <- 2
    } else {
        kappa_cases <- pars$kappa_cases
    }
    if (is.null(pars$p_NC)){
        p_NC <- 0.002
    } else {
        p_NC <- pars$p_NC
    }
    if (is.null(pars$rho_tests)){
        rho_tests <- 0.01
    } else {
        rho_tests <- pars$rho_tests
    }
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
    
    # Modelled hospitalisations
    model_hosps <- state["hosps", ] # for locations without age-stratified data
    model_hosps_0_39 <- state["hosps_0_39", ]
    model_hosps_40_49 <- state["hosps_40_49", ]
    model_hosps_50_59 <- state["hosps_50_59", ]
    model_hosps_60_69 <- state["hosps_60_69", ]
    model_hosps_70_plus <- state["hosps_70_plus", ]
    # model_hosps_70_79 <- state["hosps_70_79", ]
    # model_hosps_80_plus <- state["hosps_80_plus", ]
    
    # Modelled deaths
    model_deaths <- state["deaths", ] # for locations without age-stratified data
    model_deaths_0_39 <- state["deaths_0_39", ]
    model_deaths_40_49 <- state["deaths_40_49", ]
    model_deaths_50_59 <- state["deaths_50_59", ]
    model_deaths_60_69 <- state["deaths_60_69", ]
    model_deaths_70_plus <- state["deaths_70_plus", ]
    # model_deaths_70_79 <- state["deaths_70_79", ]
    # model_deaths_80_plus <- state["deaths_80_plus", ]
    
    # Modelled cases
    model_cases_0_9 <- state["cases_0_9", ]
    model_cases_10_19 <- state["cases_10_19", ]
    model_cases_20_29 <- state["cases_20_29", ]
    model_cases_30_39 <- state["cases_30_39", ]
    model_cases_40_49 <- state["cases_40_49", ]
    model_cases_50_59 <- state["cases_50_59", ]
    model_cases_60_69 <- state["cases_60_69", ]
    model_cases_70_plus <- state["cases_70_plus", ]
    
    model_confirmed_cases_0_9 <- phi_cases * model_cases_0_9
    model_confirmed_cases_10_19 <- phi_cases * model_cases_10_19
    model_confirmed_cases_20_29 <- phi_cases * model_cases_20_29
    model_confirmed_cases_30_39 <- phi_cases * model_cases_30_39
    model_confirmed_cases_40_49 <- phi_cases * model_cases_40_49
    model_confirmed_cases_50_59 <- phi_cases * model_cases_50_59
    model_confirmed_cases_60_69 <- phi_cases * model_cases_60_69
    model_confirmed_cases_70_plus <- phi_cases * model_cases_70_plus
    
    # Modelled seropositives
    # model_sero_pos_1 can go above N_tot (I think due to reinfection) so cap it to avoid probabilities > 1
    model_sero_pos_1 <- state["sero_pos_1", ]
    model_sero_pos_1_capped <- pmin(model_sero_pos_1, sum(pars$N_tot))
    model_sero_prob_pos_1 <- test_prob_pos(model_sero_pos_1_capped,
                                           sum(pars$N_tot) - model_sero_pos_1_capped,
                                           pars$sero_sensitivity_1,
                                           pars$sero_specificity_1,
                                           exp_noise)
    # model_sero_pos_1_X can go above N_tot_X (I think due to reinfection) so cap to avoid probabilities > 1
    model_sero_pos_1_20_29 <- state["sero_pos_1_20_29", ]
    model_sero_pos_1_20_29_capped <- pmin(model_sero_pos_1_20_29, pars$N_tot[3])
    model_sero_prob_pos_1_20_29 <- test_prob_pos(model_sero_pos_1_20_29_capped,
                                                 pars$N_tot[3] - model_sero_pos_1_20_29_capped,
                                                 pars$sero_sensitivity_1,
                                                 pars$sero_specificity_1,
                                                 exp_noise)
    model_sero_pos_1_30_39 <- state["sero_pos_1_30_39", ]
    model_sero_pos_1_30_39_capped <- pmin(model_sero_pos_1_30_39, pars$N_tot[4])
    model_sero_prob_pos_1_30_39 <- test_prob_pos(model_sero_pos_1_30_39_capped,
                                                 pars$N_tot[4] - model_sero_pos_1_30_39_capped,
                                                 pars$sero_sensitivity_1,
                                                 pars$sero_specificity_1,
                                                 exp_noise)
    model_sero_pos_1_40_49 <- state["sero_pos_1_40_49", ]
    model_sero_pos_1_40_49_capped <- pmin(model_sero_pos_1_40_49, pars$N_tot[5])
    model_sero_prob_pos_1_40_49 <- test_prob_pos(model_sero_pos_1_40_49_capped,
                                                 pars$N_tot[5] - model_sero_pos_1_40_49_capped,
                                                 pars$sero_sensitivity_1,
                                                 pars$sero_specificity_1,
                                                 exp_noise)
    model_sero_pos_1_50_59 <- state["sero_pos_1_50_59", ]
    model_sero_pos_1_50_59_capped <- pmin(model_sero_pos_1_50_59, pars$N_tot[6])
    model_sero_prob_pos_1_50_59 <- test_prob_pos(model_sero_pos_1_50_59_capped,
                                                 pars$N_tot[6] - model_sero_pos_1_50_59_capped,
                                                 pars$sero_sensitivity_1,
                                                 pars$sero_specificity_1,
                                                 exp_noise)
    model_sero_pos_1_60_69 <- state["sero_pos_1_60_69", ]
    model_sero_pos_1_60_69_capped <- pmin(model_sero_pos_1_60_69, pars$N_tot[7])
    model_sero_prob_pos_1_60_69 <- test_prob_pos(model_sero_pos_1_60_69_capped,
                                                 pars$N_tot[7] - model_sero_pos_1_60_69_capped,
                                                 pars$sero_sensitivity_1,
                                                 pars$sero_specificity_1,
                                                 exp_noise)
    model_sero_pos_1_70_plus <- state["sero_pos_1_70_plus", ]
    model_sero_pos_1_70_plus_capped <- pmin(model_sero_pos_1_70_plus, pars$N_tot[8])
    model_sero_prob_pos_1_70_plus <- test_prob_pos(model_sero_pos_1_70_plus_capped,
                                                   pars$N_tot[8] - model_sero_pos_1_70_plus_capped,
                                                   pars$sero_sensitivity_1,
                                                   pars$sero_specificity_1,
                                                   exp_noise)
    
    # Modelled cases and "non-variant" cases (1st strain)
    model_cases <- state["cases", ]
    model_cases_non_variant <- state["cases_non_variant", ]
    model_strain_prob_pos <- test_prob_pos(
        model_cases_non_variant,
        model_cases - model_cases_non_variant,
        1, 1, exp_noise)
    
    # Weekly modelled cases
    # model_cases_week <- state["cases_week", ]
    # test_negs <- pars$p_NC * (pars$N_tot - model_cases_week)
    test_negs <- p_NC * (sum(pars$N_tot) - model_cases)
    model_test_prob_pos <- test_prob_pos(model_cases,
                                         test_negs,
                                         pars$test_sensitivity,
                                         pars$test_specificity,
                                         exp_noise)
    
    # Log-likelihoods for deaths
    # ll_hosps <- ll_nbinom(observed$hosps,model_hosps,kappa_hosp,exp_noise)
    ll_hosps_0_39 <- ll_nbinom(observed$hosps_0_39,model_hosps_0_39,kappa_hosp,exp_noise)
    ll_hosps_40_49 <- ll_nbinom(observed$hosps_40_49,model_hosps_40_49,kappa_hosp,exp_noise)
    ll_hosps_50_59 <- ll_nbinom(observed$hosps_50_59,model_hosps_50_59,kappa_hosp,exp_noise)
    ll_hosps_60_69 <- ll_nbinom(observed$hosps_60_69,model_hosps_60_69,kappa_hosp,exp_noise)
    ll_hosps_70_plus <- ll_nbinom(observed$hosps_70_plus,model_hosps_70_plus,kappa_hosp,exp_noise)
    # ll_hosps_70_79 <- ll_nbinom(observed$hosps_70_79,model_hosps_70_79,kappa_hosp,exp_noise)
    # ll_hosps_80_plus <- ll_nbinom(observed$hosps_80_plus,model_hosps_80_plus,kappa_hosp,exp_noise)
    
    # hosps_by_age <- matrix(c(observed$hosps_0_39,observed$hosps_40_49,observed$hosps_50_59,observed$hosps_60_69,observed$hosps_70_plus),nrow = 1)
    # model_hosps_by_age <- cbind(model_hosps_0_39,model_hosps_40_49,model_hosps_50_59,model_hosps_60_69,model_hosps_70_plus)
    # ll_hosps <- ll_dirmnom(hosps_by_age,model_hosps_by_age,size_hosp,exp_noise)
    
    # Log-likelihoods for deaths
    # ll_deaths <- ll_nbinom(observed$deaths,model_deaths,kappa_death,exp_noise)
    ll_deaths_0_39 <- ll_nbinom(observed$deaths_0_39,model_deaths_0_39,kappa_death,exp_noise)
    ll_deaths_40_49 <- ll_nbinom(observed$deaths_40_49,model_deaths_40_49,kappa_death,exp_noise)
    ll_deaths_50_59 <- ll_nbinom(observed$deaths_50_59,model_deaths_50_59,kappa_death,exp_noise)
    ll_deaths_60_69 <- ll_nbinom(observed$deaths_60_69,model_deaths_60_69,kappa_death,exp_noise)
    ll_deaths_70_plus <- ll_nbinom(observed$deaths_70_plus,model_deaths_70_plus,kappa_death,exp_noise)
    # # ll_deaths_70_79 <- ll_nbinom(observed$deaths_70_79,model_deaths_70_79,kappa_death,exp_noise)
    # # ll_deaths_80_plus <- ll_nbinom(observed$deaths_80_plus,model_deaths_80_plus,kappa_death,exp_noise)
    
    # deaths_by_age <- matrix(c(observed$deaths_0_39,observed$deaths_40_49,observed$deaths_50_59,observed$deaths_60_69,observed$deaths_70_plus),nrow = 1)
    # model_deaths_by_age <- cbind(model_deaths_0_39,model_deaths_40_49,model_deaths_50_59,model_deaths_60_69,model_deaths_70_plus)
    # ll_deaths <- ll_dirmnom(deaths_by_age,model_deaths_by_age,size_death,exp_noise)
    
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
    
    # Log-likelihoods for cases
    ll_cases_0_9 <- ll_nbinom(observed$cases_0_9,model_confirmed_cases_0_9,kappa_cases,exp_noise)
    ll_cases_10_19 <- ll_nbinom(observed$cases_10_19,model_confirmed_cases_10_19,kappa_cases,exp_noise)
    ll_cases_20_29 <- ll_nbinom(observed$cases_20_29,model_confirmed_cases_20_29,kappa_cases,exp_noise)
    ll_cases_30_39 <- ll_nbinom(observed$cases_30_39,model_confirmed_cases_30_39,kappa_cases,exp_noise)
    ll_cases_40_49 <- ll_nbinom(observed$cases_40_49,model_confirmed_cases_40_49,kappa_cases,exp_noise)
    ll_cases_50_59 <- ll_nbinom(observed$cases_50_59,model_confirmed_cases_50_59,kappa_cases,exp_noise)
    ll_cases_60_69 <- ll_nbinom(observed$cases_60_69,model_confirmed_cases_60_69,kappa_cases,exp_noise)
    ll_cases_70_plus <- ll_nbinom(observed$cases_70_plus,model_confirmed_cases_70_plus,kappa_cases,exp_noise)
    
    # Log-likelihoods for seroprevalence
    # ll_sero_pos_1 <- ll_binom(observed$sero_pos_1,observed$sero_tot_1,model_sero_prob_pos_1)
    ll_sero_pos_1_20_29 <- ll_binom(observed$sero_pos_1_20_29,observed$sero_tot_1_20_29,model_sero_prob_pos_1_20_29)
    ll_sero_pos_1_30_39 <- ll_binom(observed$sero_pos_1_30_39,observed$sero_tot_1_30_39,model_sero_prob_pos_1_30_39)
    ll_sero_pos_1_40_49 <- ll_binom(observed$sero_pos_1_40_49,observed$sero_tot_1_40_49,model_sero_prob_pos_1_40_49)
    ll_sero_pos_1_50_59 <- ll_binom(observed$sero_pos_1_50_59,observed$sero_tot_1_50_59,model_sero_prob_pos_1_50_59)
    ll_sero_pos_1_60_69 <- ll_binom(observed$sero_pos_1_60_69,observed$sero_tot_1_60_69,model_sero_prob_pos_1_60_69)
    ll_sero_pos_1_70_plus <- ll_binom(observed$sero_pos_1_70_plus,observed$sero_tot_1_70_plus,model_sero_prob_pos_1_70_plus)

    # Log-likelihood for variant proportion
    ll_strain <- ll_binom(observed$strain_non_variant,observed$strain_tot,model_strain_prob_pos)
    
    # Log-likelihood for positive tests
    ll_tests <- ll_betabinom(observed$pos,
                             observed$tot,
                             model_test_prob_pos,
                             rho_tests)
    
    # Calculate total log-likelihood
    # ll_hosps + ll_deaths
    # ll_hosps + ll_hosps_0_39 + ll_hosps_40_49 + ll_hosps_50_59 + ll_hosps_60_69 + ll_hosps_70_plus + #ll_hosps_70_79 + ll_hosps_80_plus +
    # ll_deaths + ll_deaths_0_39 + ll_deaths_40_49 + ll_deaths_50_59 + ll_deaths_60_69 + ll_deaths_70_plus #+ ll_deaths_70_79 + ll_deaths_80_plus
    # ll_hosps_70_plus + ll_deaths_70_plus
    # ll_hosps + ll_deaths + ll_sero_pos_1 + 
    ll_strain + ll_tests +
    ll_cases_0_9 + ll_cases_10_19 + ll_cases_20_29 + ll_cases_30_39 + ll_cases_40_49 + ll_cases_50_59 + ll_cases_60_69 + ll_cases_70_plus +
        ll_hosps_0_39 + ll_hosps_40_49 + ll_hosps_50_59 + ll_hosps_60_69 + ll_hosps_70_plus +
        ll_deaths_0_39 + ll_deaths_40_49 + ll_deaths_50_59 + ll_deaths_60_69 + ll_deaths_70_plus +
        ll_sero_pos_1_20_29 + ll_sero_pos_1_30_39 + ll_sero_pos_1_40_49 + ll_sero_pos_1_50_59 + ll_sero_pos_1_60_69 + ll_sero_pos_1_70_plus
}


# Function for plotting fitted trajectories
plot_particle_filter <- function(history, true_history, times, idx, obs_end = NULL) {
    if (is.null(obs_end)) {
        obs_end = max(times)
    }
    
    # par(mfrow = c(2,4), oma=c(2,3,0,0))
    par(mfrow = c(1,1), oma=c(2,3,0,0))
    idx_E <- idx$state[grep("E",names(idx$state))]
    for (i in 1:n_age){
        par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
        cols <- c(S = "#8c8cd9", E = "#ffff00", I_P = "#cc0044", I_A = "green", I_C = "blue", R = "#999966", D = "#000000")
        # matplot(times, t(history[i, ,-1]), type = "l", # Offset to access numbers in age compartment
        #         xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
        #         col = alpha(cols[["S"]],0.1), lty = 1, ylim=range(history))
        matplot(times, t(history[i + n_age*n_vax + n_age, ,-1]), type = "l", # Offset to access numbers in age compartment
                xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
                col = alpha(cols[["E"]],0.1), lty = 1, ylim=range(true_history[idx_E - 2, ,-1]))
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
        points(times[1:obs_end], true_history[idx_E[1] - 2 + n_age + i, ,-1], pch = 19,
               col = cols[2])
        # legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
        legend("right", lwd = 1, col = cols[2:5], legend = names(cols[2:5]), bty = "n")
        axis(2, las = 2)
    }
}


# Probably delete these transform functions as they're out of date now
covid_multi_strain_transform <- function(dt,
                                         n_age,
                                         n_vax,
                                         m,
                                         beta_date,
                                         beta_type,
                                         gamma_E,
                                         gamma_P,
                                         gamma_A,
                                         gamma_C,
                                         gamma_H,
                                         gamma_G,
                                         gamma_pre_1,
                                         gamma_P_1,
                                         theta_A,
                                         p_C,
                                         p_H,
                                         p_G,
                                         p_D,
                                         p_P_1,
                                         population,
                                         # start_date,
                                         initial_seed_size,
                                         initial_seed_pattern,
                                         # strain_transmission,
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
                                         vacc_skip_progression_rate,
                                         vacc_skip_to,
                                         vacc_skip_weight,
                                         waning_rate,
                                         cross_immunity,
                                         sero_sensitivity_1,
                                         sero_specificity_1){
    
    function(pars){
        # beta <- pars[["beta"]]
        beta_value <- unname(pars[paste0("beta", seq_along(beta_date))])
        # ngm <- t(t(transmission)*(p_C*(1/pars[["gamma"]]+1/pars[["gamma"]])+(1-p_C)*theta_A/gamma_A)*population)
        # beta <- pars[["R0"]]/eigen(ngm)$value[1]
        start_date <- pars[["start_date"]]
        rel_strain_transmission <- pars[["rel_strain_transmission"]]
        strain_seed_date <- pars[["strain_seed_date"]]
        p_H_max <- pars[["p_H_max"]]
        p_D_max <- pars[["p_D_max"]]
        phi_cases <- pars[["phi_cases"]]
        kappa_cases <- 1/pars[["alpha_cases"]]
        
        # Parameters for 1st epoch
        p <- parameters(dt,
                        n_age,
                        n_vax,
                        m,
                        beta_date,
                        beta_value = beta_value,
                        beta_type,
                        gamma_E,
                        gamma_P,
                        gamma_A,
                        gamma_C,
                        gamma_H,
                        gamma_G,
                        gamma_pre_1,
                        gamma_P_1,
                        theta_A,
                        p_C,
                        # p_H,
                        p_H_max*p_H,
                        p_G,
                        # p_D,
                        p_D_max*p_D,
                        p_P_1,
                        population,
                        start_date = start_date,
                        initial_seed_size,
                        initial_seed_pattern,
                        strain_transmission = c(1,rel_strain_transmission),
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
                        vaccine_index_booster,
                        vaccine_catchup_fraction,
                        n_doses,
                        vacc_skip_progression_rate,
                        vacc_skip_to,
                        vacc_skip_weight,
                        waning_rate,
                        cross_immunity,
                        phi_cases,
                        kappa_cases,
                        sero_sensitivity_1,
                        sero_specificity_1)
        p
    }
}


covid_multi_strain_transform_multistage <- function(dt,
                                                    n_age,
                                                    n_vax,
                                                    m,
                                                    beta_date,
                                                    beta_type,
                                                    gamma_E,
                                                    gamma_P,
                                                    gamma_A,
                                                    gamma_C,
                                                    gamma_H,
                                                    gamma_G,
                                                    gamma_pre_1,
                                                    gamma_P_1,
                                                    theta_A,
                                                    p_C,
                                                    p_H,
                                                    p_G,
                                                    p_D,
                                                    p_P_1,
                                                    population,
                                                    # start_date,
                                                    initial_seed_size,
                                                    initial_seed_pattern,
                                                    # strain_transmission,
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
                                                    vacc_skip_progression_rate,
                                                    vacc_skip_to,
                                                    vacc_skip_weight,
                                                    waning_rate,
                                                    cross_immunity,
                                                    start_date1,
                                                    strain_rel_p_sympt1,
                                                    strain_rel_p_hosp_if_sympt1,
                                                    strain_rel_p_death1,
                                                    rel_susceptibility1,
                                                    rel_p_sympt1,
                                                    rel_p_hosp_if_sympt1,
                                                    rel_p_death1,
                                                    rel_infectivity1,
                                                    cross_immunity1,
                                                    sero_sensitivity_1,
                                                    sero_specificity_1){
    
    function(pars){
        # beta <- pars[["beta"]]
        beta_value <- unname(pars[paste0("beta", seq_along(beta_date))])
        # ngm <- t(t(transmission)*(p_C*(1/pars[["gamma"]]+1/pars[["gamma"]])+(1-p_C)*theta_A/gamma_A)*population)
        # beta <- pars[["R0"]]/eigen(ngm)$value[1]
        start_date <- pars[["start_date"]]
        rel_strain_transmission <- pars[["rel_strain_transmission"]]
        strain_seed_date <- pars[["strain_seed_date"]]
        p_H_max <- pars[["p_H_max"]]
        p_D_max <- pars[["p_D_max"]]
        rel_strain_transmission1 <- pars[["rel_strain_transmission1"]]
        strain_seed_date1 <- pars[["strain_seed_date1"]]
        phi_cases <- pars[["phi_cases"]]
        kappa_cases <- 1/pars[["alpha_cases"]]
        
        # Parameters for 1st epoch
        p <- parameters(dt,
                        n_age,
                        n_vax,
                        m,
                        beta_date,
                        beta_value = beta_value,
                        beta_type,
                        gamma_E,
                        gamma_P,
                        gamma_A,
                        gamma_C,
                        gamma_H,
                        gamma_G,
                        gamma_pre_1,
                        gamma_P_1,
                        theta_A,
                        p_C,
                        # p_H,
                        p_H_max*p_H,
                        p_G,
                        # p_D,
                        p_D_max*p_D,
                        p_P_1,
                        population,
                        start_date = start_date,
                        initial_seed_size,
                        initial_seed_pattern,
                        strain_transmission = c(1,rel_strain_transmission),
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
                        vaccine_index_booster,
                        vaccine_catchup_fraction,
                        n_doses,
                        vacc_skip_progression_rate,
                        vacc_skip_to,
                        vacc_skip_weight,
                        waning_rate,
                        cross_immunity,
                        phi_cases,
                        kappa_cases,
                        sero_sensitivity_1,
                        sero_specificity_1)
        
        # Parameters for 2nd epoch
        p1 <- parameters(dt,
                         n_age,
                         n_vax,
                         m,
                         beta_date,
                         beta_value = beta_value,
                         beta_type,
                         gamma_E,
                         gamma_P,
                         gamma_A,
                         gamma_C,
                         gamma_H,
                         gamma_G,
                         gamma_pre_1,
                         gamma_P_1,
                         theta_A,
                         p_C,
                         # p_H,
                         p_H_max*p_H,
                         p_G,
                         # p_D,
                         p_D_max*p_D,
                         p_P_1,
                         population,
                         start_date = start_date1,
                         initial_seed_size = 0,
                         initial_seed_pattern,
                         strain_transmission = c(rel_strain_transmission,rel_strain_transmission1),
                         strain_seed_date = strain_seed_date1,
                         strain_seed_size,
                         strain_seed_pattern,
                         strain_rel_p_sympt1,
                         strain_rel_p_hosp_if_sympt1,
                         strain_rel_p_death1,
                         rel_susceptibility1,
                         rel_p_sympt1,
                         rel_p_hosp_if_sympt1,
                         rel_p_death1,
                         rel_infectivity1,
                         vaccine_progression_rate,
                         schedule,
                         vaccine_index_dose2,
                         vaccine_index_booster,
                         vaccine_catchup_fraction,
                         n_doses,
                         vacc_skip_progression_rate,
                         vacc_skip_to,
                         vacc_skip_weight,
                         waning_rate,
                         cross_immunity1,
                         phi_cases,
                         kappa_cases,
                         sero_sensitivity_1,
                         sero_specificity_1)
        
        epochs <- list(
            multistage_epoch(start_date1, pars = p1, transform_state = transform_state)
        )
        p_multistage <- multistage_parameters(p, epochs)
        p_multistage
    }
}


plot_outcome <- function(incidence_modelled, incidence_observed, vrble, phi = NULL, 
                         by_age = FALSE, moving_avg = FALSE, pred_intvl = T, alpha = NULL){
    # Extract modelled incidence for outcome vrble from output
    nms <- dimnames(incidence_modelled$state)[[1]]
    if (by_age){
        idx <- grep(paste0(vrble,"_[0-9]+"),nms)
    } else {
        idx <- which(nms == vrble)
    }
    # Select variables to plot and drop initial conditions
    state <- incidence_modelled$state[idx,,-1,drop = F]
    # Convert to data table
    inc_dt <- as.data.table(state)
    names(inc_dt) <- c("age_group","sample","date","value")
    inc_dt[,date := covid_multi_strain_date_as_date(incidence_modelled$date[date])]
    inc_dt[,age_group := sub("_","-",sub(paste0(vrble,"_"),"",age_group))]
    if (vrble == "cases"){
        # Convert reporting fraction vector to data table with sample number
        phi_dt <- data.table(sample = seq_along(phi),phi = phi)
        # Add phi to inc_dt
        inc_dt[,phi := phi_dt[match(inc_dt[,sample],sample),phi]]
        # Multiply incidence by reporting rate
        inc_dt[,value := phi * value]
    }
    if (pred_intvl){
        # Convert overdispersion parameter vector to data table with sample number
        alpha_dt <- data.table(sample = seq_along(alpha),alpha = alpha)
        # Add alpha to inc_dt
        inc_dt[,alpha := alpha_dt[match(inc_dt[,sample],sample),alpha]]
        # Sample observed value
        inc_dt[,value_obs := rnbinom(.N,size = 1/alpha,mu = value)]
        # Calculate median and 95% quantiles
        q_inc_dt <- inc_dt[,.(med = quantile(value,0.5),
                              q95l = quantile(value,0.025),
                              q95u = quantile(value,0.975),
                              q95l_obs = quantile(value_obs,0.025),
                              q95u_obs = quantile(value_obs,0.975)),
                           by = .(age_group,date)]
    } else {
        # Calculate median and 95% quantiles
        q_inc_dt <- inc_dt[,.(med = quantile(value,0.5),
                              q95l = quantile(value,0.025),
                              q95u = quantile(value,0.975)),
                           by = .(age_group,date)]
    }
    
    # Extract observed incidence for outcome vrble from data frame
    if (by_age){
        idx_obs <- names(incidence_observed)[grep(paste0(vrble,"_"),names(incidence_observed))]
    } else {
        idx_obs <- names(incidence_observed)[names(incidence_observed) == vrble]    
    }
    inc_obs_dt <- as.data.table(incidence_observed[c("day_end",idx_obs)])
    inc_obs_dt <- melt(inc_obs_dt,measure.vars = patterns(vrble),variable.name = "age_group")
    inc_obs_dt[,date := covid_multi_strain_date_as_date(day_end)]
    inc_obs_dt[,age_group := sub("_","-",sub(paste0(vrble,"_"),"",age_group))]
    if (moving_avg){
        inc_obs_dt[,value := as.numeric(value)]
        inc_obs_dt[,value := frollmean(value,7,align = "center"),by = .(age_group)]
    }
    
    # Plot
    ylbl <- sub("(.)", ("\\U\\1"), tolower(vrble), pe=TRUE)
    if (ylbl == "Hosps"){
        ylbl <- "Hospitalisations"  
    } else if (ylbl == "Deaths"){
        ylbl <- "Hospital deaths"
    }
    p <- ggplot()
    if (moving_avg){
        p <- p + geom_line(aes(x = date, y = value), inc_obs_dt, size = 0.5, color = "red")
    } else {
        p <- p + geom_point(aes(x = date, y = value), inc_obs_dt, size = 0.5, color = "red")
    }
    p <- p + geom_line(aes(x = date, y = med),q_inc_dt) + 
        geom_ribbon(aes(x = date, ymin = q95l, ymax = q95u),q_inc_dt,alpha = 0.5) + 
        labs(x = "Date", y = ylbl) + 
        theme_cowplot(font_size = 12)
    if (pred_intvl){
        p <- p + geom_ribbon(aes(x = date, ymin = q95l_obs, ymax = q95u_obs),q_inc_dt,alpha = 0.2)
    }
    if (by_age){
        if (vrble == "cases"){ # 2 columns for cases as there are more age groups in the data 
            p <- p + facet_wrap(~age_group, ncol = 2) + theme(strip.background = element_blank())
        } else {
            p <- p + facet_wrap(~age_group, ncol = 1) + theme(strip.background = element_blank())
        }
    }
    return(p)
}

    
plot_tests <- function(trajectories, data, p_NC, population, moving_avg = F){
    state <- trajectories$state
    pos <- state["cases",,-1]
    test_negs <- p_NC * (population - pos)
    model_prob_pos <- test_prob_pos(pos,test_negs,1,1,1e6)
    
    q_model_prob_pos <- as.data.table(t(apply(model_prob_pos, 2, function(x) quantile(x,probs = c(0.025,0.5,0.975)))))
    q_model_prob_pos[,date := covid_multi_strain_date_as_date(trajectories$date)]
    
    data_dt <- as.data.table(data[c("day_end","pos","tot")])
    data_dt[,date := covid_multi_strain_date_as_date(day_end)]
    data_dt[,prop_pos := pos/tot]
    if (moving_avg){
        data_dt[,prop_pos := frollmean(prop_pos,7,align = "center")]
    }
    
    p <- ggplot()
    if (moving_avg){
        p <- p + geom_line(aes(x = date,y = prop_pos),data_dt,color = "red")
    } else {
        p <- p + geom_point(aes(x = date,y = prop_pos),data_dt,color = "red")
    }
    p <- p + geom_line(aes(x = date,y = `50%`),q_model_prob_pos) + 
        geom_ribbon(aes(x = date,ymin = `2.5%`,ymax = `97.5%`),q_model_prob_pos,alpha = 0.5) + 
        labs(x = "Date",y = "Proportion positive") + 
        theme_cowplot(font_size = 12)
    return(p)
}


plot_sero <- function(seroprev_modelled, seroprev_observed, population, by_age = FALSE){
    # Extract modelled seroprevalence from output
    nms <- dimnames(seroprev_modelled$state)[[1]]
    if (by_age){
        idx_sero <- grep("sero_pos_1_",nms) 
    } else {
        idx_sero <- which(nms == "sero_pos_1")
    }
    # Select variables to plot and drop initial conditions
    state <- seroprev_modelled$state[idx_sero,,-1,drop = F]
    # Convert to data table
    sero_dt <- as.data.table(state)
    names(sero_dt) <- c("age_group","sample","date","value")
    sero_dt[,date := covid_multi_strain_date_as_date(seroprev_modelled$date[date])]
    sero_dt[,age_group := sub("_","-",sub(paste0("sero_pos_1_"),"",age_group))]
    q_sero_dt <- sero_dt[,.(med = quantile(value,0.5),
                           q95l = quantile(value,0.025),
                           q95u = quantile(value,0.975)),
                        by = .(age_group,date)]
    if (by_age){
        pop_dt <- copy(population)
        pop_dt[age_group == "70+",age_group := "70-plus"]
        q_sero_dt <- merge(q_sero_dt,pop_dt,by = "age_group",all.x = T)
    } else {
        q_sero_dt[,population := population[!(age_group %in% c("0-9","10-19")),sum(population)]]
    }
    # Extract observed seroprevalence
    if (by_age){
        idx_sero_obs <- names(seroprev_observed)[grep("sero_pos_1_",names(seroprev_observed))]
        idx_sero_tot_obs <- names(seroprev_observed)[grep("sero_tot_1_",names(seroprev_observed))]
    } else {
        idx_sero_obs <- names(seroprev_observed)[names(seroprev_observed) == "sero_pos_1"]
        idx_sero_tot_obs <- names(seroprev_observed)[names(seroprev_observed) == "sero_tot_1"]        
    }
    sero_obs_dt <- as.data.table(seroprev_observed[c("day_end",idx_sero_obs,idx_sero_tot_obs)])
    sero_obs_dt <- melt(sero_obs_dt, measure.vars = patterns(positive = "sero_pos_1",total = "sero_tot_1"))
    sero_obs_dt[,date := covid_multi_strain_date_as_date(day_end)]
    # sero_obs_dt[,age_group := sub("_","-",sub(paste0("sero_pos_1_"),"",age_group))]
    sero_obs_dt[,age_group := sero_dt[,unique(age_group)][variable]]
    sero_obs_dt[!is.na(positive),`:=`(ci_lb = mapply(function(n,x) binom.test(n,x)$conf.int[1],positive,total),
                                      ci_ub = mapply(function(n,x) binom.test(n,x)$conf.int[2],positive,total))]
    
    # Plot
    p <- ggplot() + 
        geom_point(aes(x = date,y = positive/total),sero_obs_dt,color = "red") + 
        geom_errorbar(aes(x = date,ymin = ci_lb,ymax = ci_ub),sero_obs_dt,color = "red") + 
        geom_line(aes(x = date,y = pmin(med/population,1)),q_sero_dt) +
        geom_ribbon(aes(x = date,ymin = pmin(q95l/population,1),ymax = pmin(q95u/population,1)),q_sero_dt,alpha = 0.5) + 
        labs(x = "Date",y = "Seroprevalence") +
        theme_cowplot(font_size = 12)
    if (by_age){
        p <- p + facet_wrap(~age_group,ncol = 1) + theme(strip.background = element_blank())
    }
    return(p)
}


plot_outcome_by_age <- function(trajectories,vrble,phi,ttls,n_smpls,seed = 1){
    # Extract modelled incidence for outcome vrble from output
    nms <- dimnames(trajectories$state)[[1]]
    idx <- grep(paste0(vrble,"_[0-9]+",collapse = "|"),nms)
    # Select variables to plot and drop initial values
    state <- trajectories$state[idx,,-1,drop = F]
    set.seed(1)
    smpl <- sample.int(ncol(state),n_smpls) 
    state <- state[,smpl,,drop = F]
    # Convert to data table
    state_dt <- as.data.table(state)
    names(state_dt) <- c("age_group","sample","date","value")
    state_dt[,date := covid_multi_strain_date_as_date(trajectories$date[date])]
    state_dt[,state := gsub("_[0-9a-z]+","",age_group)]
    state_dt[,age_group := sub("_","-",sub(".*_([0-9]+_[0-9a-z]+)","\\1",age_group))]
    # Aggregate cases in <40 year-olds
    state_dt[state == "cases" & get_min_age(age_group) < 40,age_group := "0-39"]
    state_dt <- state_dt[,.(value = sum(value)),by = .(state,age_group,sample,date)]
    if ("cases" %in% vrble){
        # Convert to data table with sample number
        phi_dt <- data.table(sample = seq_along(phi),phi = phi)
        # Add phi to inc_dt
        state_dt[,phi := phi_dt[match(state_dt[,sample],sample),phi]]
        # Multiply incidence by reporting rate
        state_dt[,value := phi * value]        
    }
    state_dt[,prop := value/sum(value),by = .(state,sample,date)]
    q_state_dt <- state_dt[,unlist(lapply(.SD,function(x) 
        list(q95l = quantile(x,probs = 0.025,na.rm = T),
             med = quantile(x,probs = 0.5,na.rm = T),
             q95u = quantile(x,probs = 0.975,na.rm = T))),
        recursive = F),.SDcols = c("value","prop"),by = .(state,age_group,date)]
    # q_state_dt[state == "cases" & get_min_age(age_group) < 40,age_group := "0-39"]
    # q_state_dt <- q_state_dt[,lapply(.SD,sum),.SDcols = c("med","q95l","q95u"),by = .(state,age_group,date)]
    q_state_dt[,state := factor(state,levels = vrble)]
    p <- ggplot(q_state_dt,aes(x = date,y = value.med,fill = age_group)) + 
        geom_area() +
        labs(x = "Date") + 
        scale_fill_discrete(name = "Age group") + 
        facet_wrap(~state,scales = "free",labeller = labeller(state = ttls)) +
        theme_cowplot(font_size = 12) + 
        theme(axis.title.y = element_blank(),legend.position = "bottom",strip.background = element_blank())
    
    p1 <- ggplot(q_state_dt,aes(x = date,y = prop.med,fill = age_group)) +
        geom_area() +
        labs(x = "Date",y = "Proportion") +
        scale_fill_discrete(name = "Age group") +
        facet_wrap(~state,scales = "free",labeller = labeller(state = ttls)) + 
        theme_cowplot(font_size = 12) + 
        theme(legend.position = "bottom", strip.background = element_blank())
    
    return(list(p = p,p1 = p1))    
}


plot_variant_proportion <- function(trajectories, data, start_date, ylbl = "Variant proportion"){
    # Select variables and drop initial conditions
    state <- trajectories$state[c("cases","cases_non_variant"),,-1,drop = F]
    # Convert to data table
    state_dt <- as.data.table(state)
    names(state_dt) <- c("state","sample","date","value")
    state_dt[,date := covid_multi_strain_date_as_date(trajectories$date[date])]
    state_dt <- dcast(state_dt,sample + date ~ state,value.var = "value")
    state_dt[,vrnt_prop := 1 - cases_non_variant/cases]
    
    q_state_dt <- state_dt[,.(med = quantile(vrnt_prop,0.5),
                            q95l = quantile(vrnt_prop,0.025),
                            q95u = quantile(vrnt_prop,0.975)),
                         by = .(date)]
    
    data_dt <- as.data.table(data[c("day_end","strain_non_variant","strain_tot")])
    data_dt[,date := covid_multi_strain_date_as_date(day_end)]
    
    p <- ggplot() + 
        geom_point(aes(x = date,y = 1 - strain_non_variant/strain_tot), data_dt[date > start_date], color = "red") + 
        geom_line(aes(x = date,y = med),q_state_dt[date > start_date,]) + 
        geom_ribbon(aes(x = date,ymin = q95l,ymax = q95u),q_state_dt[date > start_date],alpha = 0.5) + 
        labs(x = "Date", y = ylbl) + 
        theme_cowplot(font_size = 12)
    return(p)
}


plot_cases <- function(cases_modelled, cases_observed, times){
    # par(mfrow = c(2,1), oma = c(2,3,0,0))
    par(mfrow = c(1,1))
    nms <- dimnames(cases_modelled)[[1]]
    idx_cases <- which(nms %in% c("cases","cases_non_variant"))
    idx_cases_obs <- which(names(cases_observed) %in% c("cases","cases_non_variant"))
    if (ncol(cases_modelled)>1){
        idx_plot <- seq(10,ncol(cases_modelled),by=10)    
    } else {
        idx_plot <- 1
    }
    dates_plot <- seq.Date(times[1], times[length(times)], by = 30)
    for (i in 1){#1:2){#
        par(mar = c(3, 4, 2, 0.5))
        if (ncol(cases_modelled)>1){
            y <- t(cases_modelled[idx_cases[1]-1+i,idx_plot,-1])
        } else {
            y <- cases_modelled[idx_cases[1]-1+i,idx_plot,-1]
        }
        matplot(times, y,
                type = "l", col = alpha("black",0.1), xlab = "Day", ylab = "Cases", xaxt = "n"
                # , #yaxt = "n",
                # main = rownames(cases_modelled)[idx_cases[1]-1+i]
                )
        points(times, cases_observed[[idx_cases_obs[1]-1+i]], pch = 19, col = "red")
        axis(1, dates_plot, format(dates_plot,"%Y-%m-%d"))
    }
}


plot_transmission_rate <- function(pars,beta_type,beta_date,dt,end_date){
    beta_value_post <- t(apply(pars[,seq_along(beta_date)],2,
                               function(x) quantile(x,probs = c(0.5,0.025,0.975))))
    if (beta_type == "piecewise-linear"){
        beta_step <- apply(
            beta_value_post,2,
            function(x) parameters_piecewise_linear(beta_date,
                                                    x %||% 0.1, dt))
    } else if (beta_type == "piecewise-constant"){
        beta_step <- apply(
            beta_value_post,2,
            function(x) parameters_piecewise_constant(beta_date,
                                                      x %||% 0.1, dt))
    }
    beta_t <- seq(0,end_date,by = dt)
    beta_step <- as.data.table(apply(
        beta_step,2,function(x) parameters_expand_step(seq_along(beta_t),x)))
    beta_step[,date := covid_multi_strain_date_as_date(beta_t)]
    p <- ggplot() + 
        geom_line(aes(x = date,y = `50%`,color = "Fitted"),beta_step) + 
        geom_ribbon(aes(x = date,ymin = `2.5%`,ymax = `97.5%`,fill = "Fitted"),beta_step,alpha = 0.5) + 
        labs(x = "Date",y = "beta(t)") + 
        theme_cowplot(font_size = 12)
    return(list(p = p,beta_step = beta_step))
}


plot_immune_status <- function(output,pop,age_groups,n_smpls = 1000,seed = 1){
    dat <- readRDS(output)
    
    set.seed(seed)
    smpl <- sample.int(ncol(dat$samples$trajectories$state),n_smpls)
    state <- dat$samples$trajectories$state[,smpl,,drop = F]
    
    rm(dat)
    gc()
    
    index_S <- grep("S_",dimnames(state)[[1]])
    nms_S <- dimnames(state)[[1]][index_S]
    # names(index_S) <- sub("_.*","",sub("_.*_","",nms_S))
    names(index_S) <- nms_S
    
    index_R <- grep("R_",dimnames(state)[[1]])
    nms_R <- dimnames(state)[[1]][index_R]
    # names(index_R) <- gsub("_[0-9]+|_[0-9]+_S[0-9]+","",nms_R)
    names(index_R) <- gsub("_S[0-9]+","",nms_R)
    
    
    # S <- state[index_S,,,drop = F]
    # index_R1 <- index_R[1:100]
    # index_R2 <- index_R[101:length(index_R)]
    # R1 <- state[index_R1,,,drop = F]
    # R2 <- state[index_R2,,,drop = F]
    
    # S_R <- abind1(abind1(S,R1),R2)
    
    index_S_R <- c(index_S,index_R)
    S_R <- state[index_S_R,,,drop = F]
    
    rm(state)
    gc()
    
    S_R_by_vacc <- sapply(
        seq_len(dim(S_R)[3]),
        function(k) rowsum(S_R[,,k],names(index_S_R),reorder = F), 
        simplify = "array")
    x <- sapply(
        seq_len(dim(S_R_by_vacc)[3]),
        function(k) rowsum(S_R_by_vacc[,,k],sub("[A-Za-z]+_([0-9]+).*","\\1",dimnames(S_R_by_vacc)[[1]]),reorder = F),
        simplify = "array")
    population <- pop[,.(sum(total)),by = .(age_group)][,V1]
    inf <- population - x
    
    # inf1 <- array(sum(population) - apply(S_R_by_vacc,c(2,3),sum),dim = c(1,dim(S_R_by_vacc)[2:3]))
    dimnames(inf)[[1]] <- paste0("Inf_",dimnames(inf)[[1]])
    
    S_R_by_vacc_inf <- abind1(S_R_by_vacc,inf)
    
    rm(S_R,S_R_by_vacc,inf,x)
    gc()
    
    # med_S_R_by_vacc_inf <- apply(S_R_by_vacc_inf,c(1,3),median)
    # 
    # matplot(t(med_S_R_by_vacc_inf),type = "l")
    # 
    # med_S_R_by_vacc_inf1 <- data.table(t(med_S_R_by_vacc_inf))
    # med_S_R_by_vacc_inf1[, time := as.numeric(rownames(med_S_R_by_vacc_inf1))]
    # med_S_R_by_vacc_inf_long <- melt(med_S_R_by_vacc_inf1,id.vars = "time")
    # 
    # # Calculate proportion
    # med_S_R_by_vacc_inf_long[,prop := value/sum(value), by = .(time)]
    # 
    # # Set plotting order
    # med_S_R_by_vacc_inf_long[,variable := factor(variable,levels = c("S","S1","S2","S3","S4","Infected","R","R_V1","R_V2","R_V3","R_V4"))]
    # 
    # ggplot(med_S_R_by_vacc_inf_long, aes(x = time,y = prop,fill = variable)) + geom_area()
    # ggsave("output/pop_immune_status.png",height = 4,width = 6)
    # 
    # View(med_S_R_by_vacc_inf_long[time==657,])
    
    S_R_by_vacc_inf_long <- as.data.table(S_R_by_vacc_inf)
    names(S_R_by_vacc_inf_long) <- c("state","sample","date","value")
    
    # Convert dates
    S_R_by_vacc_inf_long[, date := covid_multi_strain_date_as_date(date)]
    
    # Add age and vaccination status variables
    S_R_by_vacc_inf_long[,age := sub("[A-Za-z]+_([0-9]+).*","\\1",state)]
    min_ages <- get_min_age(age_groups)
    S_R_by_vacc_inf_long[,age_group := cut(as.numeric(age),c(min_ages,Inf),labels = age_groups,right = F)]
    S_R_by_vacc_inf_long[,vacc := sub("_[0-9]+","",state)]
    
    # Calculate proportion
    S_R_by_vacc_inf_long[,prop := value/sum(value),by = .(sample,date,age_group)]
    
    # Calculate median and 95% CI
    q_S_R_by_vacc_inf <- S_R_by_vacc_inf_long[,unlist(lapply(.SD,function(x) 
        list(q95l = quantile(x,probs = 0.025),
             med = quantile(x,probs = 0.5),
             q95u = quantile(x,probs = 0.975))),
        recursive = F),.SDcols = c("value","prop"),by = .(state,vacc,date,age_group)]
    
    # Set plotting order
    q_S_R_by_vacc_inf[,vacc := factor(vacc,levels = c("S","S_1","S_2","S_3","S_4","Inf","R","R_V1","R_V2","R_V3","R_V4"))]
    
    # Plot
    lbls = c(S = "Fully susc",
             S_1 = "Susc 1 dose",
             S_2 = "Susc 2 dose",
             S_3 = "Susc waned",
             S_4 = "Susc boosted",
             "Inf" = "Inf",
             R = "Rec unvacc",
             R_V1 = "Rec 1 dose",
             R_V2 = "Rec 2 dose",
             R_V3 = "Rec waned",
             R_V4 = "Rec boosted")
    p <- ggplot(q_S_R_by_vacc_inf, aes(x = date,y = prop.med,fill = vacc)) + 
        geom_area() + 
        labs(x = "Date",y = "Proportion",fill = "State") + 
        scale_fill_viridis_d(labels = lbls) +
        # theme(axis.text.x = element_text(angle = 90)) +
        facet_wrap(~age_group,nrow = 2) + 
        theme_cowplot(font_size = 10) + 
        theme(strip.background = element_blank())
    # ggsave(paste0("output/pop_immune_status_by_age",run,".pdf"),height = 5,width = 11)
    
    pop_immune_status_by_age <- q_S_R_by_vacc_inf[date == max(date),]
    pop_immune_status_by_age <- pop_immune_status_by_age[order(age_group,vacc),]
    
    tbl <- pop_immune_status_by_age[,.(
        State = lbls[vacc],
        `Age group` = age_group,
        Number = med_and_CI(value.med,value.q95l,value.q95u,
                            d = 3,method = "signif"),
        Percentage = med_and_CI(prop.med,prop.q95l,prop.q95u,
                                f = 100,d = 3,method = "signif")
    )]
    # write.csv(tbl,paste0("output/table2_by_age_",run,".csv"),row.names = F)
    
    # Proportion of total population in each immune state
    S_R_by_vacc_inf_long1 <- S_R_by_vacc_inf_long[,.(value = sum(value)),by = .(vacc,sample,date)]
    S_R_by_vacc_inf_long1[,prop := value/sum(value),by = .(sample,date)]
    
    # Calculate median and 95% CI
    q_S_R_by_vacc_inf1 <- S_R_by_vacc_inf_long1[,unlist(lapply(.SD,function(x) 
        list(q95l = quantile(x,probs = 0.025),
             med = quantile(x,probs = 0.5),
             q95u = quantile(x,probs = 0.975))),
        recursive = F),.SDcols = c("value","prop"),by = .(vacc,date)]
    
    # Set plotting order
    q_S_R_by_vacc_inf1[,vacc := factor(vacc,levels = c("S","S_1","S_2","S_3","S_4","Inf","R","R_V1","R_V2","R_V3","R_V4"))]
    
    # Plot
    p1 <- ggplot(q_S_R_by_vacc_inf1, aes(x = date,y = prop.med,fill = vacc)) + 
        geom_area() + 
        labs(x = "Date",y = "Proportion",fill = "State") + 
        scale_fill_viridis_d(labels = lbls) + 
        theme_cowplot(font_size = 12) + 
        theme(strip.background = element_blank())
    # ggsave(paste0("output/pop_immune_status",run,".pdf"),height = 4,width = 6)
    
    pop_immune_status <- q_S_R_by_vacc_inf1[date == max(date),]
    pop_immune_status <- pop_immune_status[q_S_R_by_vacc_inf1[,order(unique(vacc))],]
    
    tbl1 <- pop_immune_status[,.(
        State = lbls[vacc],
        Number = med_and_CI(value.med,value.q95l,value.q95u,
                            d = 3,method = "signif"),
        Percentage = med_and_CI(prop.med,prop.q95l,prop.q95u,
                                f = 100,d = 3,method = "signif")
    )]
    # write.csv(tbl1,paste0("output/table2_",run,".csv"),row.names = F)
    return(list(p = p,p1 = p1,tbl = tbl, tbl1 = tbl1))
}


plot_traces <- function(pars,u){
    pars_dt <- as.data.table(pars[-1,u])
    pars_dt[,iter := seq_len(nrow(pars)-1)]
    pars_long_dt <- melt(pars_dt,id.vars = "iter")
    
    p <- ggplot(pars_long_dt,aes(x = iter,y = value)) + geom_line() +
        facet_wrap(~variable,scales = "free",ncol = 3) + 
        labs(x = "Iteration") + 
        theme_cowplot(font_size = 10) + 
        theme(axis.title.y = element_blank(),
              axis.text.x = element_text(size = 8),
              strip.background = element_blank())
    
    return(p)
    # for (i in u){
    #     print(ggplot(pars_df,aes(x = iter,y = !!sym(colnames(pars)[i]))) + geom_line())
    # }
}


plot_posteriors <- function(pars,u,priors,pars_min,pars_max){
    pars <- pars[,u]
    pars_dt <- as.data.table(pars)
    pars_dt[,iter := seq_len(nrow(pars))]
    pars_long_dt <- melt(pars_dt,id.vars = "iter")
    pars_long_dt[,variable := factor(variable,levels = unique(variable))]
    
    priors_dt <- data.table()
    for (i in seq_along(u)){
        nm <- names(priors)[u[i]]
        x <- seq(min(pars_dt[[i]]),max(pars_dt[[i]]),length.out = 100)
        priors_dt[,c(paste0("x_",nm),paste0("y_",nm)) := 
                      list(x,
                           if(all(exp(priors[[u[i]]](x))==1)){
                               rep(1/(pars_max[u[i]]-pars_min[u[i]]),length(x))
                           } else {
                               exp(priors[[u[i]]](x))
                           })]
    }
    
    priors_long_dt <- melt(priors_dt,measure.vars = patterns(x = "x_",y = "y_"))
    priors_long_dt[,variable := names(priors)[u[variable]]]
    priors_long_dt[,variable := factor(variable,levels = unique(variable))]
        
    p <- ggplot() + 
        geom_histogram(aes(x = value,y = ..density..),pars_long_dt) +
        geom_line(aes(x = x,y = y),priors_long_dt,colour = "red") +
        facet_wrap(~variable,scales = "free",ncol = 3) + 
        labs(y = "Density") +
        theme_cowplot(font_size = 10) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 6),
              strip.background = element_blank())
    return(p)
}

plot_pairwise_correlation <- function(pars,u){
    pars <- pars[,u]
    pars_dt <- as.data.table(pars)
    p <- ggpairs(pars_dt) + 
        theme_cowplot(font_size = 10) +
        theme(axis.text.x = element_text(angle = 90),
              strip.background = element_blank(),
              strip.text.x = element_text(angle = 90),
              strip.text.y = element_text(angle = 0)
              )
    return(p)
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


combine_steps_groups <- function(step, n_groups, n_time_steps, n_strains,
                                 n_vacc_classes, p_step, rel_p,
                                 strain_rel_p) {
    
    ret <- vapply(
        seq_len(n_groups),
        function(i) {
            if (!is.null(dim(p_step))){
                outer(
                    parameters_expand_step(step, p_step[, i]),
                    rel_p[i, , ] * strain_rel_p
                )                
            } else {
                outer(
                    parameters_expand_step(step, p_step),
                    rel_p[i, , ] * strain_rel_p
                ) 
            }
        },
        array(0, c(n_time_steps, n_strains, n_vacc_classes))
    )
    
    ret <- pmin(ret, 1)
    ret <- aperm(ret, c(4, 2, 3, 1))
    
    ret
}

## Calculates the index of a given state and adds a suffix corresponding to
##  how the state is disaggregated (e.g. _V1, _V2 for two vacc classes,
##  or _S1V2 for strain 1 vaccine 2).
calculate_index <- function(index, state, suffix_list, suffix0 = NULL,
                            state_name = state) {
    if (is.null(suffix0)) {
        suffixes <- list()
    } else {
        suffixes <- list(suffix0)
    }
    for (i in seq_along(suffix_list)) {
        nm <- names(suffix_list)[[i]]
        if (length(nm) == 0) {
            nm <- ""
        }
        suffixes <- c(suffixes,
                      list(c("", sprintf("_%s%s", nm,
                                         seq_len(suffix_list[[i]] - 1L)))))
    }
    suffixes <- expand.grid(suffixes)
    nms <- apply(suffixes, 1,
                 function(x) sprintf("%s%s",
                                     state_name, paste0(x, collapse = "")))
    set_names(index[[state]], nms)
}