calculate_rel_param <- function(vax_eff_arr){
    # Create matrix of 0s for "vaccine effectiveness" in unvaccinated individuals
    unvax_eff <- array(0, dim = dim(vax_eff_arr)[c(1:2,4)])
    # Bind to vaccine effectiveness array
    e <- abind(unvax_eff,vax_eff_arr,along = 3)
    # Name first layer
    dimnames(e)[[3]][1] <- "unvaccinated"
    
    # Calculate relative parameters
    rel_susceptibility <- 1 - e[,,,"infection"]
    rel_p_sympt <- (e[,,,"symptoms"] - e[,,,"infection"])/(1 - e[,,,"infection"])
    rel_p_hosp_if_sympt <- (e[,,,"hospitalisation"] - e[,,,"symptoms"])/(1 - e[,,,"symptoms"])
    rel_p_death <- (e[,,,"death"] - e[,,,"hospitalisation"])/(1 - e[,,,"hospitalisation"])
    rel_infectivity <- 1 - e[,,,"infectiousness"]
    
    # Ensure first layers for unvaccinated individuals are 1 (as parameters are relative)
    rel_p_sympt[,,1] <- 1
    rel_p_hosp_if_sympt[,,1] <- 1
    rel_p_death[,,1] <- 1
    
    return(list(rel_susceptibility = rel_susceptibility,
                rel_p_sympt = rel_p_sympt,
                rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
                rel_p_death = rel_p_death,
                rel_infectivity = rel_infectivity))
}

# build_rel_param <- function(rel_param, n_groups, n_vacc_classes, name_param) {
#     if (length(rel_param) == 1) {
#         mat_rel_param <- array(rel_param, c(n_groups, n_vacc_classes))
#     } else if (is.array(rel_param)) {
#         if (length(dim(rel_param)) != 2) {
#             stop(paste(name_param, "should be a two dimensional array with",
#                        "dimensions: age groups, vaccine classes"))
#         }
#         if (nrow(rel_param) != n_groups) {
#             stop(paste(name_param, "should have as many rows as age groups"))
#         }
#         if (dim(rel_param)[2] != n_vacc_classes) {
#             stop(paste(name_param,
#                        "should have number of vaccine classes as 3rd dimension"))
#         }
#         mat_rel_param <- rel_param
#     } else { # create array by repeating rel_param for each age group
#         mat_rel_param <-
#             array(rep(rel_param, each = n_groups),
#                   dim = c(n_groups, n_vacc_classes))
#     }
#     mat_rel_param
# }


# TODO: check this works for seiragevax
build_rel_param <- function(rel_param, n_groups, n_strains, n_vacc_classes, name_param) {
    if (length(rel_param) == 1) {
        mat_rel_param <- array(rel_param, c(n_groups, n_strains, n_vacc_classes))
    } else if (is.array(rel_param)) {
        if (length(dim(rel_param)) != 3) {
            stop(paste(name_param, "should be a three dimensional array with",
                       "dimensions: age groups, strains, vaccine classes"))
        }
        if (nrow(rel_param) != n_groups) {
            stop(paste(name_param, "should have as many rows as age groups"))
        }
        if (ncol(rel_param) != n_strains) {
            stop(paste(name_param, "should have as many columns as strains"))
        }
        if (dim(rel_param)[3] != n_vacc_classes) {
            stop(paste(name_param,
                       "should have number of vaccine classes as 3rd dimension"))
        }
        mat_rel_param <- rel_param
    } else { # create array by repeating rel_param for each age group and strain
        mat_rel_param <-
            array(rep(rel_param, each = n_groups * n_strains),
                  dim = c(n_groups, n_strains, n_vacc_classes))
    }
    mat_rel_param
}


build_vaccine_progression_rate <- function(vaccine_progression_rate, n_groups, n_vacc_classes, index_dose) {
    # if NULL, set vaccine_progression_rate to 0
    if (is.null(vaccine_progression_rate)) {
        mat_vaccine_progression_rate <- matrix(0, n_groups, n_vacc_classes)
    } else {
        if (is.matrix(vaccine_progression_rate)) {
            if (nrow(vaccine_progression_rate) != n_groups) {
                stop(
                    "'vaccine_progression_rate' must have as many rows as age groups")
            }
            if (ncol(vaccine_progression_rate) != n_vacc_classes) {
                stop(
                    "'vaccine_progression_rate' must have 'n_vacc_classes' columns")
            }
            if (any(vaccine_progression_rate < 0)) {
                stop("'vaccine_progression_rate' must have only non-negative values")
            }
            mat_vaccine_progression_rate <- vaccine_progression_rate
        } else { # vaccine_progression_rate vector of length n_vacc_classes
            if (!is.vector(vaccine_progression_rate) ||
                length(vaccine_progression_rate) != n_vacc_classes) {
                m1 <- "'vaccine_progression_rate' must be either:"
                m2 <- "a vector of length 'n_vacc_classes' or"
                m3 <- "a matrix with 'n_groups' rows and 'n_vacc_classes' columns"
                stop(paste(m1, m2, m3))
            }
            if (any(vaccine_progression_rate < 0)) {
                stop("'vaccine_progression_rate' must have only non-negative values")
            }
            # create matrix by repeating vaccine_progression_rate for each age group
            mat_vaccine_progression_rate <-
                matrix(rep(vaccine_progression_rate, each = n_groups), nrow = n_groups)
        }
    }
    for (i in seq_along(index_dose)) {
        j <- index_dose[[i]]
        if (!all(mat_vaccine_progression_rate[, j] == 0)) {
            stop(sprintf(
                "Column %d of 'vaccine_progression_rate' must be zero (dose %d)",
                j, i))
        }
    }
    mat_vaccine_progression_rate
}


##' The French Ministry of Health priority groups, in descending order are:
##'
##'  1. care home residents and staff over 50
##'  2. individuals 75+ years not in a care home
##'  3. individuals aged 65-74 years
##'  4. individuals age 50-64 years
##'  5. all individuals 18+ years 
##' @title Compute vaccination order
##'
##' @param uptake A vector of length 19 with fractional uptake per
##'   group. If a single number is given it is shared across all
##'   groups (note that this includes under-18s)
##'
##' @param prop_hcw Assumed fraction of healthcare workers in each
##'   group (length 19) - if `NULL` we use a default that is a guess
##'   with hopefully the right general shape.
##'
##' @param prop_very_vulnerable Assumed fraction "very vulnerable" in
##'   each group (length 19) - if `NULL` we use a default that is a
##'   guess with hopefully the right general shape.
##'
##' @param prop_underlying_condition Assumed fraction "underlying
##'   condition" in each group (length 19) - if `NULL` we use a
##'   default that is a guess with hopefully the right general shape.
##'
##' @return A matrix with n_groups rows (19) and columns representing
##'   priority groups, and element (i, j) is the proportion in group i
##'   who should be vaccinated as part of priority group j, accounting
##'   for uptake so that the sum over the rows corresponds to the
##'   total fractional uptake in that group. For
##'   `vaccine_priority_population`, the total number of
##'   individuals replaces the proportion (based on the demography
##'   used by sircovid).
##'
##' @rdname vaccine_priority
##' @export
vaccine_priority_proportion <- function(n_groups,
                                        uptake,
                                        prop_hcw = NULL,
                                        prop_very_vulnerable = NULL,
                                        prop_underlying_condition = NULL) {
    
    if (is.null(uptake)) {
        uptake <- rep(1, n_groups)
    } else {
        uptake <- recycle(uptake, n_groups)
    }
    
    prop_hcw <- prop_hcw %||%
        c(rep(0, 4), rep(0.1, 10), rep(0, 5))
    prop_very_vulnerable <- prop_very_vulnerable %||%
        c(rep(0, 4), rep(0.05, 5), rep(0.1, 5), rep(0.15, 5))
    prop_underlying_condition <- prop_underlying_condition %||%
        c(rep(0, 4), rep(0.05, 5), rep(0.1, 5), rep(0.15, 5))
    
    n_priority_groups <- 4
    p <- matrix(0, n_groups, n_priority_groups)
    
    ## Age-based priority list
    priority <- list(
        ## the age groups targeted in each priority group (see comments above)
        8, 7, 6, 1:5)
    
    # ## 1. Start with non age-based priority:
    # ## helper function
    # add_prop_to_vacc <- function(j, idx, prop_to_vaccinate, p) {
    #     p[idx, j] <- prop_to_vaccinate[idx]
    #     p
    # }
    # 
    # ## Group 2 includes frontline health and social care workers
    # p <- add_prop_to_vacc(j = 2,
    #                       idx = seq_len(priority[[2]] - 1),
    #                       prop_to_vaccinate = prop_hcw,
    #                       p)
    # 
    # ## Group 4 includes clinically extremely vulnerable individuals
    # p <- add_prop_to_vacc(j = 4,
    #                       idx = seq_len(priority[[4]] - 1),
    #                       prop_to_vaccinate = (1 - prop_hcw) *
    #                           prop_very_vulnerable,
    #                       p)
    # 
    # ## Group 6 includes all individuals aged 16 years to 64 years with
    # ## underlying health conditions which put them at higher risk of
    # ## serious disease and mortality
    # p <- add_prop_to_vacc(j = 6,
    #                       idx = 4:13,
    #                       prop_to_vaccinate = (1 - prop_hcw) *
    #                           prop_underlying_condition,
    #                       p)
    
    ## 2. Add age-based priority
    for (j in seq_along(priority)) {
        if (!is.null(priority[[j]])) {
            ## discount those already vaccinated as part of non age-based priority
            p[priority[[j]], j] <- 1 - rowSums(p)[priority[[j]]]
        }
    }
    
    ## 3. Account for uptake
    uptake_mat <- matrix(rep(uptake, n_priority_groups),
                         nrow = n_groups)
    
    p * uptake_mat
}


##' @param pop Vector of total population numbers in groups
##'
##' @rdname vaccine_priority
##' @export
vaccine_priority_population <- function(pop,
                                        uptake,
                                        prop_hcw = NULL,
                                        prop_very_vulnerable = NULL,
                                        prop_underlying_condition = NULL) {
    n_groups <- length(pop)
    p <- vaccine_priority_proportion(n_groups,
                                     uptake,
                                     prop_hcw,
                                     prop_very_vulnerable,
                                     prop_underlying_condition)
    pop_mat <- matrix(rep(pop, ncol(p)), nrow = nrow(p))
    round(p * pop_mat)
}


##' Create future vaccination schedule from a projected number of
##' daily doses.
##'
##' @title Create vaccination schedule
##'
##' @param start Either a [sircovid_date] object corresponding to the
##'   first date in daily_doses_value, or a [vaccine_schedule] object
##'   corresponding to previously carried out vaccination.
##'
##' @param daily_doses_value A vector of doses per day.
##'
##' @param mean_days_between_doses Assumed mean days between doses one
##'   and two
##'
##' @param priority_population Output from
##'   [vaccine_priority_population], giving the number of people
##'   to vaccinate in each age (row) and priority group (column)
##'
##' @param lag_groups Row indices, corresponding to age
##'   groups in which a lag should be added to the start time of the dose
##'   schedule returned by [vaccine_schedule], if NULL then no lag is added.
##'   Ignored if `lag_groups` is NULL.
##'
##' @param lag_days If `lag_groups` is not NULL then specifies the number of
##'  days to add the start of the dose schedule for the given groups. Ignored
##'  if `lag_groups` is NULL.
##'
##' @param booster_daily_doses_value A vector of booster doses per day.
##'
##' @param booster_proportion Proportion of the groups in
##'  `priority_population` to boost, default is all groups; ignored if
##'  `booster_daily_doses_value` is NULL.
##'
##' @export
vaccine_schedule_future <- function(start,
                                    daily_doses_value,
                                    mean_days_between_doses,
                                    priority_population,
                                    lag_groups = NULL,
                                    lag_days = NULL,
                                    booster_daily_doses_value = NULL,
                                    booster_proportion = rep(1L, 8)) {
    
    has_booster <- !is.null(booster_daily_doses_value)
    
    n_groups <- nrow(priority_population)
    n_priority_groups <- ncol(priority_population)
    
    if (has_booster) {
        n_doses <- 3L
        n_days <- max(length(daily_doses_value), length(booster_daily_doses_value))
    } else {
        n_doses <- 2L
        n_days <- length(daily_doses_value)
    }
    
    population_to_vaccinate_mat <-
        array(0, c(n_groups, n_priority_groups, n_doses, n_days))
    
    population_left <- array(rep(c(priority_population), n_doses),
                             c(n_groups, n_priority_groups, n_doses))
    
    if (inherits(start, "vaccine_schedule")) {
        for (dose in seq_len(min(n_doses, ncol(start$doses)))) {
            n <- rowSums(start$doses[, dose, ])
            for (i in seq_len(n_priority_groups)) {
                m <- pmin(n, population_left[, i, dose])
                n <- n - m
                population_left[, i, dose] <- population_left[, i, dose] - m
            }
        }
        
        daily_doses_prev <- apply(start$doses, c(2, 3), sum)
        n_prev <- ncol(daily_doses_prev)
        daily_doses_date <- start$date
    } else {
        daily_doses_prev <- matrix(0, n_doses, 0)
        n_prev <- 0L
        daily_doses_date <- start
    }
    
    if (nrow(daily_doses_prev) < n_doses) {
        daily_doses_prev <-
            rbind(daily_doses_prev,
                  matrix(0, n_doses - nrow(daily_doses_prev),
                         ncol(daily_doses_prev)))
    }
    daily_doses_tt <- cbind(daily_doses_prev, matrix(0, n_doses, n_days))
    
    population_to_vaccinate_mat <- vaccination_schedule_exec(
        daily_doses_tt, daily_doses_value, population_left,
        population_to_vaccinate_mat, mean_days_between_doses, n_prev, 1:2)
    if (has_booster) {
        booster_population_left <- round(population_left * booster_proportion)
        
        population_to_vaccinate_mat <- vaccination_schedule_exec(
            daily_doses_tt, booster_daily_doses_value, booster_population_left,
            population_to_vaccinate_mat, Inf, n_prev, 3)
    }
    
    doses <- apply(population_to_vaccinate_mat, c(1, 3, 4), sum)
    
    if (inherits(start, "vaccine_schedule")) {
        if (ncol(start$doses) < n_doses) {
            sdoses <- abind2(start$doses,
                             array(0, c(nrow(start$doses),
                                        n_doses - ncol(start$doses),
                                        nlayer(start$doses))))
        } else {
            sdoses <- start$doses
        }
        doses <- mcstate::array_bind(sdoses, doses)
        dimnames(doses) <- NULL
    }
    
    schedule <- vaccine_schedule(daily_doses_date, doses, n_doses)
    
    if (!is.null(lag_groups) || !is.null(lag_days)) {
        if (has_booster) {
            stop("Someone should think about how boost and lag interact")
        }
        if (is.null(lag_days) || is.null(lag_groups)) {
            stop("'lag_days' must be non-NULL iff 'lag_groups' is non_NULL")
        }
        
        ## we could do this in apply but loop is more readable
        for (i in lag_groups) {
            ## original schedule for group i
            old_schedule_group_i <- schedule$doses[i, , ]
            nr <- nrow(old_schedule_group_i)
            nc <- ncol(old_schedule_group_i)
            
            ## get start of dose schedule for the group
            start <- which(old_schedule_group_i[1, ] > 0)[1]
            ## catch case when groups are ineligible
            
            if (!is.na(start)) {
                ## initialize 0 array with same dimensions
                new_schedule_group_i <- array(0, dim(old_schedule_group_i))
                ## add dose schedule to new array after adding the given lag (also
                ##  truncate end by lag amount)
                new_schedule_group_i[, seq.int(start + lag_days, nc)] <-
                    old_schedule_group_i[, seq.int(start, nc - lag_days)]
                ## save new schedule
                schedule$doses[i, , ] <- new_schedule_group_i
            }
        }
    }
    
    schedule
}


##' Create a vaccine schedule
##'
##' @title Create vaccine schedule
##'
##' @param date A single date, representing the first day that
##'   vaccines will be given
##'
##' @param doses A 3d array of doses representing (1) the model group
##'   (19 rows for the lancelot model), (2) the dose (must be length
##'   2 at present) and (3) time (can be anything nonzero). The values
##'   represent the number of vaccine doses in that group for that
##'   dose for that day. So for `doses[i, j, k]` then it is for the
##'   ith group, the number of jth doses on day `(k - 1) + date`
##'
##' @param n_doses The number of doses in the schedule. Typically (and
##'   by default) this will be 2, but if using booster doses 3 (and in
##'   future we may extend further).
##'
##' @return A `vaccine_schedule` object
##' @export
vaccine_schedule <- function(date, doses, n_doses = 2L) {
    n_groups <- nrow(doses)
    
    if (length(dim(doses)) != 3L) {
        stop("Expected a 3d array for 'doses'")
    }
    if (nrow(doses) != n_groups) {
        stop(sprintf("'doses' must have %d rows", n_groups))
    }
    if (ncol(doses) != n_doses) {
        stop(sprintf("'doses' must have %d columns", n_doses))
    }
    if (dim(doses)[[3]] == 0) {
        stop("'doses' must have at least one element in the 3rd dimension")
    }
    if (any(is.na(doses))) {
        stop("'doses' must all be non-NA")
    }
    if (any(doses < 0)) {
        stop("'doses' must all be non-negative")
    }
    
    ret <- list(date = date, doses = doses, n_doses = n_doses)
    class(ret) <- "vaccine_schedule"
    ret
}


vaccination_schedule_exec <- function(daily_doses_tt, daily_doses_value,
                                      population_left,
                                      population_to_vaccinate_mat,
                                      mean_days_between_doses,
                                      n_prev, dose_index) {
    
    n_priority_groups <- ncol(population_left)
    i1 <- dose_index[[1]]
    i2 <- if (length(dose_index) == 2) dose_index[[2]] else dose_index[[1]]
    for (t in seq_along(daily_doses_value)) {
        tt <- t + n_prev
        tt_dose_1 <- tt - mean_days_between_doses
        if (tt_dose_1 >= 1) {
            ## If we have promised more 2nd doses than we can deliver, we
            ## move our debt forward in time by one day. If doses fluctuate
            ## this will eventually be paid off.
            if (daily_doses_tt[i1, tt_dose_1] > daily_doses_value[t]) {
                daily_doses_tt[i1, tt_dose_1 + 1] <-
                    daily_doses_tt[i1, tt_dose_1 + 1] +
                    (daily_doses_tt[i1, tt_dose_1] - daily_doses_value[t])
            }
            daily_doses_tt[i2, tt] <- min(daily_doses_value[t],
                                          daily_doses_tt[i1, tt_dose_1])
            daily_doses_tt[i1, tt] <- daily_doses_value[t] - daily_doses_tt[i2, tt]
        } else {
            ## Only distribute first doses
            daily_doses_tt[i2, tt] <- 0
            daily_doses_tt[i1, tt] <- daily_doses_value[t]
        }
        daily_doses_today <- daily_doses_tt[, tt]
        
        for (dose in dose_index) {
            eligible <- colSums(population_left[, , dose, drop = FALSE])
            ## Vaccinate the entire of the top priority groups
            n_full_vacc <- findInterval(daily_doses_today[dose], cumsum(eligible))
            if (n_full_vacc > 0) {
                i_full_vacc <- seq_len(n_full_vacc)
                population_to_vaccinate_mat[, i_full_vacc, dose, t] <-
                    population_left[, i_full_vacc, dose]
            }
            
            ## Then partially vaccinate the next priority group, if possible
            if (n_full_vacc < n_priority_groups) {
                if (n_full_vacc == 0) {
                    remaining_eligible <- daily_doses_today[dose]
                } else {
                    remaining_eligible <- daily_doses_today[dose] -
                        cumsum(eligible)[n_full_vacc]
                }
                i_vacc <- n_full_vacc + 1L
                
                ## Split remaining doses according to age
                population_to_vaccinate_mat[, i_vacc, dose, t] <-
                    round(remaining_eligible * population_left[, i_vacc, dose] /
                              sum(population_left[, i_vacc, dose]))
            }
            
            population_left[, , dose] <- population_left[, , dose] -
                population_to_vaccinate_mat[, , dose, t]
        }
    }
    
    population_to_vaccinate_mat
}