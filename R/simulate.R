simulate <- function(gen_mod, p, n_steps, deterministic = FALSE, 
                     keep_all_states = TRUE, min_ages = seq(0,70,by = 10), 
                     Rt = FALSE, p1 = NULL, n_steps1 = NULL, transform = NULL){
    # Create instance of model
    mod <- gen_mod$new(p,step = 0,n_particles = 1,n_threads = 1,seed = 1,
                       deterministic = deterministic)
    
    info <- mod$info()
    initial_state <- initial(info, NULL, p)
    
    mod$update_state(state = initial_state)
    
    # Create an array to contain outputs after looping the model.
    x <- array(NA, dim = c(info$len, 1, n_steps+1))
    
    # For loop to run the model iteratively
    x[ , ,1] <- mod$state()
    for (t in seq_len(n_steps)) {
        x[ , ,t+1] <- mod$run(t)
    }
    
    if (!is.null(p1)){
        # Apply transform function to model state
        state <- mod$state()
        state1 <- transform(state,info)
        
        # Update model parameters and state
        mod$update_state(pars = p1,state = state1)
        
        # Create data to be fitted to
        x1 <- array(NA, dim = c(info$len, 1, n_steps1-n_steps))
        
        # For loop to run the model iteratively
        # x1[ , ,1] <- mod$state()
        for (t in seq_len(n_steps1-n_steps)) {
            # print(t)
            x1[ , ,t] <- mod$run(n_steps+t)
        }
        
        # Join with first epoch
        x <- array_bind(x,x1)        
    }
    
    if (!keep_all_states){
        idx <- index(info, min_ages = min_ages, Rt = Rt)
        x <- x[idx$state, , ,drop = FALSE]
    }
    
    return(list(x = x, info = info))
}


plot_trajectories <- function(time,x,n_age,n_strains,n_vax){
    # Plot trajectories
    for (k in 1:n_vax){
        for (j in 1:n_strains) {
            par(mfrow = c(2,4), oma=c(2,3,0,0))
            for (i in 1:n_age) {
                par(mar = c(3, 4, 2, 0.5))
                cols <- c(S = "#8c8cd9", E = "#ffff00", I_P = "#cc0044", I_A = "green", I_C = "blue", R = "#999966", D = "#000000")
                matplot(time, x[29 + i + (k-1)*n_age, ,-1], type = "l", # Offset to access numbers in age compartment
                        xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
                        col = cols[["S"]], lty = 1, ylim=range(x[-(1:29),,]))
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains, ,-1], col = cols[["E"]], lty = 1)
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains + n_age*n_strains*n_vax, ,-1], col = cols[["I_P"]], lty = 1)
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 2*n_age*n_strains*n_vax, ,-1], col = cols[["I_A"]], lty = 1)
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 3*n_age*n_strains*n_vax, ,-1], col = cols[["I_C"]], lty = 1)
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 4*n_age*n_strains*n_vax, ,-1], col = cols[["R"]], lty = 1)
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 7*n_age*n_strains*n_vax, ,-1], col = cols[["D"]], lty = 1)
                legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
                axis(2, las = 2)
            }
            mtext("Number of individuals", side=2, line=1, outer=T, las=0)
            mtext("Time", side = 1, line = 0, outer =T)    
        }    
    }
}


change_beta_date <- function(p,pars,beta_type,beta_date_cntfctl){
    beta_value <- pars[sprintf("beta%s", seq_along(beta_date_cntfctl))]
    # Construct time-varying transmission coefficient
    if (beta_type == "piecewise-linear") {
        beta_step <- parameters_piecewise_linear(beta_date_cntfctl, 
                                                 beta_value, p$dt)
    } else {
        beta_step <- parameters_piecewise_constant(beta_date_cntfctl, 
                                                   beta_value, p$dt)
    }
    p$beta_step <- beta_step
    return(p)
}


change_vaccine_schedule <- function(p,schedule_cntfctl){
    # Construct vaccination parameters
    vaccination <- parameters_vaccination(p$N_tot,
                                          p$dt,
                                          p$rel_susceptibility,
                                          p$rel_p_sympt,
                                          p$rel_p_hosp_if_sympt,
                                          p$rel_p_death,
                                          p$rel_infectivity,
                                          p$vaccine_progression_rate_base,
                                          vaccine_schedule = schedule_cntfctl,
                                          p$index_dose[2],
                                          p$index_dose[3],
                                          p$n_strains, # now up to 2 strains
                                          p$vaccine_catchup_fraction,
                                          p$n_doses)
    
    # Overwrite vaccine_dose_step with value for counterfactual schedule
    p$vaccine_dose_step <- vaccination$vaccine_dose_step
    
    return(p)
}


change_booster_timing <- function(schedule, days_earlier){
    # date_booster_start <- min(which(schedule$doses[,3,] != 0,arr.ind = T)[,2])
    schedule_cntfctl <- schedule
    doses <- schedule$doses
    if (days_earlier >= 0){
        schedule_cntfctl$doses[,3,] <- doses[,3,c((days_earlier+1):nlayer(doses),rep(nlayer(doses),days_earlier))]    
    } else {
        days_later <- -days_earlier
        # FOR NOW: Copy first day's doses for days_later days (N.B. assumes 0 booster doses on first day)
        schedule_cntfctl$doses[,3,] <- doses[,3,c(rep(1,days_later),1:(nlayer(doses)-days_later))]
    }
    return(schedule_cntfctl)
}


## Run counterfactual simulations
simulate_counterfactual <- function(output,n_smpls,beta_date_cntfctl,schedule_cntfctl,burnin = NULL,seed = 1){
    # Load MCMC output
    load(output)
    
    if (is.null(burnin)){
        burnin <- round((n_iters/thinning)/10)
    }
    # Remove burn-in
    pars <- res$pars[-(1:burnin),]
        
    out <- vector("list", n_smpls)
    set.seed(seed)
    smpl <- sample.int(nrow(pars),n_smpls)
    
    # Extract posterior samples of trajectories
    states <- res$trajectories$state[,burnin + smpl,]
    
    # Remove res object as it is very large
    rm(res)
    
    for (i in seq_len(n_smpls)){
        j <- smpl[i] #sample.int(nrow(pars),1)
        pars_i <- pars[j,]
        if (is.null(names(pars_i))){
            names(pars_i) <- names(init_pars)
        }
        transform_pars <- transform(pars_i)
        if (class(transform_pars) == "multistage_parameters"){
            # TODO: Make this work with one list for p rather than separate 
            # objects for each epoch, so it works for an arbitrary number of epochs
            p_i <- transform_pars[[1]]$pars
            p_i <- change_beta_date(p_i, pars_i, beta_type, beta_date_cntfctl)
            p_i <- change_vaccine_schedule(p_i, schedule_cntfctl)
            p1_i <- transform_pars[[2]]$pars
            p1_i <- change_beta_date(p1_i, pars_i, beta_type, beta_date_cntfctl)
            p1_i <- change_vaccine_schedule(p1_i, schedule_cntfctl)
        } else {
            p_i <- transform_pars
            p_i <- change_beta_date(p_i, pars_i, beta_type, beta_date_cntfctl)
            p_i <- change_vaccine_schedule(p_i, schedule_cntfctl)
            p1_i <- NULL
            n_steps1 <- NULL
            transform_state <- NULL
        }
        out[[i]] <- simulate(covid_multi_strain, p_i, n_steps, 
                             deterministic, keep_all_states = F,
                             min_ages = min_ages, Rt = Rt,
                             p1 = p1_i, n_steps1 = n_steps1, 
                             transform = transform_state)
    }
    
    info <- out[[1]]$info
    states_cntfctl <- lapply(out,"[[",1)
    
    # Extract daily values
    states_cntfctl <- lapply(states_cntfctl,function(y) y[,,seq(1,nlayer(y),by = 1/dt),drop = F])
    
    states_cntfctl <- abind(states_cntfctl,along = 2)
    
    return(list(states_cntfctl = states_cntfctl,states = states,smpl = smpl,info = info))
}

# calculate_quantiles <- function(x, probs = c(0.025,0.5,0.975)){
#     q <- apply(x,1,function(y) quantile(y, probs = probs))
#     return(q)
# }

calculate_outcome_quantiles <- function(x,info){
    # If input is a 3-D array convert to a data table
    if (length(dim(x)) == 3){
        x <- arr_to_dt(x,info)
    }
    
    out <- x[,.(med = median(value),
                q95l = quantile(value,0.025),
                q95u = quantile(value,0.975)), by = .(state,day)]
    return(out)
}

calculate_outcomes_by_wave <- function(x,wave_date,info,min_ages = seq(0,70,by = 10),Rt = FALSE){
    cols <- names(index(info,min_ages = min_ages,Rt = Rt)$state)
    dimnames(x)[[1]] <- cols
    n_waves <- length(wave_date) - 1
    
    tmp <- vector("list",n_waves)
    for (j in 1:n_waves){
        times <- wave_date[j]:wave_date[j+1]
        # Calculate total outcomes averted over wave
        tmp[[j]] <- as.data.table(t(apply(x[,,times], c(1,2), sum)))
        tmp[[j]][,smpl := .I]
    }
    total_x_waves <- rbindlist(tmp,idcol = "wave")
    total_x_both_waves <- total_x_waves[,lapply(.SD,sum),.SDcols = cols,by = .(smpl)]
    total_x_both_waves[,wave := "Total"]
    total_x <- rbind(total_x_waves,total_x_both_waves)
    
    # Calculate quantiles
    q_total_x <- total_x[,unlist(
        lapply(.SD,function(x) list(q95l = quantile(x,probs = 0.025),
                                    med = quantile(x,probs = 0.5),
                                    q95u = quantile(x,probs = 0.975))),
        recursive = F),.SDcols = cols,by = .(wave)]
    
    return(q_total_x)
}


plot_counterfactuals <- function(q_outcomes,q_outcomes_cntfctl,outcome,ylbl,ttls){
    clrs = c("Fitted" = "black","Counterfactual" = "darkgreen")
    p <- ggplot() + 
        geom_line(aes(x = date,y = med,color = "Fitted"),q_outcomes[state == outcome]) +
        geom_ribbon(aes(x = date,ymin = q95l,ymax = q95u,fill = "Fitted"),q_outcomes[state == outcome],alpha = 0.5) +
        geom_line(aes(x = date,y = med,color = "Counterfactual"),q_outcomes_cntfctl[state == outcome]) + 
        geom_ribbon(aes(x = date,ymin = q95l,ymax = q95u,fill = "Counterfactual"),q_outcomes_cntfctl[state == outcome],alpha = 0.5) +
        labs(x = "Day",y = ylbl) +
        scale_color_manual(name = "",values = clrs) +
        scale_fill_manual(name = "",values = clrs) +
        facet_wrap(~cntfctl,nrow = 2,labeller = labeller(cntfctl = ttls),dir = "v") +
        theme(legend.position = "bottom")
    return(p)
}


simulate_prepare <- function(onward,n_smpls,seed = 1){
    info <- onward$info
    
    # Take random sample of posterior parameters without replacement
    n_smpls_mcmc <- nrow(onward$pars)
    if (n_smpls > n_smpls_mcmc){
        message(sprintf(
            "Reducing n_smpls from %d to %d as too few available in MCMC output",
            n_smpls, n_smpls_mcmc))
    }
    set.seed(seed)
    idx <- sort(sample.int(n_smpls_mcmc, n_smpls, replace = F))
    pars_mcmc <- onward$pars[idx,,drop = F]
    state <- onward$state[,idx,drop = F]
    
    # Apply transform function to matrix of MCMC parameters to get list of full parameter sets
    pars <- apply(pars_mcmc, 1, onward$transform)
    
    ret <- onward[c("step","date","dt","vaccine")]
    ret$pars <- pars
    ret$state <- state
    ret$info <- info
    ret$idx <- idx
    ret
}


simulate_future_scenario <- function(args, onward, new_strain = FALSE){
    multistrain <- onward$info$multistrain
    if (multistrain){
        n_strain <- 4
    } else {
        n_strain <- 1
    }
    
    step_start <- onward$step
    print(step_start)
    date_start <- sircovid_date(onward$date)
    print(date_start)
    end_date <- sircovid_date(args$end_date)
    print(end_date)
    dates <- seq(date_start,end_date)
    print(dates)
    steps <- as.integer(dates/onward$dt)
    step_end <- last(steps)
    
    info <- onward$info$info
    index <- simulate_index(info, args$output_keep,
                            args$output_vaccination,
                            multistrain)
    
    # Matrix of starting values for model states
    state_start <- onward$state
    
    # Set up future vaccination schedule
    pars <- simulate_pars_vaccination(onward$pars[[1]]$N_tot, args, onward, n_strain)
    
    # pars <- unlist(pars, F, F)
    
    pars <- setup_future_betas(pars, step_start, step_end, dt)
    
    if (!is.null(args$strain_transmission) && new_strain){
        state_start <- rotate_strains(state_start, info)
    }
    
    # Create instance of model
    mod <- covid_multi_strain$new(pars,step_start,n_particles = NULL,pars_multi = T,
                                  n_threads = args$n_threads,seed = args$seed)
    # Set starting values for model states for simulations
    mod$update_state(state = state_start)
    mod$set_index(index$run)
    # Simulate
    state <- mod$simulate(steps)
    
    ret <- list(date = dates)
    
    if (args$output_time_series){
       ret$state <- state[args$output_keep,,] 
    }
    
    ret
}


simulate_pars_vaccination <- function(pop, args, onward, n_strain) {
    priority_population <- vaccine_priority_population(
        pop,
        uptake = args$vaccine_uptake * args$vaccine_eligibility)

    pars <- onward$pars
    vaccine <- onward$vaccine

    vaccine_index_dose2 <- pars[[1]]$index_dose[2]
    vaccine_progression_rate <- pars[[1]]$vaccine_progression_rate_base
    N_tot <- pars[[1]]$N_tot

    ## TODO: in the validation, if booster doses is non-empty, we should
    ## check that we have a model with boosters
    if (!is.null(args$vaccine_booster_daily_doses)) {
        # args$vaccine_efficacy <-
        #     Map(cbind, args$vaccine_efficacy, args$vaccine_booster_efficacy)
        # args$strain_vaccine_efficacy <-
        #     Map(cbind, args$strain_vaccine_efficacy,
        #         args$strain_vaccine_booster_efficacy)
        ## TODO: index_dose for boosters is assumed to be fourth column, when in
        ## fact it's now the fifth!! This needs fixing ASAP
        vaccine_index_booster <- pars[[1]]$index_dose[3]
    } else {
        vaccine_index_booster <- NULL
    }

    mean_days_between_doses <- round(vaccine$mean_days_between_doses *
                                         args$vaccine_delay_multiplier)

    ## TODO: potentially a big problem here! We cannot fully disentangle booster
    ## eligibility from general vaccine eligibility, cna we? Especially, it seems
    ## the wya the latter is inputed (i.e. eligibility * uptake = calculation of
    ## priority population) seems to heavily affect actual boosters uptake! So
    ## we are potentially overestimating by quite a lot the effect of a
    ## boster vaccination programme!
    vaccine_schedule <- vaccine_schedule_scenario(
        ## schedule_past got renamed to schedule_real in new params task
        schedule_past = vaccine$schedule,
        doses_future = args$vaccine_daily_doses,
        end_date = args$end_date,
        mean_days_between_doses = mean_days_between_doses,
        priority_population = priority_population,
        lag_groups = args$vaccine_lag_groups,
        lag_days = args$vaccine_lag_days,
        boosters_future = args$vaccine_booster_daily_doses,
        boosters_prepend_zero = TRUE,
        booster_proportion = args$vaccine_booster_eligibility)

    ## check boosters
    # > par(mfrow = c(3, 1))
    # > image(t(vaccine_schedule$doses[, 1, ]))
    # > image(t(vaccine_schedule$doses[, 2, ]))
    # > image(t(vaccine_schedule$doses[, 3, ]))

    ## TODO: potential placeholder for manually setting
    ## vaccine_schedule$doses[, 3, ] to zero if we don't want to give boosters
    ## to all age groups
    ## but probably better done inside sircovid::vaccine_schedule_scenario
    # > vaccine_schedule$doses[1:17, 3, 0] <- 0

    # rel_list <- modify_severity(
    #     args$vaccine_efficacy,
    #     args$strain_vaccine_efficacy,
    #     args$strain_severity_modifier)

    extra <- vaccination_parameters(
        N_tot,
        pars[[1]]$dt,
        # rel_susceptibility = rel_list$rel_susceptibility,
        # rel_p_sympt = rel_list$rel_p_sympt,
        # rel_p_hosp_if_sympt = rel_list$rel_p_hosp_if_sympt,
        # rel_p_death = rel_list$rel_p_death,
        # rel_infectivity = rel_list$rel_infectivity,
        rel_susceptibility = pars[[1]]$rel_susceptibility,
        rel_p_sympt = pars[[1]]$rel_p_sympt,
        rel_p_hosp_if_sympt = pars[[1]]$rel_p_hosp_if_sympt,
        rel_p_death = pars[[1]]$rel_p_death,
        rel_infectivity = pars[[1]]$rel_infectivity,
        vaccine_schedule = vaccine_schedule,
        vaccine_index_dose2 = vaccine_index_dose2,
        vaccine_index_booster = vaccine_index_booster,
        vaccine_progression_rate = vaccine_progression_rate,
        n_strains = n_strain,
        n_doses = pars[[1]]$n_doses)
    
    n_group <- nrow(pars[[1]]$rel_p_sympt)
    n_vacc_class <- extra$n_vacc_classes

    # rel_severity <- build_rel_param(extra$rel_p_death, n_strain,
    #                                 n_vacc_class, "rel_p_death")
    # 
    # extra$rel_p_ICU <- extra$rel_p_R <-
    #     array(1, c(n_group, n_strain, n_vacc_class))
    # 
    # extra$rel_p_ICU_D <- extra$rel_p_H_D <- extra$rel_p_W_D <-
    #     extra$rel_p_G_D <- rel_severity

    if (!is.null(args$strain_transmission)) {
        strain_params <- parameters_strain(
            args$strain_transmission,
            sircovid_date(args$strain_seed_date),
            args$strain_seed_size,
            args$strain_seed_pattern,
            pars[[1]]$dt)
        strain_params$cross_immunity <- args$strain_cross_immunity
        strain_params$waning_rate <- rep(args$waning_rate, 8)
        extra <- c(extra, strain_params)
    }

    for (i in seq_along(pars)) {
        pars[[i]][names(extra)] <- extra
    }
    
    pars
}


simulate_index <- function(info, keep, calculate_vaccination, multistrain) {
    
    index_S <- info$index$S
    names(index_S) <- paste0("S_", seq_along(index_S))
    
    index_D <- info$index$D
    names(index_D) <- paste0("D_", seq_along(index_D))
    
    # index_I <- info$index$cum_infections_disag
    # names(index_I) <- paste0("I_", seq_along(index_I))
    # 
    # index_A <- info$index$diagnoses_admitted
    # names(index_A) <- paste0("A_", seq_along(index_A))
    
    # if (calculate_vaccination || multistrain) {
    #     index_n_vaccinated <- set_names(
    #         info$index$cum_n_vaccinated,
    #         paste0("n_vaccinated_", seq_along(index_S))
    #     )
    #     index_R <- info$index$R
    #     names(index_R) <- paste0("R_", seq_along(index_R))
    # } else {
        index_n_vaccinated <- NULL
        index_R <- NULL
    # }
    
    if (multistrain) {
        index_prob_strain <- info$index$prob_strain
        names(index_prob_strain) <- paste0(
            "prob_strain_", seq_along(index_prob_strain))
    } else {
        index_prob_strain <- NULL
    }
    
    index <- c(index(info)$state[keep],
               index_S, index_D, #index_I, index_A, 
               index_n_vaccinated,
               index_R, index_prob_strain)
    
    list(run = index,
         n_vaccinated = index_n_vaccinated,
         S = index_S,
         # I = index_I,
         # A = index_A,
         R = index_R,
         D = index_D,
         prob_strain = index_prob_strain)
}


# FOR NOW: This just repeats the estimated most recent value of beta up to the 
# end of the simulation
# TODO: Implement with Rt?
# setup_future_betas <- function(pars, rt_future, S, rt_type,
#                                step_current, step_end, dt, seasonality, 
#                                R, prob_strain) {
setup_future_betas <- function(pars, step_current, step_end, dt) {
    # ## For seasonality, we assume a sine wave with a trough at day 228 = 15th Aug
    # ## (and a peak 6 months earlier on day 46 = 15th Feb)
    # seasonality_date_peak <- sircovid::sircovid_date("2020-02-15")
    
    ## Past betas, as inferred from the pmcmc
    beta <- t(vapply(pars, "[[", numeric(length(pars[[1]]$beta_step)),
                     "beta_step"))
    
    n_beta_add <- step_end - ncol(beta)
    beta <- cbind(beta, matrix(rep(beta[, ncol(beta)], n_beta_add), nrow(beta)))
    # beta <- array(beta, c(dim(pars), ncol(beta)))
    
    # rt <- vnapply(seq_along(pars), function(i)
    #     sircovid::lancelot_Rt(step_current, S[, i, drop = FALSE], pars[[i]],
    #                           type = rt_type, R = R[, i, drop = FALSE],
    #                           prob_strain = prob_strain[, i, drop = FALSE],
    #                           weight_Rt = FALSE)[[rt_type]][1])
    # beta_rt_ratio <- beta[, , step_current] / rt
    # 
    # for (region_index in seq_len(ncol(pars))) {
    #     r <- colnames(pars)[[region_index]]
    #     rt_future_r <- rt_future[rt_future$region == r, ]
    #     rt_future_r$step_start <- sircovid::sircovid_date(rt_future_r$date) / dt
    #     rt_future_r$step_end <- c(rt_future_r$step_start[-1L] - 1L, step_end)
    #     
    #     for (i in seq_rows(rt_future_r)) {
    #         j <- seq(rt_future_r$step_start[[i]], rt_future_r$step_end[[i]])
    #         
    #         if (rt_future_r$Rt_sd[[i]] > 0) {
    #             dpars <- data.frame(mean = rt_future_r$Rt[[i]],
    #                                 sd = rt_future_r$Rt_sd[[i]])
    #             
    #             # order by increasing Rt to preserve direction of Rt scaling
    #             rt_i <- sort(distr6::dstr("Lognormal",
    #                                       mean = rt_future_r$Rt[[i]],
    #                                       sd = rt_future_r$Rt_sd[[i]])$rand(nrow(beta)))
    #         } else {
    #             rt_i <- rt_future_r$Rt[[i]]
    #         }
    #         
    #         beta[, region_index, j] <- rt_i * beta_rt_ratio[, region_index]
    #     }
    # }
    # 
    # beta <- mcstate::array_flatten(beta, 1:2)
    # 
    # ## Apply seasonality to the new betas:
    # 
    # ## Back-calculate dates here to deal with any truncating of
    # ## step_current:step_end caused by simulating with new data.
    # date <- seq(to = step_end, length.out = n_beta_add) * dt
    # ## Relative distance between us and the peak date (on [0..1])
    # beta_mult <- spim_calc_seasonality(date, seasonality_date_peak, seasonality)
    # 
    # i <- seq(to = ncol(beta), length.out = n_beta_add)
    # beta[, i] <- beta[, i] * rep(beta_mult, each = nrow(beta))

    for (i in seq_along(pars)) {
        pars[[i]]$beta_step <- beta[i, ]
    }
    
    pars
}


arr_to_dt <- function(x,info,min_ages = seq(0,70,by = 10),Rt = FALSE){
    x_long <- as.data.table(x)
    
    names(x_long)[1:3] <- c("state","smpl","day")
    
    x_long[,state := names(index(info,min_ages = min_ages,Rt = Rt)$state[state])]
    
    return(x_long)
}


med_and_CI = function(x,l,u,f=1,d=1,method="round"){
    if (method=="signif"){
        paste0(signif(f*x,d)," (",signif(f*l,d),"-",signif(f*u,d),")")
    } else if (method=="round"){
        paste0(round(f*x,d)," (",round(f*l,d),"-",round(f*u,d),")")
    }
}
