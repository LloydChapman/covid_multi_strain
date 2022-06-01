# Process MCMC output from fit_covid_multi_strain()
process_fit <- function(samples,init_pars,idx,n_steps,dt,transform,index,filter,date,
                        beta_date,mean_days_between_doses,schedule,data,
                        burnin = NULL, simulate_object = TRUE){
    # # Temporary: add names to parameters matrix from MCMC
    # colnames(samples$pars) <- names(init_pars)
    # # TODO: add parameter names in MCMC code
    # # Temporary: add dimnames to state array
    # dimnames(samples$trajectories$state) <- list(names(idx$state))
    # # TODO: add dimnames in MCMC code
    # 
    # samples <- list()
    # samples$pars <- res$pars
    # samples$probabilities <- cbind(log_prior = res$lprior,
    #                                log_likelihood = res$ll,
    #                                log_posterior = res$lpost)
    # samples$state <- rbind(rep(n_steps * dt,ncol(res$states)),res$states[,,nlayer(res$states)])
    # samples$trajectories <- list(step = as.integer((0:(nlayer(res$states)-1))/dt),
    #                              rate = as.integer(1/dt),
    #                              state = res$states,
    #                              predicted = rep(F,nlayer(res$states)))
    
    info <- covid_multi_strain$new(transform(samples$pars[1,]),step = 0,
                                   n_particles = 1)$info() #filter$model$new(samples$pars[1,],0,1)$info()
    
    samples$predict <- list(transform = transform,
                            index = index(info)$state,
                            rate = as.integer(1/dt),
                            filter = filter)
    
    samples$info <- list(info = info,
                         date = date,
                         multistrain = info$dim$prob_strain > 1,
                         beta_date = beta_date,
                         pars = names(init_pars))
    
    if (is.null(burnin)){
        burnin <- round((n_iters/thinning)/10)
    }
    # Remove burn-in
    idx_drop <- -(1:(burnin+1))
    samples$pars <- samples$pars[idx_drop, ]
    samples$probabilities <- samples$probabilities[idx_drop,]
    samples$state <- samples$state[,idx_drop]
    samples$trajectories$state <- samples$trajectories$state[,idx_drop,]
    
    samples$trajectories$date <- samples$trajectories$step/samples$trajectories$rate
    
    samples$vaccine <- list(mean_days_between_doses = mean_days_between_doses,
                            schedule = schedule)
    
    if (simulate_object){
        start_date_sim <- "2021-11-21"
        simulate <- create_simulate_object(samples, start_date_sim, samples$info$date)
    }
    
    list(samples = samples,
         simulate = simulate,
         data = data)
}


create_simulate_object <- function(samples,start_date_sim,date){
    start_date_sim <- sircovid_date(start_date_sim)
    fit_dates <- samples$trajectories$date
    idx_dates <- (fit_dates >= start_date_sim) & 
        (fit_dates <= sircovid_date(date))
    date <- fit_dates[idx_dates]
    
    state_keep <- c("hosps", "deaths")
    state_full <- samples$trajectories$state
    
    state <- state_full[state_keep, ,idx_dates]
    
    list(date = date,state = state)
}


load_fit <- function(filename){
    dat <- readRDS(filename)
    
    ret <- list()
    ret$info <- dat$samples$info[c("date","multistrain","beta_date")]
    ret$info$date <- as.Date(ret$info$date)
    
    ret$samples <- dat$samples
    
    ret$onward <- create_onward(ret)
    
    ret
}


create_onward <- function(dat){
    date <- dat$info$date
    dt <- 1/dat$samples$trajectories$rate
    list(date = date,
         step = sircovid_date(date)/dt,
         dt = dt,
         pars = dat$samples$pars,
         # base = dat$parameters$base,
         state = dat$samples$state,
         data = dat$data,
         transform = dat$samples$predict$transform,
         info = dat$samples$info,
         vaccine = dat$samples$vaccine,
         simulate = dat$simulate)
}
