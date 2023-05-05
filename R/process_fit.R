# Process MCMC output from fit_covid_multi_strain()
process_fit <- function(samples,pars,data_raw,burnin = NULL, simulate_object = TRUE){
    
    covid_multi_strain <- odin_dust("inst/odin/covid_multi_strain.R")
    info <- covid_multi_strain$new(pars$transform(samples$pars[1,])[[1]]$pars,step = 0,
                                   n_particles = 1)$info() #filter$model$new(samples$pars[1,],0,1)$info()
    
    # Extract baseline parameters
    base <- pars$base
    
    # Convert raw data to required format for particle filter
    data <- particle_filter_data(data_raw,"day",1/base$dt)
    
    # Create particle filter object
    if (deterministic){
        n_particles <- 1
        filter <- particle_deterministic$new(
            data, covid_multi_strain, compare, 
            index = function(info) 
                index(info, min_ages = min_ages, Rt = Rt), 
            initial = initial)
    } else {
        n_particles <- 200
        filter <- particle_filter$new(
            data, covid_multi_strain, n_particles, compare, 
            index = function(info) 
                index(info, min_ages = min_ages, Rt = Rt), 
            initial = initial
        )
    }
    
    samples$predict <- list(transform = transform,
                            index = index(info)$state,
                            rate = as.integer(1/base$dt),
                            filter = filter)
    
    samples$info <- list(info = info,
                         date = covid_multi_strain_date_as_date(dim(samples$trajectories$state)[3] - 1), # SORT THIS OUT, E.G. ADD end_date TO BASELINE PARS
                         multistrain = info$dim$prob_strain > 1,
                         beta_date = base$beta_date,
                         pars = pars$info$name)
    
    if (is.null(burnin)){
        burnin <- round((n_iters/thinning)/10)
    }
    # # Remove burn-in
    # idx_drop <- -(1:(burnin+1))
    # samples$pars <- samples$pars[idx_drop, ]
    # samples$probabilities <- samples$probabilities[idx_drop,]
    # samples$state <- samples$state[,idx_drop]
    # samples$trajectories$state <- samples$trajectories$state[,idx_drop,]
    
    samples$trajectories$date <- samples$trajectories$step/samples$trajectories$rate
    
    # samples$vaccine <- list(mean_days_between_doses = mean_days_between_doses,
    #                         schedule = schedule)
    
    if (simulate_object){
        start_date_sim <- "2021-11-21"
        simulate <- create_simulate_object(samples, start_date_sim, samples$info$date)
    }
    
    list(samples = samples,
         simulate = simulate,
         data = data)
}


create_simulate_object <- function(samples,start_date_sim,date){
    start_date_sim <- covid_multi_strain_date(start_date_sim)
    fit_dates <- samples$trajectories$date
    idx_dates <- (fit_dates >= start_date_sim) & 
        (fit_dates <= covid_multi_strain_date(date))
    date <- fit_dates[idx_dates]
    
    state_keep <- c("hosps", "deaths")
    state_full <- samples$trajectories$state
    
    state <- state_full[state_keep, ,idx_dates]
    
    list(date = date,state = state)
}


calculate_parameter_quantiles <- function(dat, digits = 3, burnin = NULL){
    if (is.null(burnin)){
        burnin <- round(nrow(dat$samples$pars)/10)
    }
    # Remove burn-in
    idx_drop <- -(1:(burnin+1))
    pars <- dat$samples$pars[idx_drop, ]
    
    tmp <- signif(apply(pars, 2, calculate_median_and_ci),digits = digits)
    tmp[,c("start_date","strain_seed_date","strain_seed_date1")] <- 
        apply(tmp[,c("start_date","strain_seed_date","strain_seed_date1")],2,
              function(x) as.character(covid_multi_strain_date_as_date(x)))
    tmp <- as.data.table(t(tmp),keep.rownames = T)
    
    data.table(Parameter = tmp[,rn],
               `Median (95% CI)` = sprintf("%s (%s, %s)", tmp[,V1], tmp[,`2.5%`], tmp[,`97.5%`]))
}


calculate_median_and_ci <- function(x){
    c(median(x), quantile(x,c(0.025, 0.975)))
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
         step = covid_multi_strain_date(date)/dt,
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
