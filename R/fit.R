covid_multi_strain_particle_filter <- function(data, pars, deterministic = TRUE, Rt = FALSE, initial_date = 0, n_particles = 200){
    base <- pars$base
    min_ages <- get_min_age(base$age_groups)
    initial_date <- as_covid_multi_strain_date(initial_date)
    
    # Compile model
    covid_multi_strain <- odin_dust("inst/odin/covid_multi_strain.R")
    
    # Convert raw data to required format for particle filter
    data <- particle_filter_data(data,"day",1/base$dt,initial_date)
    
    # Create particle filter object
    if (deterministic){
        filter <- particle_deterministic$new(
            data, covid_multi_strain, compare, 
            index = function(info) 
                index(info, min_ages = min_ages, Rt = Rt), 
            initial = initial)
    } else {
        filter <- particle_filter$new(
            data, covid_multi_strain, n_particles, compare, 
            index = function(info) 
                index(info, min_ages = min_ages, Rt = Rt), 
            initial = initial
        )
    }
    return(filter)
}


fit_run <- function(pars,filter,u,n_iters,deterministic = TRUE,fixed_initial = T,Rt = FALSE,thinning = 1,rerun_every = Inf){
    if (n_iters < 100){
        stop("n_iters must be at least 100")
    }
    
    # Construct parameters object for 1st epoch
    p <- pars$mcmc$model(pars$mcmc$initial())
    # Get model info
    info <- filter$model$new(p[[1]]$pars,step = 0,n_particles = 1)$info()
    
    min_ages <- get_min_age(pars$base$age_groups)
    
    
    #### Fit to multiple age-stratified data streams ####
    
    if (deterministic){
        n_particles <- 1
    } else {
        n_particles <- filter$n_particles
    }
    
    # Extract objects required for MCMC from pars
    transform <- pars$mcmc$model
    pars_init <- pars$mcmc$initial()
    pars_info <- pars$info
    if (fixed_initial){
        pars_init[u] <- pars$mcmc$propose(pars$mcmc$initial(),1000)[u]
    } else {
        pars_init[u] <- mapply(sample_prior,
                            split(pars$prior,pars$prior$name)[unique(pars$prior$name)],
                            split(pars_info,pars_info$name)[unique(pars_info$name)],
                            SIMPLIFY = T)[u]
    }
    priors <- lapply(pars$mcmc$.__enclos_env__$private$parameters,"[[","prior")
    pars_min <- pars_info$min
    pars_max <- pars_info$max
    proposal <- pars$proposal
    scaling_factor_start <- 1
    iter0 <- 100
    discrete <- pars_info$integer
    
    # Run MCMC
    tstart <- Sys.time()
    res <- mcmc(transform,filter,pars_init,priors,n_particles,n_iters,index(info,min_ages,Rt),
                scaling_factor_start,proposal,pars_min,pars_max,iter0,discrete,
                u,thinning,rerun_every)
    tend <- Sys.time()
    print(tend - tstart)
    
    return(res)
    
}