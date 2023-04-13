fit_covid_multi_strain <- function(data_raw,pars,u,n_iters,run,deterministic = TRUE,Rt = FALSE,thinning = 1){
    #### Set up model ####
    
    if (n_iters < 100){
        stop("n_iters must be at least 100")
    }
    
    # Compile model
    covid_multi_strain <- odin_dust("inst/odin/covid_multi_strain.R")
    
    # Construct parameters object for 1st epoch
    base <- pars$base
    init <- pars$mcmc$initial()
    p <- parameters(base$dt,
                    base$n_age,
                    base$n_vax,
                    base$m,
                    base$beta_date,
                    beta_value = init[grep("beta",pars$mcmc$names())],
                    base$beta_type,
                    base$gamma_E,
                    base$gamma_P,
                    base$gamma_A,
                    base$gamma_C,
                    base$gamma_H,
                    base$gamma_G,
                    base$gamma_pre_1,
                    base$gamma_P_1,
                    base$theta_A,
                    base$p_C,
                    p_H = init["p_H_max"]*base$p_H,
                    base$p_G,
                    p_D = init["p_D_max"]*base$p_D,
                    base$p_P_1,
                    base$population,
                    start_date = init["start_date"],
                    base$initial_seed_size,
                    base$initial_seed_pattern,
                    strain_transmission = c(1,init["rel_strain_transmission"]),
                    strain_seed_date = init["strain_seed_date"],
                    base$strain_seed_size,
                    base$strain_seed_pattern,
                    base$strain_rel_p_sympt,
                    base$strain_rel_p_hosp_if_sympt,
                    base$strain_rel_p_death,
                    base$rel_susceptibility,
                    base$rel_p_sympt,
                    base$rel_p_hosp_if_sympt,
                    base$rel_p_death,
                    base$rel_infectivity,
                    base$vaccine_progression_rate,
                    base$vaccine_schedule,
                    base$vaccine_index_dose2,
                    base$vaccine_index_booster,
                    base$vaccine_catchup_fraction,
                    base$n_doses,
                    base$vacc_skip_progression_rate,
                    base$vacc_skip_to,
                    base$vacc_skip_weight,
                    base$waning_rate,
                    base$cross_immunity,
                    phi_cases = init["phi_cases"],
                    kappa_cases = 1/init["alpha_cases"],
                    base$sero_sensitivity_1,
                    base$sero_specificity_1)
    
    mod <- covid_multi_strain$new(p,step = 0,n_particles = 1,seed = 1L,
                                  deterministic = deterministic)
    
    info <- mod$info()
    min_ages <- get_min_age(base$age_groups)
    
    
    #### Fit to multiple age-stratified data streams ####
    
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
    
    # Extract objects required for MCMC from pars
    transform <- pars$transform
    pars_info <- pars$info
    pars_init <- pars$mcmc$initial()
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
                u,thinning)
    tend <- Sys.time()
    print(tend - tstart)
    
    # Save output
    save(list = ls(all.names = T), file = paste0("output/MCMCoutput",run,".RData"), envir = environment())
    
}