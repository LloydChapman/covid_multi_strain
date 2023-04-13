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
    
    # Plot p.w. constant beta to check it looks right
    dt <- base$dt
    beta_date <- base$beta_date
    beta_t <- seq(0, beta_date[length(beta_date)], by = dt)
    pdf(paste0("output/plots",run,".pdf")) # save all plots into one pdf
    par(mfrow = c(1,1))
    
    n_smpls <- round(n_iters/thinning)
    beta_value_post <- apply(res$pars[seq(round(n_smpls/10),n_smpls,by = min(round(n_smpls/10),10)),seq_along(beta_date)],2,median)
    if (base$beta_type == "piecewise-linear") {
        beta_step <- parameters_piecewise_linear(beta_date, 
                                                 beta_value_post %||% 0.1, dt)
    } else if (base$beta_type == "piecewise-constant") {
        beta_step <- parameters_piecewise_constant(beta_date, 
                                                   beta_value_post %||% 0.1, dt)
    }
    plot(covid_multi_strain_date_as_date(beta_t), beta_step, type = "o", cex = 0.25, xaxt = "n", xlab = "Date", ylab = "beta")
    dates_plot <- seq.Date(covid_multi_strain_date_as_date(beta_t[1]),covid_multi_strain_date_as_date(beta_t[length(beta_t)]),by = 30)
    axis(1, dates_plot, format(dates_plot,"%Y-%m-%d"))
    
    # Trace plots
    par(mar = c(2.5,4,2,1), mfrow = c(ceiling(length(u)/2),2))
    for (i in u){
        plot(res$pars[,i],type = "l",xlab = "Iteration",ylab = names(pars_init)[i])
    }
    
    # Scaling factor
    par(mfrow = c(1,1))
    plot(res$scaling_factor,type="l")
    
    par(mfrow = c(1,1))
    # pairs(res$pars[seq(round(n_smpls/10),n_smpls,by=10), ],labels = names(pars_init))
    pairs(res$pars[seq(round(n_smpls/10),n_smpls,by = min(round(n_smpls/10),10)),u],labels = names(pars_init)[u])
    
    # Plot fitted hospitalisations and deaths against data
    print(plot_outcome(res$trajectories$state,data,"hosps",by_age = T))
    print(plot_outcome(res$trajectories$state,data,"deaths",by_age = T))
    print(plot_outcome(res$trajectories$state,data,"cases",res$pars[,"phi_cases"],by_age = T))
    print(plot_sero(res$trajectories$state,data,pop[,.(population = sum(total)),by = .(age_group)],by_age = T))
    print(plot_outcome(res$trajectories$state,data,"hosps"))
    print(plot_outcome(res$trajectories$state,data,"deaths"))
    print(plot_outcome(res$trajectories$state,data,"cases",res$pars[,"phi_cases"]))
    
    dev.off()
    
}