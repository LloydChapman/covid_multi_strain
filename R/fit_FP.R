run_fitting <- function(run,assumptions,u,n_iters,n_chains){
    # Load data
    data_raw <- fread("data/data_cases_hosps_deaths_serology.csv")
    # data_raw <- data_raw[day < covid_multi_strain_date(as.Date("2021-06-01")),]
    # cols <- setdiff(names(data_raw),"day")
    # data_raw <- data_raw[!(day <= 213 | day >= covid_multi_strain_date(as.Date("2020-09-15"))),(cols):= NA]
    vax <- fread("data/data_vaccination.csv", colClasses = c(number = "numeric"))
    pop <- fread("data/population.csv")
    
    # Set assumption for booster waning rate
    # assumptions <- "central" #-log(67.7/82.8)/(105-25) # (Stowe Nat Comm 2022 Table S11)
    # assumptions <- "pessimistic"
    # assumptions <- "optimistic" # -log(0.923)/140 (Barnard Nat Com 2022 Table S4)
    
    ## Load parameters
    # Output pars is a list containing:
    # info - the loaded info.csv
    # prior - the loaded prior.csv
    # proposal - the loaded proposal.csv
    # transform - the functions in transform.R
    # raw - exact output of first 3 csvs, without any treatment
    # base - the base.rds object, which contains the fixed parameters
    # mcmc - initialisation object built from the above to pass to the mcmc
    pars <- fit_pars_load("parameters",assumptions)
    
    # Plot vaccination coverage by age
    # vaccination_coverage_plot(pars$base$vaccine_schedule,pars$base$age_groups,vax,pop)
    dir.create("output")
    # ggsave("output/vax_cov_by_dose.pdf",width = 9,height = 3)
    
    # Fit covid_multi_strain to FP data
    # u <- c(1:5,7:9,10:12,14,15:19) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, reporting rate for confirmed cases
    # # u <- c(1:6,8:10,11:13,15,16,17:21) # beta parameters, seed date, strain seed date, IHR scaling, 2nd strain seed date, reporting rate for confirmed cases
    # # u <- c(1:3,8,10,13,16:19)
    # # u <- c(2:3,8,10,13,16:19)
    # n_iters <- 100 #5e4 #3e4 #1e4 #3e4 #
    # # Change run number for different assumption on booster waning rate
    # # run <- 77
    # # run <- 78
    # # run <- 130 #118 #117 #
    # n_chains <- 4 #2 #1 #
    deterministic <- T # flag for whether to use deterministic model or not
    fixed_initial <- T #F # flag for whether to use fixed initial values for MCMC chains or not
    Rt <- T #F # flag for whether to return variables needed for calculating Rt in "state" object
    initial_date <- pars$info$min[pars$info$name == "start_date"] - 1
    
    # Construct particle filter
    filter <- covid_multi_strain_particle_filter(data_raw,pars,deterministic,Rt,initial_date)
    
    # Run fitting
    thinning <- 10
    results <- as.list(paste0("output/MCMCsamples",run,"_",seq_len(n_chains),".RDS"))
    # for (i in seq_len(n_chains)) {
    #     set.seed(i)
    #     tmp <- fit_run(pars,filter,u,n_iters,deterministic,fixed_initial,Rt,thinning)
    #     saveRDS(tmp,results[[i]])
    #     plot_traces(tmp$pars,u)
    #     ggsave(paste0("output/par_traces",run,"_",i,".pdf"),width = 6,height = 6)
    #     rm(tmp)
    #     gc()
    # }
    
    tstart <- Sys.time()
    mclapply(seq_len(n_chains), function(i){
        set.seed(i)
        tmp <- fit_run(pars,filter,u,n_iters,deterministic,fixed_initial,Rt,thinning)
        saveRDS(tmp,results[[i]])
        plot_traces(tmp$pars,u)
        ggsave(paste0("output/par_traces",run,"_",i,".pdf"),width = 6,height = 6)
        rm(tmp)
        gc()
    }, mc.cores = 4L)
    tend <- Sys.time()
    print(tend - tstart)
    
    ## Post processing
    # Set burn-in
    burnin <- 4000 #500 #3000 #2000 #
    
    # Get results
    samples_list <- vector("list",n_chains)
    pars_list <- vector("list",n_chains)
    for (i in seq_len(n_chains)){
        tmp <- readRDS(results[[i]])
        pars_list[[i]] <- as.data.table(tmp$pars[-1,u])
        pars_list[[i]][,iter := seq_len(nrow(tmp$pars)-1)]
        tmp <- remove_burnin(tmp,burnin)
        samples_list[[i]] <- tmp
        rm(tmp)
        gc()
    }
    
    # samples <- lapply(results, readRDS)
    # 
    # # Remove burn-in
    # for (i in seq_len(n_chains)){
    #     samples[[i]] <- remove_burnin(samples[[i]],burnin)
    # }
    # samples <- lapply(samples, function(x) remove_burnin(x,burnin))
    
    # Combine chains
    samples <- chains_combine(samples = samples_list)
    
    # Process MCMC output
    dat <- fit_process(samples,pars,filter,simulate_object = T)
    rm(samples_list,samples)
    gc()
    
    # Save output
    saveRDS(dat,paste0("output/MCMCoutput",run,".RDS"))
    
    # Set whether to plot moving average of data
    moving_avg <- F
    
    # Set whether to plot prediction interval from model
    pred_intvl <- T
    
    # Set number of posterior samples for age-decomposition plots
    n_smpls <- 1000 #500 #
    
    # Plot output
    lbls <- c(paste0("$\\beta_",seq_along(pars$base$beta_date),"$"),
              "$\\sigma_{D elta}$","$t_0$","$t_{D elta}$","${p_H}_{max}$",
              "${p_D}_{max,1}$","${p_D}_{max,2}$","${p_D}_{max,3}$",
              "$\\sigma_{O micron}$","$t_{O micron}$",
              "${\\pi_H}_{D elta/Wildtype}$","$\\phi_{cases}$",
              "$\\alpha_{cases}$","$\\alpha_{hosp}$","$\\alpha_{death}$")
    
    # Combine parameter samples
    pars_samples <- rbindlist(pars_list,idcol = "chain")
    
    # Plot traces from all chains together
    plot_traces_all(pars_samples,u,lbls)
    ggsave(paste0("output/par_traces",run,".pdf"),width = 8,height = 9)
    
    # Calculate Gelman-Rubin diagnostic for convergence
    mcmc_list <- mcmc.list(lapply(pars_list,function(x) coda::mcmc(x[(burnin+1):nrow(x),!"iter"],thin = thinning)))
    gr_diag <- gelman.diag(mcmc_list,autoburnin = F)
    gr_diag$max <- max(gr_diag$psrf[,"Point est."])
    
    # Calculate effective sample size
    ESS <- effectiveSize(mcmc_list)
    ESS <- list(ESS = ESS, min = min(ESS))
    
    # Save convergence diagnostics
    cnvgnce_diag <- list(gr_diag = gr_diag,ESS = ESS)
    saveRDS(cnvgnce_diag,paste0("output/cnvgnce_diag",run,".RDS"))
    
    # Plot model fit
    plot_fit(dat,pars,run,pop,u,lbls,moving_avg,pred_intvl,n_smpls)
    
    # Process fit
    pars_qntls <- calculate_parameter_quantiles(dat)
    write.csv(pars_qntls,paste0("output/parameter_quantiles",run,".csv"), row.names = F)
}
