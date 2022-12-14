fit_covid_multi_strain <- function(u,n_iters,run,deterministic = TRUE,Rt = FALSE,thinning = 1){
    #### Set up model and parameters ####
    
    
    # Age-structured SEIR model with deaths and vaccination
    # gen_seirhd_age_vax_multistrain_sero_time_dep_beta <- odin_dust("test_examples/seirhdagevaxmultistrainserotimedepbeta.R")
    covid_multi_strain <- odin_dust("inst/odin/covid_multi_strain.R")
    
    # Get age-dependent symptomatic fraction, IHR and IFR
    sympt_frac <- qread("data/2-linelist_both_fit_fIa0.5-rbzvih.qs")
    p_C <- unname(colMeans(sympt_frac[,13:20]))
    prob_death_given_hosp <- fread("data/prob_death_given_hosp_salje.csv")
    p_D <- prob_death_given_hosp[,median_perc_mean/100]
    IHR <- fread("data/ihr_salje.csv")
    ihr <- IHR[,median_perc_mean/100]
    IFR_salje <- fread("data/ifr_salje.csv")
    ifr_salje <- IFR_salje[,median_perc_mean/100]
    # O'Driscoll IFR
    IFR <- fread("data/ifr_odriscoll.csv")
    IFR[,age_low:=get_min_age(IFR$Age_group)]
    IFR[,age_group:=cut(age_low,c(min_ages,Inf),labels=age_groups,right=F)]
    agg_IFR <- IFR[,lapply(.SD,mean),.SDcols=setdiff(names(IFR),c("Age_group","age_low","age_group")),by=.(age_group)]
    ifr_odriscoll <- agg_IFR[,Median_perc_mean/100]
    # # For the time being just assume that IFR including non-hospital deaths is 10% higher
    # ifr <- 1.1*ifr_salje
    # p_G <- (ifr - p_D*ihr)/((1-p_D)*ihr + ifr)
    # Set probability of death outside hospital from observed proportion of non-hospital deaths
    p_G <- rep(weekly_dt[,sum(`Nombre de décès à domicile`,na.rm = T)/
                             sum(`Nombre total de nouvelles hospitalisations tous hôpitaux`,na.rm = T)],length(p_C))
    p_H <- ihr/(p_C*(1-p_G))
    p_P_1 <- 0.85
    # ifr <- p_C*p_H*(p_G + (1-p_G)*p_D)
    # # Population-weighted average IFR
    # sum(population*ifr)/sum(population) # 0.01009238
    
    # Get max values
    p_D_max0 <- max(p_D)
    p_H_max0 <- max(p_H)

    # Normalise
    p_D <- p_D/max(p_D)
    p_H <- p_H/max(p_H)
    
    # IFR <- readRDS("~/UCSF/COVIDVaccineModelling/Data/IFR_by_age_ODriscoll.RDS")
    # p_death <- IFR$median_perc[19:27]/100
    # p_death[length(p_death)-1] <- (p_death[length(p_death)-1]+p_death[length(p_death)])/2
    # p_death <- p_death[1:(length(p_death)-1)]
    
    vax_eff <- fread("data/vax_eff.csv")
    
    # # Set VE against Omicron infection in waned stratum to be the same as 2nd dose
    # vax_eff[outcome == "infection" & vaccine == "PF/MD" & dose == "waned",omicron := 44.1]
    
    # Melt to long format
    vax_eff_long <- melt(vax_eff,measure.vars = c("alpha","delta","omicron"),variable.name = "variant")
    
    # FOR NOW: Drop Omicron vax effectiveness
    vax_eff_long1 <- vax_eff_long[variant != "alpha"]
    vax_eff_long <- vax_eff_long[variant != "omicron"]
    
    # # Convert vaccine efficacy to proportion
    # vax_eff_long[,value := value/100]
    # # vax_eff_long[,class := fcase(dose == 1, 1,
    # #                              dose == 2, 2,
    # #                              dose == "Waned",3,
    # #                              dose == "Booster",4)]
    # # vax_eff_long[,strain := fcase(variant == "alpha", 1,
    # #                               variant == "delta", 2)]
    # # 
    # # # Cast to wide format
    # # vax_eff_wide <- dcast(vax_eff_long[outcome=="infection" & vaccine=="AZ"],strain ~ class,value.var = "value")
    # 
    # outcomes <- vax_eff_long[,unique(outcome)]
    # vaccines <- vax_eff_long[,unique(vaccine)]
    # doses <- vax_eff_long[,unique(dose)]
    # variants <- vax_eff_long[,unique(variant)]
    # 
    # vax_eff_by_age <- CJ(age_group = age_groups, outcome = outcomes, vaccine = vaccines, dose = doses, variant = variants, sorted = F)
    # vax_eff_by_age <- merge(vax_eff_by_age, vax_eff_long, by = c("outcome","vaccine","dose","variant"), all.x = T, sort = F)
    # 
    # prop_vax_type <- CJ(vaccine = vaccines,age_group = age_groups)
    # # FOR NOW: Assume all Pfizer based on Mai's comment that most vaccinations were
    # # Pfizer with some Janssen
    # prop_vax_type[,prop := fcase(vaccine == "AZ",0,
    #                              vaccine == "PF/MD",1)]
    # 
    # vax_eff_by_age <- merge(vax_eff_by_age, prop_vax_type, by = c("vaccine","age_group"), all.x = T, sort = F)
    # # Average effectiveness according to proportions of vaccine types by age
    # vax_eff_by_age <- vax_eff_by_age[,.(value = sum(value * prop)), by = .(age_group, outcome, dose, variant)]
    # 
    # vax_eff_arr <- array(vax_eff_by_age[,value], dim = c(length(variants),length(doses),length(outcomes),length(age_groups)),
    #                      dimnames = list(variant = variants, dose = doses, outcome = outcomes, age_group = age_groups))
    # # vax_eff_arr <- array(vax_eff_long$value, dim = vax_eff_long[,sapply(.SD,function(x) length(unique(x))),.SDcols = names(vax_eff_long)[names(vax_eff_long)!="value"]],
    # #                      dimnames = lapply(vax_eff_long[,.SD,.SDcols = names(vax_eff_long)[names(vax_eff_long)!="value"]],function(x) unique(x)))
    # 
    # 
    # vax_eff_arr <- aperm(vax_eff_arr, c(4,1,2,3))
    # 
    # # Calculate relative susceptibility and infectiousness, and conditional 
    # # probabilities of symptoms, hospitalisation and death in different vaccination 
    # # groups according to average vaccine effectiveness for different age groups and 
    # # variants
    # rel_params <- calculate_rel_param(vax_eff_arr)
    # 
    # # Mirror parameters for pseudo-strains
    # # strain 3: strain 1 -> strain 2
    # # strain 4: strain 2 -> strain 1
    # rel_params <- lapply(rel_params, mirror_strain)
    
    # Calculate relative susceptibility and infectiousness, and conditional 
    # probabilities of symptoms, hospitalisation and death in different vaccination 
    # groups according to average vaccine effectiveness for different age groups and 
    # variants
    # Delta relative to WT/Alpha
    rel_params <- convert_eff_to_rel_param(vax_eff_long,age_groups)
    # Omicron relative to Delta
    rel_params1 <- convert_eff_to_rel_param(vax_eff_long1,age_groups)
    
    # Extract individual parameters
    rel_susceptibility <- rel_params$rel_susceptibility
    rel_p_sympt <- rel_params$rel_p_sympt
    rel_p_hosp_if_sympt <- rel_params$rel_p_hosp_if_sympt
    rel_p_death <- rel_params$rel_p_death
    rel_infectivity <- rel_params$rel_infectivity
    
    rel_susceptibility1 <- rel_params1$rel_susceptibility
    rel_p_sympt1 <- rel_params1$rel_p_sympt
    rel_p_hosp_if_sympt1 <- rel_params1$rel_p_hosp_if_sympt
    rel_p_death1 <- rel_params1$rel_p_death
    rel_infectivity1 <- rel_params1$rel_infectivity
    
    # Generate model
    dt <- 0.25
    n_age <- length(age_groups)
    n_vax <- 5
    
    # Transmission and natural history parameters
    # intvtn_date <- as.Date(c(strt_date-1,"2020-08-27","2020-10-24","2021-06-01","2021-08-12"))
    # intvtn_date <- as.Date(c(strt_date-1,"2020-10-24","2021-06-01","2021-08-12"))
    # intvtn_date <- as.Date(c(strt_date-1,"2020-10-24","2021-06-30","2021-08-12","2021-12-31","2022-02-15"))
    # intvtn_date <- as.Date(c(strt_date-1,"2020-10-24","2021-06-30","2021-08-12","2021-12-31"))
    # intvtn_date <- as.Date(c(strt_date-1,"2020-08-01","2020-08-27","2020-10-10","2021-01-19","2021-06-01","2021-07-26","2021-08-02","2021-09-20","2021-11-15"))
    intvtn_date <- as.Date(c("2020-08-27","2020-10-24","2021-06-01","2021-08-02","2021-11-15"))
    # intvtn_date <- as.Date(c("2020-08-27","2020-10-10","2020-11-01","2021-06-01","2021-07-26","2021-08-02","2021-11-15"))
    beta_date <- sircovid_date(intvtn_date) #as.integer(intvtn_date - min(intvtn_date))
    # beta_value_sim <- c(0.035,0.025,0.02,0.04,0.02) #7/8*
    # beta_value_sim <- c(0.025,0.02,0.025,0.02,0.025,0.02) #7/8*
    # beta_value_sim <- c(0.025,0.02,0.025,0.02,0.025) #7/8*
    # beta_value_sim <- c(0.025,0.024,0.022,0.02,0.022,0.024,0.022,0.02,0.022,0.024) #7/8*
    beta_value_sim <- c(0.025,0.02,0.024,0.02,0.024) #7/8*
    # beta_value_sim <- c(0.025,0.023,0.02,0.024,0.023,0.02,0.024) #7/8*
    # beta_type <- "piecewise-constant"
    beta_type <- "piecewise-linear"
    gamma_E <- 0.5
    gamma_P <- 0.4
    gamma_A <- 0.2
    gamma_C <- 0.4
    gamma_H <- 0.1
    gamma_G <- 1/3
    gamma_pre_1 <- 1/13
    gamma_P_1 <- -log(0.86)/365
    theta_A <- 0.5
    
    # relative transmissibilities of 1st and 2nd strains
    # strain_transmission <- c(1,2)
    
    # Seeding parameters
    # 1st strain
    # start_date <- 1L
    initial_seed_size <- 10
    initial_seed_pattern <- 1
    
    # 2nd strain
    # strain_seed_date <- 70
    strain_seed_size <- 10
    strain_seed_pattern <- 1
    
    # Vaccination parameters
    # vaccine_progression_rate <- c(0,0,1/(26*7),0,-log(67.7/82.8)/(105-25)) # (Stowe Nat Comm 2022 Table S11)
    # # vaccine_progression_rate <- c(0,0,1/(26*7),0,-log(66/89)/((5-1)*30) # (Ferdinands BMJ 2022))
    vaccine_progression_rate <- c(0,0,1/(26*7),0,-log(0.923)/140) # (Barnard Nat Comm 2022 Table S4)
    
    vaccine_index_dose2 <- 2L
    vaccine_index_booster <- 4L
    vaccine_catchup_fraction <- 1
    n_doses <- 3L #2L
    
    vacc_skip_progression_rate <- rep(0, n_vax)
    vacc_skip_to <- c(0L,0L,5L,0L,0L) #integer(n_vax) #
    vacc_skip_weight <- c(0,0,1,0,0) #integer(n_vax) #
    
    # Relative probabilities of symptoms, hospitalisation and death for different strains
    strain_rel_p_sympt <- 1
    strain_rel_p_hosp_if_sympt <- c(1,1.6*1.85) #1 #
    strain_rel_p_death <- 1
    
    strain_rel_p_sympt1 <- 1
    strain_rel_p_hosp_if_sympt1 <- c(1,strain_rel_p_hosp_if_sympt[2]*0.3) #c(1,1.85*0.3) #
    strain_rel_p_death1 <- 1
    
    # # Parameters for impact of vaccination on susceptibility and infectiousness
    # rel_susceptibility <- c(1,0.8,0.5,0.8,0.5) # relative susceptibility to infection in each vaccine stratum
    # rel_p_sympt <- c(1,0.6,0.3,0.6,0.3) # relative risk of symptoms in each vaccine stratum
    # rel_p_hosp_if_sympt <- c(1,0.95,0.95,0.95,0.95) # relative risk of hospitalisation given infection in each vaccine stratum
    # rel_p_death <- c(1,0.9,0.9,0.9,0.9) # relative risk of death in each vaccine stratum
    # rel_infectivity <- c(1,0.5,0.5,0.5,0.5) # relative infectiousness of infected individuals in each vaccine stratum
    
    # Waning parameters
    waning_rate <- 1/(6*365)
    
    # Cross immunity parameters
    cross_immunity <- c(0.95,1) # 1 #
    
    # Sensitivity and specificity of serological tests
    sero_sensitivity_1 <- 1 #0.9
    sero_specificity_1 <- 0.99
    
    # Construct parameters object for 1st epoch
    p <- parameters(dt,
                    n_age,
                    n_vax,
                    transmission,
                    beta_date,
                    beta_value = beta_value_sim,
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
                    p_H_max0*p_H,
                    p_G,
                    # p_D,
                    p_D_max0*p_D,
                    p_P_1,
                    population,
                    start_date = 1L,
                    initial_seed_size,
                    initial_seed_pattern,
                    strain_transmission = c(1,2),
                    strain_seed_date = 330L,
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
                    phi_cases = 1,
                    kappa_cases = 2,
                    sero_sensitivity_1,
                    sero_specificity_1)
    
    # Number of steps in 1st epoch
    n_steps <- 500/dt
    
    # Parameters for 2nd epoch (in which new strain is introduced)
    # relative transmissibility of strain 3 is 4 times that of strain 1 
    # (which is combined with strain 2 in model in 2nd epoch as strain 1, and 
    # strain 3 becomes strain 2)
    # rel_strain_transmission1 <- 4
    # lower cross-immunity from previous infection with strains 1 or 2, but can't be
    # infected with strains 1 or 2 after being infected with strain 3
    cross_immunity1 <- c(0.55,1) #c(0.2,1) #
    p1 <- parameters(dt,
                     n_age,
                     n_vax,
                     transmission,
                     beta_date,
                     beta_value = beta_value_sim,
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
                     p_H_max0*p_H,
                     p_G,
                     # p_D,
                     p_D_max0*p_D,
                     p_P_1,
                     population,
                     start_date = n_steps*dt,
                     initial_seed_size = 0,
                     initial_seed_pattern,
                     c(2,4),
                     strain_seed_date = 510L,
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
                     phi_cases = 1,
                     kappa_cases = 2,
                     sero_sensitivity_1,
                     sero_specificity_1)

    n_steps1 <- max(data_raw$day)/dt # number of steps to run model up to
    
    
    
    #### Try simulating some data to see if the parameters give the right general pattern ####
    
    
    # out <- simulate(gen_seirhd_age_vax_multistrain_sero_time_dep_beta, p, n_steps)
    out <- simulate(covid_multi_strain, p, n_steps, deterministic,
                    p1 = p1, n_steps1 = n_steps1, transform = rotate_strains)
    
    # # Drop time row
    # x <- out$x
    # time <- x[1,1,-1]
    # x <- x[-1, , ,drop=FALSE]
    # 
    # n_strains <- 4
    # 
    # # # Plot trajectories
    # # plot_trajectories(time,x,n_age,n_strains,n_vax)
    # 
    # # Extract true history of model states
    # true_history <- x[ , ,seq(0,n_steps1+1,by=1/dt)+1,drop=F]
    
    # Add noise to simulated data
    info <- out$info
    idx <- index(info, min_ages, Rt)
    # hosps <- true_history[idx$state[grep("hosps_",names(idx$state))]-1, ,-1]
    # deaths <- true_history[idx$state[grep("deaths_",names(idx$state))]-1, ,-1]
    # sero_pos <- true_history[idx$state[grep("sero_pos_1_",names(idx$state))]-1, ,-1]
    # cases <- matrix(true_history[idx$state[match("cases",names(idx$state))]-1, ,-1],nrow = 1)
    # cases_non_variant <- matrix(true_history[idx$state[grep("cases_non_variant",names(idx$state))]-1, ,-1],nrow = 1)
    # 
    # par(mfrow = c(1,1))
    # days <- seq(1,n_steps1*dt)
    # matplot(days,t(hosps),type="l",xlab="Day",ylab="Hospitalisations")
    # matplot(days,t(deaths),type="l",xlab="Day",ylab="Deaths")
    # matplot(days,t(sero_pos),type="l",xlab="Day",ylab="Seropositive")
    # matplot(days,t(cases),type="l",xlab="Day",ylab="Cases")
    # matplot(days,t(cases_non_variant),type="l",xlab="Day",ylab="Non-variant cases")
    
    
    
    #### Fit to multiple age-stratified data streams ####
    
    
    # Convert raw data to required format for particle filter
    data <- particle_filter_data(data_raw,"day",1/dt)
    
    # data_raw1 <- data_raw
    # data1 <- data
    # data_raw <- data_raw[1:250,]
    # data <- data[1:250,]
    
    # Create multistage parameters object
    epochs <- list(
        multistage_epoch(n_steps*dt, pars = p1, transform_state = transform_state))
    pars <- multistage_parameters(p,epochs)
    
    # Create particle filter object
    # filter <- particle_filter$new(data, gen_seirhd_age_vax_multistrain_sero_time_dep_beta, n_particles, 
    #                               compare, index, initial)
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
    
    # tstart <- Sys.time()
    # filter$run(
    #     pars = pars,
    #     save_history = TRUE)
    # tend <- Sys.time()
    # print(tend - tstart)
    
    # # Plot filtered trajectories
    # dates <- sircovid_date_as_date(data_raw$day)
    # # plot_particle_filter(filter$history(),true_history,data_raw$day,idx)
    # plot_hosps_age(filter$history(),data,dates,n_age,n_vax,n_strains)
    # plot_deaths_age(filter$history(),data,dates,n_age,n_vax,n_strains)
    # plot_sero(filter$history(),data,dates,population[3:length(population)])
    # plot_cases(filter$history(),data,dates)
    
    # Infer parameters by pMCMC
    # beta_value_list <- list(name = c("beta1","beta2","beta3","beta4","beta5"), initial = c(0.035,0.025,0.02,0.04,0.02),
    #                         min = rep(0,5), max = rep(Inf,5), discrete = rep(F,5),
    #                         prior = replicate(5,function(x) dgamma(x, shape = 1, scale = 1, log = TRUE)))
    # beta_value_list <- list(
    #     name = c("beta1","beta2","beta3","beta4","beta5","beta6"), initial = c(0.025,0.02,0.025,0.02,0.025,0.02),
    #     min = rep(0,6), max = rep(Inf,6), discrete = rep(F,6),
    #     prior = #replicate(4,function(x) dgamma(x, shape = 1, scale = 1, log = TRUE))
    #         list(function(x) dgamma(x, shape = 4, scale = 0.025/4, log = TRUE),
    #              function(x) dgamma(x, shape = 4, scale = 0.02/4, log = TRUE),
    #              function(x) dgamma(x, shape = 4, scale = 0.025/4, log = TRUE),
    #              function(x) dgamma(x, shape = 4, scale = 0.02/4, log = TRUE),
    #              function(x) dgamma(x, shape = 4, scale = 0.025/4, log = TRUE),
    #              function(x) dgamma(x, shape = 4, scale = 0.02/4, log = TRUE)
    #         )
    # )
    # beta_value_list <- list(
    #     name = c("beta1","beta2","beta3","beta4","beta5"), initial = c(0.025,0.02,0.025,0.02,0.025),
    #     min = rep(0,5), max = rep(Inf,5), discrete = rep(F,5),
    #     prior = #replicate(4,function(x) dgamma(x, shape = 1, scale = 1, log = TRUE))
    #         list(function(x) dgamma(x, shape = 4, scale = 0.025/4, log = TRUE),
    #              function(x) dgamma(x, shape = 4, scale = 0.02/4, log = TRUE),
    #              function(x) dgamma(x, shape = 4, scale = 0.025/4, log = TRUE),
    #              function(x) dgamma(x, shape = 4, scale = 0.02/4, log = TRUE),
    #              function(x) dgamma(x, shape = 4, scale = 0.025/4, log = TRUE)
    #         )
    # )
    #
    n_betas <- length(beta_date)
    beta_names <- paste0("beta", seq_len(n_betas))
    # beta_initial <- c(0.025,0.024,0.022,0.02,0.022,0.024,0.022,0.02,0.022,0.024)
    beta_initial <- c(0.025,0.02,0.024,0.02,0.024)
    # beta_initial <- c(0.025,0.023,0.02,0.024,0.023,0.02,0.024)
    beta_value_list <- list(
        name = beta_names, initial = beta_initial,
        min = rep(0,n_betas), max = rep(Inf,n_betas), discrete = rep(F,n_betas),
        prior = replicate(n_betas,function(x) dgamma(x, shape = 4, scale = 0.02/4, log = TRUE))
    )
    beta_value <- Map(pmcmc_parameter,beta_value_list$name,beta_value_list$initial,
                      beta_value_list$min,beta_value_list$max,
                      beta_value_list$discrete,beta_value_list$prior)
    # R0 <- pmcmc_parameter("R0",2,min = 0)
    # gamma <- pmcmc_parameter("gamma",0.4,min = 0,
    #                          prior = function(x) dgamma(x,shape = 1,scale = 1,log = TRUE))
    rel_strain_transmission <- pmcmc_parameter("rel_strain_transmission",2.8,min = 0.25, max = 4)
    start_date <- pmcmc_parameter("start_date",1,min = 1,max = 10)
    strain_seed_date <- pmcmc_parameter("strain_seed_date",330,min = 305,max = 345)
    p_H_max <- pmcmc_parameter("p_H_max",p_H_max0/2,min = 0,max = 1,
                               prior = function(x) dbeta(x, 1, 1, log = TRUE)
                               # prior = function(x) dbeta(x, 10, 30, log = TRUE),
                               # prior = function(x) dbeta(x, 6, 30, log = TRUE)
                               )
    p_D_max <- pmcmc_parameter("p_D_max",p_D_max0,min = 0,max = 1,
                               prior = function(x) dbeta(x, 1, 1, log = TRUE))
    rel_strain_transmission1 <- pmcmc_parameter("rel_strain_transmission1",3.5,min = 2, max = 6)
    strain_seed_date1 <- pmcmc_parameter("strain_seed_date1",505,min = 500,max = 512)
    phi_cases <- pmcmc_parameter("phi_cases",0.5,min = 0,max = 1,
                                 prior = function(x) dbeta(x, 1, 1, log = TRUE))
    alpha_cases <- pmcmc_parameter("alpha_cases",0.5,min = 0,max = 1,
                                   prior = function(x) dbeta(x, 1, 1, log = TRUE))
    # proposal <- matrix(c(0.01^2,0,0,0.01^2),nrow = 2,ncol = 2,byrow = TRUE)
    # proposal <- matrix(c(0.01^2,0,0,0,0.01^2,0,0,0,2),nrow = 3,ncol = 3,byrow = TRUE)
    # proposal <- matrix(c(0.01^2,0,0,0,
    #                      0,0.01^2,0,0,
    #                      0,0,2,0,
    #                      0,0,0,2),nrow = 4,ncol = 4,byrow = TRUE)
    # proposal <- matrix(c(0.01^2,0,0,0,0,
    #                      0,0.01^2,0,0,0,
    #                      0,0,2,0,0,
    #                      0,0,0,2,0,
    #                      0,0,0,0,2),nrow = 5,ncol = 5,byrow = TRUE)
    # proposal <- diag(c(rep(1e-7,length(beta_value)),0.1^2,rep(2^2,2)))
    proposal <- 0.1*diag(c(rep(1e-7,length(beta_value)),0.1^2,rep(2^2,2),rep(1e-5,2),0.1^2,2^2,1e-5,1e-5))
    # transform <- make_transform(dt,
    #                             n_age,
    #                             n_vax,
    #                             transmission,
    #                             beta_date,
    #                             beta_type,
    #                             gamma_E,
    #                             gamma_P,
    #                             gamma_A,
    #                             gamma_C,
    #                             gamma_H,
    #                             gamma_G,
    #                             gamma_pre_1,
    #                             gamma_P_1,
    #                             theta_A,
    #                             p_C,
    #                             p_H,
    #                             p_G,
    #                             p_D,
    #                             p_P_1,
    #                             population,
    #                             # start_date,
    #                             initial_seed_size,
    #                             initial_seed_pattern,
    #                             # strain_transmission,
    #                             # strain_seed_date,
    #                             strain_seed_size,
    #                             strain_seed_pattern,
    #                             strain_rel_p_sympt,
    #                             strain_rel_p_hosp_if_sympt,
    #                             strain_rel_p_death,
    #                             rel_susceptibility,
    #                             rel_p_sympt,
    #                             rel_p_hosp_if_sympt,
    #                             rel_p_death,
    #                             rel_infectivity,
    #                             vaccine_progression_rate,
    #                             schedule,
    #                             vaccine_index_dose2,
    #                             vaccine_index_booster,
    #                             vaccine_catchup_fraction,
    #                             n_doses,
    #                             vacc_skip_progression_rate,
    #                             vacc_skip_to,
    #                             vacc_skip_weight,
    #                             waning_rate,
    #                             cross_immunity,
    #                             sero_sensitivity_1,
    #                             sero_specificity_1)
    # pars_mcmc <- c(beta_value,
    #                list(rel_strain_transmission = rel_strain_transmission,
    #                     start_date = start_date,
    #                     strain_seed_date = strain_seed_date))
    transform <- make_transform_multistage(dt,
                                           n_age,
                                           n_vax,
                                           transmission,
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
                                           start_date1 = n_steps * dt,
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
                                           sero_specificity_1)
    # # pars_mcmc <- c(beta_value,
    # #                list(gamma = gamma,
    # #                     start_date = start_date,
    # #                     strain_seed_date = strain_seed_date,
    # #                     strain_seed_date1 = strain_seed_date1))
    # mcmc_pars <- pmcmc_parameters$new(pars_mcmc,proposal,transform = transform)
    # # mcmc_pars <- pmcmc_parameters$new(list(R0 = R0,gamma = gamma),
    # #                                   proposal,transform = transform)
    # 
    # # Run MCMC
    # # Using mcstate
    # n_samples <- 100 #1000
    # control <- pmcmc_control(
    #     n_samples,
    #     save_state = TRUE,
    #     save_trajectories = TRUE,
    #     progress = TRUE)
    # pmcmc_run <- pmcmc(mcmc_pars,filter,control = control)
    # 
    # # Check variable names in saved state trajectories
    # dimnames(pmcmc_run$trajectories$state)
    # 
    # plot_particle_filter(pmcmc_run$trajectories$state,true_history,data_raw$day,idx)
    # plot_hosps_and_deaths_age(pmcmc_run$trajectories$state,data,data_raw$day,n_age,n_vax,n_strains)
    # plot_sero(pmcmc_run$trajectories$state,data,data_raw$day)
    # plot_cases(pmcmc_run$trajectories$state,data,data_raw$day)
    # 
    # # Plot MCMC output
    # mcmc_out <- as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars))
    # summary(mcmc_out)
    # plot(mcmc_out)
    
    # Using custom accelerated adaptive MCMC algorithm
    # init_pars <- c(beta_value_list$initial,rel_strain_transmission$initial,start_date$initial,strain_seed_date$initial)
    # priors <- c(beta_value_list$prior,rel_strain_transmission$prior,replicate(2,function(x) 0))
    init_pars <- c(beta_value_list$initial,
                   rel_strain_transmission$initial,
                   start_date$initial,strain_seed_date$initial,
                   p_H_max$initial,p_D_max$initial,
                   rel_strain_transmission1$initial,
                   strain_seed_date1$initial,
                   phi_cases$initial,
                   alpha_cases$initial)
    priors <- c(beta_value_list$prior,
                rel_strain_transmission$prior,
                start_date$prior,
                strain_seed_date$prior,
                p_H_max$prior,
                p_D_max$prior,
                rel_strain_transmission1$prior,
                strain_seed_date1$prior,
                phi_cases$prior,
                alpha_cases$prior)
    # n_iters <- 100 #100
    scaling_factor_start <- 1
    # pars_min <- c(beta_value_list$min,rel_strain_transmission$min,start_date$min,strain_seed_date$min)
    # pars_max <- c(beta_value_list$max,rel_strain_transmission$max,start_date$max,strain_seed_date$max)
    pars_min <- c(beta_value_list$min,
                  rel_strain_transmission$min,
                  start_date$min,
                  strain_seed_date$min,
                  p_H_max$min,
                  p_D_max$min,
                  rel_strain_transmission1$min,
                  strain_seed_date1$min,
                  phi_cases$min,
                  alpha_cases$min)
    pars_max <- c(beta_value_list$max,
                  rel_strain_transmission$max,
                  start_date$max,
                  strain_seed_date$max,
                  p_H_max$max,
                  p_D_max$max,
                  rel_strain_transmission1$max,
                  strain_seed_date1$max,
                  phi_cases$max,
                  alpha_cases$max)
    iter0 <- 100 #10
    # discrete <- c(rep(F,length(beta_value)+1),rep(T,2))
    # discrete <- c(rep(F,length(beta_value)+1),rep(T,2),rep(F,2))
    discrete <- c(rep(F,length(init_pars)))
    # names(init_pars) <- names(priors) <- names(pars_min) <- names(pars_max) <- names(discrete) <- c(beta_value_list$name,"rel_strain_transmission","start_date","strain_seed_date")
    names(init_pars) <- names(priors) <- names(pars_min) <- names(pars_max) <- 
        names(discrete) <- 
        c(beta_value_list$name,
          "rel_strain_transmission",
          "start_date",
          "strain_seed_date",
          "p_H_max",
          "p_D_max",
          "rel_strain_transmission1",
          "strain_seed_date1",
          "phi_cases",
          "alpha_cases")
    
    tstart <- Sys.time()
    # u <- 1:8 # all parameters 1:5 # only update beta parameters
    # u <- 1:7 # all parameters 1:5 # only update beta parameters
    # u <- c(1,6:7) #c(1:4,6:8) #c(1:4,6:9) #1:9 # all parameters 1:5 # only update beta parameters
    res <- mcmc(transform,filter,init_pars,priors,n_particles,n_iters,idx,
                scaling_factor_start,proposal,pars_min,pars_max,iter0,discrete,
                u,thinning)
    tend <- Sys.time()
    print(tend - tstart)
    ## Time difference of 26.56045 mins (R 4.0.5)
    # Time difference of 15.13305 mins (R 4.1.0 for Mac M1)
    
    save(list = ls(all.names = T), file = paste0("output/MCMCoutput",run,".RData"), envir = environment())
    
    # Plot p.w. constant beta to check it looks right
    beta_t <- seq(0, beta_date[length(beta_date)], by = dt)
    pdf(paste0("output/plots",run,".pdf")) # save all plots into one pdf
    par(mfrow = c(1,1))
    # plot(beta_t, p$beta_step, type="o", cex = 0.25)
    # points(beta_date, beta_value_sim, pch = 19, col = "red")
    
    n_smpls <- round(n_iters/thinning)
    beta_value_post <- apply(res$pars[seq(round(n_smpls/10),n_smpls,by=10),1:n_betas],2,median)
    if (beta_type == "piecewise-linear") {
        beta_step <- parameters_piecewise_linear(beta_date, 
                                                 beta_value_post %||% 0.1, dt)
    } else if (beta_type == "piecewise-constant") {
        beta_step <- parameters_piecewise_constant(beta_date, 
                                                   beta_value_post %||% 0.1, dt)
    }
    plot(sircovid_date_as_date(beta_t), beta_step, type = "o", cex = 0.25, xaxt = "n", xlab = "Date", ylab = "beta")
    dates_plot <- seq.Date(sircovid_date_as_date(beta_t[1]),sircovid_date_as_date(beta_t[length(beta_t)]),by = 30)
    axis(1, dates_plot, format(dates_plot,"%Y-%m-%d"))
    
    # Trace plots
    # par(mfrow = c(ceiling(length(init_pars)/2),2))
    # for (i in seq_along(init_pars)){
    #     plot(res$pars[,i],type = "l",xlab = "Iteration",ylab = names(init_pars)[i])
    # }
    par(mar = c(2.5,4,2,1), mfrow = c(ceiling(length(u)/2),2))
    for (i in u){
        plot(res$pars[,i],type = "l",xlab = "Iteration",ylab = names(init_pars)[i])
    }
    
    # Scaling factor
    par(mfrow = c(1,1))
    plot(res$scaling_factor,type="l")
    
    par(mfrow = c(1,1))
    # pairs(res$pars[seq(round(n_smpls/10),n_smpls,by=10), ],labels = names(init_pars))
    pairs(res$pars[seq(round(n_smpls/10),n_smpls,by=10),u],labels = names(init_pars)[u])
    
    # Plot fitted hospitalisations and deaths against data
    print(plot_outcome(res$trajectories$state,data,"hosps",by_age = T))
    print(plot_outcome(res$trajectories$state,data,"deaths",by_age = T))
    print(plot_outcome(res$trajectories$state,data,"cases",res$pars[,"phi_cases"],by_age = T))
    print(plot_sero(res$trajectories$state,data,agg_pop,by_age = T))
    print(plot_outcome(res$trajectories$state,data,"hosps"))
    print(plot_outcome(res$trajectories$state,data,"deaths"))
    print(plot_outcome(res$trajectories$state,data,"cases",res$pars[,"phi_cases"]))
    
    dev.off()
    
}