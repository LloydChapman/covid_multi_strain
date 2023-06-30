create_baseline <- function(model_type,epoch_dates){
    # Load vaccination and population data
    vax <- fread("data/data_vaccination.csv", colClasses = c(number = "numeric"))
    pop <- fread("data/population.csv")
    
    # Set end date for data
    end_date <- as.Date("2022-05-06") # last death date in data files
    
    # Generate model
    dt <- 0.25
    
    # Set age groups
    age_groups <- c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+")
    n_age <- length(age_groups)
    n_vax <- 5
    
    # Transmission matrix
    # FOR NOW: use contact matrix for France (is Fiji an alternative as a Pacific island?)
    m <- transmission_matrix("FJI", pop, age_groups)
    
    # Transmission and natural history parameters
    intvtn_date <- as.Date(c("2020-08-27","2020-10-24","2021-06-01","2021-08-02","2021-11-15"))
    beta_date <- covid_multi_strain_date(intvtn_date)
    beta_names <- sprintf("beta%d",seq_along(beta_date))
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
    
    sympt_frac <- qread("data/2-linelist_both_fit_fIa0.5-rbzvih.qs")
    p_C <- unname(colMeans(sympt_frac[,13:20]))
    prob_death_given_hosp <- fread("data/prob_death_given_hosp_salje.csv")
    p_D <- prob_death_given_hosp[,median_perc_mean/100]
    IHR <- fread("data/ihr_salje.csv")
    ihr <- IHR[,median_perc_mean/100]
    # Set probability of death outside hospital from observed proportion of non-hospital deaths
    p_G <- rep(read.csv("data/prob_death_community.csv")[1,1],n_age)
    p_H <- ihr/(p_C*(1-p_G))
    p_P_1 <- 0.85
    
    # Normalise
    p_D <- p_D/max(p_D)
    p_H <- p_H/max(p_H)
    
    # Combine severity parameters into list
    severity <- list(p_C = p_C,
                     p_H = p_H,
                     p_G = p_G,
                     p_D = p_D,
                     p_P_1 = p_P_1)
    
    # Get population
    population <- pop[,.(sum(total)),by = .(age_group)][,V1]
    
    # Vaccine efficacy
    vax_eff <- fread("data/vax_eff.csv",colClasses = c(alpha_central = "numeric"))
    
    # Melt to long format
    vax_eff_long <- melt(vax_eff,measure.vars = patterns(central = "central",pessimistic = "pessimistic",optimistic = "optimistic"),variable.name = "variant")
    variants <- c("alpha","delta","omicron","omicron_ba2")
    vax_eff_long[, variant := variants[variant]]
    
    # Make vaccine efficacy data.tables for two-strain setup
    vax_eff_long2 <- vax_eff_long[variant %in% c("omicron","omicron_ba2")]
    vax_eff_long1 <- vax_eff_long[variant %in% c("delta","omicron")]
    vax_eff_long <- vax_eff_long[variant %in% c("alpha","delta")]
        
    # Calculate relative susceptibility and infectiousness, and conditional 
    # probabilities of symptoms, hospitalisation and death in different vaccination 
    # groups according to average vaccine effectiveness for different age groups and 
    # variants
    # Delta relative to WT/Alpha
    rel_params <- convert_eff_to_rel_param(vax_eff_long,age_groups)
    rel_params <- reverse_list_structure(rel_params)
    # Omicron BA.1 relative to Delta
    rel_params1 <- convert_eff_to_rel_param(vax_eff_long1,age_groups)
    rel_params1 <- reverse_list_structure(rel_params1)
    names(rel_params1) <- paste0(names(rel_params1),"1")
    # Omicron BA.2 relative to Omicron BA.1
    rel_params2 <- convert_eff_to_rel_param(vax_eff_long2,age_groups)
    rel_params2 <- reverse_list_structure(rel_params2)
    names(rel_params2) <- paste0(names(rel_params2),"2")
    
    
    # # Extract individual parameters
    # rel_susceptibility <- rel_params$rel_susceptibility
    # rel_p_sympt <- rel_params$rel_p_sympt
    # rel_p_hosp_if_sympt <- rel_params$rel_p_hosp_if_sympt
    # rel_p_death <- rel_params$rel_p_death
    # rel_infectivity <- rel_params$rel_infectivity
    # 
    # rel_susceptibility1 <- rel_params1$rel_susceptibility
    # rel_p_sympt1 <- rel_params1$rel_p_sympt
    # rel_p_hosp_if_sympt1 <- rel_params1$rel_p_hosp_if_sympt
    # rel_p_death1 <- rel_params1$rel_p_death
    # rel_infectivity1 <- rel_params1$rel_infectivity
    
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
    vaccine_progression_rate <- list(
        central = c(0,0,1/(26*7),0,-log(67.7/82.8)/(105-25)), # (Stowe Nat Comm 2022 Table S11)
        pessimistic = c(0,0,1/(26*7),0,-log(67.7/82.8)/(105-25)), # (Stowe Nat Comm 2022 Table S11)
        optimistic = c(0,0,1/(26*7),0,-log(0.923)/140) # (Barnard Nat Comm 2022 Table S4)
    )
    
    # Create vaccination schedule
    # Set delays for immune response to different vaccine doses
    delay_dose1 <- 28
    delay_dose2 <- 14
    
    # Get vaccination age groups
    age_groups_vax <- vax[,unique(age_group)]
    
    # Matrix of uptake rates (age group x dose)
    uptake <- matrix(1,nrow = length(age_groups),ncol = vax[,length(unique(dose))])
    
    # Make vaccine schedule
    vaccine_schedule <- vaccination_data(vax,delay_dose1,delay_dose2,pop,age_groups_vax,
                                         age_groups,end_date,uptake)
    
    vaccine_index_dose2 <- 2L
    vaccine_index_booster <- 4L
    vaccine_catchup_fraction <- 1
    n_doses <- vaccine_schedule$n_doses #3L #2L
    
    vacc_skip_progression_rate <- rep(0, n_vax)
    vacc_skip_to <- c(0L,0L,5L,0L,0L) #integer(n_vax) #
    vacc_skip_weight <- c(0,0,1,0,0) #integer(n_vax) #
    
    # Relative probabilities of symptoms, hospitalisation and death for different strains
    # Wildtype/Delta
    strain_rel_p_sympt <- 1
    strain_rel_p_hosp_if_sympt <- 1.6*1.85 #c(1,1.6*1.85) #1 #
    strain_rel_p_death <- 1
    
    # Delta/Omicron BA.1
    strain_rel_p_sympt1 <- 1
    strain_rel_p_hosp_if_sympt1 <- strain_rel_p_hosp_if_sympt*0.3 #c(1,strain_rel_p_hosp_if_sympt[2]*0.3) #c(1,1.85*0.3) #
    strain_rel_p_death1 <- 1
    
    # Omicron BA.1/Omicron BA.2
    strain_rel_p_sympt2 <- 1
    strain_rel_p_hosp_if_sympt2 <- 1
    strain_rel_p_death2 <- 1
    
    # Relative duration of serial interval (compared to Wildtype)
    # We assume Delta is 13% and Omicron BA.1 and BA.2 are 25% shorter than Wildtype (Perez-Guzman 2023)
    # Wildtype
    rel_si_wildtype <- 1
    # Delta
    rel_si_delta <- 0.87
    # Omicron BA.1
    rel_si_omicronba1 <- 0.75
    # Omicron BA.2
    rel_si_omicronba2 <- 0.75
    
    # Serial interval of variants
    # Wildtype/Delta
    rel_gamma_wildtype_delta <- list(E = c(1/rel_si_wildtype,1/rel_si_delta),
                                     P = c(1/rel_si_wildtype,1/rel_si_delta),
                                     C = c(1/rel_si_wildtype,1/rel_si_delta),
                                     A = c(1/rel_si_wildtype,1/rel_si_delta))
    
    # Delta/Omicron BA.1
    rel_gamma_delta_omicronba1 <- list(E = c(1/rel_si_delta,1/rel_si_omicronba1),
                                       P = c(1/rel_si_delta,1/rel_si_omicronba1),
                                       C = c(1/rel_si_delta,1/rel_si_omicronba1),
                                       A = c(1/rel_si_delta,1/rel_si_omicronba1))
    
    # Omicron BA.1/Omicron BA.2
    rel_gamma_omicronba1_omicronba2 <- list(E = c(1/rel_si_omicronba1,1/rel_si_omicronba2),
                                            P = c(1/rel_si_omicronba1,1/rel_si_omicronba2),
                                            C = c(1/rel_si_omicronba1,1/rel_si_omicronba2),
                                            A = c(1/rel_si_omicronba1,1/rel_si_omicronba2))
    
    # Waning parameters
    waning_rate <- list(central = 1/(6*365),
                        pessimistic = 1/(3*365),
                        optimistic = 1/(6*365))
    
    # Cross immunity parameters
    cross_immunity <- list(central = c(0.95,1),
                           pessimistic = c(0.75,1),
                           optimistic = c(1,1))# 1 #
    cross_immunity1 <- list(central = c(0.55,1),
                            pessimistic = c(0.25,1),
                            optimistic = c(0.55,1))
    cross_immunity2 <- list(central = c(0.5,0.8),
                            pessimistic = c(0.3,0.7),
                            optimistic = c(0.75,1))
    
    severity_cross_multiplier_delta <- 0.85
    severity_cross_multiplier_omicronba1 <- list(
        rel_p_hosp_if_sympt = 0.55,
        rel_p_death = 0.18
    )
    
    # Sensitivity and specificity of serological tests
    sero_sensitivity_1 <- 1 #0.9
    sero_specificity_1 <- 0.99
    
    # Sensitivity and specificity of PCR and rapid Ag tests
    test_sensitivity <- 1
    test_specificity <- 1
    
    baseline <- list(
        model_type = model_type,
        epoch_dates = covid_multi_strain_date(epoch_dates),
        dt = dt,
        age_groups = age_groups,
        n_age  = n_age,
        n_vax = n_vax,
        m = m,
        beta_date = beta_date,
        beta_names = beta_names,
        beta_type = beta_type,
        gamma_E = gamma_E,
        gamma_P = gamma_P,
        gamma_A = gamma_A,
        gamma_C = gamma_C,
        gamma_H = gamma_H,
        gamma_G = gamma_G,
        gamma_pre_1 = gamma_pre_1,
        gamma_P_1 = gamma_P_1,
        theta_A = theta_A,
        severity = severity,
        population = population,
        initial_seed_size = initial_seed_size,
        initial_seed_pattern = initial_seed_pattern,
        strain_seed_size = strain_seed_size,
        strain_seed_pattern = strain_seed_pattern,
        vaccine_progression_rate = vaccine_progression_rate,
        vaccine_schedule = vaccine_schedule,
        vaccine_index_dose2 = vaccine_index_dose2,
        vaccine_index_booster = vaccine_index_booster,
        vaccine_catchup_fraction = vaccine_catchup_fraction,
        n_doses = n_doses,
        vacc_skip_progression_rate = vacc_skip_progression_rate,
        vacc_skip_to = vacc_skip_to,
        vacc_skip_weight = vacc_skip_weight,
        strain_rel_p_sympt = strain_rel_p_sympt,
        strain_rel_p_hosp_if_sympt = strain_rel_p_hosp_if_sympt,
        strain_rel_p_death = strain_rel_p_death,
        strain_rel_p_sympt1 = strain_rel_p_sympt1,
        strain_rel_p_hosp_if_sympt1 = strain_rel_p_hosp_if_sympt1,
        strain_rel_p_death1 = strain_rel_p_death1,
        strain_rel_p_sympt2 = strain_rel_p_sympt2,
        strain_rel_p_hosp_if_sympt2 = strain_rel_p_hosp_if_sympt2,
        strain_rel_p_death2 = strain_rel_p_death2,
        rel_gamma_wildtype_delta = rel_gamma_wildtype_delta,
        rel_gamma_delta_omicronba1 = rel_gamma_delta_omicronba1,
        rel_gamma_omicronba1_omicronba2 = rel_gamma_omicronba1_omicronba2,
        waning_rate = waning_rate,
        cross_immunity = cross_immunity,
        cross_immunity1 = cross_immunity1,
        cross_immunity2 = cross_immunity2,
        severity_cross_multiplier_delta = severity_cross_multiplier_delta,
        severity_cross_multiplier_omicronba1 = severity_cross_multiplier_omicronba1,
        sero_sensitivity_1 = sero_sensitivity_1,
        sero_specificity_1 = sero_specificity_1,
        test_sensitivity =  test_sensitivity,
        test_specificity = test_specificity
    )
    baseline <- c(baseline, rel_params, rel_params1, rel_params2)
    
}