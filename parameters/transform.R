compute_severity <- function(pars,severity,dt){
    expected <- c("p_D","p_D_2","p_D_3")
    stopifnot(all(expected %in% names(pars)))
    
    p_D <- pars["p_D"]
    p_D_2 <- pars["p_D_2"]
    p_D_3 <- pars["p_D_3"]
    
    # Probability of death given severe disease
    # p_D_date <- covid_multi_strain_date(c("2019-12-31","2021-07-01","2021-12-01"))
    p_D_date <- covid_multi_strain_date(c("2020-07-01","2021-07-01","2021-12-01"))
    p_D_value <- c(p_D,p_D_2,p_D_3)
    # p_D_value <- p_D
    
    # Probability of hospitalisation
    p_H_max = pars["p_H_max"]
    severity$p_H <- p_H_max * severity$p_H
    
    severity <- parameters_severity(
        dt,
        severity,
        p_D = list(value = p_D_value,date = p_D_date)
        # p_D = list(value = p_D_value)
        )
}


compute_observation <- function(pars,model_type){
    if (model_type == "NB"){
        expected_model_specific <- c("phi_cases","alpha_cases")
    } else if (model_type == "BB"){
        expected_model_specific <- c("p_NC","rho_tests")
    }
    expected <- c("alpha_hosp","alpha_death",expected_model_specific)
    stopifnot(all(expected %in% names(pars)))
    
    observation <- list()
    if (model_type == "NB"){
        observation$phi_cases <- pars[["phi_cases"]]
        observation$kappa_cases <- 1/pars[["alpha_cases"]]    
    } else if (model_type == "BB"){
        observation$p_NC <- pars[["p_NC"]]
        observation$rho_tests <- pars[["rho_tests"]]
    }
    observation$kappa_hosp <- 1/pars[["alpha_hosp"]]
    observation$kappa_death <- 1/pars[["alpha_death"]]
    
    observation
}


apply_assumptions <- function(baseline, assumptions){
    stopifnot(assumptions %in% names(baseline$vaccine_progression_rate))
    baseline$vaccine_progression_rate <-
        baseline$vaccine_progression_rate[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$waning_rate))
    baseline$waning_rate <-
        baseline$waning_rate[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$cross_immunity))
    baseline$cross_immunity <-
        baseline$cross_immunity[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$cross_immunity1))
    baseline$cross_immunity1 <-
        baseline$cross_immunity1[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$cross_immunity2))
    baseline$cross_immunity2 <-
        baseline$cross_immunity2[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_susceptibility))
    baseline$rel_susceptibility <- 
        baseline$rel_susceptibility[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_p_sympt))
    baseline$rel_p_sympt <- 
        baseline$rel_p_sympt[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_p_hosp_if_sympt))
    baseline$rel_p_hosp_if_sympt <- 
        baseline$rel_p_hosp_if_sympt[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_p_death))
    baseline$rel_p_death <- 
        baseline$rel_p_death[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_infectivity))
    baseline$rel_infectivity <- 
        baseline$rel_infectivity[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_susceptibility1))
    baseline$rel_susceptibility1 <- 
        baseline$rel_susceptibility1[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_p_sympt1))
    baseline$rel_p_sympt1 <- 
        baseline$rel_p_sympt1[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_p_hosp_if_sympt1))
    baseline$rel_p_hosp_if_sympt1 <- 
        baseline$rel_p_hosp_if_sympt1[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_p_death1))
    baseline$rel_p_death1 <- 
        baseline$rel_p_death1[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_infectivity1))
    baseline$rel_infectivity1 <- 
        baseline$rel_infectivity1[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_susceptibility2))
    baseline$rel_susceptibility2 <- 
        baseline$rel_susceptibility2[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_p_sympt2))
    baseline$rel_p_sympt2 <- 
        baseline$rel_p_sympt2[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_p_hosp_if_sympt2))
    baseline$rel_p_hosp_if_sympt2 <- 
        baseline$rel_p_hosp_if_sympt2[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_p_death2))
    baseline$rel_p_death2 <- 
        baseline$rel_p_death2[[assumptions]]
    
    stopifnot(assumptions %in% names(baseline$rel_infectivity2))
    baseline$rel_infectivity2 <- 
        baseline$rel_infectivity2[[assumptions]]
    
    baseline
}


make_transform <- function(baseline){
    
    # Expected fixed parameters
    expected <- c("model_type",
                  "epoch_dates",
                  "dt",
                  "age_groups",
                  "n_age",
                  "n_vax",
                  "m",
                  "beta_date",
                  "beta_names",
                  "beta_type",
                  "gamma_E",
                  "gamma_P",
                  "gamma_A",
                  "gamma_C",
                  "gamma_H",
                  "gamma_G",
                  "gamma_pre_1",
                  "gamma_P_1",
                  "theta_A",
                  "severity",
                  "population",
                  "initial_seed_size",
                  "initial_seed_pattern",
                  "strain_seed_size",
                  "strain_seed_pattern",
                  "strain_rel_p_sympt",
                  "strain_rel_p_hosp_if_sympt",
                  "strain_rel_p_death",
                  "rel_susceptibility",
                  "rel_p_sympt",
                  "rel_p_hosp_if_sympt",
                  "rel_p_death",
                  "rel_infectivity",
                  "vaccine_progression_rate",
                  "vaccine_schedule",
                  "vaccine_index_dose2",
                  "vaccine_index_booster",
                  "vaccine_catchup_fraction",
                  "n_doses",
                  "vacc_skip_progression_rate",
                  "vacc_skip_to",
                  "vacc_skip_weight",
                  "waning_rate",
                  "rel_gamma_wildtype_delta",
                  "cross_immunity",
                  "strain_rel_p_sympt1",
                  "strain_rel_p_hosp_if_sympt1",
                  "strain_rel_p_death1",
                  "rel_susceptibility1",
                  "rel_p_sympt1",
                  "rel_p_hosp_if_sympt1",
                  "rel_p_death1",
                  "rel_infectivity1",
                  "rel_gamma_delta_omicronba1",
                  "cross_immunity1",
                  "strain_rel_p_sympt2",
                  "strain_rel_p_hosp_if_sympt2",
                  "strain_rel_p_death2",
                  "rel_susceptibility2",
                  "rel_p_sympt2",
                  "rel_p_hosp_if_sympt2",
                  "rel_p_death2",
                  "rel_infectivity2",
                  "rel_gamma_omicronba1_omicronba2",
                  "cross_immunity2",
                  "sero_sensitivity_1",
                  "sero_specificity_1",
                  "test_sensitivity",
                  "test_specificity")
    stopifnot(setequal(expected, names(baseline)))
    
    epoch_dates <- baseline$epoch_dates
    
    # Expected parameters for fitting
    expected <- c(baseline$beta_names,"start_date","rel_strain_transmission",
                  "strain_seed_date","p_H_max","p_D","p_D_2","p_D_3",
                  "rel_strain_transmission1","strain_seed_date1",
                  # "rel_strain_transmission2","strain_seed_date2",
                  if (baseline$model_type == "NB"){
                      c("phi_cases","alpha_cases")
                  } else if (baseline$model_type == "BB") {
                      c("p_NC","rho_tests")
                  },
                  "alpha_hosp","alpha_death")
    
    function(pars){
        stopifnot(setequal(expected, names(pars)))
        
        severity <- compute_severity(pars,baseline$severity,baseline$dt)
        observation <- compute_observation(pars,baseline$model_type)
        
        beta_value <- unname(pars[baseline$beta_names])
        start_date <- pars[["start_date"]]
        rel_strain_transmission <- pars[["rel_strain_transmission"]]
        strain_seed_date <- pars[["strain_seed_date"]]
        rel_strain_transmission1 <- pars[["rel_strain_transmission1"]]
        strain_seed_date1 <- pars[["strain_seed_date1"]]
        # rel_strain_transmission2 <- pars[["rel_strain_transmission2"]]
        # strain_seed_date2 <- pars[["strain_seed_date2"]]
        
        stage_parameters <- function(strains){
            
            if (strains == "Wildtype_Delta"){
                strain_transmission <- c(1,rel_strain_transmission)
                strain_seed_date <- strain_seed_date
                strain_rel_p_sympt <- baseline$strain_rel_p_sympt
                strain_rel_p_hosp_if_sympt <- baseline$strain_rel_p_hosp_if_sympt
                strain_rel_p_death <- baseline$strain_rel_p_death
                rel_susceptibility <- baseline$rel_susceptibility
                rel_p_sympt <- baseline$rel_p_sympt
                rel_p_hosp_if_sympt <- baseline$rel_p_hosp_if_sympt
                rel_p_death <- baseline$rel_p_death
                rel_infectivity <- baseline$rel_infectivity
                strain_rel_gamma = baseline$rel_gamma_wildtype_delta
                cross_immunity <- baseline$cross_immunity
            } else if (strains == "Delta_Omicron"){
                strain_transmission <- c(rel_strain_transmission,rel_strain_transmission1)
                strain_seed_date <- strain_seed_date1
                strain_rel_p_sympt <- baseline$strain_rel_p_sympt1
                strain_rel_p_hosp_if_sympt <- baseline$strain_rel_p_hosp_if_sympt1
                strain_rel_p_death <- baseline$strain_rel_p_death1
                rel_susceptibility <- baseline$rel_susceptibility1
                rel_p_sympt <- baseline$rel_p_sympt1
                rel_p_hosp_if_sympt <- baseline$rel_p_hosp_if_sympt1
                rel_p_death <- baseline$rel_p_death1
                rel_infectivity <- baseline$rel_infectivity1
                strain_rel_gamma = baseline$rel_gamma_delta_omicronba1
                cross_immunity <- baseline$cross_immunity1
            }
            
            # } else if (strains == "Delta_OmicronBA1"){
            #     strain_transmission <- c(rel_strain_transmission,rel_strain_transmission1)
            #     strain_seed_date <- strain_seed_date1
            #     strain_rel_p_sympt <- baseline$strain_rel_p_sympt1
            #     strain_rel_p_hosp_if_sympt <- baseline$strain_rel_p_hosp_if_sympt1
            #     strain_rel_p_death <- baseline$strain_rel_p_death1
            #     rel_susceptibility <- baseline$rel_susceptibility1
            #     rel_p_sympt <- baseline$rel_p_sympt1
            #     rel_p_hosp_if_sympt <- baseline$rel_p_hosp_if_sympt1
            #     rel_p_death <- baseline$rel_p_death1
            #     rel_infectivity <- baseline$rel_infectivity1
            #     strain_rel_gamma = baseline$rel_gamma_delta_omicronba1
            #     cross_immunity <- baseline$cross_immunity1
            # } else if (strains == "OmicronBA1_OmicronBA2"){
            #     strain_transmission <- c(rel_strain_transmission1,rel_strain_transmission2)
            #     strain_seed_date <- strain_seed_date2
            #     strain_rel_p_sympt <- baseline$strain_rel_p_sympt2
            #     strain_rel_p_hosp_if_sympt <- baseline$strain_rel_p_hosp_if_sympt2
            #     strain_rel_p_death <- baseline$strain_rel_p_death2
            #     rel_susceptibility <- baseline$rel_susceptibility2
            #     rel_p_sympt <- baseline$rel_p_sympt2
            #     rel_p_hosp_if_sympt <- baseline$rel_p_hosp_if_sympt2
            #     rel_p_death <- baseline$rel_p_death2
            #     rel_infectivity <- baseline$rel_infectivity2
            #     strain_rel_gamma = baseline$rel_gamma_omicronba1_omicronba2
            #     cross_immunity <- baseline$cross_immunity2
            # }
            
            parameters(baseline$dt,
                       baseline$n_age,
                       baseline$n_vax,
                       baseline$m,
                       baseline$beta_date,
                       beta_value = beta_value,
                       baseline$beta_type,
                       baseline$gamma_E,
                       baseline$gamma_P,
                       baseline$gamma_A,
                       baseline$gamma_C,
                       baseline$gamma_H,
                       baseline$gamma_G,
                       baseline$gamma_pre_1,
                       baseline$gamma_P_1,
                       baseline$theta_A,
                       severity,
                       baseline$population,
                       start_date,
                       baseline$initial_seed_size,
                       baseline$initial_seed_pattern,
                       strain_transmission = strain_transmission,
                       strain_seed_date = strain_seed_date,
                       baseline$strain_seed_size,
                       baseline$strain_seed_pattern,
                       strain_rel_p_sympt = strain_rel_p_sympt,
                       strain_rel_p_hosp_if_sympt = strain_rel_p_hosp_if_sympt,
                       strain_rel_p_death = strain_rel_p_death,
                       strain_rel_gamma_E = strain_rel_gamma$E,
                       strain_rel_gamma_P = strain_rel_gamma$P,
                       strain_rel_gamma_C = strain_rel_gamma$C,
                       strain_rel_gamma_A = strain_rel_gamma$A,
                       rel_susceptibility = rel_susceptibility,
                       rel_p_sympt = rel_p_sympt,
                       rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
                       rel_p_death = rel_p_death,
                       rel_infectivity = rel_infectivity,
                       baseline$vaccine_progression_rate,
                       baseline$vaccine_schedule,
                       baseline$vaccine_index_dose2,
                       baseline$vaccine_index_booster,
                       baseline$vaccine_catchup_fraction,
                       baseline$n_doses,
                       baseline$vacc_skip_progression_rate,
                       baseline$vacc_skip_to,
                       baseline$vacc_skip_weight,
                       baseline$waning_rate,
                       cross_immunity = cross_immunity,
                       observation,
                       baseline$sero_sensitivity_1,
                       baseline$sero_specificity_1,
                       baseline$test_sensitivity,
                       baseline$test_specificity)
        }
        
        # Parameters for 1st epoch
        p <- stage_parameters("Wildtype_Delta")
        
        # Parameters for 2nd epoch
        p1 <- stage_parameters("Delta_Omicron")
        
        epochs <- list(
            multistage_epoch(epoch_dates[1], pars = p1, transform_state = transform_state)
        )
        # # Parameters for 2nd epoch
        # p1 <- stage_parameters("Delta_OmicronBA1")
        # 
        # # Parameters for 3rd epoch
        # p2 <- stage_parameters("OmicronBA1_OmicronBA2")
        # 
        # epochs <- list(
        #     multistage_epoch(epoch_dates[1], pars = p1, transform_state = transform_state),
        #     multistage_epoch(epoch_dates[2], pars = p2, transform_state = transform_state)
        # )
        p_multistage <- multistage_parameters(p, epochs)
        p_multistage
    }
}