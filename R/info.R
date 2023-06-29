create_info <- function(beta_date,model_type){
    n_beta <- length(beta_date)
    
    start_date <- list(init = "2020-07-15",
                       min = "2020-07-01",
                       max = "2020-07-31")
    start_date <- lapply(start_date,covid_multi_strain_date)
    
    strain_seed_date <- list(init = "2021-06-14", 
                             min = "2021-05-20",
                             max = "2021-06-29")
    strain_seed_date <- lapply(strain_seed_date,covid_multi_strain_date)
    
    strain_seed_date1 <- list(init = "2021-12-01", 
                              min = "2021-11-21",
                              max = "2021-12-13")
    strain_seed_date1 <- lapply(strain_seed_date1,covid_multi_strain_date)
    
    strain_seed_date2 <- list(init = "2021-12-31", 
                              min = "2021-12-16",
                              max = "2022-01-10")
    strain_seed_date2 <- lapply(strain_seed_date2,covid_multi_strain_date)
    
    pars_info <- data.frame(
        name = c(paste0(rep("beta",n_beta),seq_len(n_beta)),
                 "rel_strain_transmission","start_date","strain_seed_date",
                 "p_H_max","p_D","p_D_2","p_D_3","rel_strain_transmission1","strain_seed_date1",
                 # "rel_strain_transmission2","strain_seed_date2",
                 if (model_type == "NB"){
                     c("phi_cases","alpha_cases")
                 } else if (model_type == "BB"){
                     c("p_NC","rho_tests")
                 },
                 "alpha_hosp","alpha_death"),
        initial = c(0.025,0.02,0.024,0.02,0.024,
                    2.7,start_date$init,strain_seed_date$init,
                    0.400626705/2,0.316,0.316,0.316,3.2,strain_seed_date1$init,
                    # 3.5,strain_seed_date2$init,
                    if (model_type == "NB"){
                        c(0.5,0.5)
                    } else if (model_type == "BB"){
                        c(0.01,0.01)
                    },
                    0.5,0.5),
        min = c(rep(0,n_beta),
                0.25,start_date$min,strain_seed_date$min,0,0,0,0,2,strain_seed_date1$min,
                # 3,strain_seed_date2$min,
                0,0,0,0),
        max = c(rep(1,n_beta),
                4,start_date$max,strain_seed_date$max,1,1,1,1,6,strain_seed_date1$max,
                # 9,strain_seed_date2$max,
                1,1,1,1),
        integer = c(rep(F,n_beta),
                    F,F,F,F,F,F,F,F,F,
                    # F,F,
                    F,F,F,F)
    )
    pars_info
}