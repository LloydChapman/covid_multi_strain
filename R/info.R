create_info <- function(beta_date){
    n_beta <- length(beta_date)
    pars_info <- data.frame(
        name = c(paste0(rep("beta",n_beta),seq_len(n_beta)),
                 "rel_strain_transmission","start_date","strain_seed_date",
                 "p_H_max","p_D_max","rel_strain_transmission1","strain_seed_date1",
                 "phi_cases","alpha_cases","alpha_hosp","alpha_death"),
        initial = c(0.025,0.02,0.024,0.02,0.024,2.8,1,330,0.400626705/2,0.316,3.5,505,0.5,0.5,0.5,0.5),
        min = c(rep(0,n_beta),0.25,1,305,0,0,2,500,0,0,0,0),
        max = c(rep(1,n_beta),4,10,345,1,1,6,512,1,1,1,1),
        integer = c(rep(F,n_beta),F,F,F,F,F,F,F,F,F,F,F)
    )
    pars_info
}