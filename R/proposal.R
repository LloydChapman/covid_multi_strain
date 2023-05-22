create_proposal <- function(beta_date){
    n_beta <- length(beta_date)
    nms <- c(paste0(rep("beta",n_beta),seq_len(n_beta)),
             "rel_strain_transmission","start_date","strain_seed_date",
             "p_H_max","p_D_max","rel_strain_transmission1","strain_seed_date1",
             "rel_strain_transmission2","strain_seed_date2",
             "phi_cases","alpha_cases","alpha_hosp","alpha_death")
    proposal_matrix <- 0.1*diag(c(rep(1e-7,n_beta),0.1^2,rep(2^2,2),rep(1e-5,2),0.1^2,2^2,0.1^2,2^2,1e-5,1e-5,1e-5,1e-5))
    colnames(proposal_matrix) <- nms
    proposal <- cbind(names = nms, proposal_matrix)
}