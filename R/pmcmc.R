# 
# 
# pmcmc <- function(log_likelihood = NULL) {
#     if (is.null(log_likelihood)) {
#         log_likelihood <- calc_loglikelihood
#     }
# }
# 

mcmc <- function(transform,filter,curr_pars,priors,n_particles,n_iters,proposal_matrix,pars_min,pars_max){
    # Apply transform function to bind fixed parameter values
    transform_pars <- transform(curr_pars)
    
    # Calculate initial log-likelihood
    curr_ll <- filter$run(pars = transform_pars,save_history = T)
    
    # Calculate initial log prior density
    curr_lprior <- calc_lprior(pars = curr_pars,priors = priors)
    
    # Calculate log-posterior density
    curr_lpost <- curr_ll + curr_lprior
    
    # Get sample state
    particle_idx <- sample.int(n_particles,1)
    curr_ss <- filter$history(particle_idx)[, , ,drop=T]
    
    # Create arrays for storing output
    pars <- matrix(nrow = n_iters + 1,ncol = length(curr_pars))
    ll <- numeric(n_iters + 1)
    lprior <- numeric(n_iters + 1)
    lpost <- numeric(n_iters + 1)
    states <- array(dim = c(nrow(curr_ss),n_iters + 1,ncol(curr_ss)))
    
    # Record initial parameter values and log-likelihood 
    pars[1, ] <- curr_pars
    ll[1] <- curr_ll
    lprior[1] <- curr_lprior
    lpost[1] <- curr_lpost
    states[,1, ] <- curr_ss
    
    # Initialise acceptance counter
    acc <- 0
        
    for (iter in seq_len(n_iters) + 1){
        # Propose new parameter values
        # prop_pars <- rmvn(n = 1, mu = curr_pars, sigma = proposal_matrix)
        prop_pars <- mvrnorm(n = 1, mu = curr_pars, Sigma = proposal_matrix)
        if (all(prop_pars > pars_min & prop_pars < pars_max)){
            transform_pars <- transform(prop_pars)
            prop_ll <- filter$run(pars = transform_pars,save_history = T)
            prop_lprior <- calc_lprior(pars = prop_pars,priors = priors)
            prop_lpost <- prop_ll + prop_lprior
            particle_idx <- sample.int(n_particles,1)
            prop_ss <- filter$history(particle_idx)[, , ,drop=T]
            
            acc_prob <- prop_lpost - curr_lpost
            if (log(runif(1)) < acc_prob){
                # Update parameter values and log-likelihood
                curr_pars <- prop_pars
                curr_lprior <- prop_lprior
                curr_ll <- prop_ll
                curr_lpost <- prop_lpost
                curr_ss <- prop_ss
                acc <- acc + 1
            }
        }
        
        # Save samples
        pars[iter,] <- curr_pars
        ll[iter] <- curr_ll
        lprior[iter] <- curr_lprior
        lpost[iter] <- curr_lpost
        states[,iter, ] <- curr_ss
    }
    
    res <- list(pars = pars,
                ll = ll,
                lprior = lprior,
                lpost = lpost,
                states = states,
                acc = acc)
    return(res)
    
}

# placeholder function
calc_lprior <- function(pars,priors) {
    lp <- Map(function(p,value) priors[[p]](value),names(priors),pars)
    lprior <- sum(unlist(lp))
    return(lprior)
}

# propose_parameters <- function(pars, proposal_kernel, pars_discrete, pars_min, pars_max) {
#     ## proposed jumps are normal with mean pars and sd as input for parameter
#     proposed <- pars + rmvnorm(n = 1,  sigma = proposal_kernel)
#     for_chain <- reflect_proposal(x = proposed,
#                                   floor = pars_min,
#                                   cap = pars_max)
#     
#     # discretise if necessary
#     for_eval <- proposed
#     #for_eval[pars_discrete] <- round(for_eval[pars_discrete])
#     for_eval <- reflect_proposal(x = for_eval,
#                                  floor = pars_min,
#                                  cap = pars_max)
#     return(list(for_eval = for_eval, for_chain = for_chain))
# }