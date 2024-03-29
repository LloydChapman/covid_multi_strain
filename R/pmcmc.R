mcmc <- function(transform,filter,curr_pars,priors,n_particles,n_iters,idx,
                 scaling_factor_start,proposal_matrix,pars_min,pars_max,iter0,
                 discrete = rep(FALSE,length(curr_pars)),u = seq_along(curr_pars),
                 thinning = 1,rerun_every = Inf){
    # Get number of parameters being fitted
    n_pars <- length(curr_pars)
    n_u <- length(u)
    
    # Apply transform function to bind fixed parameter values
    transform_pars <- transform(curr_pars)
    
    # Calculate initial log-likelihood
    curr_ll <- filter$run(pars = transform_pars,save_history = T)
    
    # Calculate initial log prior density
    curr_lprior <- calc_lprior(pars = curr_pars,priors = priors)
    
    # Calculate log-posterior density
    curr_lpost <- curr_ll + curr_lprior
    
    # Get sample state
    idx_particle <- sample.int(n_particles,1)
    curr_state <- filter$state()[,idx_particle]
    curr_trajectories <- filter$history(idx_particle)[, , ,drop=T]
    
    # Create arrays for storing output
    n_smpls <- round(n_iters/thinning)
    pars_full <- matrix(nrow = n_iters + 1,ncol = n_pars)
    pars <- matrix(nrow = n_smpls + 1,ncol = n_pars)
    colnames(pars_full) <- colnames(pars) <- names(curr_pars)
    ll <- numeric(n_smpls + 1)
    lprior <- numeric(n_smpls + 1)
    lpost <- numeric(n_smpls + 1)
    state <- matrix(nrow = length(curr_state),ncol = n_smpls + 1)
    trajectories_state <- array(dim = c(nrow(curr_trajectories),n_smpls + 1,ncol(curr_trajectories)))
    dimnames(trajectories_state) <- list(names(idx$state))
    scaling_factor <- numeric(n_iters + 1)
    
    # Record initial parameter values and log-likelihood 
    pars_full[1, ] <- curr_pars
    pars[1, ] <- curr_pars
    ll[1] <- curr_ll
    lprior[1] <- curr_lprior
    lpost[1] <- curr_lpost
    state[,1] <- curr_state
    trajectories_state[,1, ] <- curr_trajectories
    
    mean_pars <- curr_pars
    
    # Initialise acceptance counter
    acc <- 0
    
    # Initialise scaling factor for proposal matrix
    scaling_factor[1] <- scaling_factor_start
    
    # Target acceptance rate
    a <- 0.234
    iter_start <- 5/(a*(1 - a))
    A <- -qnorm(a/2)
    delta <- (1 - 1/n_pars)*sqrt(2*pi)*exp(A^2/2)/(2*A) + 1/(n_pars*(a*(1-a)))
    
    prop_pars <- curr_pars
        
    for (iter in seq_len(n_iters)){
        # Rerun particle filter on current accepted state every rerun_every iterations
        if (iter %% rerun_every == 0){
            transform_pars <- transform(curr_pars)
            curr_ll <- filter$run(transform_pars,save_history = T)
            curr_lpost <- curr_ll + curr_lprior
            idx_particle <- sample.int(n_particles,1)
            curr_state <- filter$state()[,idx_particle]
            curr_trajectories <- filter$history(idx_particle)[, , ,drop=T]
        }
        
        # Propose new parameter values
        # prop_pars <- rmvn(n = 1, mu = curr_pars, sigma = proposal_matrix)
        prop_pars[u] <- mvrnorm(n = 1, mu = curr_pars[u], Sigma = 2.38^2/n_u*scaling_factor[iter]^2*proposal_matrix[u,u])
        prop_pars[discrete] <- round(prop_pars[discrete])
        if (all(prop_pars >= pars_min & prop_pars <= pars_max)){
            transform_pars <- transform(prop_pars)
            prop_ll <- filter$run(transform_pars,save_history = T)
            prop_lprior <- calc_lprior(prop_pars,priors)
            prop_lpost <- prop_ll + prop_lprior
            idx_particle <- sample.int(n_particles,1)
            prop_state <- filter$state()[,idx_particle]
            prop_trajectories <- filter$history(idx_particle)[, , ,drop=T]
            
            log_acc_prob <- min(0,prop_lpost - curr_lpost)
            if (log(runif(1)) < log_acc_prob){
                # Update parameter values and log-likelihood
                curr_pars <- prop_pars
                curr_lprior <- prop_lprior
                curr_ll <- prop_ll
                curr_lpost <- prop_lpost
                curr_state <- prop_state
                curr_trajectories <- prop_trajectories
                acc <- acc + 1
            }
        } else {
            log_acc_prob <- -Inf
        }
        scaling_factor[iter + 1] <- max(1,scaling_factor[iter] * exp(delta*(exp(log_acc_prob) - a)/(iter_start + iter))) #1 #
        if (abs(log(scaling_factor[iter + 1]) - log(scaling_factor_start)) > log(3)){
            scaling_factor_start <- scaling_factor[iter + 1]
            iter_start <- 5/(a*(1 - a)) - iter
        }
        
        pars_full[iter+1,] <- curr_pars
        
        # Save samples
        if (iter %% thinning == 0){
            smpl = iter/thinning
            pars[smpl+1,] <- curr_pars
            ll[smpl+1] <- curr_ll
            lprior[smpl+1] <- curr_lprior
            lpost[smpl+1] <- curr_lpost
            state[,smpl+1] <- curr_state
            trajectories_state[,smpl+1, ] <- curr_trajectories            
        }
        
        # Update empirical covariance matrix
        tmp <- update_mean_and_cov_Spencer(mean_pars,proposal_matrix,pars_full,iter,iter0,ncol(pars))
        mean_pars <- tmp$mean_new
        proposal_matrix <- tmp$cov_new
        
        if (iter %% 100 == 0){
            print(iter)
        }
    }
    
    probabilities <- cbind(log_prior = lprior,
                           log_likelihood = ll,
                           log_posterior = lpost)
    
    if (class(transform_pars) == "multistage_parameters"){
        dt <- transform_pars[[1]]$pars$dt
    } else {
        dt <- transform_pars$dt
    }
    trajectories <- list(
        step = filter$inputs()$data$step_end,
        rate = as.integer(1/dt),
        state = trajectories_state,
        predict = rep(F, nlayer(trajectories_state))
    )
    
    return(list(pars_full = pars_full,
                pars = pars,
                probabilities = probabilities,
                state = state,
                trajectories = trajectories,
                acc = acc,
                proposal_matrix = proposal_matrix,
                scaling_factor = scaling_factor))
}

# placeholder function
calc_lprior <- function(pars,priors) {
    lp <- Map(function(p,value) priors[[p]](value),names(priors),pars)
    lprior <- sum(unlist(lp))
    return(lprior)
}

update_mean_and_cov_Spencer <- function(mean_old,cov_old,x,i,i0,n_pars){
    # Add 1 to x row index as can't index from 0
    fi <- floor(i/2)
    if (i==1){ # if i==1
        mean_new <- (x[1,] + x[2,])/2
        cov_new <- ((i0+n_pars+1)*cov_old + x[1,]%*%t(x[1,]) + x[2,]%*%t(x[2,]) - 2*mean_new%*%t(mean_new))/(i0+n_pars+3)
    }
    else if (fi==floor((i-1)/2)+1){ # if i is even, replace the oldest observation in the estimates of the mean and covariance with the new observation
        mean_new <- mean_old + (x[i+1,]-x[fi,])/(i-fi+1)
        cov_new <- cov_old + (x[i+1,]%*%t(x[i+1,]) - x[fi,]%*%t(x[fi,]) + (i-fi+1)*(mean_old%*%t(mean_old) - mean_new%*%t(mean_new)))/(i-fi+i0+n_pars+2)
    } else { # if i is odd, add the new observation to the estimates of the mean and covariance
        mean_new <- ((i-fi)*mean_old + x[i+1,])/(i-fi+1)
        cov_new <- ((i-fi+i0+n_pars+1)*cov_old + x[i+1,]%*%t(x[i+1,])+(i-fi)*(mean_old%*%t(mean_old)) - (i-fi+1)*(mean_new%*%t(mean_new)))/(i-fi+i0+n_pars+2)
    }
    return(list(mean_new=mean_new,cov_new=cov_new))  
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