simulate <- function(gen_mod, p, n_steps, deterministic = FALSE, 
                     keep_all_states = TRUE, p1 = NULL, 
                     n_steps1 = NULL, transform = NULL){
    # Create instance of model
    mod <- gen_mod$new(p,step = 0,n_particles = 1,n_threads = 1,seed = 1,
                       deterministic = deterministic)
    
    info <- mod$info()
    initial_state <- initial(info, NULL, p)
    
    mod$update_state(state = initial_state)
    
    # Create an array to contain outputs after looping the model.
    x <- array(NA, dim = c(info$len, 1, n_steps+1))
    
    # For loop to run the model iteratively
    x[ , ,1] <- mod$state()
    for (t in seq_len(n_steps)) {
        x[ , ,t+1] <- mod$run(t)
    }
    
    if (!is.null(p1)){
        # Apply transform function to model state
        state <- mod$state()
        state1 <- transform(state,info)
        
        # Update model parameters and state
        mod$update_state(pars = p1,state = state1)
        
        # Create data to be fitted to
        x1 <- array(NA, dim = c(info$len, 1, n_steps1-n_steps+1))
        
        # For loop to run the model iteratively
        x1[ , ,1] <- mod$state()
        for (t in seq_len(n_steps1-n_steps)) {
            # print(t)
            x1[ , ,t+1] <- mod$run(n_steps+t)
        }
        
        # Join with first epoch
        x <- array_bind(x,x1)        
    }
    
    if (!keep_all_states){
        idx <- index(info)
        x <- x[idx$state, , ,drop = FALSE]
    }
    
    return(list(x = x, info = info))
}


plot_trajectories <- function(time,x,n_age,n_strains,n_vax){
    # Plot trajectories
    for (k in 1:n_vax){
        for (j in 1:n_strains) {
            par(mfrow = c(2,4), oma=c(2,3,0,0))
            for (i in 1:n_age) {
                par(mar = c(3, 4, 2, 0.5))
                cols <- c(S = "#8c8cd9", E = "#ffff00", I_P = "#cc0044", I_A = "green", I_C = "blue", R = "#999966", D = "#000000")
                matplot(time, x[29 + i + (k-1)*n_age, ,-1], type = "l", # Offset to access numbers in age compartment
                        xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
                        col = cols[["S"]], lty = 1, ylim=range(x[-(1:29),,]))
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains, ,-1], col = cols[["E"]], lty = 1)
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains + n_age*n_strains*n_vax, ,-1], col = cols[["I_P"]], lty = 1)
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 2*n_age*n_strains*n_vax, ,-1], col = cols[["I_A"]], lty = 1)
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 3*n_age*n_strains*n_vax, ,-1], col = cols[["I_C"]], lty = 1)
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 4*n_age*n_strains*n_vax, ,-1], col = cols[["R"]], lty = 1)
                matlines(time, x[93 + i + (j-1)*n_age + (k-1)*n_age*n_strains + 7*n_age*n_strains*n_vax, ,-1], col = cols[["D"]], lty = 1)
                legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
                axis(2, las = 2)
            }
            mtext("Number of individuals", side=2, line=1, outer=T, las=0)
            mtext("Time", side = 1, line = 0, outer =T)    
        }    
    }
}

## Run counterfactual simulations
simulate_counterfactual <- function(output,n_smpls,beta_date_cntfctl,schedule_cntfctl,burnin = NULL,seed = 1){
    # Load MCMC output
    load(output)
    
    if (is.null(burnin)){
        burnin <- round((n_iters/thinning)/10)
    }
    # Remove burn-in
    pars <- res$pars[-(1:burnin),]
    
    out <- vector("list", n_smpls)
    set.seed(seed)
    smpl <- sample.int(nrow(pars),n_smpls)
    for (i in seq_len(n_smpls)){
        j <- smpl[i] #sample.int(nrow(pars),1)
        pars_i <- pars[j,]
        p_i <- parameters(dt,
                          n_age,
                          n_vax,
                          transmission,
                          beta_date = beta_date_cntfctl,
                          beta_value = pars_i[1:4],
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
                          pars_i[8]*p_H,
                          p_G,
                          pars_i[9]*p_D,
                          p_P_1,
                          population,
                          start_date = pars_i[6],
                          initial_seed_size,
                          initial_seed_pattern,
                          strain_transmission = c(1,pars_i[5]),
                          strain_seed_date = pars_i[7],
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
                          vaccine_schedule = schedule_cntfctl,
                          vaccine_index_dose2,
                          vaccine_index_booster,
                          vaccine_catchup_fraction,
                          n_doses,
                          waning_rate,
                          cross_immunity,
                          sero_sensitivity_1,
                          sero_specificity_1)
        out[[i]] <- simulate(covid_multi_strain, p_i, n_steps, 
                             deterministic, keep_all_states = F)
        smpl[i] <- j
    }
    
    info <- out[[1]]$info
    states_cntfctl <- lapply(out,"[[",1)
    
    # Extract daily values
    states_cntfctl <- lapply(states_cntfctl,function(y) y[,,seq(1,nlayer(y),by = 1/dt),drop = F])
    
    states_cntfctl <- abind(states_cntfctl,along = 2)
    
    # Extract corresponding posterior samples of trajectories
    states <- res$states[,burnin + smpl,]
    
    return(list(states_cntfctl = states_cntfctl,states = states,smpl = smpl,info = info))
}

# calculate_quantiles <- function(x, probs = c(0.025,0.5,0.975)){
#     q <- apply(x,1,function(y) quantile(y, probs = probs))
#     return(q)
# }

calculate_outcome_quantiles <- function(x,info){
    # If input is a 3-D array convert to a data table
    if (length(dim(x)) == 3){
        x <- arr_to_dt(x,info)
    }
    
    out <- x[,.(med = median(value),
                q95l = quantile(value,0.025),
                q95u = quantile(value,0.975)), by = .(state,day)]
    return(out)
}

calculate_outcomes_by_wave <- function(x,wave_date,info){
    cols <- names(index(info)$state)
    dimnames(x)[[1]] <- cols
    n_waves <- length(wave_date) - 1
    
    tmp <- vector("list",n_waves)
    for (j in 1:n_waves){
        times <- wave_date[j]:wave_date[j+1]
        # Calculate total outcomes averted over wave
        tmp[[j]] <- as.data.table(t(apply(x[,,times], c(1,2), sum)))
        tmp[[j]][,smpl := .I]
    }
    total_x_waves <- rbindlist(tmp,idcol = "wave")
    total_x_both_waves <- total_x_waves[,lapply(.SD,sum),.SDcols = cols,by = .(smpl)]
    total_x_both_waves[,wave := "Total"]
    total_x <- rbind(total_x_waves,total_x_both_waves)
    
    # Calculate quantiles
    q_total_x <- total_x[,unlist(
        lapply(.SD,function(x) list(q95l = quantile(x,probs = 0.025),
                                    med = quantile(x,probs = 0.5),
                                    q95u = quantile(x,probs = 0.975))),
        recursive = F),.SDcols = cols,by = .(wave)]
    
    return(q_total_x)
}


plot_counterfactuals <- function(q_outcomes,q_outcomes_cntfctl,outcome,ylbl,ttls){
    clrs = c("Fitted" = "black","Counterfactual" = "darkgreen")
    p <- ggplot() + 
        geom_line(aes(x = date,y = med,color = "Fitted"),q_outcomes[state == outcome]) +
        geom_ribbon(aes(x = date,ymin = q95l,ymax = q95u,fill = "Fitted"),q_outcomes[state == outcome],alpha = 0.5) +
        geom_line(aes(x = date,y = med,color = "Counterfactual"),q_outcomes_cntfctl[state == outcome]) + 
        geom_ribbon(aes(x = date,ymin = q95l,ymax = q95u,fill = "Counterfactual"),q_outcomes_cntfctl[state == outcome],alpha = 0.5) +
        labs(x = "Day",y = ylbl) +
        scale_color_manual(name = "",values = clrs) +
        scale_fill_manual(name = "",values = clrs) +
        facet_wrap(~cntfctl,nrow = 2,labeller = labeller(cntfctl = ttls),dir = "v") +
        theme(legend.position = "bottom")
    return(p)
}


arr_to_dt <- function(x,info){
    x_long <- as.data.table(x)
    
    names(x_long)[1:3] <- c("state","smpl","day")
    
    x_long[,state := names(index(info)$state[state])]
    
    return(x_long)
}


med_and_CI = function(x,l,u,f=1,d=1,method="round"){
    if (method=="signif"){
        paste0(signif(f*x,d)," (",signif(f*l,d),"-",signif(f*u,d),")")
    } else if (method=="round"){
        paste0(round(f*x,d)," (",round(f*l,d),"-",round(f*u,d),")")
    }
}
