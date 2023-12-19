run_simulations <- function(run,sim_run,assumptions,deterministic,Rt,n_smpls){
    ## Set MCMC output
    output <- paste0("output/MCMCoutput",run,".RDS")
    
    ## Run counterfactual simulations
    # 1: No lockdowns
    # 2: 1st lockdown 2 weeks earlier
    # 3: 1st lockdown 2 weeks later
    # 4: 2nd lockdown 2 weeks earlier
    # 5: 2nd lockdown 2 weeks later
    # 6: Both lockdowns 2 weeks earlier
    # 7: Both lockdowns 2 weeks later
    # 8: No vaccination
    # 9: Boosters 1 month earlier
    # 10: No boosters
    
    # Set random number seed
    seed <- 1L
    
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
    
    # Get minimum ages of age groups
    min_ages <- get_min_age(pars$base$age_groups)
    
    initial_date <- pars$info$min[pars$info$name == "start_date"] - 1
    
    # Set probabilities for quantiles for outcomes
    probs <- c(0.025,0.5,0.975)
    
    ## Lockdown counterfactuals
    beta_date <- pars$base$beta_date
    beta_idx_list <- replicate(10,seq_along(beta_date),F)
    beta_date_cntfctl_list <- replicate(10,beta_date,F) #replicate(10,beta_date,F)
    
    # Remove lockdowns
    beta_idx_list[[1]] <- c(1,3,5)
    
    # Move first lockdown two weeks earlier
    beta_date_cntfctl_list[[2]][1:2] <- beta_date[1:2] - 14
    # Move first lockdown two weeks later
    beta_date_cntfctl_list[[3]][1:2] <- beta_date[1:2] + 14
    # Move second lockdown two weeks earlier
    beta_date_cntfctl_list[[4]][3:4] <- beta_date[3:4] - 14
    # Move second lockdown two weeks later
    beta_date_cntfctl_list[[5]][3:4] <- beta_date[3:4] + 14
    # Move both lockdowns earlier
    beta_date_cntfctl_list[[6]][1:2] <- beta_date[1:2] - 14
    beta_date_cntfctl_list[[6]][3:4] <- beta_date[3:4] - 14
    # Move both lockdowns later
    beta_date_cntfctl_list[[7]][1:2] <- beta_date[1:2] + 14
    beta_date_cntfctl_list[[7]][3:4] <- beta_date[3:4] + 14
    # Use actual dates for vaccination counterfactual simulations
    
    ## Vaccination counterfactuals
    # Set counterfactual vaccine schedule
    schedule <- pars$base$vaccine_schedule
    schedule_cntfctl_list <- replicate(10,schedule,F) #replicate(10,schedule,F)
    # No vaccination in counterfactual
    schedule_cntfctl_list[[8]]$doses <- array(0,dim = dim(schedule$doses))
    # Boosters starting 1 month earlier
    schedule_cntfctl_list[[9]] <- change_booster_timing(schedule, 30)
    # # Boosters starting 2 weeks later
    # schedule_cntfctl_list[[9]] <- change_booster_timing(schedule, -180)
    # No boosters
    schedule_cntfctl_list[[10]]$doses[,3,] <- matrix(0, nrow = nrow(schedule$doses), 
                                                     ncol = nlayer(schedule$doses))
    
    # Create objects for storing quantiles of simulation output
    q_outcomes_list <- vector("list",length(beta_date_cntfctl_list))
    q_outcomes_cntfctl_list <- vector("list",length(beta_date_cntfctl_list))
    q_outcomes_averted_list <- vector("list",length(beta_date_cntfctl_list))
    q_prop_outcomes_averted_list <- vector("list",length(beta_date_cntfctl_list))
    q_total_outcomes_list <- vector("list",length(beta_date_cntfctl_list))
    q_total_outcomes_cntfctl_list <- vector("list",length(beta_date_cntfctl_list))
    q_total_outcomes_averted_list <- vector("list",length(beta_date_cntfctl_list))
    q_prop_total_outcomes_averted_list <- vector("list",length(beta_date_cntfctl_list))
    
    # new_infections_cntfctl_by_vacc_list <- vector("list",length(beta_date_cntfctl_list))
    # new_reinfections_cntfctl_by_vacc_list <- vector("list",length(beta_date_cntfctl_list))
    
    # Run simulations
    for (i in seq_along(beta_date_cntfctl_list)){ #9:10){ # 
        out <- simulate_counterfactual(output,n_smpls,beta_date_cntfctl_list[[i]],beta_idx_list[[i]],
                                       schedule_cntfctl_list[[i]],initial_date,deterministic,seed = seed,min_ages = min_ages, Rt = Rt)
        # Extract state counts for counterfactual and actual scenarios
        # Drop initial values (as they do not correspond to the time step before the first time step in the data)
        states_cntfctl <- out$states_cntfctl[,,-1]
        smpl <- out$smpl
        info <- out$info
        states <- out$states[,,-1]
        dates <- out$dates
        
        # Add row names to states_cntfctl array
        dimnames(states_cntfctl) <- list(names(index(info,min_ages,Rt)$state))
        
        # Calculate pairwise differences in simulated outcomes across time
        outcomes_averted <- states - states_cntfctl
        prop_outcomes_averted <- outcomes_averted/states_cntfctl
        # dimnames(outcomes_averted[[i]])[[1]] <- names(index(info)$state)
        
        ## Calculate total differences in outcomes over each wave
        # Set wave dates
        wave_date <- c(covid_multi_strain_date_as_date(min(dates)),
                       as.Date(c("2021-06-11","2021-11-21")),
                       covid_multi_strain_date_as_date(max(dates)))
        wave_date <- covid_multi_strain_date(wave_date) - min(dates) + 1
        
        q_outcomes_list[[i]] <- calculate_outcome_quantiles(states,dates,info,min_ages,Rt)
        q_outcomes_cntfctl_list[[i]] <- calculate_outcome_quantiles(states_cntfctl,dates,info,min_ages,Rt)
        q_outcomes_averted_list[[i]] <- calculate_outcome_quantiles(outcomes_averted,dates,info,min_ages,Rt)
        q_prop_outcomes_averted_list[[i]] <- calculate_outcome_quantiles(prop_outcomes_averted,dates,info,min_ages,Rt)
        
        q_total_outcomes_list[[i]] <- calculate_outcomes_by_wave(states,wave_date,info,min_ages,Rt)$q_total
        tmp1 <- calculate_outcomes_by_wave(states_cntfctl,wave_date,info,min_ages,Rt)
        q_total_outcomes_cntfctl_list[[i]] <- tmp1$q_total
        tmp2 <- calculate_outcomes_by_wave(outcomes_averted,wave_date,info,min_ages,Rt)
        q_total_outcomes_averted_list[[i]] <- tmp2$q_total
        tmp3 <- cbind(tmp1$total[,1],tmp2$total[,-1][, Map(`/`, .SD, tmp1$total[,-1])])
        cols <- names(index(info,min_ages = min_ages,Rt = Rt)$state)
        q_prop_total_outcomes_averted_list[[i]] <- tmp3[,unlist(
            lapply(.SD,function(x) list(q95l = quantile(x,probs = 0.025,na.rm = T),
                                        med = quantile(x,probs = 0.5,na.rm = T),
                                        q95u = quantile(x,probs = 0.975,na.rm = T))),
            recursive = F),.SDcols = cols,by = .(wave)]
        
        rm(out,states,states_cntfctl,outcomes_averted,prop_outcomes_averted)
        gc()
        # res <- calculate_new_infections_by_vacc(states_cntfctl)
        # new_infections_cntfctl_by_vacc_list[[i]] <- res$new_infections_by_vacc
        # new_reinfections_cntfctl_by_vacc_list[[i]] <- res$new_reinfections_by_vacc
    }
    
    q_outcomes <- rbindlist(q_outcomes_list, idcol = "cntfctl")
    q_outcomes_cntfctl <- rbindlist(q_outcomes_cntfctl_list, idcol = "cntfctl")
    q_outcomes_averted <- rbindlist(q_outcomes_averted_list, idcol = "cntfctl")
    q_prop_outcomes_averted <- rbindlist(q_prop_outcomes_averted_list, idcol = "cntfctl")
    q_outcomes[,date := covid_multi_strain_date_as_date(day)]
    q_outcomes_cntfctl[,date := covid_multi_strain_date_as_date(day)]
    q_outcomes_averted[,date := covid_multi_strain_date_as_date(day)]
    q_prop_outcomes_averted[,date := covid_multi_strain_date_as_date(day)]
    
    q_total_outcomes <- rbindlist(q_total_outcomes_list, idcol = "cntfctl")
    q_total_outcomes_cntfctl <- rbindlist(q_total_outcomes_cntfctl_list, idcol = "cntfctl")
    q_total_outcomes_averted <- rbindlist(q_total_outcomes_averted_list, idcol = "cntfctl")
    q_prop_total_outcomes_averted <- rbindlist(q_prop_total_outcomes_averted_list, idcol = "cntfctl")
    
    # new_infections_cntfctl_by_vacc <- rbindlist(new_infections_cntfctl_by_vacc_list, idcol = "cntfctl")
    # new_reinfections_cntfctl_by_vacc <- rbindlist(new_reinfections_cntfctl_by_vacc_list, idcol = "cntfctl")
    
    ttls <- c("No change in lockdown dates","No lockdowns",
              "1st lockdown 2 weeks earlier","1st lockdown 2 weeks later",
              "2nd lockdown 2 weeks earlier","2nd lockdown 2 weeks later",
              "Both lockdowns 2 weeks earlier","Both lockdowns 2 weeks later",
              "No vaccination", "Boosters 1 month earlier","No boosters")#,"Boosters")
    names(ttls) <- as.character(seq_along(ttls)-1)
    
    # pdf("output/infctns_by_vacc_status.pdf",)
    # for (i in 9:10){
    #     idx_wave3 = wave_date[3]:wave_date[4]
    #     dates <- seq.Date(covid_multi_strain_date_as_date(wave_date[3]),covid_multi_strain_date_as_date(wave_date[4]),by = 1)
    #     dates_plot <- seq.Date(covid_multi_strain_date_as_date(wave_date[3]),covid_multi_strain_date_as_date(wave_date[4]),by = 30)
    #     matplot(dates,t(apply(new_infections_cntfctl_by_vacc_list[[i]][,,idx_wave3],c(1,3),median)),type = "l",xaxt = "n",xlab = "Date",ylab = "Infections", ylim = c(0,550), main = ttls[names(ttls) == i])
    #     nc <- nrow(new_infections_cntfctl_by_vacc_list[[i]])
    #     ltyp <- rep(1:5, times=nc/5, each=1) # the line types matplot uses
    #     cols <- rep(1:6, times=round(nc/6), each=1) # the cols matplot uses 
    #     legend("topright", legend = c("Unvacc","1st dose","2nd dose","Waned","Booster"), lty = ltyp, col = cols)
    #     axis(1,dates_plot,format(dates_plot,"%Y-%m-%d"))
    # }
    # 
    # for (i in 9:10){    
    #     matplot(dates,t(apply(new_reinfections_cntfctl_by_vacc_list[[i]][,,idx_wave3],c(1,3),median)),type = "l",xaxt = "n",xlab = "Date",ylab = "Reinfections", ylim = c(0,550), main = ttls[names(ttls) == i])
    #     nc <- nrow(new_reinfections_cntfctl_by_vacc_list[[i]])
    #     ltyp <- rep(1:5, times=nc/5, each=1) # the line types matplot uses
    #     cols <- rep(1:6, times=round(nc/6), each=1) # the cols matplot uses 
    #     legend("topright", legend = c("Unvacc","1st dose","2nd dose","Waned","Booster"), lty = ltyp, col = cols) 
    #     axis(1,dates_plot,format(dates_plot,"%Y-%m-%d"))
    # }
    # dev.off()
    # 
    # pdf("output/cntfctl_infctns_booster6.pdf")
    # plot(dates,colSums(apply(new_infections_cntfctl_by_vacc_list[[9]][,,idx_wave3],c(1,3),median)),type = "l",col = "red",xaxt = "n",xlab = "Date",ylab = "Infections",ylim = c(0,1000))
    # lines(dates,colSums(apply(new_infections_cntfctl_by_vacc_list[[10]][,,idx_wave3],c(1,3),median)))
    # legend("topright",legend = c("Boosters", "No boosters"),lty = 1,col = c("black", "red"))
    # axis(1,dates_plot,format(dates_plot,"%Y-%m-%d"))
    # 
    # plot(dates,colSums(apply(new_reinfections_cntfctl_by_vacc_list[[9]][,,idx_wave3],c(1,3),median)),type = "l",col = "red",xaxt = "n",xlab = "Date",ylab = "Reinfections",ylim = c(0,1000))
    # lines(dates,colSums(apply(new_reinfections_cntfctl_by_vacc_list[[10]][,,idx_wave3],c(1,3),median)))
    # legend("topright",legend = c("Boosters", "No boosters"),lty = 1, col = c("black", "red"))
    # axis(1,dates_plot,format(dates_plot,"%Y-%m-%d"))
    # dev.off()
    
    
    # Bind non-counterfactual and counterfactual total outcomes
    tmp <- q_total_outcomes[cntfctl == 1]
    tmp[,cntfctl := 0]
    total_outcomes <- rbind(tmp,q_total_outcomes_cntfctl)
    
    # Create table of total outcomes by wave and overall for different counterfactuals
    tbl <- write_table(total_outcomes,ttls,f_cases = 0.001)
    write.csv(tbl[c(1,2,9,11),],paste0("output/table1_",sim_run,".csv"),row.names = F)
    
    write.csv(tbl[c(1,9,11),],paste0("output/table1_vax_",sim_run,".csv"),row.names = F)
    write.csv(tbl[c(1,3:8),],paste0("output/table1_lockdown_",sim_run,".csv"),row.names = F)
    
    # Create table of total outcomes averted by wave and overall for different counterfactuals
    tbl1 <- write_table(q_total_outcomes_averted,ttls,f_cases = 0.001,d_cases = 1,d_hosps = 0,d_deaths = 0,method = "round")
    write.csv(tbl1, paste0("output/total_outcomes_averted_",sim_run,".csv"),row.names = F)
    
    # Create table of proportion of total outcomes averted by wave and overall for different counterfactuals
    tbl2 <- write_table(q_prop_total_outcomes_averted,ttls,method = "round")
    write.csv(tbl2,paste0("output/prop_total_outcomes_averted_",sim_run,".csv"),row.names = F)
    
    # Plot counterfactuals
    plot_simulations(q_outcomes,q_outcomes_cntfctl,q_total_outcomes_averted,ttls,sim_run,cols,dates)
    
    # Save output
    res <- list(q_outcomes = q_outcomes,
                q_outcomes_cntfctl = q_outcomes_cntfctl,
                q_outcomes_averted = q_outcomes_averted,
                q_prop_outcomes_averted = q_prop_outcomes_averted,
                q_total_outcomes = q_total_outcomes,
                q_total_outcomes_cntfctl = q_total_outcomes_cntfctl,
                q_total_outcomes_averted = q_total_outcomes_averted,
                q_prop_total_outcomes_averted = q_prop_total_outcomes_averted)
    saveRDS(res,paste0("output/cntfctl_output",sim_run,".RDS"))
}
