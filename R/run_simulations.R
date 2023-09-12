run_simulations <- function(run,sim_run,assumptions){
    ## Set MCMC output
    # Change run number for different assumption on booster waning rate
    # run <- 130 #77
    # run <- 78
    # run <- 80
    output <- paste0("output/MCMCoutput",run,".RDS")
    
    ## Run counterfactual simulations
    # Set run number
    # sim_run <- 28 #18
    # sim_run <- 19
    # sim_run <- 20
    # Set whether model is deterministic or stochastic
    deterministic <- T 
    # Set whether states required to calculate Rt have been output 
    Rt <- T
    # Set number of parameter samples and burn-in to remove
    n_smpls <- 500
    # burnin <- 1500
    seed <- 1L
    
    # Set assumption for booster waning rate
    # assumptions <- "central" #-log(67.7/82.8)/(105-25) # (Stowe Nat Comm 2022 Table S11)
    # assumptions <- "optimistic" # -log(0.923)/140 (Barnard Nat Com 2022 Table S4)
    # assumptions <- "pessimistic"
    
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
    
    # Move first lockdown a week earlier
    beta_date_cntfctl_list[[2]][1:2] <- beta_date[1:2] - 14
    # Move first lockdown a week later
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
    
    
    
    tbl <- total_outcomes[,.(Counterfactual = ttls[match(cntfctl,names(ttls))],
                             Wave = wave,
                             Cases = med_and_CI(cases.med,cases.q95l,cases.q95u,f = 0.001,d = 3,method = "signif"),
                             Hospitalisations = med_and_CI(hosps.med,hosps.q95l,hosps.q95u,d = 3,method = "signif"),
                             `Hospital deaths` = med_and_CI(deaths.med,deaths.q95l,deaths.q95u,d = 3,method = "signif"))]
    tbl[,Counterfactual := factor(Counterfactual, levels = unique(Counterfactual))]
    tbl <- dcast(tbl, Counterfactual ~ Wave, value.var = c("Cases","Hospitalisations","Hospital deaths"))
    write.csv(tbl[c(1,2,9,11),],paste0("output/table1_",sim_run,".csv"),row.names = F)
    
    write.csv(tbl[c(1,9,11),],paste0("output/table1_vax_",sim_run,".csv"),row.names = F)
    write.csv(tbl[c(1,3:8),],paste0("output/table1_lockdown_",sim_run,".csv"),row.names = F)
    
    tbl1 <- q_total_outcomes_averted[,.(Counterfactual = ttls[match(cntfctl,names(ttls))],
                                        Wave = wave,
                                        Cases = med_and_CI(cases.med,cases.q95l,cases.q95u,f = 0.001,d = 1,method = "round"),
                                        Hospitalisations = med_and_CI(hosps.med,hosps.q95l,hosps.q95u,d = 0,method = "round"),
                                        `Hospital deaths` = med_and_CI(deaths.med,deaths.q95l,deaths.q95u,d = 0,method = "round"))]
    tbl1[,Counterfactual := factor(Counterfactual, levels = unique(Counterfactual))]
    tbl1 <- dcast(tbl1, Counterfactual ~ Wave, value.var = c("Cases","Hospitalisations","Hospital deaths"))
    write.csv(tbl1, paste0("output/total_outcomes_averted_",sim_run,".csv"),row.names = F)
    
    tbl2 <- q_prop_total_outcomes_averted[,.(Counterfactual = ttls[match(cntfctl,names(ttls))],
                                             Wave = wave,
                                             Cases = med_and_CI(cases.med,cases.q95l,cases.q95u,d = 3,method = "round"),
                                             Hospitalisations = med_and_CI(hosps.med,hosps.q95l,hosps.q95u,d = 3,method = "round"),
                                             `Hospital deaths` = med_and_CI(deaths.med,deaths.q95l,deaths.q95u,d = 3,method = "round"))]
    tbl2[,Counterfactual := factor(Counterfactual, levels = unique(Counterfactual))]
    tbl2 <- dcast(tbl2, Counterfactual ~ Wave, value.var = c("Cases","Hospitalisations","Hospital deaths"))
    write.csv(tbl2,paste0("output/prop_total_outcomes_averted_",sim_run,".csv"),row.names = F)
    
    # Plot counterfactuals
    
    # No lockdown counterfactual
    idx <- 1
    # Cases
    p_cases_lockdown <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"cases","Cases",ttls[names(ttls) %in% idx])
    # Hospitalisations
    p_hosps_lockdown <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"hosps","Hospitalisations",ttls[names(ttls) %in% idx])
    # Deaths in hospital
    p_deaths_lockdown <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"deaths","Hospital deaths",ttls[names(ttls) %in% idx])
    pp0 <- plot_grid(p_cases_lockdown + theme(legend.position = "none"),
                     p_hosps_lockdown + theme(legend.position = "none"),
                     p_deaths_lockdown + theme(legend.position = "none"), nrow = 1)
    l0 <- get_legend(p_deaths_lockdown)
    ggsave(paste0("output/cntfctl_cases_hosps_and_deaths_lockdown",sim_run,".pdf"),plot_grid(pp0,l0,nrow = 2,rel_heights = c(1,0.1)),width = 9,height = 3)
    
    # Alternative lockdown timing counterfactuals
    idx <- 2:7
    # Cases
    p_cases <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"cases","Cases",ttls[names(ttls) %in% idx])
    ggsave(paste0("output/cntfctl_cases",sim_run,".pdf"),p_cases,width = 8,height = 5.2)
    # Hospitalisations
    p_hosps <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"hosps","Hospitalisations",ttls[names(ttls) %in% idx])
    ggsave(paste0("output/cntfctl_hosps",sim_run,".pdf"),p_hosps,width = 8,height = 5.2)
    # Deaths in hospital
    p_deaths <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"deaths","Hospital deaths",ttls[names(ttls) %in% idx])
    ggsave(paste0("output/cntfctl_deaths",sim_run,".pdf"),p_deaths,width = 8,height = 5.2)
    
    # No vaccination counterfactual
    idx <- 8
    p_cases_vax <- plot_counterfactuals(q_outcomes[cntfctl == idx],q_outcomes_cntfctl[cntfctl == idx],"cases","Cases",ttls[names(ttls) %in% idx])
    p_hosps_vax <- plot_counterfactuals(q_outcomes[cntfctl == idx],q_outcomes_cntfctl[cntfctl == idx],"hosps","Hospitalisations",ttls[names(ttls) %in% idx])
    # ggsave("output/cntfctl_hosps_vax.pdf",p_hosps_vax,width = 8,height = 6)
    p_deaths_vax <- plot_counterfactuals(q_outcomes[cntfctl == idx],q_outcomes_cntfctl[cntfctl == idx],"deaths","Hospital deaths",ttls[names(ttls) %in% idx])
    pp <- plot_grid(p_cases_vax + theme(legend.position = "none"),
                    p_hosps_vax + theme(legend.position = "none"),
                    p_deaths_vax + theme(legend.position = "none"), nrow = 1)
    l <- get_legend(p_deaths_vax)
    ggsave(paste0("output/cntfctl_hosps_and_deaths_vax",sim_run,".pdf"),plot_grid(pp,l,nrow = 2,rel_heights = c(1,0.1)),width = 9,height = 3)
    
    # No boosters counterfactual
    idx <- 10 #8:9
    p_cases_booster <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"cases","Cases",ttls[names(ttls) %in% idx])
    p_hosps_booster <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"hosps","Hospitalisations",ttls[names(ttls) %in% idx])
    # ggsave("output/cntfctl_hosps_vax.pdf",p_hosps_vax,width = 8,height = 6)
    p_deaths_booster <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"deaths","Hospital deaths",ttls[names(ttls) %in% idx])
    pp1 <- plot_grid(p_cases_booster + theme(legend.position = "none"),
                     p_hosps_booster + theme(legend.position = "none"),
                     p_deaths_booster + theme(legend.position = "none"), nrow = 1)
    l1 <- get_legend(p_deaths_booster)
    ggsave(paste0("output/cntfctl_hosps_and_deaths_booster",sim_run,".pdf"),plot_grid(pp1,l1,nrow = 2,rel_heights = c(1,0.1)),width = 9,height = 3)
    
    # pp2 <- plot_grid(pp0,pp,pp1,l1,nrow = 4,ncol = 1,rel_heights = c(1,1,1,0.1),width = 9, height = 10)
    outcome <- c("cases","hosps","deaths")
    ttls1 <- c("Cases","Hospitalisations","Hospital deaths")
    names(ttls1) <- outcome
    p_cntfctl <- plot_counterfactuals_together(q_outcomes[cntfctl %in% c(1,8,10)],q_outcomes_cntfctl[cntfctl %in% c(1,8,10)],outcome,ttls1,ttls[names(ttls) %in% c(1,8,10)])
    ggsave(paste0("output/cntfctl_cases_hosps_and_deaths",sim_run,".pdf"),p_cntfctl,width = 9,height = 3)
    
    # Plot numbers of cases, hospitalisations and deaths averted in each wave and overall
    p_dt <- melt(q_total_outcomes_averted,measure.vars = patterns(c("med","q95l","q95u")),value.name = c("med","q95l","q95u"))
    p_dt[,variable := cols[variable]]
    p_dt <- p_dt[(variable %in% outcome) & (cntfctl %in% c(1,8,10))]
    p_dt[,variable := factor(variable,levels = outcome)]
    lbls <- c(`1` = "Lockdowns",`8` = "Vaccination",`10` = "Boosters")
    p_dt[,cntfctl := lbls[match(cntfctl,names(lbls))]]
    p_dt[,cntfctl := factor(cntfctl,levels = unique(cntfctl))]
    ggplot(p_dt,aes(x = wave,fill = factor(cntfctl))) +
        geom_bar(aes(y = -med),position = position_dodge(),stat = "identity") + 
        geom_errorbar(aes(ymin = -q95l,ymax = -q95u),width = 0.2,position = position_dodge(0.9)) +
        labs(x = "Wave",y = "Number averted",fill = "Counterfactual") +
        facet_wrap(~variable,scales = "free",labeller = labeller(variable = ttls1)) +
        theme_cowplot(font_size = 12) +
        theme(legend.position = "bottom",
              strip.background = element_blank())
    ggsave(paste0("output/cases_hosps_and_deaths_averted",sim_run,".pdf"),width = 8, height = 3.5)
    
    # Same plot but with date rather than wave number on the x-axis
    wave_mid_date <- as.Date(c("2020-11-01","2021-08-15","2022-02-01","2022-07-01"))
    names(wave_mid_date) <- c("1","2","3","Total")
    p_dt[,date := wave_mid_date[match(wave,names(wave_mid_date))]]
    p <- ggplot(p_dt,aes(x = date,fill = factor(cntfctl))) +
        geom_bar(aes(y = -med),position = position_dodge(width = 70),width = 70,stat = "identity") + 
        geom_errorbar(aes(ymin = -q95l,ymax = -q95u),position = position_dodge(width = 70),width = 50) +
        labs(x = "Wave",y = "Change",fill = "Counterfactual") +
        xlim(covid_multi_strain_date_as_date(min(dates)),covid_multi_strain_date_as_date(max(dates)+100)) +
        facet_wrap(~variable,scales = "free",labeller = labeller(variable = ttls1)) +
        theme_cowplot(font_size = 12) +
        theme(legend.position = "bottom",
              strip.background = element_blank())
    ggsave(paste0("output/cases_hosps_and_deaths_averted_date",sim_run,".pdf"),p,width = 10, height = 3.5)
    
    # Combine numbers averted plot with counterfactual epidemic waves plot
    align_labels = function(plotsList){
        # Extract labels of y-axis
        # Note: Don't use the as.character on the maximum limit, 
        #       as decimal places in labels may increase the character count 
        y.labels <- lapply(plotsList, function(y){lapply(ggplot_build(y)$layout$panel_params, function(x){x$y$get_labels()})})
        
        # Calculate the maximum number of characters for each plot's labels
        maxChars <- lapply(y.labels, function(y){sapply(y, function(x){max(nchar(x),na.rm = T)})})
        
        # Define a function that would space-pad the labels and apply
        format.labels = function(label){str_pad(label, max(sapply(1:length(maxChars[[1]]),function(i) max(sapply(maxChars,"[[",i)))), pad = " ")}
        return(lapply(plotsList, function(x){return(x + scale_y_continuous(labels = format.labels))}))
    }
    
    p_list <- list(p + theme(axis.line.x = element_blank(),axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),axis.title.x = element_blank(),
                            legend.position = "none"),
                  p_cntfctl + xlim(NA_Date_,as.Date("2022-05-06")+100) + theme(strip.text = element_blank()))
    p_list <- align_labels(p_list)
    plot_grid(plotlist = p_list,nrow = 2,rel_heights = c(0.6,1))
    ggsave(paste0("output/cntfctl_and_averted_cases_hosps_and_deaths",sim_run,".pdf"),width = 10,height = 4)
    
    # Save workspace
    save(list = ls(all.names=T),file = paste0("output/cntfctl_output",sim_run,".RData"),envir = environment())
}
