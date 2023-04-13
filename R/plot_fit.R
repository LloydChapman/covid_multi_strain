plot_fit <- function(output,run,pop,burnin = NULL,moving_avg = FALSE,n_smpls = 1000){
    load(output)
    
    if (is.null(burnin)){
        burnin <- round(nrow(res$pars)/10)
    }
    
    plot_outcome(res$trajectories$state,data,"cases",res$pars[,"phi_cases"],by_age = T,burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/cases_by_age_fit",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 8, height = 6.4)
    plot_outcome(res$trajectories$state,data,"hosps",by_age = T,burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/hosps_by_age_fit",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 4, height = 8)
    plot_outcome(res$trajectories$state,data,"deaths",by_age = T,burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/deaths_by_age_fit",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 4, height = 8)
    plot_sero(res$trajectories$state,data,pop[,.(population = sum(total)),by = .(age_group)],by_age = T,burnin = burnin)
    ggsave(paste0("output/sero_by_age_fit",run,".pdf"),width = 4, height = 8)
    plot_outcome(res$trajectories$state,data,"cases",res$pars[,"phi_cases"],burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/cases",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 4, height = 2.7)
    plot_outcome(res$trajectories$state,data,"hosps",burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/hosps",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 4, height = 2.7)
    plot_outcome(res$trajectories$state,data,"deaths",burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/deaths",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 4, height = 2.7)
    
    vrble <- c("cases","hosps","deaths")
    ttls <- c("Cases","Hospitalisations","Hospital deaths")
    names(ttls) <- vrble
    p_list <- plot_outcome_by_age(res$trajectories$state,vrble,res$pars[,"phi_cases"],ttls,n_smpls,burnin = burnin)
    ggsave(paste0("output/outcomes_by_age",run,".pdf"),p_list$p,width = 9,height = 3)
    ggsave(paste0("output/prop_outcomes_by_age",run,".pdf"),p_list$p1,width = 9,height = 3)
    
    base <- pars$base
    res1 <- plot_transmission_rate(res$pars,base$beta_type,base$beta_date,base$dt,max(data$day_end),burnin = burnin)
    p <- res1$p
    res2 <- plot_transmission_rate(res$pars[,c(1,3,5:14)],base$beta_type,base$beta_date[c(1,3,5)],base$dt,max(data$day_end),burnin = burnin)
    beta_step <- res2$beta_step
    p1 <- p + geom_line(aes(x = date,y = `50%`,color = "No lockdowns"),beta_step) + 
        geom_ribbon(aes(x = date,ymin = `2.5%`,ymax = `97.5%`,fill = "No lockdowns"),beta_step,alpha = 0.5) + 
        scale_color_manual(name = "",values = c("Fitted" = "black","No lockdowns" = "blue")) +
        scale_fill_manual(name = "",values = c("Fitted" = "black","No lockdowns" = "blue")) + 
        theme(legend.position = "bottom")
    ggsave(paste0("output/transmission_rate",run,".pdf"),p1,width = 6,height = 6)
    
    plot_traces(res$pars,u)
    ggsave(paste0("output/par_traces",run,".pdf"),width = 6,height = 6)
    
    plot_posteriors(res$pars,u,priors,pars_min,pars_max,burnin = burnin)
    ggsave(paste0("output/par_posteriors",run,".pdf"),width = 6,height = 6)
    
    p2 <- plot_pairwise_correlation(res$pars,u,burnin = burnin)
    ggsave(paste0("output/par_pairwise_corr",run,".pdf"),p2,width = 10,height = 10)
}
