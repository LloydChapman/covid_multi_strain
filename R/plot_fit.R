plot_fit <- function(output,run,burnin = NULL,moving_avg = FALSE){
    load(output)
    
    if (is.null(burnin)){
        burnin <- round(nrow(res$pars)/10)
    }
    
    plot_outcome(res$trajectories$state,data,"cases",res$pars[,"phi_cases"],by_age = T,burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/cases_by_age_fit",run,".pdf"),width = 4, height = 12)
    plot_outcome(res$trajectories$state,data,"hosps",by_age = T,burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/hosps_by_age_fit",run,".pdf"),width = 4, height = 8)
    plot_outcome(res$trajectories$state,data,"deaths",by_age = T,burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/deaths_by_age_fit",run,".pdf"),width = 4, height = 8)
    plot_sero(res$trajectories$state,data,agg_pop,by_age = T,burnin = burnin)
    ggsave(paste0("output/sero_by_age_fit",run,".pdf"),width = 4, height = 8)
    plot_outcome(res$trajectories$state,data,"cases",res$pars[,"phi_cases"],burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/cases",run,".pdf"),width = 4, height = 2.7)
    plot_outcome(res$trajectories$state,data,"hosps",burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/hosps",run,".pdf"),width = 4, height = 2.7)
    plot_outcome(res$trajectories$state,data,"deaths",burnin = burnin,moving_avg = moving_avg)
    ggsave(paste0("output/deaths",run,".pdf"),width = 4, height = 2.7)
    
    plot_transmission_rate(res$pars,beta_date,n_betas,max(data$day_end),burnin = burnin)
    ggsave(paste0("output/transmission_rate",run,".pdf"),width = 6,height = 6)
    
    plot_traces(res$pars,u)
    ggsave(paste0("output/par_traces",run,".pdf"),width = 6,height = 6)
    
    plot_posteriors(res$pars,u,priors,pars_min,pars_max,burnin = burnin)
    ggsave(paste0("output/par_posteriors",run,".pdf"),width = 6,height = 6)
    
    p <- plot_pairwise_correlation(res$pars,u,burnin = burnin)
    ggsave(paste0("output/par_pairwise_corr",run,".pdf"),p,width = 10,height = 10)
}
