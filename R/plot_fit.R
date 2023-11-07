plot_fit <- function(dat,pars,run,pop,u,lbls,moving_avg = FALSE,pred_intvl = FALSE,n_smpls = 1000){
    model_type <- pars$base$model_type
    if (model_type == "NB"){
        plot_outcome(dat$samples$trajectories,dat$data,"cases",phi = dat$samples$pars[,"phi_cases"],by_age = T,moving_avg = moving_avg,
                     pred_intvl = pred_intvl,alpha = dat$samples$pars[,"alpha_cases"])
        ggsave(paste0("output/cases_by_age_fit",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 8, height = 6.4)
    }
    plot_outcome(dat$samples$trajectories,dat$data,"hosps",by_age = T,moving_avg = moving_avg,
                 pred_intvl = pred_intvl,alpha = dat$samples$pars[,"alpha_hosp"])
    ggsave(paste0("output/hosps_by_age_fit",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 4, height = 8)
    plot_outcome(dat$samples$trajectories,dat$data,"deaths",by_age = T,moving_avg = moving_avg,
                 pred_intvl = pred_intvl,alpha = dat$samples$pars[,"alpha_death"])
    ggsave(paste0("output/deaths_by_age_fit",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 4, height = 8)
    plot_sero(dat$samples$trajectories,dat$data,pop[,.(population = sum(total)),by = .(age_group)],by_age = T)
    ggsave(paste0("output/sero_by_age_fit",run,".pdf"),width = 4, height = 10)
    if (model_type == "NB"){
        plot_outcome(dat$samples$trajectories,dat$data,"cases",phi = dat$samples$pars[,"phi_cases"],by_age = F,moving_avg = moving_avg,
                     pred_intvl = pred_intvl,alpha = dat$samples$pars[,"alpha_cases"])
        ggsave(paste0("output/cases",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 4, height = 2.7)
    } else if (model_type == "BB"){
        plot_tests(dat$samples$trajectories,dat$data,dat$samples$pars[,"p_NC"],pop[,sum(total)],moving_avg = moving_avg)
        ggsave(paste0("output/tests",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 4, height = 2.7)
    }
    plot_outcome(dat$samples$trajectories,dat$data,"hosps",moving_avg = moving_avg,
                 pred_intvl = pred_intvl,alpha = dat$samples$pars[,"alpha_hosp"])
    ggsave(paste0("output/hosps",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 4, height = 2.7)
    plot_outcome(dat$samples$trajectories,dat$data,"deaths",moving_avg = moving_avg,
                 pred_intvl = pred_intvl,alpha = dat$samples$pars[,"alpha_death"])
    ggsave(paste0("output/deaths",ifelse(moving_avg,"_moving_avg",""),run,".pdf"),width = 4, height = 2.7)
    # plot_variant_proportion(dat$samples$trajectories,dat$data,as.Date("2022-01-01"),"Omicron BA.2 proportion")
    # ggsave(paste0("output/variant_proportion",run,".pdf"),width = 4,height = 2.7)
    
    if (model_type == "NB"){
        vrble <- c("cases","hosps","deaths")
        ttls <- c("Cases","Hospitalisations","Hospital deaths")
    } else if (model_type == "BB"){
        vrble <- c("hosps","deaths")
        ttls <- c("Hospitalisations","Hospital deaths")
    }
    names(ttls) <- vrble
    p_list <- plot_outcome_by_age(dat$samples$trajectories,vrble,dat$samples$pars[,"phi_cases"],ttls,n_smpls)
    ggsave(paste0("output/outcomes_by_age",run,".pdf"),p_list$p,width = 9,height = 3)
    ggsave(paste0("output/prop_outcomes_by_age",run,".pdf"),p_list$p1,width = 9,height = 3)
    
    base <- pars$base
    res1 <- plot_transmission_rate(dat$samples$pars,base$beta_type,base$beta_date,base$dt,max(dat$data$day_end))
    p <- res1$p
    res2 <- plot_transmission_rate(dat$samples$pars[,c(1,3,5:14)],base$beta_type,base$beta_date[c(1,3,5)],base$dt,max(dat$data$day_end))
    beta_step <- res2$beta_step
    p1 <- p + geom_line(aes(x = date,y = `50%`,color = "No lockdowns"),beta_step) +
        geom_ribbon(aes(x = date,ymin = `2.5%`,ymax = `97.5%`,fill = "No lockdowns"),beta_step,alpha = 0.5) +
        scale_color_manual(name = "",values = c("Fitted" = "black","No lockdowns" = "blue")) +
        scale_fill_manual(name = "",values = c("Fitted" = "black","No lockdowns" = "blue")) +
        theme(legend.position = "bottom")
    ggsave(paste0("output/transmission_rate",run,".pdf"),p1,width = 6,height = 6)
    
    # plot_traces(dat$samples$pars,u)
    # ggsave(paste0("output/par_traces",run,".pdf"),width = 6,height = 6)
    
    priors <- lapply(pars$mcmc$.__enclos_env__$private$parameters,"[[","prior")
    pars_min <- pars$info$min
    pars_max <- pars$info$max
    plot_posteriors(dat$samples$pars,u,priors,pars_min,pars_max,lbls)
    ggsave(paste0("output/par_posteriors",run,".pdf"),width = 6,height = 6)
    
    p2 <- plot_pairwise_correlation(dat$samples$pars,u,lbls)
    ggsave(paste0("output/par_pairwise_corr",run,".pdf"),p2,width = 14,height = 14)
}
