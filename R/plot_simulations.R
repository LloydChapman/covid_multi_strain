plot_simulations <- function(q_outcomes,q_outcomes_cntfctl,q_total_outcomes_averted,ttls,sim_run,cols,dates){
    # Plot counterfactuals
    # No lockdown counterfactual
    idx <- 1
    p_lockdown <- plot_cases_hosps_deaths_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],ttls[names(ttls) %in% idx])
    ggsave(paste0("output/cntfctl_cases_hosps_and_deaths_lockdown",sim_run,".pdf"),p_lockdown,width = 9,height = 3)
    
    # Alternative lockdown timing counterfactuals
    idx <- 2:7
    p_alt_lockdown_timings <- plot_cases_hosps_deaths_counterfactuals1(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],ttls[names(ttls) %in% idx])
    ggsave(paste0("output/cntfctl_cases",sim_run,".pdf"),p_alt_lockdown_timings$p_cases,width = 8,height = 5.2)
    ggsave(paste0("output/cntfctl_hosps",sim_run,".pdf"),p_alt_lockdown_timings$p_hosps,width = 8,height = 5.2)
    ggsave(paste0("output/cntfctl_deaths",sim_run,".pdf"),p_alt_lockdown_timings$p_deaths,width = 8,height = 5.2)
    
    # No vaccination counterfactual
    idx <- 8
    p_vax <- plot_cases_hosps_deaths_counterfactuals(q_outcomes[cntfctl == idx],q_outcomes_cntfctl[cntfctl == idx],ttls[names(ttls) %in% idx])
    ggsave(paste0("output/cntfctl_hosps_and_deaths_vax",sim_run,".pdf"),p_vax,width = 9,height = 3)
    
    # No boosters counterfactual
    idx <- 10
    p_booster <- plot_cases_hosps_deaths_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],ttls[names(ttls) %in% idx])
    ggsave(paste0("output/cntfctl_hosps_and_deaths_booster",sim_run,".pdf"),p_booster,width = 9,height = 3)
    
    # Plot no lockdown, no vaccination, and no booster counterfactuals together
    outcome <- c("cases","hosps","deaths")
    ttls1 <- c("Cases","Hospitalisations","Hospital deaths")
    names(ttls1) <- outcome
    p_cntfctl <- plot_counterfactuals_together(q_outcomes[cntfctl %in% c(1,8,10)],q_outcomes_cntfctl[cntfctl %in% c(1,8,10)],outcome,ttls1,ttls[names(ttls) %in% c(1,8,10)])
    ggsave(paste0("output/cntfctl_cases_hosps_and_deaths",sim_run,".pdf"),p_cntfctl,width = 9,height = 3)
    
    # Plot numbers of cases, hospitalisations and deaths averted in each wave and overall
    lbls <- c(`1` = "Lockdowns",`8` = "Vaccination",`10` = "Boosters")
    p_outcomes_averted <- plot_outcomes_averted(q_total_outcomes_averted[cntfctl %in% c(1,8,10)],cols,outcome,lbls,ttls1,by_wave = T)
    ggsave(paste0("output/cases_hosps_and_deaths_averted",sim_run,".pdf"),p_outcomes_averted,width = 8, height = 3.5)
    
    # Same plot but with date rather than wave number on the x-axis
    p_outcomes_averted_date <- plot_outcomes_averted(q_total_outcomes_averted[cntfctl %in% c(1,8,10)],cols,outcome,lbls,ttls1,by_wave = F,dates = dates)
    ggsave(paste0("output/cases_hosps_and_deaths_averted_date",sim_run,".pdf"),p_outcomes_averted_date,width = 10, height = 3.5)
    
    # Combine numbers averted plot with counterfactual epidemic waves plot
    p_list <- list(p_outcomes_averted_date + theme(axis.line.x = element_blank(),axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),axis.title.x = element_blank(),
                             legend.position = "none"),
                   p_cntfctl + xlim(NA_Date_,as.Date("2022-05-06")+100) + theme(strip.text = element_blank()))
    p_list <- align_labels(p_list)
    plot_grid(plotlist = p_list,nrow = 2,rel_heights = c(0.6,1))
    ggsave(paste0("output/cntfctl_and_averted_cases_hosps_and_deaths",sim_run,".pdf"),width = 10,height = 4)
}


plot_cases_hosps_deaths_counterfactuals <- function(q_outcomes,q_outcomes_cntfctl,ttls){
    p <- plot_cases_hosps_deaths_counterfactuals1(q_outcomes,q_outcomes_cntfctl,ttls)
    p1 <- plot_grid(p$p_cases + theme(legend.position = "none"),
                   p$p_hosps + theme(legend.position = "none"),
                   p$p_deaths + theme(legend.position = "none"), nrow = 1)
    l <- get_legend(p$p_deaths)
    
    p2 <- plot_grid(p1,l,nrow = 2,rel_heights = c(1,0.1))
    return(p2)
}


plot_cases_hosps_deaths_counterfactuals1 <- function(q_outcomes,q_outcomes_cntfctl,ttls){
    # Cases
    p_cases <- plot_counterfactuals(q_outcomes,q_outcomes_cntfctl,"cases","Cases",ttls)
    # Hospitalisations
    p_hosps <- plot_counterfactuals(q_outcomes,q_outcomes_cntfctl,"hosps","Hospitalisations",ttls)
    # Deaths in hospital
    p_deaths <- plot_counterfactuals(q_outcomes,q_outcomes_cntfctl,"deaths","Hospital deaths",ttls)
    
    p <- list(p_cases = p_cases,p_hosps = p_hosps,p_deaths = p_deaths)
    return(p)
}


plot_outcomes_averted <- function(q_total_outcomes_averted,cols,outcome,lbls,ttls1,by_wave = T,dates = NULL){
    p_dt <- melt(q_total_outcomes_averted,measure.vars = patterns(c("med","q95l","q95u")),value.name = c("med","q95l","q95u"))
    p_dt[,variable := cols[variable]]
    p_dt <- p_dt[variable %in% outcome]
    p_dt[,variable := factor(variable,levels = outcome)]
    p_dt[,cntfctl := lbls[match(cntfctl,names(lbls))]]
    p_dt[,cntfctl := factor(cntfctl,levels = unique(cntfctl))]
    if (by_wave){
        p <- ggplot(p_dt,aes(x = wave,fill = factor(cntfctl))) +
            geom_bar(aes(y = -med),position = position_dodge(),stat = "identity") + 
            geom_errorbar(aes(ymin = -q95l,ymax = -q95u),width = 0.2,position = position_dodge(0.9)) +
            labs(x = "Wave",y = "Number averted",fill = "Counterfactual") +
            facet_wrap(~variable,scales = "free",labeller = labeller(variable = ttls1)) +
            theme_cowplot(font_size = 12) +
            theme(legend.position = "bottom",
                  strip.background = element_blank())
    } else {
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
    }
    return(p)
}


align_labels = function(plot_list){
    # Extract labels of y-axis
    # Note: Don't use the as.character on the maximum limit, 
    #       as decimal places in labels may increase the character count 
    ylbls <- lapply(plot_list, function(y){lapply(ggplot_build(y)$layout$panel_params, function(x){x$y$get_labels()})})
    
    # Calculate the maximum number of characters for each plot's labels
    maxChars <- lapply(ylbls, function(y){sapply(y, function(x){max(nchar(x),na.rm = T)})})
    
    # Define a function that would space-pad the labels and apply
    format_labels = function(label){str_pad(label, max(sapply(1:length(maxChars[[1]]),function(i) max(sapply(maxChars,"[[",i)))), pad = " ")}
    return(lapply(plot_list, function(x){return(x + scale_y_continuous(labels = format_labels))}))
}
