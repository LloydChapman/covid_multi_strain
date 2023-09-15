plot_sensitivity_analysis <- function(sim_runs,assumption){
    q_total_outcomes_averted_SA_list <- vector("list",length(sim_runs))
    q_outcomes_SA_list <- vector("list",length(sim_runs))
    q_outcomes_cntfctl_SA_list <- vector("list",length(sim_runs))
    cntfctls <- c(1,8,10)
    outcomes <- c("cases","hosps","deaths")
    for (j in seq_along(sim_runs)){
        # load(paste0("output/cntfctl_output",sim_runs[j],".RData"))
        res <- readRDS(paste0("output/cntfctl_output",sim_runs[j],".RDS"))
        tmp <- res$q_total_outcomes_averted[cntfctl %in% cntfctls & wave == "Total",c("cntfctl",paste(rep(outcomes,each = 3),c("med","q95l","q95u"),sep = ".")),with = F]
        q_total_outcomes_averted_SA_list[[j]] <- tmp
        tmp1 <- res$q_outcomes[cntfctl %in% cntfctls & state %in% outcomes]
        tmp2 <- res$q_outcomes_cntfctl[cntfctl %in% cntfctls & state %in% outcomes]
        q_outcomes_SA_list[[j]] <- tmp1
        q_outcomes_cntfctl_SA_list[[j]] <- tmp2
    }
    
    q_total_outcomes_averted_SA <- rbindlist(q_total_outcomes_averted_SA_list,idcol = "assumptions")
    q_outcomes_SA <- rbindlist(q_outcomes_SA_list,idcol = "assumptions")
    q_outcomes_SA[,state := factor(state,levels = outcomes)]
    q_outcomes_cntfctl_SA <- rbindlist(q_outcomes_cntfctl_SA_list,idcol = "assumptions")
    q_outcomes_cntfctl_SA[,state := factor(state,levels = outcomes)]
    
    q_total_outcomes_averted_SA_long <- melt(q_total_outcomes_averted_SA,id.vars = c("assumptions","cntfctl")) #measure.vars = patterns("cases","hosps","deaths"),value.name = c("med","q95l","q95u"))
    q_total_outcomes_averted_SA_long[,assumptions := assumption[assumptions]]
    q_total_outcomes_averted_SA_long[,assumptions := factor(assumptions,levels = unique(assumptions))]
    q_total_outcomes_averted_SA_long[,quantile := sub(".*\\.","",variable)]
    q_total_outcomes_averted_SA_long[,variable := sub("\\..*","",variable)]
    q_total_outcomes_averted_SA_long[,variable := factor(variable,levels = unique(variable))]
    
    q_total_outcomes_averted_SA_long <- dcast(q_total_outcomes_averted_SA_long,assumptions + cntfctl + variable ~ quantile,value.var = "value")
    lbls <- c(`1` = "Lockdowns",`8` = "Vaccination",`10` = "Boosters") 
    q_total_outcomes_averted_SA_long[,cntfctl := lbls[match(cntfctl,names(lbls))]]
    q_total_outcomes_averted_SA_long[,cntfctl := factor(cntfctl,levels = unique(cntfctl))]
    
    ttls1 <- c("Cases","Hospitalisations","Hospital deaths")
    names(ttls1) <- outcomes
    ggplot(q_total_outcomes_averted_SA_long,aes(x = cntfctl,fill = assumptions)) + 
        geom_bar(aes(y = -med),position = position_dodge(),stat = "identity") + 
        geom_errorbar(aes(ymin = -q95l,ymax = -q95u),width = 0.2,position = position_dodge(0.9)) + 
        labs(x = "Intervention",y = "Number averted",fill = "Parameter assumptions") +
        facet_wrap(~variable,scales = "free",labeller = labeller(variable = ttls1)) + 
        theme_cowplot(font_size = 11) +
        theme(legend.position = "bottom",strip.background = element_blank())
    ggsave("output/cases_hosps_deaths_averted_SA.pdf",width = 8,height = 3.5)
    
    q_outcomes_SA1 <- q_outcomes_SA[,.(med = median(med),q95l = min(q95l),q95u = max(q95u)),by = .(cntfctl,state,day,date)]
    q_outcomes_cntfctl_SA[,assumptions := assumption[assumptions]]
    q_outcomes_cntfctl_SA[,assumptions := factor(assumptions,levels = unique(assumptions))]
    
    ttls <- c("No change in lockdown dates","No lockdowns",
              "1st lockdown 2 weeks earlier","1st lockdown 2 weeks later",
              "2nd lockdown 2 weeks earlier","2nd lockdown 2 weeks later",
              "Both lockdowns 2 weeks earlier","Both lockdowns 2 weeks later",
              "No vaccination", "Boosters 1 month earlier","No boosters")#,"Boosters")
    names(ttls) <- as.character(seq_along(ttls)-1)
    
    ggplot() + 
        geom_line(aes(x = date,y = med,group = assumptions,color = assumptions),q_outcomes_cntfctl_SA) + 
        geom_ribbon(aes(x = date,ymin = q95l,ymax = q95u,fill = assumptions),q_outcomes_cntfctl_SA,alpha = 0.5) + 
        geom_line(aes(x = date,y = med),q_outcomes_SA1,color = "black",linetype = "dashed") + 
        geom_ribbon(aes(x = date,ymin = q95l,ymax = q95u),q_outcomes_SA1,fill = "black",alpha = 0.2) + 
        labs(x = "Date",color = "Parameter assumptions",fill = "Parameter assumptions") + 
        facet_grid2(cntfctl ~ state,scales = "free",independent = "y",labeller = labeller(state = ttls1,cntfctl = ttls)) +
        theme_cowplot(font_size = 11) +
        theme(axis.title.y = element_blank(),
              legend.position = "bottom",
              strip.background = element_blank())
    ggsave("output/cntfctl_cases_hosps_deaths_SA.pdf",width = 8,height = 6)
}
