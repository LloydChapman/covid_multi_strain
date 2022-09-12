plot_fit <- function(output,run){
    load(output)
    
    plot_outcome(res$trajectories$state,data,"cases",phi_cases$initial,by_age = T)
    ggsave(paste0("output/cases_by_age_fit",run,".pdf"),width = 4, height = 12)
    plot_outcome(res$trajectories$state,data,"hosps",by_age = T)
    ggsave(paste0("output/hosps_by_age_fit",run,".pdf"),width = 4, height = 8)
    plot_outcome(res$trajectories$state,data,"deaths",by_age = T)
    ggsave(paste0("output/deaths_by_age_fit",run,".pdf"),width = 4, height = 8)
    plot_sero(res$trajectories$state,data,agg_pop,by_age = T)
    ggsave(paste0("output/sero_by_age_fit",run,".pdf"),width = 4, height = 8)
    plot_outcome(res$trajectories$state,data,"cases",phi_cases$initial)
    ggsave(paste0("output/cases",run,".pdf"),width = 4, height = 2.7)
    plot_outcome(res$trajectories$state,data,"hosps")
    ggsave(paste0("output/hosps",run,".pdf"),width = 4, height = 2.7)
    plot_outcome(res$trajectories$state,data,"deaths")
    ggsave(paste0("output/deaths",run,".pdf"),width = 4, height = 2.7)
    
}
