plot_fit <- function(output,run){
    load(output)
    
    dimnames(res$trajectories$state) <- list(names(idx$state)) #dimnames(pmcmc_run$trajectories$state)
    
    pdf(paste0("output/hosps_by_age_fit",run,".pdf"),width = 4, heigh = 8)
    plot_hosps_age(res$trajectories$state,data,strt_date+data_raw$day-1,n_age,n_vax,n_strains)
    dev.off()
    pdf(paste0("output/deaths_by_age_fit",run,".pdf"),width = 4, heigh = 8)
    plot_deaths_age(res$trajectories$state,data,strt_date+data_raw$day-1,n_age,n_vax,n_strains)
    dev.off()
    pdf(paste0("output/sero_by_age_fit",run,".pdf"),width = 4, heigh = 8)
    plot_sero(res$trajectories$state,data,strt_date+data_raw$day-1,population[3:length(population)])
    dev.off()
    pdf(paste0("output/cases",run,".pdf"),width = 4, height = 2.7)
    plot_cases(res$trajectories$state,data,strt_date+data_raw$day-1)
    dev.off()
}
