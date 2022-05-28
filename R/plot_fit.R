plot_fit <- function(output){
    load(output)
    
    dimnames(res$states) <- list(names(idx$state)) #dimnames(pmcmc_run$trajectories$state)
    
    # pdf("output/hosps_by_age_fit.pdf",width = 4, heigh = 8)
    plot_hosps_age(res$states,data,strt_date+data_raw$day-1,n_age,n_vax,n_strains)
    # dev.off()
    # pdf("output/deaths_by_age_fit.pdf",width = 4, heigh = 8)
    plot_deaths_age(res$states,data,strt_date+data_raw$day-1,n_age,n_vax,n_strains)
    # dev.off()
    # pdf("output/sero_by_age_fit.pdf",width = 4, heigh = 8)
    plot_sero(res$states,data,strt_date+data_raw$day-1,population[3:length(population)])
    # dev.off()
    # pdf("output/cases.pdf",width = 4, height = 2.7)
    plot_cases(res$states,data,strt_date+data_raw$day-1)
    # dev.off()
}
