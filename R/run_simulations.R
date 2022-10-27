library(odin)
library(odin.dust)
library(dust)
library(mcstate)
library(socialmixr)
library(qs)
library(data.table)
library(readxl)
library(lubridate)
library(ISOweek)
library(ggplot2)
library(cowplot)
library(abind)

source("R/utils.R")
source("R/date.R")
source("R/vaccination.R")
source("R/covid_multi_strain.R")
source("R/simulate.R")

## Process FP data
source("R/process_FP_data.R")

# Create dust model generator
covid_multi_strain <- odin_dust("inst/odin/covid_multi_strain.R")

## Set MCMC output
mcmc_run <- 70
output <- paste0("output/MCMCoutput",mcmc_run,".RData")

## Run counterfactual simulations
# Set run number
run <- 9
# Set whether states required to calculate Rt have been output 
Rt <- T
# Set number of parameter samples and burn-in to remove
n_smpls <- 500
burnin <- NULL #2000
seed <- 1

# Set probabilities for quantiles for outcomes
probs <- c(0.025,0.5,0.975)

## Lockdown counterfactuals
# Make a list of alternatives for lockdown dates
# intvtn_date <- as.Date(c(strt_date-1,"2020-10-24","2021-06-30","2021-08-12","2021-12-31","2022-02-15"))
# intvtn_date <- as.Date(c(strt_date-1,"2020-10-24","2021-06-30","2021-08-12","2021-12-31"))
intvtn_date <- as.Date(c("2020-08-27","2020-10-24","2021-06-01","2021-08-02","2021-11-15"))
beta_date <- sircovid_date(intvtn_date) #as.integer(intvtn_date - min(intvtn_date))
beta_date_cntfctl_list <- replicate(9,beta_date,F) #replicate(10,beta_date,F)

# Move first lockdown a week earlier
beta_date_cntfctl_list[[1]][1:2] <- beta_date[1:2] - 14
# Move first lockdown a week later
beta_date_cntfctl_list[[2]][1:2] <- beta_date[1:2] + 14
# Move second lockdown two weeks earlier
beta_date_cntfctl_list[[3]][3:4] <- beta_date[3:4] - 14
# Move second lockdown two weeks later
beta_date_cntfctl_list[[4]][3:4] <- beta_date[3:4] + 14
# Move both lockdowns earlier
beta_date_cntfctl_list[[5]][1:2] <- beta_date[1:2] - 14
beta_date_cntfctl_list[[5]][3:4] <- beta_date[3:4] - 14
# Move both lockdowns later
beta_date_cntfctl_list[[6]][1:2] <- beta_date[1:2] + 14
beta_date_cntfctl_list[[6]][3:4] <- beta_date[3:4] + 14
# Use actual dates for vaccination counterfactual simulations

## Vaccination counterfactuals
# Set counterfactual vaccine schedule
schedule_cntfctl_list <- replicate(9,schedule,F) #replicate(10,schedule,F)
# No vaccination in counterfactual
schedule_cntfctl_list[[7]]$doses <- array(0,dim = dim(schedule$doses))
# Boosters starting 1 month earlier
schedule_cntfctl_list[[8]] <- change_booster_timing(schedule, 30)
# # Boosters starting 2 weeks later
# schedule_cntfctl_list[[9]] <- change_booster_timing(schedule, -180)
# No boosters
schedule_cntfctl_list[[9]]$doses[,3,] <- matrix(0, nrow = nrow(schedule$doses), 
                                              ncol = nlayer(schedule$doses))

# Create objects for storing quantiles of simulation output
q_outcomes_list <- vector("list",length(beta_date_cntfctl_list))
q_outcomes_cntfctl_list <- vector("list",length(beta_date_cntfctl_list))
q_outcomes_averted_list <- vector("list",length(beta_date_cntfctl_list))
q_total_outcomes_list <- vector("list",length(beta_date_cntfctl_list))
q_total_outcomes_cntfctl_list <- vector("list",length(beta_date_cntfctl_list))
q_total_outcomes_averted_list <- vector("list",length(beta_date_cntfctl_list))

new_infections_cntfctl_by_vacc_list <- vector("list",length(beta_date_cntfctl_list))
new_reinfections_cntfctl_by_vacc_list <- vector("list",length(beta_date_cntfctl_list))
# Run simulations
for (i in seq_along(beta_date_cntfctl_list)){ #9:10){ # 
    out <- simulate_counterfactual(output,n_smpls,beta_date_cntfctl_list[[i]],
                                   schedule_cntfctl_list[[i]],burnin = burnin,seed = seed)
    states_cntfctl <- out$states_cntfctl
    smpl <- out$smpl
    info <- out$info
    states <- out$states
    
    # Add row names to states_cntfctl array
    dimnames(states_cntfctl) <- list(names(index(info,min_ages,Rt)$state))
    
    # Calculate pairwise differences in simulated outcomes across time
    outcomes_averted <- states - states_cntfctl
    # dimnames(outcomes_averted[[i]])[[1]] <- names(index(info)$state)

    ## Calculate total differences in outcomes over each wave
    # Set wave dates
    wave_date <- c(1,300,500,nlayer(states_cntfctl))

    q_outcomes_list[[i]] <- calculate_outcome_quantiles(states,info,min_ages,Rt)
    q_outcomes_cntfctl_list[[i]] <- calculate_outcome_quantiles(states_cntfctl,info,min_ages,Rt)
    q_outcomes_averted_list[[i]] <- calculate_outcome_quantiles(outcomes_averted,info,min_ages,Rt)

    q_total_outcomes_list[[i]] <- calculate_outcomes_by_wave(states,wave_date,info,min_ages,Rt)
    q_total_outcomes_cntfctl_list[[i]] <- calculate_outcomes_by_wave(states_cntfctl,wave_date,info,min_ages,Rt)
    q_total_outcomes_averted_list[[i]] <- calculate_outcomes_by_wave(outcomes_averted,wave_date,info,min_ages,Rt)

    # res <- calculate_new_infections_by_vacc(states_cntfctl)
    # new_infections_cntfctl_by_vacc_list[[i]] <- res$new_infections_by_vacc
    # new_reinfections_cntfctl_by_vacc_list[[i]] <- res$new_reinfections_by_vacc
}

q_outcomes <- rbindlist(q_outcomes_list, idcol = "cntfctl")
q_outcomes_cntfctl <- rbindlist(q_outcomes_cntfctl_list, idcol = "cntfctl")
q_outcomes_averted <- rbindlist(q_outcomes_averted_list, idcol = "cntfctl")
q_outcomes[,date := strt_date + day - 1]
q_outcomes_cntfctl[,date := strt_date + day - 1]
q_outcomes_averted[,date := strt_date + day - 1]

q_total_outcomes <- rbindlist(q_total_outcomes_list, idcol = "cntfctl")
q_total_outcomes_cntfctl <- rbindlist(q_total_outcomes_cntfctl_list, idcol = "cntfctl")
q_total_outcomes_averted <- rbindlist(q_total_outcomes_averted_list, idcol = "cntfctl")

# new_infections_cntfctl_by_vacc <- rbindlist(new_infections_cntfctl_by_vacc_list, idcol = "cntfctl")
# new_reinfections_cntfctl_by_vacc <- rbindlist(new_reinfections_cntfctl_by_vacc_list, idcol = "cntfctl")

ttls <- c("No change in lockdown dates",
          "1st lockdown 2 weeks earlier","1st lockdown 2 weeks later",
          "2nd lockdown 2 weeks earlier","2nd lockdown 2 weeks later",
          "Both lockdowns 2 weeks earlier","Both lockdowns 2 weeks later",
          "No vaccination", "Boosters 1 month earlier","No boosters")#,"Boosters")
names(ttls) <- as.character(seq_along(ttls)-1)

# pdf("output/infctns_by_vacc_status.pdf",)
# for (i in 9:10){
#     idx_wave3 = wave_date[3]:wave_date[4]
#     dates <- seq.Date(sircovid_date_as_date(wave_date[3]),sircovid_date_as_date(wave_date[4]),by = 1)
#     dates_plot <- seq.Date(sircovid_date_as_date(wave_date[3]),sircovid_date_as_date(wave_date[4]),by = 30)
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
                         Deaths = med_and_CI(deaths.med,deaths.q95l,deaths.q95u,d = 3,method = "signif"))]
tbl[,Counterfactual := factor(Counterfactual, levels = unique(Counterfactual))]
tbl <- dcast(tbl, Counterfactual ~ Wave, value.var = c("Cases","Hospitalisations","Deaths"))
write.csv(tbl,paste0("output/table1_",run,".csv"),row.names = F)

write.csv(tbl[c(1,8,10),],paste0("output/table1_vax_",run,".csv"),row.names = F)

# Plot counterfactuals

# Hospitalisations
idx <- 1:6
p_cases <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"cases","Cases",ttls[names(ttls) %in% idx])
ggsave(paste0("output/cntfctl_cases",run,".pdf"),p_cases,width = 8,height = 5.2)
p_hosps <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"hosps","Hospitalisations",ttls[names(ttls) %in% idx])
ggsave(paste0("output/cntfctl_hosps",run,".pdf"),p_hosps,width = 8,height = 5.2)
# Deaths in hospital
p_deaths <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"deaths","Deaths",ttls[names(ttls) %in% idx])
ggsave(paste0("output/cntfctl_deaths",run,".pdf"),p_deaths,width = 8,height = 5.2)

idx <- 7
p_cases_vax <- plot_counterfactuals(q_outcomes[cntfctl == idx],q_outcomes_cntfctl[cntfctl == idx],"cases","Cases",ttls[names(ttls) %in% idx])
p_hosps_vax <- plot_counterfactuals(q_outcomes[cntfctl == idx],q_outcomes_cntfctl[cntfctl == idx],"hosps","Hospitalisations",ttls[names(ttls) %in% idx])
# ggsave("output/cntfctl_hosps_vax.pdf",p_hosps_vax,width = 8,height = 6)
# Deaths in hospital
p_deaths_vax <- plot_counterfactuals(q_outcomes[cntfctl == idx],q_outcomes_cntfctl[cntfctl == idx],"deaths","Deaths",ttls[names(ttls) %in% idx])
pp <- plot_grid(p_cases_vax + theme(legend.position = "none"),
                p_hosps_vax + theme(legend.position = "none"),
                p_deaths_vax + theme(legend.position = "none"), nrow = 1)
l <- get_legend(p_deaths_vax)
ggsave(paste0("output/cntfctl_hosps_and_deaths_vax",run,".pdf"),plot_grid(pp,l,nrow = 2,rel_heights = c(1,0.1)),width = 9,height = 3)

idx <- 9 #8:9
p_cases_booster <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"cases","Cases",ttls[names(ttls) %in% idx])
p_hosps_booster <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"hosps","Hospitalisations",ttls[names(ttls) %in% idx])
# ggsave("output/cntfctl_hosps_vax.pdf",p_hosps_vax,width = 8,height = 6)
# Deaths in hospital
p_deaths_booster <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"deaths","Deaths",ttls[names(ttls) %in% idx])
pp1 <- plot_grid(p_cases_booster + theme(legend.position = "none"),
                 p_hosps_booster + theme(legend.position = "none"),
                 p_deaths_booster + theme(legend.position = "none"), nrow = 1)
l1 <- get_legend(p_deaths_booster)
ggsave(paste0("output/cntfctl_hosps_and_deaths_booster",run,".pdf"),plot_grid(pp1,l1,nrow = 2,rel_heights = c(1,0.1)),width = 9,height = 3)

save.image(paste0("output/cntfctl_output",mcmc_run,".RData"))

