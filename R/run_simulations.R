library(odin)
library(odin.dust)
library(dust)
library(socialmixr)
library(qs)
library(data.table)
library(readxl)
library(ggplot2)
library(cowplot)
library(abind)

source("R/utils.R")
source("R/vaccination.R")
source("R/covid_multi_strain.R")
source("R/simulate.R")

## Process FP data
source("R/process_FP_data.R")

## Set MCMC output
output <- "output/MCMCoutput40.RData"

## Run counterfactual simulations
# Set number of parameter samples and burn-in to remove
n_smpls <- 500
burnin <- 2000
seed <- 1

# Set probabilities for quantiles for outcomes
probs <- c(0.025,0.5,0.975)

# Make a list of alternatives for lockdown dates
intvtn_date <- as.Date(c(strt_date,"2020-10-24","2021-06-30","2021-08-12"))
beta_date <- as.integer(intvtn_date - min(intvtn_date))
# beta_date_cntfctl_list <- replicate(2,beta_date,F)
beta_date_cntfctl_list <- replicate(7,beta_date,F)

# Move first lockdown a week earlier
beta_date_cntfctl_list[[1]][2] <- beta_date[2] - 7
# Move first lockdown a week later
beta_date_cntfctl_list[[2]][2] <- beta_date[2] + 7
# Move second lockdown a week earlier
beta_date_cntfctl_list[[3]][4] <- beta_date[4] - 7
# Move second lockdown a week later
beta_date_cntfctl_list[[4]][4] <- beta_date[4] + 7
# Move both lockdowns a week earlier
beta_date_cntfctl_list[[5]][c(2,4)] <- beta_date[c(2,4)] - 7
# Move both lockdowns a week later
beta_date_cntfctl_list[[6]][c(2,4)] <- beta_date[c(2,4)] + 7
# Use actual dates for vaccine counterfactual simulation

# Set counterfactual vaccine schedule
schedule_cntfctl_list <- replicate(7,schedule,F)

# No vaccination in counterfactual
schedule_cntfctl_list[[7]]$doses <- array(0,dim = dim(schedule$doses))

q_outcomes_list <- vector("list",length(beta_date_cntfctl_list))
q_outcomes_cntfctl_list <- vector("list",length(beta_date_cntfctl_list))
q_outcomes_averted_list <- vector("list",length(beta_date_cntfctl_list))
q_total_outcomes_list <- vector("list",length(beta_date_cntfctl_list))
q_total_outcomes_cntfctl_list <- vector("list",length(beta_date_cntfctl_list))
q_total_outcomes_averted_list <- vector("list",length(beta_date_cntfctl_list))
for (i in seq_along(beta_date_cntfctl_list)){
    out <- simulate_counterfactual(output,n_smpls,beta_date_cntfctl_list[[i]],schedule_cntfctl_list[[i]],burnin,seed)
    states_cntfctl <- out$states_cntfctl
    smpl <- out$smpl
    info <- out$info
    states <- out$states
    
    # Calculate pairwise differences in simulated outcomes across time
    outcomes_averted <- states - states_cntfctl
    # dimnames(outcomes_averted[[i]])[[1]] <- names(index(info)$state)
    
    ## Calculate total differences in outcomes over each wave
    # Set wave dates
    wave_date <- c(1,300,nlayer(states_cntfctl))

    q_outcomes_list[[i]] <- calculate_outcome_quantiles(states,info)
    q_outcomes_cntfctl_list[[i]] <- calculate_outcome_quantiles(states_cntfctl,info)
    q_outcomes_averted_list[[i]] <- calculate_outcome_quantiles(outcomes_averted,info)
        
    q_total_outcomes_list[[i]] <- calculate_outcomes_by_wave(states,wave_date,info)
    q_total_outcomes_cntfctl_list[[i]] <- calculate_outcomes_by_wave(states_cntfctl,wave_date,info)
    q_total_outcomes_averted_list[[i]] <- calculate_outcomes_by_wave(outcomes_averted,wave_date,info)
    
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

# Bind out
tmp <- q_total_outcomes[cntfctl == 1]
tmp[,cntfctl := 0]
total_outcomes <- rbind(tmp,q_total_outcomes_cntfctl)

ttls <- c("No change in lockdown dates",
          "1st lockdown 1 week earlier","1st lockdown 1 week later",
          "2nd lockdown 1 week earlier","2nd lockdown 1 week later",
          "Both lockdowns 1 week earlier","Both lockdowns 1 week later",
          "No vaccination")
names(ttls) <- as.character(seq_along(ttls)-1)

tbl <- total_outcomes[,.(Counterfactual = ttls[match(cntfctl,names(ttls))],
                         Wave = wave,
                         Hospitalisations = med_and_CI(hosps.med,hosps.q95l,hosps.q95u,d = 3,method = "signif"),
                         Deaths = med_and_CI(deaths.med,deaths.q95l,deaths.q95u,d = 3,method = "signif"))]
tbl[,Counterfactual := factor(Counterfactual, levels = unique(Counterfactual))]
tbl <- dcast(tbl, Counterfactual ~ Wave, value.var = c("Hospitalisations","Deaths"))
write.csv(tbl,"output/table1_1.csv",row.names = F)

# Plot counterfactuals

# Hospitalisations
idx <- 1:6
p_hosps <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"hosps","Hospitalisations",ttls[names(ttls) %in% idx])
ggsave("output/cntfctl_hosps1.pdf",p_hosps,width = 8,height = 7)
# Deaths in hospital
p_deaths <- plot_counterfactuals(q_outcomes[cntfctl %in% idx],q_outcomes_cntfctl[cntfctl %in% idx],"deaths","Deaths",ttls[names(ttls) %in% idx])
ggsave("output/cntfctl_deaths1.pdf",p_deaths,width = 8,height = 7)

idx <- 7
p_hosps_vax <- plot_counterfactuals(q_outcomes[cntfctl == idx],q_outcomes_cntfctl[cntfctl == idx],"hosps","Hospitalisations",ttls[names(ttls) %in% idx])
# ggsave("output/cntfctl_hosps_vax.pdf",p_hosps_vax,width = 8,height = 6)
# Deaths in hospital
p_deaths_vax <- plot_counterfactuals(q_outcomes[cntfctl == idx],q_outcomes_cntfctl[cntfctl == idx],"deaths","Deaths",ttls[names(ttls) %in% idx])
pp <- plot_grid(p_hosps_vax + theme(legend.position = "none"),
                p_deaths_vax + theme(legend.position = "none"))
l <- get_legend(p_deaths_vax)
ggsave("output/cntfctl_hosps_and_deaths_vax1.pdf",plot_grid(pp,l,nrow = 2,rel_heights = c(1,0.1)),width = 6,height = 4)

# save.image("../covid_multistrain_wip3.RData")

