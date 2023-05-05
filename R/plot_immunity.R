library(data.table)
library(ggplot2)
library(cowplot)

source("R/utils.R")
source("R/date.R")
source("R/covid_multi_strain.R")
source("R/simulate.R")

# Load data
pop <- fread("data/population.csv")

# Set age groups
age_groups <- c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+")

run <- 77 #78

# Set output to use
output <- paste0("output/MCMCoutput",run,".RDS")

# Set burn-in
burnin <- 1500 #1 #

# Set number of posterior samples for age-decomposition plots
n_smpls <- 1000 #10 #

# Set seed for immune status plot
seed <- 1L

# Plot breakdown of immunity in the population over time
res <- plot_immune_status(output,pop,age_groups,burnin,n_smpls,seed)
p <- res$p
ggsave(paste0("output/pop_immune_status_by_age",run,".pdf"),p,height = 5,width = 11)
p1 <- res$p1
ggsave(paste0("output/pop_immune_status",run,".pdf"),p1,height = 4,width = 6)
tbl <- res$tbl
write.csv(tbl,paste0("output/table2_by_age_",run,".csv"),row.names = F)
tbl1 <- res$tbl1
write.csv(tbl1,paste0("output/table2_",run,".csv"),row.names = F)
