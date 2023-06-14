source("R/utils.R")
source("R/date.R")
source("R/baseline.R")
source("R/info.R")
source("R/priors.R")
source("R/proposal.R")
source("R/vaccination.R")
source("R/covid_multi_strain.R")

# model_type <- "NB"
model_type <- "BB"

# Set epoch date
epoch_dates <- c("2021-11-21","2021-12-16")

baseline <- create_baseline(model_type,epoch_dates)
saveRDS(baseline,"parameters/base.rds")

beta_date <- baseline$beta_date

info <- create_info(beta_date,model_type)
write.csv(info,"parameters/info.csv",row.names = F)

prior <- create_priors(beta_date,info)
write.csv(prior,"parameters/prior.csv",row.names = F)

proposal <- create_proposal(beta_date,model_type)
write.csv(proposal,"parameters/proposal.csv",row.names = F)

