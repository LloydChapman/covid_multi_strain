library(data.table)
library(qs)
library(abind)

source("R/utils.R")
source("R/date.R")
source("R/baseline.R")
source("R/info.R")
source("R/priors.R")
source("R/proposal.R")
source("R/vaccination.R")
source("R/covid_multi_strain.R")

# Set model type
# "NB" for fitting to cases with negative binomial likelihood
# "BB" for fitting to testing data with beta-binomial likelihood
model_type <- "NB"
# model_type <- "BB"

# Set epoch date
epoch_dates <- "2021-11-21" #c("2021-11-21","2021-12-16")

# Set assumptions
# "central" for default central values
# "alt_contact_matrix" for using alternative contact matrix for Fiji
# assumptions <- "central"
assumptions <- "alt_contact_matrix"

baseline <- create_baseline(model_type,epoch_dates,assumptions)
saveRDS(baseline,"parameters/base.rds")

beta_date <- baseline$beta_date

info <- create_info(beta_date,model_type)
write.csv(info,"parameters/info.csv",row.names = F)

prior <- create_priors(beta_date,info)
write.csv(prior,"parameters/prior.csv",row.names = F)

proposal <- create_proposal(beta_date,model_type)
write.csv(proposal,"parameters/proposal.csv",row.names = F)

