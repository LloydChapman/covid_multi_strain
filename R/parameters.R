fit_pars_load <- function(path, assumptions){
    parameters <- pars_pmcmc_load(path)
    info <- parameters$info
    prior <- parameters$prior
    proposal <- pars_proposal(info, parameters$proposal)
    
    dat <- load_transform(path, assumptions)
    base <- dat$base
    transform <- dat$transform
    
    mcmc <- pars_mcmc(info, prior, proposal, transform)
    
    list(info = info,
         prior = prior,
         proposal = proposal,
         transform = transform,
         raw = parameters,
         base = base,
         mcmc = mcmc)
}

load_transform <- function(path, assumptions){
    e <- new.env()
    sys.source(file.path(path, "transform.R"), e)
    stopifnot(is.function(e$make_transform),
              is.function(e$apply_assumptions))
    make_transform <- e$make_transform
    base <- e$apply_assumptions(
        readRDS(file.path(path, "base.rds")), assumptions)
    transform <- make_transform(base)
    list(base = base, transform = transform)
}