pars_pmcmc_load <- function(path, info = "info.csv", prior = "prior.csv",
                            proposal = "proposal.csv"){
    list(info = read_csv(file.path(path, info)),
         prior = read_csv(file.path(path, prior)),
         proposal = read_csv(file.path(path, proposal)))
}

make_prior <- function(d) {
    if (d$type == "gamma") {
        ## TODO: as_duration was droppd from here as never used, but if it
        ## is, then we'd transform p to 1/p
        shape <- d$gamma_shape
        scale <- d$gamma_scale
        function(p) {
            dgamma(p, shape = shape, scale = scale, log = TRUE)
        }
    } else if (d$type == "beta") {
        shape1 <- d$beta_shape1
        shape2 <- d$beta_shape2
        function(p) {
            dbeta(p, shape1 = shape1, shape2 = shape2, log = TRUE)
        }
    } else if (d$type == "null") {
        NULL
    } else {
        stop("Unknown prior type")
    }
}

sample_prior <- function(d,i){
    if (d$type == "gamma"){
        shape = d$gamma_shape
        scale = d$gamma_scale
        rgamma(1, shape = shape, scale = scale)
    } else if (d$type == "beta"){
        shape1 = d$beta_shape1
        shape2 = d$beta_shape2
        rbeta(1, shape1 = shape1, shape2 = shape2)
    } else if (d$type == "null"){
        if (i$integer){
            sample(i$min:i$max, 1)
        } else {
            runif(1, min = i$min, max = i$max)    
        }
    } else {
        stop("Unknown prior type")
    }
}

pars_mcmc <- function(info, prior, proposal, transform) {
    pars_mcmc <- Map(
        pmcmc_parameter,
        name = info$name,
        initial = info$initial,
        min = info$min,
        max = info$max,
        integer = info$integer,
        prior = lapply(split(prior, prior$name)[unique(prior$name)], make_prior))
    
    ret <- pmcmc_parameters$new(pars_mcmc, proposal, transform)
    
    ## Try and transform a single case and see if it works:
    ret$model(ret$initial())
    
    ret
}

pars_proposal <- function(info, proposal){
    proposal <- as.matrix(proposal[match(info$name, proposal$name), info$name])
    rownames(proposal) <- info$name
    proposal
}
