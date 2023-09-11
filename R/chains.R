chains_combine <- function(..., samples = list(...)) {
    pars <- lapply(samples, "[[", "pars")
    probabilities <- lapply(samples, "[[", "probabilities")
    state <- lapply(samples, "[[", "state")
    trajectories <- lapply(samples, "[[", "trajectories")
    
    chain <- rep(seq_along(samples), each = nrow(samples[[1]]$pars))
    pars <- array_bind(arrays = pars, dimension = 1)
    probabilities <- array_bind(arrays = probabilities, dimension = 1)
    
    if (is.null(state[[1]])) {
        state <- NULL
    } else {
        state <- array_bind(arrays = state)
    }
    
    if (is.null(trajectories[[1]])) {
        trajectories <- NULL
    } else {
        trajectories <- combine_state(trajectories)
    }
    
    list(pars = pars, 
         probabilities = probabilities, 
         state = state, 
         trajectories = trajectories,
         chain = chain)
}

combine_state <- function(x) {
    base <- lapply(x, function(el) el[names(el) != "state"])
    if (length(unique(base)) != 1L) {
        stop(sprintf("%s data is inconsistent", deparse(substitute(x))))
    }
    
    state <- lapply(x, "[[", "state")
    state <- array_bind(arrays = state, dimension = length(dim(state[[1]])) - 1)
    
    ret <- x[[1]]
    ret$state <- state
    ret
}