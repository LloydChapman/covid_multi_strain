create_priors <- function(beta_date,pars_info) {
    n_beta <- length(beta_date)
    beta_hps <- data.frame(type = rep("gamma",n_beta),
                           name = paste0("beta",seq_len(n_beta)),
                           scale = rep(0.02/4,n_beta),
                           shape = rep(4,n_beta),
                           shape1 = NA_real_,
                           shape2 = NA_real_)
    
    p_hps <- data.frame(type = "beta",
                        name = c("p_H_max","p_D_max"),
                        scale = NA_real_,
                        shape = NA_real_,
                        shape1 = 1, 
                        shape2 = 1)
    
    strain_hps <- data.frame(
        type = "null",
        name = c("rel_strain_transmission","start_date","strain_seed_date",
                 "rel_strain_transmission1","strain_seed_date1",
                 "rel_strain_transmission2","strain_seed_date2"),
        scale = NA_real_,
        shape = NA_real_,
        shape1 = NA_real_,
        shape2 = NA_real_
    )
    
    case_obs_hps <- data.frame(
        type = "beta",
        name = c("phi_cases","alpha_cases"),
        scale = NA_real_,
        shape = NA_real_,
        shape1 = 1,
        shape2 = 1
    )
    
    hosp_obs_hps <- data.frame(
        type = "beta",
        name = "alpha_hosp",
        scale = NA_real_,
        shape = NA_real_,
        shape1 = 1,
        shape2 = 1
    )
    
    death_obs_hps <- data.frame(
        type = "beta",
        name = "alpha_death",
        scale = NA_real_,
        shape = NA_real_,
        shape1 = 1,
        shape2 = 1
    )
    
    ret <- rbind(
        beta_hps,
        strain_hps[strain_hps$name %in% 
                       c("rel_strain_transmission","start_date","strain_seed_date"),],
        p_hps,
        strain_hps[strain_hps$name %in%
                       c("rel_strain_transmission1","strain_seed_date1"),],
        strain_hps[strain_hps$name %in% 
                       c("rel_strain_transmission2","strain_seed_date2"),],
        case_obs_hps,
        hosp_obs_hps,
        death_obs_hps)
    names(ret)[match(c("scale","shape","shape1","shape2"),names(ret))] <- 
        c("gamma_scale","gamma_shape","beta_shape1","beta_shape2")
    
    nms_expected <- unique(pars_info$name)
    nms_found <- unique(ret$name)
    msg <- setdiff(nms_expected, nms_found)
    if (length(msg) > 0) {
        stop(sprintf("Missing parameters, update priors (missing %s)",
                     paste(msg, collapse = ", ")))
    }
    extra <- setdiff(nms_found, nms_expected)
    if (length(extra)) {
        message(sprintf("Dropping %d unused priors: %s",
                        length(extra), paste(extra, collapse = ", ")))
        ret <- ret[ret$name %in% nms_expected, ]
    }
    rownames(ret) <- NULL
    
    ret
}
