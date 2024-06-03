##### Auxiliary functions specific to prediction #####

create_predict_json <- function(model, JSON_filename, t, ws, individuals,
                                CI_method, sig, MC_B, numint_nodes){

  # Draw sample from asymptotic distribution (if Monte Carlo method)
  if (CI_method == 'MC'){
    thetas <- asympt_draw(model, n = MC_B)
  } else {
    thetas <- NULL
  }

  # Initialize list
  JSON_list <- list()

  for (ind in individuals){

    # Check formatting and add defaults
    ind_f <- plot_ind_formatting(ind, CI_method)

    # Consider only history prior to t
    hist_prev <- ind_f$hist[ind_f$hist <= t]

    # Get predicted risk
    probs_l<-lapply(ws, predict, obj=model, t=t,
                    covs = ind_f$covs,
                    hist = hist_prev,
                    type = ind_f$type,
                    CI = ind_f$ci_f,
                    sig=sig, MC_B=MC_B, numint_nodes=numint_nodes,
                    par_sample=thetas
    )

    # Predictions
    probs <- unlist(probs_l)
    lower <- sapply(probs_l, attr, which='CI_lower')
    upper <- sapply(probs_l, attr, which='CI_upper')

    # Set to NA if they are null
    if (all(sapply(lower, is.null))) lower <- jsonlite::unbox(NA)
    if (all(sapply(upper, is.null))) upper <- jsonlite::unbox(NA)

    #JSON style output
    JSON_list <- rlist::list.append(JSON_list,
                                    list(rec_times = hist_prev,
                                         pred_point = jsonlite::unbox(t),
                                         pred_times = t + ws,
                                         pred = probs,
                                         CI_lower = lower,
                                         CI_upper = upper,
                                         display_name = jsonlite::unbox(ind$name)
                                    )
    )
  }

  # Write the temporal json
  write(jsonlite::toJSON(JSON_list, na = 'null', auto_unbox = FALSE, pretty = TRUE),
        file = JSON_filename)

  return(TRUE)
}

plot_ind_formatting <- function(individual, CI_method){

  # Covariates (mandarory)
  if (! 'covs' %in% names(individual)){
    stop('Must include individual covariates. Specify individual$covs.',
         call. = FALSE)
  }

  # Must include terminal and recurrent covariates
  if (! all('terminal' %in% names(individual$covs),
            'recurrent' %in% names(individual$covs)) ){
    stop('Must include terminal and recurrent covariates. Specify
         individual$covs$terminal and individual$covs$recurrent.',
         call. = FALSE)
  }

  # Default value for history
  if (! 'hist' %in% names(individual)){
    individual$hist <- c()
  }

  # Default value for type
  if (! 'type' %in% names(individual)){
    individual$type <- 'unconditional'
  }

  # Default value for ci
  if ('ci' %in% names(individual)){

    if (individual$ci){
      individual$ci_f <- CI_method
    } else {
      individual$ci_f <- 'No'
    }

  } else {
    individual$ci_f <- 'No'
  }

  # Default value for name
  if (! 'name' %in% names(individual)){
    individual$name <- NA
  }

  return(individual)
}

pred_den <- function(St_d, S_r, J, alpha, parbeta, covs, pdfinfo,
                     intcontrol, par_pos, type, gradient){

  #Cox index for terminal and exponentiate it
  C_d <- sum(as.vector(covs$terminal) * parbeta$terminal)
  C_d<-exp(C_d)

  #Cox index for recurrent model and exponentiate it
  C_r <- sum(as.vector(covs$recurrent) * parbeta$recurrent)
  C_r<-exp(C_r)

  #Probability not conditioning on history
  if (type=='unconditional'){

    #Denominator: Sexp1_d^(u^alpha)
    DEN<-surv_intfrailty(S_d1=St_d^C_d, S_d2=0, S_r=1,
                         expo_u=0, alpha=alpha,
                         pdfinfo=pdfinfo, intcontrol=intcontrol)

    # Conditional on recurrent event history
  } else if (type=='conditional'){

    #Denominator: u^J * Sexp1_d^(u^alpha) * Sexp_r^u
    DEN<-surv_intfrailty(S_d1=St_d^C_d, S_d2=0, S_r=S_r^C_r,
                         expo_u=J, alpha=alpha,
                         pdfinfo=pdfinfo, intcontrol=intcontrol)
  }

  # Compute gradient if needed
  if (gradient){

    # List Cox indexes
    coxind <- list(terminal=C_d, recurrent=C_r)

    grad <- pred_den_grad(St_d, S_r, J, alpha, coxind, covs,
                          pdfinfo, intcontrol, par_pos, type)

    attr(DEN, 'gradient') <- grad
  }

  return(DEN)
}

#.............................................................................#
#.............................................................................#

pred_den_grad <- function(St_d, S_r, J, alpha, coxind, covs, pdfinfo,
                          intcontrol, par_pos, type){

  #Initialize the derivative vector
  Deriv<-rep(NA, length = max(sapply(par_pos, max)))

  # Cox indexes
  C_d <- coxind$terminal
  C_r <- coxind$recurrent

  # Gradients of survival functions
  St_d_grad <- attr(St_d, 'gradient')
  S_r_grad <- attr(S_r, 'gradient')

  # Gradient of integral
  if (type == 'unconditional'){

    int_grad <- surv_intfrailty_gradient(S_d1=St_d^C_d, S_d2=0, S_r=1,
                                         expo_u=0, alpha=alpha,
                                         pdfinfo=pdfinfo, intcontrol=intcontrol,
                                         terms=c('S_d1', 'alpha', 'pdfparam'))

  } else if (type=='conditional'){

    int_grad <- surv_intfrailty_gradient(S_d1=St_d^C_d, S_d2=0, S_r=S_r^C_r,
                                         expo_u=J, alpha=alpha,
                                         pdfinfo=pdfinfo, intcontrol=intcontrol,
                                         terms=c('S_d1', 'S_r', 'alpha', 'pdfparam'))
  }

  ## Derivative wrt the Terminal event model coefficients

  # Exponentiate
  St_d_exp <- St_d^C_d

  # What multiplies each covariate
  multiplier_d <- St_d_exp  * log(St_d_exp) * int_grad$S_d1

  # Deriv wrt beta_d
  Deriv_bd <- multiplier_d * as.vector(covs$terminal)

  # If result is NaN but covariates are zero, replace by zero
  Deriv_bd[covs$terminal == 0 & is.nan(Deriv_bd)] <- 0

  # Derivatives
  Deriv[par_pos$beta_d] <- Deriv_bd

  ## Derivative wrt the Recurrent event model coefficients
  if (type == 'unconditional'){

    Deriv[par_pos$beta_r] <- 0

  } else if (type == 'conditional'){

    S_r_exp <- S_r^C_r # Exponentiate

    # What multiplies each covariate
    multiplier_r <- S_r_exp  * log(S_r_exp) * int_grad$S_r

    # Deriv wrt beta_r
    Deriv_br <- multiplier_r * as.vector(covs$recurrent)

    # If result is NaN but covariates are zero, replace by zero
    Deriv_br[covs$recurrent == 0 & is.nan(Deriv_br)] <- 0

    # Derivatives
    Deriv[par_pos$beta_r] <- Deriv_br
  }

  ## Derivative wrt its Terminal Hazard parameters
  Deriv[par_pos$a_d] <- C_d * St_d^(C_d-1) * St_d_grad * int_grad$S_d1

  ## Derivative wrt its Recurrent Hazard parameters
  if (type == 'unconditional'){

    Deriv[par_pos$a_r] <- 0

  } else if (type == 'conditional'){

    Deriv[par_pos$a_r] <- C_r * S_r^(C_r-1) * S_r_grad * int_grad$S_r
  }

  ## Derivative wrt alpha
  Deriv[par_pos$alpha] <- int_grad$alpha

  ## Derivative wrt the params of frailty distribution
  Deriv[par_pos$sig] <- int_grad$pdfparam

  return(Deriv)
}

#.............................................................................#
#.............................................................................#

pred_num <- function(St_d, Sw_d, S_r, J, alpha, parbeta, covs, pdfinfo,
                     intcontrol, par_pos, type, gradient){

  #Cox index for terminal and exponentiate it
  C_d <- sum(as.vector(covs$terminal) * parbeta$terminal)
  C_d<-exp(C_d)

  #Cox index for recurrent model and exponentiate it
  C_r <- sum(as.vector(covs$recurrent) * parbeta$recurrent)
  C_r<-exp(C_r)

  #Probability not conditioning on history
  if (type=='unconditional'){

    #Numerator: Sexp1_d^(u^alpha) - Sexp2_d^(u^alpha)
    NUM<-surv_intfrailty(S_d1=St_d^C_d, S_d2=Sw_d^C_d, S_r=1,
                         expo_u=0, alpha=alpha,
                         pdfinfo=pdfinfo, intcontrol=intcontrol)

  # Conditional on recurrent event history
  } else if (type=='conditional'){

    #Numerator: u^J * ( Sexp1_d^(u^alpha) - Sexp2_d^(u^alpha) ) * Sexp_r^u
    NUM<-surv_intfrailty(S_d1=St_d^C_d, S_d2=Sw_d^C_d, S_r=S_r^C_r,
                         expo_u=J, alpha=alpha,
                         pdfinfo=pdfinfo, intcontrol=intcontrol)

  }

  # Compute gradient if needed
  if (gradient){

    # List Cox indexes
    coxind <- list(terminal=C_d, recurrent=C_r)

    grad <- pred_num_grad(St_d, Sw_d, S_r, J, alpha, coxind, covs,
                          pdfinfo, intcontrol, par_pos, type)

    attr(NUM, 'gradient') <- grad
  }

  return(NUM)
}

#.............................................................................#
#.............................................................................#

pred_num_grad <- function(St_d, Sw_d, S_r, J, alpha, coxind, covs, pdfinfo,
                          intcontrol, par_pos, type){

  #Initialize the derivative vector
  Deriv<-rep(NA, length = max(sapply(par_pos, max)))

  # Cox indexes
  C_d <- coxind$terminal
  C_r <- coxind$recurrent

  # Gradients of survival functions
  St_d_grad <- attr(St_d, 'gradient')
  Sw_d_grad <- attr(Sw_d, 'gradient')
  S_r_grad <- attr(S_r, 'gradient')

  # Gradient of integral
  if (type == 'unconditional'){

    int_grad <- surv_intfrailty_gradient(S_d1=St_d^C_d, S_d2=Sw_d^C_d, S_r=1,
                                         expo_u=0, alpha=alpha,
                                         pdfinfo=pdfinfo, intcontrol=intcontrol,
                                         terms=c('S_d1', 'S_d2', 'alpha', 'pdfparam'))

  } else if (type=='conditional'){

    int_grad <- surv_intfrailty_gradient(S_d1=St_d^C_d, S_d2=Sw_d^C_d, S_r=S_r^C_r,
                                         expo_u=J, alpha=alpha,
                                         pdfinfo=pdfinfo, intcontrol=intcontrol,
                                         terms=c('S_d1', 'S_d2', 'S_r', 'alpha', 'pdfparam'))
  }

  ## Derivative wrt the Terminal event model coefficients

  # Exponentiate
  St_d_exp <- St_d^C_d
  Sw_d_exp <- Sw_d^C_d

  # What multiplies each covariate
  multiplier_d <- St_d_exp  * log(St_d_exp) * int_grad$S_d1 +
                  Sw_d_exp  * log(Sw_d_exp) * int_grad$S_d2

  # Deriv wrt beta_d
  Deriv_bd <- multiplier_d * as.vector(covs$terminal)

  # If result is NaN but covariates are zero, replace by zero
  Deriv_bd[covs$terminal == 0 & is.nan(Deriv_bd)] <- 0

  # Derivatives
  Deriv[par_pos$beta_d] <- Deriv_bd

  ## Derivative wrt the Recurrent event model coefficients
  if (type == 'unconditional'){

    Deriv[par_pos$beta_r] <- 0

  } else if (type == 'conditional'){

    S_r_exp <- S_r^C_r # Exponentiate

    # What multiplies each covariate
    multiplier_r <- S_r_exp  * log(S_r_exp) * int_grad$S_r

    # Deriv wrt beta_r
    Deriv_br <- multiplier_r * as.vector(covs$recurrent)

    # If result is NaN but covariates are zero, replace by zero
    Deriv_br[covs$recurrent == 0 & is.nan(Deriv_br)] <- 0

    # Derivatives
    Deriv[par_pos$beta_r] <- Deriv_br
  }

  ## Derivative wrt its Terminal Hazard parameters
  Deriv[par_pos$a_d] <- C_d * St_d^(C_d-1) * St_d_grad * int_grad$S_d1 +
                        C_d * Sw_d^(C_d-1) * Sw_d_grad * int_grad$S_d2

  ## Derivative wrt its Recurrent Hazard parameters
  if (type == 'unconditional'){

    Deriv[par_pos$a_r] <- 0

  } else if (type == 'conditional'){

    Deriv[par_pos$a_r] <- C_r * S_r^(C_r-1) * S_r_grad * int_grad$S_r
  }

  ## Derivative wrt alpha
  Deriv[par_pos$alpha] <- int_grad$alpha

  ## Derivative wrt the params of frailty distribution
  Deriv[par_pos$sig] <- int_grad$pdfparam

  return(Deriv)
}

#.............................................................................#
#.............................................................................#
