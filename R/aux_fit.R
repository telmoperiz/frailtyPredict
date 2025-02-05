##### Auxiliary functions specific to estimation #####

#' Valid positivity constraint
#'
#' Get the positivity constraints for each function in the model.
#'
#' @param par_pos list with indexes of subparameters, see [SharedModel()].
#' @param model_funs list with function shapes for the model, see [SharedModel()].
#'
#' @return numeric vector with constraint parameter positions
#' @export
valid_posit_constraint <- function(par_pos, model_funs){

  # Initialize it to all FALSE (not needed constaint)
  L <- max(sapply(par_pos, max)) #length of parameter vector
  posit_cons <- vector(mode = 'logical', length = L)

  # Call the positivity constraint functions
  posit_cons[par_pos$sig]<-do.call(paste(model_funs$pdf, 'posit', sep = '_'),
                                   list())
  posit_cons[par_pos$a_d]<-do.call(paste(model_funs$hazard_d, 'posit', sep = '_'),
                                   list())
  posit_cons[par_pos$a_r]<-do.call(paste(model_funs$hazard_r, 'posit', sep = '_'),
                                   list())

  # Positivity constraint for piecewise Cs
  posit_cons[par_pos$pw_C] <- TRUE

  return(posit_cons)
}

#.............................................................................#
#.............................................................................#

#' Part of the loglikelihood that integrates w.r.t frailty (u)
#'
#' @param obs Observation unit information. See [individual_loglikelihood()].
#' @param model [SharedModel()] object with model information.
#'
#' @return numeric Integral part of the log-likelihood
#' @export
iloglik_intfrailty <- function(obs, model){

  ##Survival functions exponentiated to the corresponding cox index
  #Get Cox-part parameters
  beta<-param_coef(model)

  #Get the coefficients of the Cox indexes
  C_d <- as.vector(obs$ter_covs %*% beta$terminal)
  C_r <- as.vector(obs$rec_covs %*% beta$recurrent)

  #Exponentiate them
  C_d<-exp(C_d)
  C_r<-exp(C_r)

  #Get the Cox coeffients (for recurrent event) corresponding to T^d
  C_r_Td<-C_r[!uncensored(obs$rec_times)]

  # Times and number of recurrent events
  Tjs<-get_rec_times(obs$rec_times)

  #Get baseline survival functions at T^d
  S_d<-base_survival(model, t=times(obs$ter_time), process = 'terminal')
  S_r<-surv_given_history(model, t=times(obs$ter_time), Tjs=Tjs)

  #Exponentiate the baseline survival functions
  Sexp_d<-S_d^C_d
  Sexp_r<-S_r^C_r_Td #C_r corresponding to covariate values at T^d

  ## Other integration info
  #Check if the terminal event was observed or censored
  uncens<-as.integer(uncensored(obs$ter_time))

  #Number of recurrent events
  J<-length(Tjs)

  #Frailty parameters: alpha and pdf shape parameters
  parfrail<-param_frailty(model)

  #Get the integral part of the likelihood
  likInt <- surv_intfrailty(S_d1=Sexp_d, S_d2=0, S_r=Sexp_r,
                            expo_u= uncens * parfrail['alpha'] + J, alpha=parfrail['alpha'],
                            pdfinfo=list(func=model$funtion_shapes$pdf,
                                         param = parfrail[grepl('shape', names(parfrail))]),
                            intcontrol=model$optim_control)

  #Take logs
  loglikInt<-log(likInt)

  return(loglikInt)
}

#.............................................................................#
#.............................................................................#

#' Derivative of the part of the loglikelihood with the integral wrt frailty
#'
#' @param obs Observation unit information. See [individual_loglikelihood()].
#' @param model [SharedModel()] object with model information.
#'
#' @return numeric vector Derivative of the integral part of the log-likelihood
#' @export
iloglik_intfrailty_gradient <- function(obs, model){

  #Initialize the derivative vector
  Deriv<-rep(NA, length = length(model$par))

  ## Common computations:

  #Cox indexes
  #Get Cox-part parameters
  beta<-param_coef(model)

  #Get the coefficients of the Cox indexes
  C_d <- as.vector(obs$ter_covs %*% beta$terminal)
  C_r <- as.vector(obs$rec_covs %*% beta$recurrent)

  #Exponentiate them
  C_d<-exp(C_d)
  C_r<-exp(C_r)

  # Save only the one corresponding to T_d
  C_r<-C_r[!uncensored(obs$rec_times)]

  #Times and number of recurrent events
  Tjs<-get_rec_times(obs$rec_times)
  J<-length(Tjs)

  # Compute the value of the integral
  log_intfrailty<-iloglik_intfrailty(obs, model)
  intfrailty<-exp(log_intfrailty)

  #Survival functions and gradients at T^d
  S_d<-base_survival(model, t=times(obs$ter_time), process = 'terminal',
                     gradient=TRUE)
  grad_d<-attr(S_d, "gradient") #gradient as attribute

  S_r<-surv_given_history(model, t=times(obs$ter_time), Tjs=Tjs, gradient=TRUE)
  grad_r<-attr(S_r, "gradient") #gradient as attribute

  # Gradient wrt piecewise constants
  if (model$rec_timescale == 'piecewise-renewal'){
    grad_r_C <- attr(S_r, "gradient_C")
  }

  #Check if the terminal event was observed or censored
  uncens<-as.integer(uncensored(obs$ter_time))

  #Frailty parameters: alpha and pdf shape parameters
  parfrail<-param_frailty(model)

  #Derivative of the integral
  dlikInt <- surv_intfrailty_gradient(S_d1=S_d^C_d, S_d2=0, S_r=S_r^C_r,
                                      expo_u= uncens * parfrail['alpha'] + J, alpha=parfrail['alpha'],
                                      pdfinfo=list(func=model$funtion_shapes$pdf,
                                                   param = parfrail[grepl('shape', names(parfrail))]),
                                      intcontrol=model$optim_control,
                                      terms=c('S_d1', 'S_r', 'expo_u', 'alpha', 'pdfparam'))

  ## Derivative wrt the Terminal event model coefficients

  S_d_exp <- S_d^C_d # Exponentiate

  # What multiplies each covariate
  multiplier_d <- S_d_exp  * log(S_d_exp) * dlikInt$S_d1 / intfrailty

  # Score wrt beta_d
  Score_d <- multiplier_d * as.vector(obs$ter_covs)

  # If result is NaN but covariates are zero, replace by zero
  Score_d[obs$ter_covs == 0 & is.nan(Score_d)] <- 0

  # Derivatives
  Deriv[model$par_pos$beta_d] <- Score_d

  ## Derivative wrt the Recurrent event model coefficients

  S_r_exp <- S_r^C_r # Exponentiate

  # What multiplies each covariate
  multiplier_r <- S_r_exp  * log(S_r_exp) * dlikInt$S_r / intfrailty

  # Covariates at T_d
  rec_covs_Td <- obs$rec_covs[!uncensored(obs$rec_times),]

  # Score wrt beta_r
  Score_r <- multiplier_r * rec_covs_Td

  # If result is NaN but covariates are zero, replace by zero
  Score_r[rec_covs_Td == 0 & is.nan(Score_r)] <- 0

  # Derivative
  Deriv[model$par_pos$beta_r] <- Score_r

  ## Derivative wrt its Terminal Hazard parameters
  Deriv[model$par_pos$a_d] <- C_d * S_d^(C_d-1) * grad_d * dlikInt$S_d1 / intfrailty

  ## Derivative wrt its Recurrent Hazard parameters
  Deriv[model$par_pos$a_r] <- C_r * S_r^(C_r-1) * grad_r * dlikInt$S_r / intfrailty

  ## Derivative wrt Recurrent Hazard piecewise constants
  if (model$rec_timescale == 'piecewise-renewal'){
    Deriv[model$par_pos$pw_C] <- C_r * S_r^(C_r-1) * grad_r_C * dlikInt$S_r / intfrailty
  }

  ## Derivative wrt alpha
  Deriv[model$par_pos$alpha]<-( uncens * dlikInt$expo_u + dlikInt$alpha) / intfrailty

  ## Derivative wrt the params of frailty distribution
  Deriv[model$par_pos$sig] <- dlikInt$pdfparam / intfrailty

  #Return the gradient vector
  return(Deriv)
}

#.............................................................................#
#.............................................................................#

#' Part of the log-likelihood that does not depend on the frailty
#'
#' @param obs Observation unit information. See [individual_loglikelihood()].
#' @param model [SharedModel()] object with model information.
#'
#' @return Part of the log-likelihood that does not depend on the frailty
#' @export
iloglik_nofrailty <- function(obs, model){

  #Parameters of the Cox part of the model
  beta<-param_coef(model)

  #Get the Cox indexes
  C_d <- as.vector(obs$ter_covs %*% beta$terminal)
  C_r <- as.vector(obs$rec_covs %*% beta$recurrent)

  #Get only the Cox indexes corresponding to T_j^r's
  #This removes the last which corresponds to T^d or censoring
  C_r_js<-C_r[uncensored(obs$rec_times)]

  # Get recurrent event times
  Tjs<-get_rec_times(obs$rec_times)

  #Check if the terminal event was observed or censored
  uncens<-as.integer(uncensored(obs$ter_time))

  #Get the hazards given history a T^d (terminal) and T_j^r's
  h_d<-base_hazard(model, t=times(obs$ter_time), process = 'terminal')
  h_r<-hazard_given_history(model, t=Tjs, Tjs=Tjs)

  #Sum of log hazards w/ the corresponding Cox indexes
  hC_r<- log(h_r) + C_r_js

  # Hazard part of the likelihood (w/out frailty)
  haz_part<- uncens * (C_d + log(h_d)) +  sum(hC_r)
  return(haz_part)
}

#.............................................................................#
#.............................................................................#

#' Derivative of the part of the loglikelihood independent of frailty
#'
#' @param obs Observation unit information. See [individual_loglikelihood()].
#' @param model [SharedModel()] object with model information.
#'
#' @return numeric vector Derivative of the part of the loglikelihood
#' independent of frailty
#' @export
iloglik_nofrailty_gradient <- function(obs, model){

  #Initialize the derivative vector
  Deriv<-rep(NA, length = length(model$par))

  ## Derivative wrt the Terminal event model coefficients
  #Check if the terminal event was observed or censored
  uncens<-as.integer(uncensored(obs$ter_time))

  #Derivative are the covariates (only entering for uncensored observations)
  Deriv[model$par_pos$beta_d]<-as.vector(obs$ter_covs) * uncens

  ## Derivative wrt the Recurrent event model coefficients

  #Get only the covariate matrix corresponding to recurrent event times
  # (remove the last which corresponds to T^d or censoring)
  # subset function helps keep the matrix structure
  covs_js<-subset(obs$rec_covs, uncensored(obs$rec_times))

  #Time sum of the covariates (sum for each column)
  Deriv[model$par_pos$beta_r]<-as.vector(colSums(covs_js))

  ## Derivative of log(Terminal hazard) wrt its parameters
  #Get the hazard and gradient
  h_d<-base_hazard(model, t=times(obs$ter_time), process = 'terminal',
                   gradient=TRUE)
  grad_d<-attr(h_d, "gradient") #gradient as attribute

  #Derivative of log hazard (only if not censored)
  Deriv[model$par_pos$a_d]<- (grad_d/h_d) * uncens

  ## Derivative of  log(Recurrent hazard) wrt its parameters
  #Get the recurrent event times
  Tjs<-get_rec_times(obs$rec_times)

  #If there is at least one recurrent event
  if (length(Tjs)>0){

    #Get the hazards and gradients
    h_r<-hazard_given_history(model, t=Tjs, Tjs=Tjs, gradient=TRUE)
    grad_r<-attr(h_r, "gradient") #gradient as attribute

    #For each T_j^r (in rows), get the ratio between gradient and hazard
    #replicate gives a matrix with copies of h_r in columns
    ratios_r<-grad_r/replicate(n=ncol(grad_r), h_r)

    #Get the derivative: sum across rows
    Deriv[model$par_pos$a_r]<-colSums(ratios_r)

    # Gradient wrt piecewise constants
    if (model$rec_timescale == 'piecewise-renewal'){

      grad_r_C <- attr(h_r, "gradient_C")

      ratios_r_C <- grad_r_C / replicate(n=ncol(grad_r_C), h_r)

      # This of the positions where Tj falls (by definition)
      Deriv[model$par_pos$pw_C] <- colSums(ratios_r_C)
    }

    #If there is no recurrent event, the term does not appear in the loglikelihood
  } else {

    Deriv[model$par_pos$a_r]<-0

    if (model$rec_timescale == 'piecewise-renewal'){
      Deriv[model$par_pos$pw_C] <- 0
    }
  }

  ## Derivative of the hazard part wrt alpha is zero
  Deriv[model$par_pos$alpha]<-0

  ## Derivative of the hazard part wrt the params of frailty distr. is zero
  Deriv[model$par_pos$sig]<-0

  #Return the gradient vector
  return(Deriv)
}

#.............................................................................#
#.............................................................................#

#' Log-likelihood computation for an individual
#'
#' Computes the log-likelihood for an observation unit in the sample.
#'
#' @param obs Observation unit information. List with components:
#' * `$ter_time` = Time and censoring indicator for terminal process.
#' * `$ter_covs` = Matrix with covariates in terminal process model.
#' * `$rec_times` =  Times and censoring indicators for recurrent process.
#' * `$rec_covs` = Matrix with covariates in recurrent process model.
#' @param model [SharedModel()] object with model information.
#'
#' @return numeric with log-likelihooh value.
individual_loglikelihood <- function(obs, model){

  #Part that does not depend on the frailty term
  nofrail<-iloglik_nofrailty(obs, model)

  #Part that integrates w.r.t. the frailty term
  intfrail<-iloglik_intfrailty(obs, model)

  #Sum the two terms
  return(nofrail + intfrail)
}

#.............................................................................#
#.............................................................................#

#' Score computation for and individual
#'
#' @param obs Observation unit information. See [individual_loglikelihood()].
#' @param model [SharedModel()] object with model information.
#'
#' @return numeric vector Score
#' @export
individual_score <- function(obs, model){

  #Part that does not depend on the frailty term
  nofrail<-iloglik_nofrailty_gradient(obs, model)

  #Part that integrates w.r.t. the frailty term
  intfrail<-iloglik_intfrailty_gradient(obs, model)

  #Sum the two terms
  return(nofrail + intfrail)
}

#.............................................................................#
#.............................................................................#

#' Log-likelihood function for the whole data
#'
#' Computes the log-likelihood for the whole data given a parameter guess nu.
#' Input for optimization via [stats::nlm()].
#'
#' @param nu Guess value of the parameter
#' @param OBSL List with all observations
#' @param par_pos Position of each specific parameter in nu
#' @param posit_cons Elements in nu subject to positivity constraints
#' @param model_funs List with model function shapes
#' @param rec_timescale Timescale for the recurrent events
#' @param int_mode How to perform the integral
#' @param MC_N Monte Carlo sample size
#' @param GC_nodes Gaussian Quadrature number of nodes
#' @param int_seq Additional sequence needed for integration in some modes
#' @param is_control Control for Importance Sampling: distribution and parameters
#' @param anal.grand Whether to use the analytic gradient or no
#' @param BHHH.hessian Whether to compute BHHH style hessian or not
#'
#' @return Log-likelihood as numeric with attribute 'gradient' and 'hessian'
#' @export
loglikelihood <- function(nu, OBSL, par_pos, posit_cons, model_funs,
                          rec_timescale, int_mode, MC_N, GC_nodes,
                          int_seq, is_control, anal.grad, BHHH.hessian){

  #Transform to impose positivity
  theta<-posit_constraint(nu, posit_cons, type="direct")

  # Get (Q)MC sample for integration
  if (int_mode == 'MC'){
    MC_samp <- dist.draw(n=MC_N, pdf=model_funs$pdf, param = theta[par_pos$sig])

  } else if(int_mode == 'QMC'){
    MC_samp <- dist.transform(seq=int_seq, pdf=model_funs$pdf,
                              param = theta[par_pos$sig])

  } else if (int_mode %in% c('ISMC', 'ISQMC'))  {
    MC_samp <- int_seq

  } else if (int_mode == 'GQ'){
    MC_samp <- NULL
  }

  #Build the model
  model<-SharedModel(par=theta, varmat = NULL, par_pos=par_pos,
                     var_names=list(terminal=colnames(OBSL[[1]]$ter_covs),
                                    recurrent=colnames(OBSL[[1]]$rec_covs)),
                     OBSL=OBSL, logLikelihood = NULL,
                     function_shapes=model_funs, rec_timescale=rec_timescale,
                     Calls=NULL, optim_control=list(int_mode = int_mode,
                                                    MC_N = MC_N,
                                                    MC_samp = MC_samp,
                                                    GC_nodes = GC_nodes,
                                                    is_control = is_control)
  )

  #Get the the log-likelihood
  ell <- logLik(model)

  # Return minus the likelihood to maximize
  res <- -ell

  if (anal.grad){
    #Get the Scores for each observation
    #number of observations times number of parameters matrix
    Scores<-scores(model)

    #gradient of the positivity transformation at nu
    dF<-posit_constraint(nu, posit_cons, type = "gradient")

    #Account for positivity constraint
    Scores_posit <- Scores %*% dF

    # Set gradient attribute for minimization
    attr(res, "gradient") <- - colSums(Scores_posit)
  }

  if (BHHH.hessian){
    # Get the Hessian of the model
    Hessian <- BHHH_hessian(model, S=Scores)

    # second derivatives of the positivity transformation
    HF<-posit_constraint(nu, posit_cons, type = "diag_hessian")

    # Set hessian attribute for minimization
    attr(res, "hessian") <- - dF %*% Hessian %*% dF - colSums(Scores) * HF
  }

  return(res)
}

#.............................................................................#
#.............................................................................#

#' Impose positivity constraint
#'
#' Function to impose and reverse the square transformation
#' for positivity constraint.
#'
#' @param par Parameter of the model, see [SharedModel()].
#' @param impose_ind Indexes to impose the constraint, see [SharedModel()].
#' @param type Either `"direct"` (positive-valued function: `x^2`), `"reverse"`
#' (inverse of the function), `"gradient"` (gradient of the transformation), or
#' `"diag_hessian"` (diagonal of the Hessian of the transformation)
#'
#' @return numeric vector Depending on `type`.
#' @export
posit_constraint <- function(par, impose_ind, type="direct"){

  #Transform the relevant indexes
  if (type=="direct"){

    out_par<-par
    out_par[impose_ind] <- par[impose_ind]^2

  } else if (type == "reverse"){

    out_par<-par
    out_par[impose_ind] <- sqrt(par[impose_ind])

    #Jacobia matrix of direct transformation (first derivatives in diag)
  } else if (type == "gradient"){

    out_par <- diag(length(par)) #identity matrix of len times len
    diag(out_par)[impose_ind]<-2*par[impose_ind]

    # Second derivative of each direct transformation in the diagonal
  } else if (type == "diag_hessian"){
    out_par <- diag(length(par)) #identity matrix of len times len
    diag(out_par)[impose_ind]<- 2
  }

  return(out_par)
}

#.............................................................................#
#.............................................................................#
