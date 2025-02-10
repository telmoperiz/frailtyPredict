##### Builder for SharedModel S3 class #####

#' Build a SharedModel object
#'
#' SharedModel S3 objects contain fit information and information regarding
#' control of the optimization process. See Details for a list of S3 methods.
#'
#' @details
#'
#' List of methods:
#' * [asymptvar.SharedModel()]
#' * [asympt_draw.SharedModel()]
#' * [base_hazard.SharedModel()]
#' * [base_survival.SharedModel()]
#' * [BHHH_hessian.SharedModel()]
#' * [check_intsample.SharedModel()]
#' * [dist_relevance_test.SharedModel()]
#' * [get_rec_number.SharedModel()]
#' * [hazard_given_history.SharedModel()]
#' * [logLik.SharedModel()]
#' * [max_timespan.SharedModel()]
#' * [optim_details.SharedModel()]
#' * [param_coef.SharedModel()]
#' * [param_frailty.SharedModel()]
#' * [param_hazard.SharedModel()]
#' * [param_tables.SharedModel()]
#' * [predict.SharedModel()]
#' * [predict_plot.SharedModel()]
#' * [print.SharedModel()]
#' * [rec_boxplot.SharedModel()]
#' * [scores.SharedModel()]
#' * [summary.SharedModel()]
#' * [surv_given_history.SharedModel()]
#' * [toLatex.SharedModel()]
#' * [Wald_test.SharedModel()]
#'
#' @param par Parameter of the model
#' @param varmat Asymptotic covariance matrix.
#' @param par_pos list with indexes of subparameters:
#' * `$alpha`= frailty parameter.
#' * `$sig` = parameters of the frailty distribution.
#' * `$a_d` = parameters of the baseline hazard for terminal event.
#' * `$a_r` = parameters of the baseline hazard for recurrent events.
#' * `$beta_d` = covariate parameters in the terminal event model.
#' * `$beta_r` = covariate parameters in the recurrent events model.
#' @param var_names list with covariate variable names:
#' * `$terminal` = for terminal event model.
#' * `$recurrent` = for recurrent events model.
#' @param OBSL List with observation unit information. See [individual_loglikelihood()].
#' @param logLikelihood Value of the log-likelihood at optimum.
#' @param function_shapes list with function shapes for the model:
#' * `$hazard_d` = model for baseline hazard of the terminal event.
#' * `$hazard_r` = model for baseline hazard of recurrent events.
#' * `$surv_d` = model for baseline survival function of the terminal event.
#' * `$surv_r` = model for baseline survival function of the recurrent events.
#' * `$pdf` = model for the distribution of the frailty variable.
#' @param rec_timescale Model specification for the hazard of the recurrent event.
#' @param Calls list with calls with formulae:
#' * `$terminal` = for terminal event model.
#' * `$recurrent` = for recurrent events model.
#' @param optim_control list with information for optimization
#' * `$gradtol` = Tolerance for the gradient, see [stats::nlm()].
#' * `$steptol` = Tolerance for the steps, see [stats::nlm()].
#' * `$iterlim` = Iteration limit, see [stats::nlm()].
#' * `$iter` = Number of iterations until convergence.
#' * `$nlm_code` = Optimization exit status, see [stats::nlm()].
#' * `$nlm_hessian` = Hessian returned by [stats::nlm()].
#' * `$nlm_warnings` = Number of warnings during optimization. Mostly due to
#' expressions evaluating to infinity.
#' * `$int_mode` = Numerical integration mode. See [shared_frailty_fit()].
#' * `$MC_N` = Number of points for Monte Carlo integration.
#' * `$GC_nodes` = Number of points for Gaussian Quadrature integration.
#' * `$is_control` = Importance Sampling control parameters. See [shared_frailty_fit()].
#' * `$anal.grad` = If analytic gradient was used.
#' * `$BHHH.hessian` = If BHHH Hessian was used. See [shared_frailty_fit()].
#' * `$posit_cons` = Parameters with positivity constrains (transformations).
#' * `$par_scale` = Scale of each parameter.
#' * `$MC_samp` (Only during computarions) = Monte Carlo sample for numerical integration.
#'
#' @return SharedModel object.
#'
#' @export
SharedModel <- function(par, varmat, par_pos, var_names, OBSL, logLikelihood,
                              function_shapes, rec_timescale, Calls, optim_control){

  model<-list(par=par, varmat=varmat, par_pos=par_pos, var_names=var_names,
              OBSL=OBSL, logLikelihood=logLikelihood, funtion_shapes=function_shapes,
              rec_timescale=rec_timescale, Calls=Calls, optim_control=optim_control)

  class(model)<-"SharedModel" #Define the SharedModel class

  return(model)
}

#.............................................................................#
#.............................................................................#

##### Methods for SharedModel S3 class #####

#' Asymptotic covariance matrix
#'
#' Computes the asymptotic covariance matrix of the fitted parameters.
#'
#' @details
#' The asymptotic covariance matrix is the inverse of the Information Matrix.
#'
#' `BHHH = TRUE` uses the sum of squared scores to estimate the Information Matrix.
#' `BHHH = FALSE` uses the Hessian returned by [stats::nlm()] to estimate the
#' Information Matrix. Note that the `BHHH = FALSE` method requires small
#' scores to obtain a correct approximation (due to transformations that ensure
#' positivity of certain parameters).
#'
#'
#' @param obj [SharedModel()] object
#' @param BHHH Whether to approximate the Information Matrix by the sum of
#' squared scores. Values = `TRUE` (default) or `FALSE`.
#'
#' @return numeric matrix with asymptotic covariances
#' @export
asymptvar.SharedModel <- function(obj, BHHH = TRUE){

  # Get the Hessian of the log-likelihood
  if (BHHH){
    H <- BHHH_hessian(obj)

    # Invert
    V <- solve(H)

    # Hessian from optimization algorithm (approximation)
  } else {

    # Impose positivity constraint
    nu <- posit_constraint(obj$par, obj$optim_control$posit_cons, type = "reverse")

    # Inverse of the asymptotic variance of the transformed parameter
    H <- obj$optim_control$nlm_hessian

    # Invert
    nu_var <- solve(H)

    #gradient of the positivity transformation at nu
    dF<-posit_constraint(nu, obj$optim_control$posit_cons, type = "gradient")

    # Adjust for positivity constraints
    # This approximation is not perfect in the presence of positivity constraints
    # (a term is missing), but is perfect without them.
    V<-t(dF) %*% nu_var %*% dF
  }

  return(V)
}

#.............................................................................#
#.............................................................................#

#' Draw from the asymptotic distribution of the estimate
#'
#' Draw `n` points from the asymptotic distribution of the estimated parameter
#' of the model. It takes into account positivity constraints.
#'
#' @details
#' If `enforce_posit=TRUE`, the function ensures positivity of some parameters
#' (e.g., Weibull scale), even if this was not enforced during optimization.
#' This essentially draws from the asymptotic distribution of the parameters,
#' conditonal of these being positive.
#'
#'
#' @param obj [SharedModel()] object
#' @param n Sample size.
#' @param enforce_posit `TRUE` to enforce positivity (default).
#'
#' @return n times number of parameters matrix.
#'
#' @export
asympt_draw.SharedModel <- function(obj, n = 100, enforce_posit=TRUE){

  # Exit if sample size is zero
  if (n==0){ return(NULL)}

  # Impose positivity constraint used for optimization
  nu <- posit_constraint(obj$par, obj$optim_control$posit_cons, type = "reverse")

  # Asymp. var of the transformed parameter is the inverse of the Hessian
  H <-obj$optim_control$nlm_hessian

  # Invert
  V <- solve(H)

  # Draw from the normal distribution
  nus <- mvtnorm::rmvnorm(n, mean = nu, sigma = V)

  # Revert the transformation and return
  thetas <- t(apply(nus, 1, posit_constraint,
                  impose_ind = obj$optim_control$posit_cons, type = "direct"))

  # Ensure positivity if required
  if (enforce_posit){

    # Get parameters with a positivity constraint
    val_posit_cons <- valid_posit_constraint(obj$par_pos, obj$funtion_shapes)

    # Check which rows satisfy positivity for required parameters
    keep_r <- apply(thetas, 1, function(theta) all(theta[val_posit_cons]>0))

    # Keep valid rows
    thetas <- subset(thetas, subset=keep_r)

    # Re run the algorithm until the sample is completed
    thetas <- rbind(thetas,
                    asympt_draw(obj, n = n - nrow(thetas), enforce_posit=TRUE)
                    )
  }

  return(thetas)

}

#.............................................................................#
#.............................................................................#

#' Baseline hazard functions
#'
#' Get the baseline hazard functions at evaluated at `t`.
#'
#' @param obj [SharedModel()] object
#' @param t time to evaluate the hazard
#' @param process Process for which the hazard is computed. Either `terminal`
#' (`t`) or `recurrent` (`r`).
#' @param gradient `TRUE` to compute the gradient (default is `FALSE`).
#'
#' @return numeric value of the hazard at `t`. If `gradient=TRUE`, the returned
#' value has an attribute `"gradient"` with its derivative.
#'
#' @export
base_hazard.SharedModel <- function(obj, t, process, gradient = FALSE){

  #Get the parameters
  par_both<-param_hazard(obj)

  #Choose the baseline function and parameters in terms of the process
  if (process == "terminal" | process == "d"){

    basefun<-obj$funtion_shapes$hazard_d
    param<-par_both[['terminal']]

  } else if (process=="recurrent" | process == "r"){

    basefun<-obj$funtion_shapes$hazard_r
    param<-par_both[['recurrent']]

  } else {
    stop('Process must be "terminal" or "recurrent".')
  }

  #Evaluate the baseline function at t's
  h<-do.call(basefun, list(t=t, param=param))

  #If gradient=TRUE, also compute the gradient
  if (gradient){

    #Rename the function
    gradfun <- paste(basefun, "gradient", sep = "_")

    #Evaluate function
    grad<-do.call(gradfun, list(t=t, param=param))

    #We set the gradient as an attribute of the main result (hazard)
    attr(h, "gradient")<-grad
  }

  return(h)
}

#.............................................................................#
#.............................................................................#

#' Baseline survival functions
#'
#' Get the baseline survival functions at evaluated at `t`.
#'
#' @param obj [SharedModel()] object
#' @param t time to evaluate the survival function
#' @param process Process for which the survival function is computed. Either
#'  `terminal` (`t`) or `recurrent` (`r`).
#' @param gradient `TRUE` to compute the gradient (default is `FALSE`).
#'
#' @return numeric value of the survival function at `t`. If `gradient=TRUE`,
#' the returned value has an attribute `"gradient"` with its derivative.
#'
#' @export
base_survival.SharedModel <- function(obj, t, process, gradient = FALSE){

  #If the argument has no times, just return an empty vector
  if (length(t)==0){
    return(numeric(0))
  }

  #Get the parameters
  par_both<-param_hazard(obj)

  #Choose the baseline function and parameters in terms of the process
  if (process == "terminal" | process == "d"){

    basefun<-obj$funtion_shapes$surv_d
    param<-par_both[['terminal']]

  } else if (process=="recurrent" | process == "r"){

    basefun<-obj$funtion_shapes$surv_r
    param<-par_both[['recurrent']]

  } else {
    stop('Process must be "terminal" or "recurrent".')
  }

  #Evaluate the baseline function at t's
  S<-do.call(basefun, list(t=t, param=param))

  #If gradient=TRUE, also compute the gradient
  if (gradient){

    #Rename the function
    gradfun <- paste(basefun, "gradient", sep = "_")

    #Evaluate function
    grad<-do.call(gradfun, list(t=t, param=param))

    #We set the gradient as an attribute of the main result (hazard)
    attr(S, "gradient")<-grad
  }

  return(S)
}

#.............................................................................#
#.............................................................................#

#' Berndt–Hall–Hall–Hausman Hessian
#'
#' Estimates the Asymptotic Information matrix by aggregating the squared
#' scores.
#'
#' @param obj [SharedModel()] object
#' @param S matrix with scores (individuals in rows). Defaul is `NULL`, so that
#' the function computes the matrix.
#'
#' @return Estimate of the Asymptotic Information matrix
#'
#' @export
BHHH_hessian.SharedModel <- function(obj, S = NULL){

  # Check whether there is sample for MC integration (if it is needed)
  obj <- check_intsample(obj)

  if(is.null(S)){
    S <- scores(obj)
  }

  S_squared <- lapply(1:nrow(S), function(i) as.matrix(S[i,]) %*% t(S[i,]))

  H <- Reduce('+', S_squared)

  return(H)
}

#.............................................................................#
#.............................................................................#

#' Check whether there is Monte Carlo sample is stored. If not, it computes it.
#'
#' @param obj [SharedModel()] object
#'
#' @return [SharedModel()] object with a stored numerical integration sample.
#'
#' @export
check_intsample.SharedModel <- function(obj){

  # MC sample not present
  if (!('MC_samp' %in% names(obj$optim_control))){

    # Compute MC sample depending on the method
    if (obj$optim_control$int_mode == 'MC'){
      obj$optim_control$MC_samp <- dist.draw(n=obj$optim_control$MC_N,
                                             pdf=obj$funtion_shapes$pdf,
                                             param = obj$par[obj$par_pos$sig])

    } else if (obj$optim_control$int_mode == 'QMC'){
      seq <- vipor::vanDerCorput(n = obj$optim_control$MC_N, base = 2)
      obj$optim_control$MC_samp <- dist.transform(seq=seq,
                                                  pdf=obj$funtion_shapes$pdf,
                                                  param = obj$par[obj$par_pos$sig])

    } else if (obj$optim_control$int_mode == 'ISMC'){
      obj$optim_control$MC_samp <- do.call(paste('r',
                                                 obj$optim_control$is_control$dist,
                                                 sep = ''),
                                           c(list(obj$optim_control$MC_N),
                                             obj$optim_control$is_control$param)
      )
    } else if (obj$optim_control$int_mode == 'ISQMC'){
      seq <- vipor::vanDerCorput(n = obj$optim_control$MC_N, base = 2)
      obj$optim_control$MC_samp <- do.call(paste('q',
                                                 obj$optim_control$is_control$dist,
                                                 sep = ''),
                                           c(list(seq),
                                             obj$optim_control$is_control$param)
      )
    } else {
      obj$optim_control$MC_samp <- NULL
    }
  }

  return(obj)
}

#.............................................................................#
#.............................................................................#

#' Test relevance of recurrent event temporal distribution
#'
#' !!!!!! NEED TO DOCUMENT THIS.
#'
#' @param obj [SharedModel()] object
#' @param verbose
#' @param BHHH description
#'
#' @return
#' @export
#'
#' @examples
dist_relevance_test.SharedModel <- function(obj, verbose = FALSE, BHHH = TRUE){

  # Not possible if model is Poisson
  if (obj$rec_timescale == "Poisson") {
    warning("Test cannot be performed for Poisson models.\n",
            "The distribution of recurrent event times is irrelevant.")
    return(NULL)
  }

  #Get the frailty alpha
  alpha_pos <- obj$par_pos$alpha
  alpha <- obj$par[alpha_pos]

  #Variance (use BHHH variance)
  if (BHHH){
    VAR <- obj$varmat

    # Use optimization algorithm hessian
  } else {
    VAR <- asymptvar(obj, BHHH=FALSE)
  }

  #Wald test for alpha*(sigma-1)=0 in the Weibull case
  if (obj$funtion_shapes$hazard_r=="hazard_Weibull"){

    #Get the  Weibull shape
    shape_pos <- obj$par_pos$a_r[2]
    shape <-obj$par[shape_pos]

    # Get the covariance matrix
    V <- VAR[c(alpha_pos, shape_pos),c(alpha_pos, shape_pos)]

    # Transformation and gradient
    R <- unname(alpha * (shape-1))
    delta_R <- unname(as.matrix(c(shape-1, alpha)))

    # Test statistic
    DEN <- as.numeric(t(delta_R) %*% V %*% delta_R)
    W <- R^2 / DEN

    # p-value
    p <- 1 - stats::pchisq(W, df = 1)

    # Save result as wald.test object
    RES <- structure(list(Sigma = V,
                          b = c(alpha, shape),
                          Terms = NULL,
                          H0 = 0,
                          L = delta_R,
                          result = list(chi2 = c(chi2 = W, df = 1, P = p)),
                          verbose = verbose,
                          df = NULL),
                     class = "wald.test")

    #Default
  } else {
    test<-NULL
  }

  return(RES)
}

#.............................................................................#
#.............................................................................#

#' Number of recurrent events for each individual
#'
#' @param obj [SharedModel()] object
#'
#' @return numeric vector with the number of recurrent event for each individual.
#'
#' @export
get_rec_number.SharedModel <- function(obj){
  return(sapply(obj$OBSL, function(obs) sum(uncensored(obs$rec_times))))
}

#.............................................................................#
#.............................................................................#

#' Hazard given history
#'
#' Finds the hazard of the recurrent event at `t`, given the history of
#' recurrent events prior to `t`.
#'
#' @param obj [SharedModel()] object
#' @param t time to evaluate the hazard
#' @param Tjs numeric vector with recurrent event times prior to `t`.
#' @param gradient `TRUE` to compute the gradient (default is `FALSE`).
#'
#' @return numeric value of the hazard given history at `t`. If `gradient=TRUE`,
#'  the returned value has an attribute `"gradient"` with its derivative.
#'
#' @export
hazard_given_history.SharedModel<- function(obj, t, Tjs=NULL, gradient=FALSE){

  #If the argument has no times, just return an empty vector
  if (length(t)==0){
    return(numeric(0))
  }

  #Get the time scale of the recurrent process model
  timescale<-obj$rec_timescale

  # Append a zero at the begining of the Tjs
  # (in case t is smaller than the first Tj)
  Tjs<-c(0, Tjs)

  if (timescale=="Poisson"){

    #If process is Poisson, history does not matter
    #(we evaluate baseline hazard at t)
    t_eval<-t

  } else if (timescale=="renewal"){

    #For each t, find the index of first Tjs that is smaller than t
    ind_t<-sapply(t, function(ti) max(which(Tjs<ti)))

    #Find the gaps between t and the first Tjs smaller than it
    gaps<-t-Tjs[ind_t]

    #Baseline survival is evaluated at gaps
    t_eval<-gaps

  } else if (timescale=="piecewise-renewal"){

    # Poisson: piecewise-exponential

    # thresholds for pieces (add t=0)
    taus <- c(0, obj$funtion_shapes$rec_piecewise_ts)

    # Constants (add C0 = 1 as first constant, normalized)
    pw_C <- c(1, param_hazard(obj)$piecewise_Cs)

    #For each t, find the number of taus smaller than t (index for constant)
    ind_C<-sapply(t, function(ti) sum(taus<ti))

    # Constants for each t
    hazard_Consts <- pw_C[ind_C]

    # Renewal: user defined hazard

    #For each t, find the index of first Tjs that is smaller than t
    ind_t<-sapply(t, function(ti) max(which(Tjs<ti)))

    #Find the gaps between t and the first Tjs smaller than it
    gaps<-t-Tjs[ind_t]

    # Evaluate at gaps (renewal part)
    h_ren <- base_hazard(obj, t=gaps, process = 'recurrent', gradient = gradient)

    # Hazard = Poisson * renewal
    h <- hazard_Consts * h_ren

    # Computations for gradient
    if (gradient){

      # Gradient wrt renewal part parameters
      attr(h, 'gradient') <- hazard_Consts * attr(h_ren, 'gradient')

      # Gradient wrt Poisson part parameters (constants)
      grad_C <- sapply(seq(length(pw_C) - 1), function(col){
        (ind_C == col + 1) * h_ren
      })

      attr(h, 'gradient_C') <- matrix(grad_C, ncol = length(pw_C) - 1)
    }

    # Return piecewise-renewal hazard
    return(h)

  } else {
    stop('Incorrect definition of the recurrent event process.')
  }

  #Evaluate baseline hazard and ask for the gradient, if needed
  h<-base_hazard(obj, t=t_eval, process = 'recurrent', gradient = gradient)

  return(h)
}

#.............................................................................#
#.............................................................................#

#' Log-likelihood of the model
#'
#' Computes the log-likelihood of the model.
#'
#' @param obj [SharedModel()] object
#'
#' @return numeric value of the log-likelihood
#' @export
logLik.SharedModel <- function(obj){

  # Check whether there is sample for MC integration (if it is needed)
  obj <- check_intsample(obj)

  #Find the individual likelihoods
  Ls<-sapply(obj$OBSL, individual_loglikelihood, model=obj)

  #Return the sum
  return(sum(Ls))
}

#.............................................................................#
#.............................................................................#

#' Largest time observation in the sample
#'
#' @param obj [SharedModel()] object
#'
#' @return numeric value of the largest time observation in the sample
#' @export
max_timespan.SharedModel <- function(obj){
  Td <- sapply(obj$OBSL, function(obs) times(obs$ter_time))
  return(max(Td))
}

#.............................................................................#
#.............................................................................#

#' Optimization details
#'
#' Print optimization details in console.
#'
#' @param obj [SharedModel()] object
#'
#' @return Silent. Prints in console.
#'
#' @export
optim_details.SharedModel <- function(obj, digits = 1){

  # Optimization control
  cat('CONTROL PARAMETERS\n')
  cat('Using analytic gradient? ', obj$optim_control$anal.grad, '\n', sep = '')
  cat('Using BHHH Hessian? ', obj$optim_control$BHHH.hessian, '\n', sep = '')
  cat('Parameters with positivity constraint transformation = \n')
  print(obj$optim_control$posit_cons)
  cat('Parameter scaling during optimization = \n')
  print(obj$optim_control$par_scale, digits = digits)

  cat('\n')

  # Tolerances
  cat('TOLERANCES\n')
  cat('Step tolerance: ', obj$optim_control$steptol, '\n', sep = '')
  cat('Gradient tolerance: ', obj$optim_control$gradtol, '\n', sep = '')

  cat('\n')

  # Integration method
  cat('NUMERICAL INTEGRATION\n')

  if (obj$optim_control$int_mode == 'MC'){
    cat('Method: Monte Carlo\n')
    cat('Draws: ', obj$optim_control$MC_N, '\n', sep = '')

  } else if (obj$optim_control$int_mode == 'GQ') {
    cat('Method: Gauss-Laguerre Quadrature\n')
    cat('Nodes: ', obj$optim_control$GC_nodes, '\n', sep = '')

  } else if (obj$optim_control$int_mode == 'QMC') {
    cat('Method: Quasi Monte Carlo\n')
    cat('Number of van der Corput points: ', obj$optim_control$MC_N, '\n', sep = '')

  } else if (obj$optim_control$int_mode == 'ISMC') {
    cat('Method:  Importance Sampling Quasi Monte Carlo\n')
    cat('Draws: ', obj$optim_control$MC_N, '\n', sep = '')
    cat('Importance sampling distribution: ', obj$optim_control$is_control$dist,
        '\n', sep = '')
    cat('Importance sampling parameters: ',
        paste(sapply(obj$optim_control$is_control$param, round, digits = digits),
              collapse = ', '),
        '\n', sep = '')

  } else if (obj$optim_control$int_mode == 'ISQMC') {
    cat('Method:  Importance Sampling Quasi Monte Carlo\n')
    cat('Number of van der Corput points: ', obj$optim_control$MC_N, '\n', sep = '')
    cat('Importance sampling distribution: ', obj$optim_control$is_control$dist,
        '\n', sep = '')
    cat('Importance sampling parameters: ',
        paste(sapply(obj$optim_control$is_control$param, round, digits = digits),
              collapse = ', '),
        '\n', sep = '')
  }

  cat('\n')

  # Warnings
  # !!!!! PROBABLY NEED TO IMPROVE THIS
  cat('WARNINGS\n')
  cat('There were ', obj$optim_control$nlm_warnings, ' warnings throughout',
      ' the optimization process.\n',  sep='')

}

#.............................................................................#
#.............................................................................#

#' Coefficients of covariates
#'
#' @param obj [SharedModel()] object
#'
#' @return List with the coefficients:
#' * `$terminal` = for the terminal event process.
#' * `$recurrent` = for the recurrent event process.
#'
#' @export
param_coef.SharedModel <- function(obj){

  #Get coefficients for the terminal and recurrent events
  beta_d<-obj$par[obj$par_pos$beta_d]
  beta_r<-obj$par[obj$par_pos$beta_r]

  #Get names
  names_d<-obj$var_names$terminal
  names_r<-obj$var_names$recurrent

  #Set names
  names(beta_d)<-names_d
  names(beta_r)<-names_r

  return(list(terminal=beta_d, recurrent=beta_r))
}

#.............................................................................#
#.............................................................................#

#' Parameters relative to the frailty term
#'
#' @param obj [SharedModel()] object
#'
#' @return named numeric vector with the coefficients
#'
#' @export
param_frailty.SharedModel <- function(obj){

  #Get parameters for the terminal and recurrent events
  pars<-obj$par[c(obj$par_pos$alpha, obj$par_pos$sig)]
  names_pars<-c('alpha', paste('shape', seq_along(obj$par_pos$sig), sep = '_'))

  #Set names
  names(pars)<-names_pars

  return(pars)
}

#.............................................................................#
#.............................................................................#

#' Parameters of the baseline hazards
#'
#' @param obj [SharedModel()] object
#'
#' @return List with the parameters:
#' * `$terminal` = for the terminal event process.
#' * `$recurrent` = for the recurrent event process.
#'
#' @export
param_hazard.SharedModel <- function(obj){

  #Get parameters for the terminal and recurrent events
  a_d<-obj$par[obj$par_pos$a_d]
  a_r<-obj$par[obj$par_pos$a_r]
  pw_C<-obj$par[obj$par_pos$pw_C]

  #Get names depending on the method
  names_d<-do.call(paste(obj$funtion_shapes$hazard_d, 'names', sep = '_'), list())
  names_r<-do.call(paste(obj$funtion_shapes$hazard_r, 'names', sep = '_'), list())
  names_C<-NULL
  
  if(!is.null(obj$funtion_shapes$rec_piecewise_ts)){
    names_C<-paste('C', obj$funtion_shapes$rec_piecewise_ts, '-',
                   obj$funtion_shapes$rec_piecewise_ts[-1],
                   sep = '')
    names_C[length(names_C)] <- paste('C',
                                      obj$funtion_shapes$rec_piecewise_ts[length(names_C)],
                                      sep = '>')
  }

  #Set names
  names(a_d)<-names_d
  names(a_r)<-names_r
  names(pw_C)<-names_C

  return(list(terminal=a_d, recurrent=a_r, piecewise_Cs=pw_C))
}

#.............................................................................#
#.............................................................................#

#' Parameter tables for printing and latex output
#'
#' @param obj [SharedModel()] object with fit results.
#' @param cols character vector with values to display. A subset of
#' `c('est', 'est_SE', 'hr', 'hr_SE', 'pvals')`.
#' @param col_names Names of the columns.
#' @param BHHH `TRUE` to use Score based variance for computations. `FALSE` to
#' use the variance based on the optimization routine Hessian.
#'
#' @return List with tables.
#' * `$beta_d` = Covariate coefficients for terminal event hazard.
#' * `$beta_r` = Covariate coefficients for recurrent event hazard.
#' * `$a_d` = Parameters of baseline hazard of terminal event.
#' * `$a_r` = Parameters of baseline hazard of recurrent event.
#' * `$frail` = Parameters of frailty pdf.
#'
#' @export
param_tables.SharedModel <- function(obj,
                                     cols=c('est', 'est_SE', 'hr', 'hr_SE', 'pvals'),
                                     col_names=NULL,
                                     BHHH=TRUE){

  # Check if valid columns
  valid_cols <- c('est', 'est_SE', 'est_CI', 'hr', 'hr_SE', 'hr_CI', 'pvals')
  if (! all(cols %in% valid_cols)){
    stop(paste('Invalid columns. Select values from: ',
               paste(valid_cols, collapse = ' ')),
         call. = FALSE)
  }

  # Set default column names
  if (is.null(col_names)){
    col_names <- cols
  }

  # Check if length is valid
  if (length(col_names) != length(cols)){
    stop('Invalid lenght for col_names', call. = FALSE)
  }

  #Estimates
  parbeta<-param_coef(obj)
  parhaz<-param_hazard(obj)
  parfrail<-param_frailty(obj)

  #Standard errors (use BHHH variance)
  if (BHHH){
    SE<-sqrt(diag(obj$varmat))

  # Use optimization algorithm hessian
  } else {
    SE<-sqrt(diag(asymptvar(obj, BHHH=FALSE)))
  }

  # Confidence intervals at 95%
  alpha <- 0.95
  est_CIs <- sapply(seq_along(obj$par), function(p)
    conf_intervals(obj$par[p], SE[p], level = alpha, HR=FALSE)
  )

  #z-test p values
  pvals<-sapply(seq_along(obj$par), function(p)
    z_test(obj$par[p], SE[p], one.sided = obj$optim_control$posit_cons[p])$pval
  )

  # HR Confidence intervals at 95%
  alpha <- 0.95
  hr_CIs <- sapply(seq_along(obj$par), function(p)
    conf_intervals(obj$par[p], SE[p], level = alpha, HR=TRUE)
  )

  #Get parameter values in a list
  params<-vector(mode="list", length=0) # initialize

  # Add them if asked for
  if ('est' %in% cols){
    params$beta_d$est<-unname(parbeta$terminal)
    params$beta_r$est<-unname(parbeta$recurrent)
    params$a_d$est<-unname(parhaz$terminal)
    params$a_r$est<-unname(parhaz$recurrent)
    params$frail$est<-unname(parfrail)
    
    if (obj$rec_timescale %in% c('piecewise-renewal')){
      params$piecewise_Cs$est<-unname(parhaz$piecewise_Cs)
    }
  }

  if ('est_SE' %in% cols){
    params$beta_d$est_SE<-SE[obj$par_pos$beta_d]
    params$beta_r$est_SE<-SE[obj$par_pos$beta_r]
    params$a_d$est_SE<-SE[obj$par_pos$a_d]
    params$a_r$est_SE<-SE[obj$par_pos$a_r]
    params$frail$est_SE<-SE[c(obj$par_pos$alpha, obj$par_pos$sig)]
    
    if (obj$rec_timescale %in% c('piecewise-renewal')){
      params$piecewise_Cs$est_SE<-SE[obj$par_pos$pw_C]
    }
  }

  if ('est_CI' %in% cols){
    params$beta_d$est_CI<-est_CIs[obj$par_pos$beta_d]
    params$beta_r$est_CI<-est_CIs[obj$par_pos$beta_r]
    params$a_d$est_CI<-est_CIs[obj$par_pos$a_d]
    params$a_r$est_CI<-est_CIs[obj$par_pos$a_r]
    params$frail$est_CI<-est_CIs[c(obj$par_pos$alpha, obj$par_pos$sig)]
    
    if (obj$rec_timescale %in% c('piecewise-renewal')){
      params$piecewise_Cs$est_CI<-est_CIs[obj$par_pos$pw_C]
    }
  }

  if ('pvals' %in% cols){
    params$beta_d$pvals<-pvals[obj$par_pos$beta_d]
    params$beta_r$pvals<-pvals[obj$par_pos$beta_r]
    params$a_d$pvals<-pvals[obj$par_pos$a_d]
    params$a_r$pvals<-pvals[obj$par_pos$a_r]
    params$frail$pvals<-pvals[c(obj$par_pos$alpha, obj$par_pos$sig)]
    
    if (obj$rec_timescale %in% c('piecewise-renewal')){
      params$piecewise_Cs$pvals<-pvals[obj$par_pos$pw_C]
    }
  }

  if ('hr' %in% cols){
    params$beta_d$hr <- unname(exp(parbeta$terminal))
    params$beta_r$hr <- unname(exp(parbeta$recurrent))
    params$a_d$hr <- NA
    params$a_r$hr <- NA
    params$frail$hr <- NA
    params$piecewise_Cs$hr <- NA
  }

  if ('hr_SE' %in% cols){
    params$beta_d$hr_SE<- params$beta_d$hr * params$beta_d$est_SE
    params$beta_r$hr_SE <- params$beta_r$hr * params$beta_r$est_SE
    params$a_d$hr_SE <- NA
    params$a_r$hr_SE <- NA
    params$frail$hr_SE <- NA
    params$piecewise_Cs$hr_SE <- NA
  }

  if ('hr_CI' %in% cols){
    params$beta_d$hr_CI<-hr_CIs[obj$par_pos$beta_d]
    params$beta_r$hr_CI<-hr_CIs[obj$par_pos$beta_r]
    params$a_d$hr_CI<-NA
    params$a_r$hr_CI<-NA
    params$frail$hr_CI<-NA
    params$piecewise_Cs$hr_CI <- NA
  }

  # Add names
  params$beta_d$names<-names(parbeta$terminal)
  params$beta_r$names<-names(parbeta$recurrent)
  params$a_d$names<-names(parhaz$terminal)
  params$a_r$names<-names(parhaz$recurrent)
  params$frail$names<-names(parfrail)
  
  if (obj$rec_timescale %in% c('piecewise-renewal')){
    params$piecewise_Cs$names <- names(parhaz$piecewise_Cs)
  }

  #Create a table for each combination
  tables<-vector(mode="list", length = 0)

  i<-1
  for (par in params) {

    # create new table w/ the info
    newtab<-data.frame(par[cols])
    colnames(newtab) <- col_names
    rownames(newtab) <- par$names

    #Assign
    tables[[i]]<-newtab

    #increase i
    i<-i+1
  }

  #Name each table
  if(obj$rec_timescale %in% c('piecewise-renewal')){
    names(tables)<-c('beta_d', 'beta_r', 'a_d', 'a_r', 'frail', 'piecewise_Cs')
  } else {
    names(tables)<-c('beta_d', 'beta_r', 'a_d', 'a_r', 'frail')
  }

  return(tables)
}

#.............................................................................#
#.............................................................................#

# predict.SharedModel()
#
# This method is implemented in a separate script 'shared_frailty_pred.R'
# due to its relevance.
#
# Check the script for implementation details.

#.............................................................................#
#.............................................................................#

# predict_plot.SharedModel()
#
# This method is implemented in a separate script 'shared_frailty_pred.R'
# due to its relevance.
#
# Check the script for implementation details.

#.............................................................................#
#.............................................................................#

#' Print fit results
#'
#' @param obj [SharedModel()] object with fit results.
#' @param cols character vector with values to display. A subset of
#' `c('est', 'est_SE', 'est_CI', 'hr', 'hr_SE', 'hr_CI', 'pvals')`.
#' @param col_names Names of the columns.
#' @param digits Digits to round results.
#' @param BHHH `TRUE` to use Score based variance for computations. `FALSE` to
#' use the variance based on the optimization routine Hessian.
#'
#' @return Silent. Prints information in console.
#'
#' @export
print.SharedModel <- function(obj, cols=c('est', 'hr', 'pvals'),
                              col_names=NULL,
                              digits=3, BHHH=TRUE){

  # Default names
  default_col_names <- c(est = 'Estimate', hr = 'HR', pvals='P-value',
                       est_SE = 'Est. SE', hr_SE = 'HR SE',
                       est_CI = 'Est. 95% CI', hr_CI = 'HR 95% CI')
  if (is.null(col_names)){
    col_names <- default_col_names[cols]
  }

  #Get the tables w/ parameter estimate information
  Tabs<-param_tables(obj, cols=cols, col_names=col_names, BHHH=BHHH)

  # Remove columns with only NA's
  for (i in seq_along(Tabs)){
    Tabs[[i]] <- subset(Tabs[[i]], subset=rep(TRUE, nrow(Tabs[[i]])),
                        select=!apply(Tabs[[i]], 2, function(r) all(is.na(r))))
  }

  #Output on console
  #Title
  cat('Parametric shared frailty model for terminal and recurrent events. \n\n')

  #Cox coefficients
  cat('COEFFICIENTS\n')

  #Terminal
  cat('Estimates for the terminal event:\n')
  print(round(Tabs$beta_d, digits=digits))
  cat('\n')

  #Recurrent
  cat('Estimates for the recurrent event:\n')
  print(round(Tabs$beta_r, digits=digits))
  cat('\n\n')

  #Hazard estimates
  cat('HAZARDS\n')

  #Terminal
  cat('Hazard for the terminal event = ',
      stringr::str_remove(obj$funtion_shapes$hazard_d, 'hazard_'), '.\n', sep='')
  cat('Baseline Hazard Parameters:\n')
  print(round(Tabs$a_d, digits=digits))
  cat('\n')

  #Recurrent
  cat('Model for the recurrent event = ', obj$rec_timescale, '.\n', sep='')
  cat('Hazard for the recurrent event = ',
      stringr::str_remove(obj$funtion_shapes$hazard_r, 'hazard_'), '.\n', sep='')
  cat('Baseline Hazard Parameters:\n')
  print(round(Tabs$a_r, digits=digits))
  
  if(obj$rec_timescale == c('piecewise-renewal')){
    cat('Piece constants (Poisson part):\n')
    print(round(Tabs$piecewise_Cs, digits = digits))
  }
  cat('\n\n')

  #Frailty estimates
  cat('FRAILTY\n')
  cat('Frailty distribution = ', stringr::str_remove(obj$funtion_shapes$pdf, 'pdf_'),
      '.\n', sep='')
  cat('Estimates (alpha and distribution parameters):\n')
  print(round(Tabs$frail, digits = digits))
  cat('\n\n')

  #Likelihood
  cat('Log-likelihood = ', obj$logLikelihood, '.\n', sep='')

  # Optimization result
  if (obj$optim_control$nlm_code %in% c(1,2)){
    cat('Optimization concluded succesfully in ', obj$optim_control$iter, '/',
        obj$optim_control$iterlim, ' iterations.\n', sep='')

  } else if (obj$optim_control$nlm_code == 3){
    cat('Optimization: Last global step (', obj$optim_control$iter, '/',
        obj$optim_control$iterlim, ') failed to locate a point lower than estimate.\n',
        'Either estimate is an approximate local minimum of the function',
        ' or steptol is too small.\n',  sep='')

  } else if (obj$optim_control$nlm_code == 4){
    cat('Optimization: Iteration limit (', obj$optim_control$iterlim, ') exceeded.\n',
        sep = '')

  } else if (obj$optim_control$nlm_code == 5){
    cat('Optimization: maximum step (', obj$optim_control$iter, '/',
        obj$optim_control$iterlim, ') size stepmax exceeded five consecutive times.\n',
        'Either the function is unbounded below, becomes asymptotic to a finite value',
        ' from above in some direction or stepmax is too small.\n')
  }

  # More details about optimization
  cat('For more information regarding optimization, use optim_details(object).\n')
}

#.............................................................................#
#.............................................................................#

#' Boxplot with number of recurrent events
#'
#' @param obj [SharedModel()] object
#' @param ... Additional parameters to [graphics::boxplot()]. Also, specify
#' `horizontal=TRUE` for a horizontal plot.
#'
#' @return Silent. Creates plot.
#'
#' @examples
#'
#' @export
rec_boxplot.SharedModel <- function(obj, ...){

  # Data
  recnum <- get_rec_number(obj)
  mean_num <- mean(recnum)

  # Check if horizontal was specified to adjust where the mean is positioned
  vars <- list(...)
  if ('horizontal' %in% names(vars)){

    if (vars$horizontal){
      mean_dat <- list(x=mean_num, y=1)
      mean_ax <- 1

    } else {
      mean_dat <- list(x=1, y=mean_num)
      mean_ax <- 2
    }

    # Default is false
  } else {
    mean_dat <- list(x=1, y=mean_num)
    mean_ax <- 2
  }

  graphics::boxplot(recnum, ...)
  graphics::points(mean_dat, pch=17)
  graphics::axis(mean_ax, at = round(mean_num, 2))
}

#.............................................................................#
#.............................................................................#

#' Scores of the model
#'
#' @param obj [SharedModel()] object
#'
#' @return numeric matrix with individuals in rows and scores wrt each
#' parameter in columns.
#'
#' @export
scores.SharedModel <- function(obj){

  # Check whether there is sample for MC integration (if it is needed)
  obj <- check_intsample(obj)

  #Find the individual scores
  # sample size x num par matrix
  Scores<-t(sapply(obj$OBSL, individual_score, model=obj))

  return(Scores)
}

#.............................................................................#
#.............................................................................#

#' Descriptive summary of the sample
#'
#' @param obj [SharedModel()] object
#'
#' @return Silent. Prints information in console.
#'
#' @export
summary.SharedModel <- function(obj){

  ## Computations
  #Observation list
  OBSL<-obj$OBSL

  #Sample size
  N<-length(OBSL)

  # Observed time span
  max_T <- max_timespan(obj)

  #Terminal event information
  uncensored<-sapply(OBSL, function(obs) uncensored(obs$ter_time))
  n_uncensored<-sum(uncensored)

  #Recurrent event info
  Js<-get_rec_number(obj)
  Js_total<-sum(Js)
  Js_info<-summary(Js)

  #Terminal covariate info
  ter_covs<-t(sapply(OBSL, function(obs) obs$ter_covs[1,]))
  ter_covs_info<-apply(ter_covs, 2, mean)

  #Recurrent covariate info
  rec_covs<-t(sapply(OBSL, function(obs) obs$rec_covs[1,]))
  rec_covs_info<-apply(ter_covs, 2, mean)

  ## Output
  cat('Descriptive summary of the data used to fit the model.\n')
  cat('Use print method for fit results.\n\n')

  #Sample size and terminal events
  cat('SAMPLE SIZE AND FOLLOW-UP\n')
  cat('Sample size = ', N, '.\n', sep = '')

  if (length(obj$optim_control$NA_obsvar) >0){
    cat('Removed due to NA\'s (obsvar) = ',
        paste(obj$optim_control$NA_obsvar, collapse = ', '),
        '.\n', sep = '')
  }
  cat('Largest observed time = ', max_T, '.\n', sep = '')
  cat('\n')

  cat('TERMINAL EVENT\n')
  cat('Number of terminal events = ', n_uncensored, '.\n', sep = '')
  cat('Censored observations = ', N-n_uncensored, '.\n', sep = '')
  cat('\n')

  #Recurrent events
  cat('RECURRENT EVENT\n')
  cat('Total number of recurrent events = ', Js_total ,'.\n')
  cat('Distribution of recurrent events = \n')
  print(Js_info)
  cat('\n')

  #Covariates
  cat('COVARIATES\n')
  cat('Terminal event covariate mean = \n')
  print(ter_covs_info)
  cat('\n')
  cat('Recurrent event covariate mean = \n')
  print(rec_covs_info)
  cat('\n')
}

#.............................................................................#
#.............................................................................#

#' Survival function given history
#'
#' Finds the survival function of the recurrent event at `t`, given the history
#' of recurrent events prior to `t`.
#'
#' @param obj [SharedModel()] object
#' @param t time to evaluate the survival function
#' @param Tjs numeric vector with recurrent event times prior to `t`.
#' @param gradient `TRUE` to compute the gradient (default is `FALSE`).
#'
#' @return numeric value of the survival function given history at `t`. If
#' `gradient=TRUE`, the returned value has an attribute `"gradient"` with its
#' derivative.
#'
#' @export
surv_given_history.SharedModel<- function(obj, t, Tjs=NULL, gradient = FALSE){

  #Get the time scale of the recurrent process model
  timescale<-obj$rec_timescale

  # Append a zero at the begining of the Tjs
  # (in case t is smaller than the first Tj)
  Tjs<-c(0, Tjs)

  if (timescale=="Poisson"){

    #If process is Poisson, history does not matter
    S<-base_survival(obj, t=t, process = 'recurrent', gradient = gradient)

  } else if (timescale=="renewal"){

    #Find the indexes of the Tjs that are smaller than each t
    ind_t<-lapply(t, function(ti) which(Tjs<=ti))

    #Add the time where we predict the survival conditional on history
    T_list<-lapply(1:length(t), function(i) c(Tjs[ind_t[[i]]], t[i]))

    #Find the gaps between the provided times
    gaps_l<-lapply(T_list, function(ts) diff(ts))

    #Survival function at the gaps
    Sj_l<-lapply(gaps_l, base_survival, obj=obj, process = 'recurrent',
                 gradient = gradient)

    #Product of the Sjs
    S<-sapply(Sj_l, prod)

    # Compute gradient for recurrent
    if (gradient){

      # Get gradients for each gap
      Sj_grad_l <- lapply(Sj_l, attr, which="gradient")

      S_grad <- t(sapply(1:length(t), function(ti)
        renewal_grad(Sj_l[[ti]], Sj_grad_l[[ti]])))

      # Set attribute
      attr(S, "gradient") <- S_grad
    }

  } else if (timescale == 'piecewise-renewal'){

    #Find the indexes of the Tjs that are smaller than each t
    ind_t<-lapply(t, function(ti) which(Tjs<=ti))

    #Add the time where we predict the survival conditional on history
    T_list<-lapply(1:length(t), function(i) c(Tjs[ind_t[[i]]], t[i]))

    # Find conditional survival for each t
    Sj_l <- lapply(T_list, rec_piece_Surv, model=obj, gradient=gradient)

    #Product of the Sjs
    S<-sapply(Sj_l, prod)

    if (gradient){

      # Get gradients for each gap
      Sj_grad_l <- lapply(Sj_l, attr, which="gradient")
      Sj_grad_C_l <- lapply(Sj_l, attr, which="gradient_C")

      # Compute product gradient for each t
      S_grad <- t(sapply(1:length(t), function(ti)
        renewal_grad(Sj_l[[ti]], Sj_grad_l[[ti]])))
      S_grad_C <- t(sapply(1:length(t), function(ti)
        renewal_grad(Sj_l[[ti]], Sj_grad_C_l[[ti]])))

      # Set attributes
      attr(S, "gradient") <- S_grad
      attr(S, "gradient_C") <- S_grad_C
    }


  } else {
    stop('Incorrect definition of the recurrent event process.')
  }

  return(S)
}

#.............................................................................#
#.............................................................................#

#' Latex tables
#'
#' Creates tables with fit information and exports them in latex format.
#'
#' @param obj [SharedModel()] object
#' @param filename Path where tables are saved
#' @param cols character vector with values to display. A subset of
#' `c('est', 'est_SE', 'est_CI', 'hr', 'hr_SE', 'hr_CI', 'pvals')`.
#' @param col_names Names of the columns.
#' @param merged `TRUE` (default) to build a single table. `FALSE` creates
#' separated tables for fit results.
#' @param digits Digits to round results.
#' @param horizontal `TRUE` (default) to order the results regarging the terminal
#' and recurrent events horizontally.
#' @param append `TRUE` to append results to `filename`. Default is `FALSE`.
#' @param BHHH  `TRUE` (default) to approximate the Information Matrix by the
#'  sum of squared scores. See [asymptvar.SharedModel()].
#'
#' @return Silent. Creates files with latex tables.
#'
#' @examples
#'
#' @export
toLatex.SharedModel <- function(obj, filename='tab.tex',
                                cols=c('est', 'hr', 'pvals'),
                                col_names=NULL,
                                merged=TRUE,
                                digits = 3,
                                horizontal=TRUE, append=FALSE, BHHH=TRUE){

  # Default names
  default_col_names <- c(est = 'Estimate', hr = 'HR', pvals='P-value',
                       est_SE = 'Est. SE', hr_SE = 'HR SE',
                       est_CI = 'Est. 95% CI', hr_CI = 'HR 95% CI')
  if (is.null(col_names)){
    col_names <- default_col_names[cols]
  }

  #Get the tables
  Tabs<-param_tables(obj, cols=cols, col_names=col_names, BHHH=BHHH)

  # Merge in a single table
  if (merged){

    #Get the latex string for each
    TexTabs<-lapply(Tabs, function(tab)
                            utils::toLatex(xtable::xtable(tab,
                                                          digits = digits)))

    #Get event names from the call formulas
    namesD<-get_time_event_variables(obj$Calls$terminal)
    namesR<-get_time_event_variables(obj$Calls$recurrent)

    # Get the table: horizontal format or vertical format
    if (horizontal){
      textable <- textable_horizontal(Tabs, TexTabs, namesD, namesR)
    } else {
      textable <- textable_vertical(Tabs, TexTabs, namesD, namesR)
    }

    #Delete if file already exists
    if ((!append) & file.exists(filename)) file.remove(filename)

    #Export
    writeLines(textable, con = filename)
    cat(paste('Successfully created', filename, '\n'))

  #If not, just export the tables
  } else {

    # Get table names
    tabnames<-names(Tabs)

    #Print the tables
    for (i in seq_along(Tabs)){

      # Remove columns with only NA's
      Tabs[[i]] <- subset(Tabs[[i]], subset=rep(TRUE, nrow(Tabs[[i]])),
                          select=!apply(Tabs[[i]], 2, function(r) all(is.na(r))))

      #Complete file name
      outname<-paste(dirname(filename), '/',
                     tools::file_path_sans_ext(basename(filename)), '_',
                     tabnames[i],
                     '.', tools::file_ext(filename),
                     sep = '')

      #Delete if file already exists
      if ((!append) & file.exists(outname)) file.remove(outname)

      #Write table
      writeLines(utils::toLatex(xtable::xtable(Tabs[[i]])), con = outname)
      cat(paste('Successfully created', outname, '\n'))
    }
  }

}

#.............................................................................#
#.............................................................................#

#' Wald tests for the parameters of the model
#'
#' This function is a wrapper around [aod::wald.test()].
#'
#' @param obj obj [SharedModel()] object
#' @param termpos numeric vector with positions of the terms to test. Use
#' `obj$par_pos`. See [SharedModel()].
#' @param par0 Null hypothesis to test.
#' @param BHHH `TRUE` to use Score based variance for computations. `FALSE` to
#' use the variance based on the optimization routine Hessian.
#'
#' @return An object of class [wald.test]. See [aod::wald.test()].
#'
#' @examples
#'
#' @export
Wald_test.SharedModel <- function(obj, termpos, par0 = NULL, BHHH = TRUE){

  #Set to default, which is a vector of zeros
  if (is.null(par0)) par0<-rep(0, times=length(termpos))

  #Variance (use BHHH variance)
  if (BHHH){
    VAR <- obj$varmat

  # Use optimization algorithm hessian
  } else {
    VAR <- asymptvar(obj, BHHH=FALSE)
  }

  #Call the aoe package function
  test<-aod::wald.test(b=obj$par, Sigma = VAR, Terms = termpos,
                       H0 = par0)

  return(test)
}

#.............................................................................#
#.............................................................................#
