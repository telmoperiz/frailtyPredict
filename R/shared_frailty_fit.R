#' Fit Shared Frailty Model
#'
#' `shared_frailty_fit()` fits a shared frailty model for terminal and recurrent
#'  events.
#'
#' @details
#' Model for the hazard of the terminal event (d) and recurrent event (r)
#' processes:
#' \deqn{  \lambda_i^d(t) = u_i^\alpha \cdot
#' \exp(\beta_d'Z_i) \cdot \lambda_0^d(t) }
#' \deqn{ \lambda_i^r(t) = u_i \cdot
#' \exp(\beta_r'Z_i) \cdot \lambda_0^r(t|N_i^r(t)),}
#' where \eqn{u_i} is the frailty of unit \eqn{i}, \eqn{Z_i} are covariates,
#' \eqn{h_0^d(t)} is the baseline hazard for the terminal event, and
#' \eqn{\lambda_0^r(t|N_i^r(t))} is the hazard for the recurrent event given the
#' history of recurrent events \eqn{N_i(t)}.
#'
#' If `rec_timescale = "Poisson"`, the hazard for recurrent events is
#' \eqn{\lambda_0^r(t|N_i^r(t))=h_0^r(t)}, where \eqn{h_0^r(t)} is the baseline
#' hazard for recurrent events. If `rec_timescale = "renewal"`, the hazard for
#' recurrent events is \eqn{\lambda_0^r(t|N_i^r(t))=h_0^r(t - T^r_{i, N_i(t-)})}.
#' In the preceeding equation, \eqn{T^r_{i, N_i(t-)}} represents the time of the
#' last recurrent event before time \eqn{t}, so the renewal time scale fits the
#' model in gap times.
#'
#' When `ter_hazard = "Weibull"` or `rec_hazard = "Weibull"`, the corresponding
#' baseline hazards are
#' \deqn{h_0(t) = b t^{b-1}/a^b,}
#' where \eqn{a} is a scale parameter and \eqn{b} is a shape parameter.
#'
#' For `frailty = "gamma"` the distribution of \eqn{u_i} is
#' \deqn{\frac{u^{1/\sigma-1}\exp(-u/\sigma)}{\sigma^{1/\sigma}\Gamma(1/\sigma)},}
#' where \eqn{\sigma = Var(u_i)}. Also, \eqn{E[u_i]=1}.
#'
#' Maximum likelihood estimation is performed using [stats::nlm()].
#'
#'
#' @param data data.frame with sample, in the following format. A ROW is a time
#' observations. In COLUMNS one must provide:
#' * Time of the observation.
#' * Observation unit
#' * Indicator of terminal event happening
#' * Indicator of recurrent event happening
#' * Covariate values
#' @param terminal_formula Formula for the terminal event.
#' @param recurrent_formula Formula for the recurrent event.
#' @param obsvar Varible indicating observation units (individuals)
#' @param ter_hazard Parametric specification of the baseline hazard for the
#' terminal event. Supports: `"Weibull"`.
#' @param rec_hazard Parametric specification of the baseline hazard for the
#' recurrent event. Supports: `"Weibull"`.
#' @param rec_timescale Model specification for the hazard of the recurrent event
#' Supports: `"Poisson"` and `"renewal"`.
#' @param frailty Distribution of the frailty variable. Supports: `"gamma"`.
#' @param int_mode How to compute the numerical integrals:
#' * `"GQ"` = Gaussian Quadrature.
#' * `"MC"` = Monte Carlo, with R sampler.
#' * `"QMC"` = Quasi Monte Carlo, with van der Corput sequence.
#' * `"ISMC"` = Importance Sampling MC.
#' * `"ISQMC"` = Importance Sampling QMC (the default).
#' @param MC_N  Number of points for Monte Carlo integration.
#' @param GC_nodes Number of points for Gaussian Quadrature integration.
#' @param is_control list with Importance Sampling control parameters.
#' * `$dist` = distribution to sample from. Any [stats::Distributions()] is valid.
#' Default is `"gamma"` which samples from [stats::dgamma()].
#' * `$param` = list with parameters to pass to the R sampler. Default is `list(1,1)`.
#' @param anal.grad Whether to use the analytic gradient = `TRUE` (default) or `FALSE`.
#' @param BHHH.hessian Wheter to use the Berndt-Hall-Hall-Haussman algorithm, which
#' approximates the Hessian of the loglikelihood using scores. It can be `TRUE`,
#' use BHHH Hessian in optimization, or `FALSE` (default), numerically compute the
#' Hessian.
#' @param param_scale Scale of the parameters. It may speed up convergence.
#' Values = `"auto"` (default), for automatic scalling, `"same"`, for equal scalling,
#' or `numeric(parameter length)`, for customized scalling.
#' @param positivity_cons Whether to impose a positivity constraint when maximizing.
#' Values = `"auto"` (default), for automatic selection of the constrained parameters,
#' `"none"`, for optimization without constraints, or `logical(parameter length)`,
#' for customized positivity constraints.
#' @param gradtol Tolerance for the gradient, see [stats::nlm()].
#' @param steptol Tolerance for the steps, see [stats::nlm()].
#' @param iterlim Iteration limit, see [stats::nlm()].
#' @param print_level Constrols console output. Values = `0`, `1`, or `2`;
#' where `2` indicates more verbose. See [stats::nlm()].
#'
#' @return [SharedModel()] object with fit results.
#'
#' @examples
#' @export
shared_frailty_fit <- function(data, terminal_formula, recurrent_formula, obsvar,
                              ter_hazard, rec_hazard, rec_timescale, frailty,
                              rec_piecewise_ts = NULL,
                              int_mode = "ISQMC", MC_N = 1e3, GC_nodes = 32,
                              is_control = list(dist='gamma', param=list(1, 1)),
                              anal.grad = TRUE, BHHH.hessian = FALSE,
                              param_scale = 'auto', positivity_cons = 'auto',
                              gradtol = 1e-6, steptol = 1e-6,
                              iterlim = 200, print_level = 0){

  ##### Checking correct inputs #####

  # Supported frailty distributions
  frailty_supp <- c('gamma')
  if (!(frailty %in% frailty_supp)){
    stop(paste('Invalid frailty. Supported:', paste(frailty_supp, collapse = ', ')))
  }

  # Supported hazard models
  hazard_supp <- c('Weibull')
  if (!( (ter_hazard %in% hazard_supp)  & (rec_hazard %in% hazard_supp) )){
    stop(paste('Invalid hazards. Supported:', paste(hazard_supp, collapse = ', ')))
  }

  # Supported timescales
  timescale_supp <- c('Poisson','renewal', 'piecewise-renewal')
  if (!(rec_timescale %in% timescale_supp)){
    stop(paste('Invalid timescale. Supported:', paste(timescale_supp, collapse = ', ')))
  }
  
  # Check thresholds for piecewise-renewal
  if(rec_timescale == 'piecewise-renewal' && is.null(rec_piecewise_ts)){
    stop('Must provide rec_piecewise_ts (thresholds)')
  }


  # Integral mode
  if (!int_mode %in% c('GQ', 'MC', 'QMC', 'ISMC', 'ISQMC')){
    stop(paste('int_mode must be  "GQ" (Gaussian Quadrature),',
         '"MC" (Monte Carlo), "QMC" (Quasi Monte Carlo),',
         '"ISMC" (Importance Sampling MC), or "ISQMC" (Importance Sampling QMC)'))
  }

  # Parameter scale 'auto', 'same' or numeric
  if (!is.numeric(param_scale)) {
    if (!param_scale %in% c('auto', 'same')){
      stop(paste('param_scale must be "auto", "same", or a numeric vector'))
    }
  }

  # Positivity constraint 'auto', 'none' or boolean
  if (!is.logical(positivity_cons)) {
    if (!positivity_cons %in% c('auto', 'none')){
      stop(paste('positivity_cons must be "auto", "none", or a logical vector'))
    }
  }

  # BHHH only valid with analytic gradients
  if (BHHH.hessian & !(anal.grad)){
    stop('BHHH Hessian computation only available with analytic gradients')
  }

  # Check Importance Sampling distribution
  if (int_mode %in% c('ISMC', 'ISQMC')){

    # Try to draw one point
    rpoint <- tryCatch(do.call(paste('r', is_control$dist, sep = ''),
                               c(list(1), is_control$param)),
                       error = function(cond) NULL,
                       warning = function(cond) NaN
              )

    # Try to get median
    median <- tryCatch(do.call(paste('q', is_control$dist, sep = ''),
                               c(list(0.5), is_control$param)),
                       error = function(cond) NULL,
                       warning = function(cond) NaN
              )

    # Check if it went wrong
    if (is.null(rpoint) | is.null(median)){
      stop('Invalid Importance Sampling distribution name.')
    }

    # Check correct parameters
    if ( is.nan(rpoint) | is.nan(median) | is.infinite(median)){
      stop('Invalid Importance Sampling parameter values.')
    }
  }

  ##### Specification of functional forms #####

  #List w/ functional forms
  model_funs<-list(hazard_d=paste("hazard", ter_hazard, sep = "_"),
                   hazard_r=paste("hazard", rec_hazard, sep = "_"),
                   surv_d=paste("surv", ter_hazard, sep = "_"),
                   surv_r=paste("surv", rec_hazard, sep = "_"),
                   pdf=paste("pdf", frailty, sep = "_"),
                   rec_piecewise_ts=rec_piecewise_ts
                   )

  #Functions giving default parameters
  hd_defaults<-paste("hazard", ter_hazard, "defaults", sep = "_")
  hr_defaults<-paste("hazard", rec_hazard, "defaults", sep = "_")
  pdf_defaults<-paste("pdf", frailty, "defaults", sep = "_")

  ##### Data work #####

  #Build models for terminal (D) and recurrent events (R)
  MD<-stats::model.frame(terminal_formula, data = data, na.action = NULL)
  MR<-stats::model.frame(recurrent_formula, data = data, na.action = NULL)

  #Extracting the data

  #Get the terms
  TermsD<-attr(MD, "terms")
  TermsR<-attr(MR, "terms")

  #We add a fake intercept. Why?
  #model.matrix adds all factor levels if intercept=0. we want to have a control
  #level to account for the baseline hazards
  attr(TermsD, "intercept")<-1
  attr(TermsR, "intercept")<-1

  #Construct the data matrices
  TD<-stats::model.extract(MD, 'response') #termianl event time
  ZD<-stats::model.matrix(TermsD, MD) #Covariate matrix
  TR<-stats::model.extract(MR, 'response') #Recurrent event times
  ZR<-stats::model.matrix(TermsR, MR) #Covariate matrix

  #Remove the fake intercept
  ZD<-ZD[,-1]
  ZR<-ZR[,-1]
  attr(TermsD, "intercept")<-0
  attr(TermsR, "intercept")<-0

  ##### Create list w/ times and variables for each observation #####

  #Indicator for each observation
  OBSIND <- stats::model.frame(obsvar, data=data)
  OBSIND <- OBSIND[,1]

  #Number of different observations
  OBS<-unique(OBSIND)

  #Creating list w/ information for each observation
  # For the terminal event we just take the last row
  # (the others correspond to recurrent events)
  #   $ter_time = terminal event time (an whether it is censored)
  #   $ter_covs = 1 x NCOVS vector w/ covariates
  #   $rec_times = recurrent event times (w/ censiring/termianl time at last)
  #   $rec_covs = (J+1) x NCOVS vector w/ covariates
  OBSL<-lapply(OBS, function(i) list(ter_time=utils::tail(subset(TD, OBSIND==i), n=1),
                                     ter_covs=utils::tail(subset(ZD, OBSIND==i), n=1),
                                     rec_times=subset(TR, OBSIND==i),
                                     rec_covs=subset(ZR, OBSIND==i)))

  # Remove NA's
  remove <- c()
  for (j in seq_along(OBSL)){

    obs <- OBSL[[j]]
    if ( any(is.na(obs$ter_time), is.na(obs$rec_times), is.na(obs$ter_covs),
             is.na(obs$rec_covs)) ){
        remove <- c(remove, j)
    }
  }

  if (length(remove) > 0){
    OBSL <- OBSL[-remove]
    removed <- OBS[remove]
  } else {
    removed <- c()
  }

  ##### Define the default parameters #####

  #Terminal event
  defs_d<-do.call(hd_defaults, list(data=data[!OBSIND %in% removed,],
                                    terminal_formula=terminal_formula,
                                    recurrent_formula=recurrent_formula, OBSL=OBSL,
                                    rec_timescale=rec_timescale, eventtype="terminal"))

  beta_d<-defs_d$beta
  a_d<- defs_d$a

  #Recurrent event
  defs_r<-do.call(hr_defaults, list(data=data[!OBSIND %in% removed,],
                                    terminal_formula=terminal_formula,
                                    recurrent_formula=recurrent_formula, OBSL=OBSL,
                                    rec_timescale=rec_timescale, eventtype="recurrent"))

  beta_r<-defs_r$beta
  a_r<- defs_r$a
  
  # Piecewise Constants C1, C2, ...
  if(rec_timescale == 'piecewise-renewal'){
    pw_C <- rep(1, length(rec_piecewise_ts))
  } else {
    pw_C <- c()
  }

  #Frailty term
  alpha <- 1
  sig <- do.call(pdf_defaults, list())

  #Initial guess
  theta_ini<-c(alpha, sig, a_d, a_r, beta_d, beta_r, pw_C)

  #Ordering of the parameters
  par_pos<-vector(mode='list', length = 0) #initialize
  par_pos$alpha<-1
  par_pos$sig<- seq(par_pos$alpha+1, length.out=length(sig))
  par_pos$a_d<-seq(utils::tail(par_pos$sig, n=1)+1, length.out=length(a_d))
  par_pos$a_r<-seq(utils::tail(par_pos$a_d, n=1)+1, length.out=length(a_r))
  par_pos$beta_d<-seq(utils::tail(par_pos$a_r, n=1)+1, length.out=length(beta_d))
  par_pos$beta_r<-seq(utils::tail(par_pos$beta_d, n=1)+1, length.out=length(beta_r))
  par_pos$pw_C <- seq(utils::tail(par_pos$beta_r, n=1)+1, length.out=length(pw_C))

  ##### Positivity constraint and parameter scale #####

  #Positivity constraint
  if (all(positivity_cons == 'auto')){

    # Get specified positivity constraints
    posit_cons <- valid_posit_constraint(par_pos, model_funs)

  } else if (all(positivity_cons == 'none')){
    posit_cons<-rep(FALSE, times=length(theta_ini))

  # User provided
  } else {
    if (length(positivity_cons) != length(theta_ini)) {
      stop('Provided positivity_cons must have the same length as number of parameters')
    }

    posit_cons <- positivity_cons
  }

  #Impose the positivity constraint on the initial values
  nu_ini <- posit_constraint(theta_ini, posit_cons, type = "reverse")

  # Parameter scale
  if (all(param_scale == 'auto')){
    par_scale <- abs(nu_ini) # Scale in absolute value

  } else if (all(param_scale == 'same')){
    par_scale <- rep(1, times=length(theta_ini))

    # User provided
  } else {
    if (length(param_scale) != length(theta_ini)) {
      stop('Provided param_scale must have the same length as number of parameters')
    }

    par_scale <- posit_constraint(param_scale, posit_cons, type = 'reverse')
  }

  ##### Sequence for (Quasi) Monte Carlo integration #####

  # Generate sequence for QMC integration
  if (int_mode %in% c('QMC', 'ISQMC')){

    int_seq <- vipor::vanDerCorput(MC_N, base = 2) # van der Corput sequence

    # Draw from IS distribution
    if (int_mode == 'ISQMC') {
      int_seq <- do.call(paste('q', is_control$dist, sep = ''),
                         c(list(int_seq), is_control$param))
    }

  # Importance sampling for MC
  } else if (int_mode == 'ISMC') {
    int_seq <- do.call(paste('r', is_control$dist, sep = ''),
                       c(list(MC_N), is_control$param))
  } else {
    int_seq <- NULL
  }

  ##### Call the optimization routine #####

  # clear warnings to count the optimization ones
  # !!!!!!!! THIS DOESNT WORK
  assign("last.warning", NULL)

  # Message regarding parameters, likelihood and optimization
  if (print_level > 0){
    cat('Printing optimization details provided by nlm.\n')
    cat('------------------------------------------------------------------------------\n')
    cat('Note that Parameter values may change according to the positivity constraint.\n')
    cat('Minimizing f(theta) = - loglikelihood(data, theta): Funtion Value is positive.\n')
    cat('------------------------------------------------------------------------------\n')
  }

  # Optimization
  res <- stats::nlm(loglikelihood, nu_ini,
                   OBSL = OBSL, par_pos = par_pos, posit_cons = posit_cons,
                   model_funs = model_funs, rec_timescale = rec_timescale,
                   int_mode = int_mode, MC_N = MC_N, GC_nodes = GC_nodes,
                   int_seq = int_seq, is_control = is_control,
                   anal.grad = anal.grad, BHHH.hessian = BHHH.hessian,
                   hessian = TRUE, typsize = par_scale, gradtol=gradtol,
                   steptol = steptol, iterlim = iterlim, print.level=print_level,
                   check.analyticals = FALSE)

  #Number of warnings (usually Infinity (0 inside the log))
  N_warn<-length(last.warning)

  #Recover the estimate from the transformation
  theta<-posit_constraint(res$estimate, posit_cons, type = "direct")

  ##### Build the output (SharedModel) #####
  out<-SharedModel(par=theta, varmat = NULL, par_pos=par_pos,
                    var_names=list(terminal=colnames(OBSL[[1]]$ter_covs),
                                   recurrent=colnames(OBSL[[1]]$rec_covs)),
                    OBSL=OBSL, logLikelihood = -res$minimum,
                    function_shapes=model_funs, rec_timescale=rec_timescale,
                    Calls=list(terminal=terminal_formula, recurrent=recurrent_formula),
                    optim_control=list(gradtol=gradtol, steptol = steptol, iterlim = iterlim,
                                       iter=res$iterations, nlm_code=res$code,
                                       nlm_hessian = res$hessian,
                                       nlm_warnings=N_warn, int_mode = int_mode,
                                       MC_N = MC_N, GC_nodes=GC_nodes,
                                       is_control = is_control,
                                       anal.grad = anal.grad,
                                       BHHH.hessian = BHHH.hessian,
                                       posit_cons = posit_cons,
                                       par_scale = par_scale,
                                       NA_obsvar = removed)
                    )

  # Compute asymptotic variance
  out$varmat <- asymptvar(out)

  return(out)
}
