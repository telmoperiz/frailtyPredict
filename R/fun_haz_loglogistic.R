###### Baseline hazards: Log-logistic ######

# hazard
hazard_loglogistic <- function(t, param){

  # Check time larger than 0
  if (any(t <= 0)){
    stop('Times must be larger than 0')
  }

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  #Log-logistic hazard function:
  # h(t)= (shape/scale) * (t/scale)^(shape-1) / (1 + (t/scale)^shape)
  h <- unname((shap/scal) * (t/scal)^(shap-1) / (1 + (t/scal)^shap))

  return(h)
}

#.............................................................................#
#.............................................................................#

# defaults
hazard_loglogistic_defaults<- function(data, terminal_formula, recurrent_formula,
                                       OBSL,  rec_timescale, eventtype){

  #Initialize defaults
  defs<-vector(mode="list", length = 0)

  if (eventtype == "terminal"){

    #Run parametric log-logistic model for the terminal event
    #We update the formula to always include an intercept (scale parameter)
    mod<-survival::survreg(formula = update.formula(terminal_formula, ~ . +1 ),
                           data = data, dist = "loglogistic")

  } else if (eventtype == "recurrent") {

    #Parametric survival regression for defaults (in gap or calendar times)
    #Update data to create variable t.mod (survival object in the correct timescale)
    data_mod<-data
    data_mod[,'t.mod']<-modify_timescale(OBSL, rec_timescale)

    #Update the formula to have t.mod as dependent variable (and intercept)
    mod<-survival::survreg(formula = update.formula(recurrent_formula, t.mod ~ . +1),
                           data = data_mod, dist = "loglogistic")
  }

  #Get the defaults
  #Cox coefficients
  # Proper scaling as Kalbfleisch and Prentice (p. 37)
  defs$beta <- - unname(mod$coefficients[-1]) / mod$scale

  #Weibull scale and shape
  # Note that K&P's lambda is the inverse of our scale
  defs$a <- unname(c( exp(mod$coefficients[1]), (mod$scale)^(-1) ))

  return(defs)
}

#.............................................................................#
#.............................................................................#

# gradient of hazard
hazard_loglogistic_gradient <- function(t, param){

  # Check time larger than 0
  if (any(t <= 0)){
    stop('Times must be larger than 0')
  }

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  # Denominator in derivative
  den <-scal^shap + t^shap

  #Derivative wrt scale
  dh_scal<-unname(
    - shap^2 * (t * scal)^(shap - 1) / den^2
  )

  #Derivative wrt shape
  dh_shap<- unname(
    t^(shap - 1) / den +
      shap * scal * (t * scal)^(shap - 1) * (log(t) - log(scal)) / den^2
  )

  #Return a length(t) x num_param (2) matrix with the gradient at each t
  return(matrix(c(dh_scal, dh_shap), ncol = length(param)))
}

#.............................................................................#
#.............................................................................#

# Names of the parameters
hazard_loglogistic_names<- function(){
  return(c('scale','shape'))
}

#.............................................................................#
#.............................................................................#

# Positivity constraint on the parameters
hazard_loglogistic_posit<- function(){
  return(c(TRUE, TRUE))
}

#.............................................................................#
#.............................................................................#

# survival
surv_loglogistic <- function(t, param){

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  #Loglogistic survival function: S(t)=1 / (1 + (t/scale)^shape)
  S <-unname(
    1 / (1 + (t/scal)^shap)
  )

  return(S)
}

#.............................................................................#
#.............................................................................#

# gradient of survival
surv_loglogistic_gradient <- function(t, param){

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  # Survival function
  surv <- surv_loglogistic(t, param)

  # Scaled parameter
  t_scal <- t / scal

  #Derivative wrt scale
  dh_scal<-unname(
    surv^2 * shap * t_scal^shap / scal
  )

  #Derivative wrt shape
  dh_shap<- unname( - (t_scal)^shap * log(t/scal) * surv^2 )

  #Return a length(t) x num_param (2) matrix with the gradient at each t
  return(matrix(c(dh_scal, dh_shap), ncol = length(param)))
}

#.............................................................................#
#.............................................................................#
