###### Baseline hazards: Gompertz ######

# hazard
hazard_Gompertz <- function(t, param){

  # Check time larger than 0
  if (any(t <= 0)){
    stop('Times must be larger than 0')
  }

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  #Log-logistic hazard function:
  # h(t)= (shape/scale) * exp(t/scale)
  h <- unname((shap/scal) * exp(t / scal))

  return(h)
}

#.............................................................................#
#.............................................................................#

# defaults
hazard_Gompertz_defaults<- function(data, terminal_formula, recurrent_formula,
                                       OBSL,  rec_timescale, eventtype){

  #Initialize defaults
  defs<-vector(mode="list", length = 0)

  if (eventtype == "terminal"){

    #Run parametric log-logistic model for the terminal event
    #We update the formula to always include an intercept (scale parameter)
    mod<-eha::phreg(formula = update.formula(terminal_formula, ~ . +1 ),
                    data = data, dist = 'gompertz', param='rate')

  } else if (eventtype == "recurrent") {

    #Parametric survival regression for defaults (in gap or calendar times)
    #Update data to create variable t.mod (survival object in the correct timescale)
    data_mod<-data
    data_mod[,'t.mod']<-modify_timescale(OBSL, rec_timescale)

    #Update the formula to have t.mod as dependent variable (and intercept)
    mod<-eha::phreg(formula = update.formula(recurrent_formula, t.mod ~ . +1),
                    data = data_mod, dist = 'gompertz', param='rate')
  }

  #Get the defaults
  len_fit <- length(mod$coefficients)

  #Cox coefficients (hazard parameters are at end)
  defs$beta <- unname(mod$coefficients[-c(len_fit - 1, len_fit)])

  # Scale and shape
  # The fitted hazard is h(t) = level * exp(rate * t) and log(level) is reported
  scal <- 1 / mod$coefficients['rate']
  defs$a <- unname(c(
    scal,
    scal * exp(mod$coefficients['log(level)'])
  ))

  return(defs)
}

#.............................................................................#
#.............................................................................#

# gradient of hazard
hazard_Gompertz_gradient <- function(t, param){

  # Check time larger than 0
  if (any(t <= 0)){
    stop('Times must be larger than 0')
  }

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  # time scaled
  t_scal <- t / scal

  #Derivative wrt scale
  dh_scal<-unname(
    - (shap / scal^2) * exp(t_scal) * (t_scal + 1)
  )

  #Derivative wrt shape
  dh_shap<- unname(
    exp(t_scal) / scal
  )

  #Return a length(t) x num_param (2) matrix with the gradient at each t
  return(matrix(c(dh_scal, dh_shap), ncol = length(param)))
}

#.............................................................................#
#.............................................................................#

# Names of the parameters
hazard_Gompertz_names<- function(){
  return(c('scale','shape'))
}

#.............................................................................#
#.............................................................................#

# Positivity constraint on the parameters
hazard_Gompertz_posit<- function(){
  return(c(TRUE, TRUE))
}

#.............................................................................#
#.............................................................................#

# survival
surv_Gompertz <- function(t, param){

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  #Loglogistic survival function: S(t)=exp ( -shape * (exp(t/scale) - 1))
  S <-unname(
    exp( - shap * (exp(t / scal) - 1))
  )

  return(S)
}

#.............................................................................#
#.............................................................................#

# gradient of survival
surv_Gompertz_gradient <- function(t, param){

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  # Survival function
  surv <- surv_Gompertz(t, param)

  # Scaled parameter
  t_scal <- t / scal

  #Derivative wrt scale
  dh_scal<-unname(
    surv * (shap / scal) * t_scal * exp(t_scal)
  )

  #Derivative wrt shape
  dh_shap<- unname(
    - surv *(exp(t_scal) - 1)
  )

  #Return a length(t) x num_param (2) matrix with the gradient at each t
  return(matrix(c(dh_scal, dh_shap), ncol = length(param)))
}

#.............................................................................#
#.............................................................................#
