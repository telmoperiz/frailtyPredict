###### Baseline hazards: Gompertz ######

# hazard
hazard_Gompertz <- function(t, param){

  # Check time larger than 0
  if (any(t <= 0)){
    stop('Times must be larger than 0')
  }

  # Rate parameter
  r<-param['rate']

  # Level parameter
  L<-param['level']

  # Gompertz hazard function:
  # h(t)= Level * exp(rate * t)
  h <- unname(L * exp(r * t))

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
  defs$a <- unname(c(
    mod$coefficients['rate'],
    exp(mod$coefficients['log(level)'])
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

  # Rate parameter
  r<-param['rate']

  # Level parameter
  L<-param['level']

  # Exponential of r*t
  exp_rt <- exp(r * t)

  #Derivative wrt rate
  dh_rate<-unname(
    L * t * exp_rt
  )

  #Derivative wrt level
  dh_lvl<- unname(
    exp_rt
  )

  #Return a length(t) x num_param (2) matrix with the gradient at each t
  return(matrix(c(dh_rate, dh_lvl), ncol = length(param)))
}

#.............................................................................#
#.............................................................................#

# Names of the parameters
hazard_Gompertz_names<- function(){
  return(c('rate','level'))
}

#.............................................................................#
#.............................................................................#

# Positivity constraint on the parameters
hazard_Gompertz_posit<- function(){
  return(c(FALSE, TRUE))
}

#.............................................................................#
#.............................................................................#

# survival
surv_Gompertz <- function(t, param){

  # Rate parameter
  r<-param['rate']

  # Level parameter
  L<-param['level']

  # Exponential of r*t
  exp_rt <- exp(r * t)

  # Coefficient in survival
  coef <- L / r

  # Gompertz survival function: S(t)=exp ( -level / rate * (exp(rate * t) - 1))
  S <-unname(
    exp( coef * (1 - exp_rt))
  )

  return(S)
}

#.............................................................................#
#.............................................................................#

# gradient of survival
surv_Gompertz_gradient <- function(t, param){

  # Rate parameter
  r<-param['rate']

  # Level parameter
  L<-param['level']

  # Exponential of r*t
  exp_rt <- exp(r * t)

  # Coefficient in survival
  coef <- L / r

  # Survival function
  surv <- surv_Gompertz(t, param)

  #Derivative wrt rate
  dh_rate<-unname(
    surv * coef * ((1 - r * t) * exp_rt - 1) / r
  )

  #Derivative wrt level
  dh_lvl<- unname(
    surv * (1 - exp_rt) / r
  )

  #Return a length(t) x num_param (2) matrix with the gradient at each t
  return(matrix(c(dh_rate, dh_lvl), ncol = length(param)))
}

#.............................................................................#
#.............................................................................#
