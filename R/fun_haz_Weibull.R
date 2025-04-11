###### Baseline hazards: Weibull ######

# hazard
hazard_Weibull <- function(t, param){
  
  # Check time larger than 0
  if (any(t <= 0)){
    stop('Times must be larger than 0')
  }

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  #Weibull hazard function: h(t)= (shape/scale) * (t/scale)^(shape-1)
  h <- unname((shap/scal) * (t/scal)^(shap-1))

  return(h)
}

#.............................................................................#
#.............................................................................#

# defaults
hazard_Weibull_defaults<- function(data, terminal_formula, recurrent_formula,
                                   OBSL,  rec_timescale, eventtype){

  #Initialize defaults
  defs<-vector(mode="list", length = 0)

  if (eventtype == "terminal"){

    #Run parametric Weibull model for the terminal event
    #We update the formula to always include an intercept (scale parameter)
    mod<-survival::survreg(formula = update.formula(terminal_formula, ~ . +1 ),
                           data = data, dist = "weibull")

  } else if (eventtype == "recurrent") {

    #Parametric survival regression for defaults (in gap or calendar times)
    #Update data to create variable t.mod (survival object in the correct timescale)
    data_mod<-data
    data_mod[,'t.mod']<-modify_timescale(OBSL, rec_timescale)

    #Update the formula to have t.mod as dependent variable (and intercept)
    mod<-survival::survreg(formula = update.formula(recurrent_formula, t.mod ~ . +1),
                           data = data_mod, dist = "weibull")
  }

  #Get the defaults
  #Cox coefficients
  defs$beta <- unname(mod$coefficients[-1]) #w/out intercept

  #Weibull scale and shape
  defs$a <- unname(c( exp(mod$coefficients[1]), (mod$scale)^(-1) ))

  return(defs)
}

#.............................................................................#
#.............................................................................#

# gradient of hazard
hazard_Weibull_gradient <- function(t, param){
  
  # Check time larger than 0
  if (any(t <= 0)){
    stop('Times must be larger than 0')
  }

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  #Derivative wrt scale
  dh_scal<-unname( -shap^2 * t^(shap-1) * scal^(-shap-1) )

  #Derivative wrt shape
  dh_shap<- unname( t^(shap-1) * scal^(-shap) * ( 1 + shap*log(t/scal)) )

  #Return a length(t) x num_param (2) matrix with the gradient at each t
  return(matrix(c(dh_scal, dh_shap), ncol = length(param)))
}

#.............................................................................#
#.............................................................................#

# Names of the parameters
hazard_Weibull_names<- function(){
  return(c('scale','shape'))
}

#.............................................................................#
#.............................................................................#

# Positivity constraint on the parameters
hazard_Weibull_posit<- function(){
  return(c(TRUE, TRUE))
}

#.............................................................................#
#.............................................................................#

# survival
surv_Weibull <- function(t, param){

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  #Weibull survival function: S(t)=exp[- (t/scale)^shape]
  S <-unname(exp(- (t/scal)^shap))

  return(S)
}

#.............................................................................#
#.............................................................................#

# gradient of survival
surv_Weibull_gradient <- function(t, param){

  #Scale parameter
  scal<-param['scale']

  #Shape parameter
  shap<-param['shape']

  # Survival function
  surv <- surv_Weibull(t, param)

  #Derivative wrt scale
  dh_scal<-unname( shap * ( t^shap / scal^(shap + 1) ) * surv )

  #Derivative wrt shape
  dh_shap<- unname( - (t/scal)^shap * log(t/scal) * surv )

  #Return a length(t) x num_param (2) matrix with the gradient at each t
  return(matrix(c(dh_scal, dh_shap), ncol = length(param)))
}

#.............................................................................#
#.............................................................................#
