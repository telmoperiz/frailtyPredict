##### Auxiliary functions specific to integration #####

# Draw for Monte Carlo
dist.draw<- function(n, pdf, param){

  # Select random number generation
  rfun <- paste(pdf, 'random', sep = '_')

  # Call the function
  out <- do.call(rfun, list(n = n, sig = param))

  return(out)
}

#.............................................................................#
#.............................................................................#

# Ratio between pdf and importance sampling pdf
dist.ratio <- function(seq, pdf, param, is_control){

  # Generate the data pdf
  data_pdf <- do.call(pdf, list(u = seq, sig = param))

  # Generate the Importance sampling pdf
  is_pdf <- do.call(paste('d', is_control$dist, sep = ''),
                    c(list(seq), is_control$param))

  # Return the ratio
  return(data_pdf/is_pdf)
}

#.............................................................................#
#.............................................................................#

# Transformation of van der Corput seq (uniform) to specific distribution
dist.transform <- function(seq, pdf, param){

  # Select the quantile function
  qfun <- paste(pdf, 'inv', sep = '_')

  # Call the function
  out <- do.call(qfun, list(u = seq, sig = param))

  return(out)
}

#.............................................................................#
#.............................................................................#

#Selection of the Gaussian quadrature method
gc.select <- function(n, pdf, param){

  # Select function for GQ nodes
  GQfun <- paste(pdf, 'GQ', sep = '_')

  # Call the function
  out <- do.call(GQfun, list(n = n, sig = param))

  return(out)
}

#.............................................................................#
#.............................................................................#

#' Sample and weights for numerical integral
#'
#' @param pdfinfo
#' @param intcontrol
#'
#' @return
select.intsample <- function(pdfinfo, intcontrol){

  # Nodes and weights for Gaussian Quadrature
  if (intcontrol$int_mode == 'GQ'){

    samp <- gc.select(n=intcontrol$GC_nodes, pdf=pdfinfo$func, param = pdfinfo$param)

  # For MC and QMC
  } else if (intcontrol$int_mode %in% c('MC', 'QMC')) {

    # Sample points have been previously computed
    nodes <- intcontrol$MC_samp

    # Weights to compute average (weight by sample size)
    ws <- rep(1/intcontrol$MC_N, intcontrol$MC_N)

    samp <- list(nodes=nodes, weights=ws)

  # For importance sampling
  } else if (intcontrol$int_mode %in% c('ISMC', 'ISQMC')) {

    # Sample points have been previously computed
    nodes <- intcontrol$MC_samp

    # Weights to compute average (weight by sample size)
    ws <- rep(1/intcontrol$MC_N, intcontrol$MC_N)

    # Get the pdf ratio that also goes into the weight
    pdf_ratio <- dist.ratio(nodes, pdf=pdfinfo$func, param = pdfinfo$param,
                            is_control = intcontrol$is_control)
    ws <- ws * pdf_ratio

    samp <- list(nodes=nodes, weights=ws)
  }

  return(samp)
}
#.............................................................................#
#.............................................................................#

#Computes the generic integral:
# int u^expo_u * [ S_d1^(u^alpha) - S_d2^(u^alpha) ] * S_r^u * g(u) du
# pdfinfo has in $func = shape of pdf and $param = parameters.
surv_intfrailty <- function(S_d1, S_d2=0, S_r, expo_u, alpha, pdfinfo, intcontrol){

  # Get sample and weights for numerical integral
  samp <- select.intsample(pdfinfo, intcontrol)

  # Evaluate integrand at nodes
  Integrand<-survIntegrand(samp$nodes, S_d1, S_d2, S_r, expo_u, alpha)

  #Approximate the integral
  INT<-sum(Integrand * samp$weights)

  #Return the integral
  return(INT)
}

#.............................................................................#
#.............................................................................#

#Returns the derivative of surv_intfrailty w.r.t. terms
surv_intfrailty_gradient <- function(S_d1, S_d2=0, S_r, expo_u, alpha,
                                     pdfinfo, intcontrol,
                                     terms=c('S_d1', 'S_d2', 'S_r', 'expo_u',
                                             'alpha', 'pdfparam')){

  # Get sample and weights for numerical integral
  samp <- select.intsample(pdfinfo, intcontrol)

  #Compute for each term that is asked for
  Deriv<-vector(mode = "list", length = length(terms))
  names(Deriv)<-terms

  for(term in terms){

    #Evaluate integrand in the nodes
    #Need to take a special integrand if we take derivatives of the density
    if (term=='pdfparam'){

      #Original integrand w/out taking derivatives
      RawIntegrand<-survIntegrand(samp$nodes, S_d1, S_d2, S_r, expo_u, alpha)

      #Function giving the score of the pdf wrt parameters
      scorefun<-paste(pdfinfo$func, 'score', sep='_')

      #Score of the pdf wrt the parameters
      pdfScores<-do.call(scorefun, list(u=samp$nodes, sig=pdfinfo$param))

      #Multiply the two
      Integrand<-RawIntegrand * pdfScores

    } else {
      Integrand<-survIntegrand_gradient(samp$nodes, S_d1, S_d2, S_r, expo_u, alpha,
                                        term=term)
    }

    #Approximate the integral
    INT <- sum(Integrand * samp$weights)

    #Save in list
    Deriv[term]<-INT
  }

  return(Deriv)
}

#.............................................................................#
#.............................................................................#

#generic integrand for likelihood and forecasting
# u^expo_u * [ S_d1^(u^alpha) - S_d2^(u^alpha) ] * S_r^u
survIntegrand <- function(u, S_d1, S_d2=0, S_r, expo_u, alpha){

  #frailty power parts
  u_alpha <- u^alpha
  u_power <- u^expo_u

  #Terminal power part (depending on whether there is S_d2)
  #This is to avoid cases in which 0^0=1
  if (S_d2==0){
    d_power <- S_d1^u_alpha
  } else {
    d_power <- S_d1^u_alpha - S_d2^u_alpha
  }

  #Recurrent power part
  r_power <- S_r^u

  return(u_power*d_power*r_power)
}

#.............................................................................#
#.............................................................................#

#Gradient of survIntegrand wrt term
survIntegrand_gradient <- function(u, S_d1, S_d2=0, S_r, expo_u, alpha, term){

  #Derivative of the integrand wrt S_d1
  if (term=='S_d1'){

    #Powers of u
    u_power <- u^(expo_u+alpha)
    u_alpha <- u^alpha

    #Terminal part (set infinite values to max)
    d_power<-S_d1^(u_alpha-1)
    d_power[is.infinite(d_power)] <- .Machine$double.xmax

    #Recurrent part
    r_power<-S_r^u

    #Integrand
    I<- u_power * d_power * r_power

  } else if (term=='S_d2'){

    #Powers of u
    u_power <- u^(expo_u+alpha)
    u_alpha <- u^alpha

    #Terminal part (set infinite values to max)
    d_power<- - S_d2^(u_alpha-1)
    d_power[is.infinite(d_power)] <- - .Machine$double.xmax

    #Recurrent part
    r_power<-S_r^u

    #Integrand
    I<- u_power * d_power * r_power

  } else if (term=='S_r'){

    #Powers of u
    u_power <- u^(expo_u+1)
    u_alpha <- u^alpha

    #Terminal part
    if (S_d2==0){
      d_power <- S_d1^u_alpha
    } else {
      d_power <- S_d1^u_alpha - S_d2^u_alpha
    }

    #Recurrent part (set infinite values to max)
    r_power<-S_r^(u-1)
    r_power[is.infinite(r_power)] <- .Machine$double.xmax

    #Integrand
    I<- u_power * d_power * r_power

  } else if (term=='expo_u'){

    #Powers of u
    u_power <- u^expo_u
    u_alpha <- u^alpha

    #Terminal part
    if (S_d2==0){
      d_power <- S_d1^u_alpha
    } else {
      d_power <- S_d1^u_alpha - S_d2^u_alpha
    }

    #Recurrent part
    r_power<-S_r^u

    #Integrand
    I<- log(u) * u_power * d_power * r_power

  } else if (term=='alpha'){

    #Powers of u
    u_power <- u^(expo_u+alpha)
    u_alpha <- u^alpha

    #Terminal part
    if (S_d2==0){
      d_power <- log(S_d1) * S_d1^u_alpha
    } else {
      d_power <- log(S_d1) * S_d1^u_alpha - log(S_d2) * S_d2^u_alpha
    }

    #Recurrent part
    r_power<-S_r^u

    #Integrand
    I<- log(u) * u_power * d_power * r_power

  }

  return(I)
}

#.............................................................................#
#.............................................................................#











