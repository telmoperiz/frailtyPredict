###### Frailty pdfs: gamma ######

pdf_gamma <- function(u, sig) {
  
  invpar<-1/sig #Inverse of parameter
  
  #Return the pdf
  return(stats::dgamma(u, shape = invpar, scale = sig))
}

#.............................................................................#
#.............................................................................#

#defaults
pdf_gamma_defaults <- function(){
  return(1)
}


#.............................................................................#
#.............................................................................#

pdf_gamma_GQ <- function(n, sig) {
  
  invpar<-1/sig #Inverse of parameter
  
  #Return the Gaussian Quadrature nodes and weights
  return(statmod::gauss.quad.prob(n, dist = "gamma", alpha = invpar, beta = sig))
}

#.............................................................................#
#.............................................................................#

pdf_gamma_inv <- function(u, sig) {
  
  invpar<-1/sig #Inverse of parameter
  
  #Return the quantiles
  return(stats::qgamma(u, shape = invpar, scale = sig))
}

#.............................................................................#
#.............................................................................#

# Positivity constraint on the parameters
pdf_gamma_posit<- function(){
  return(TRUE)
}

#.............................................................................#
#.............................................................................#

pdf_gamma_random <- function(n, sig) {
  
  invpar<-1/sig #Inverse of parameter
  
  #Return random numbers
  return(stats::rgamma(n, shape = invpar, scale = sig))
}

#.............................................................................#
#.............................................................................#

#derivative of log density:
# (u - log(u) + log(sig) - 1 + digamma(1/sig))/sig^2
pdf_gamma_score <- function(u, sig) {
  
  #u part of the numerator
  u_num <- u - log(u)
  
  #sig part of the numerator
  sig_num <- log(sig) - 1 + digamma(1/sig)
  
  #Return the score
  return((u_num + sig_num ) / ( sig^2 ))
}