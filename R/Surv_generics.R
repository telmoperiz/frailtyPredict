##### Generic functions for Surv not defined in survival package #####

#' Get observation time
#'
#' This is an S3 generic. See [times.Surv()] for more information.
#'
#' @export
times <- function(obj){
  UseMethod("times")
}

#' Get censoring information
#'
#' This is an S3 generic. See [uncensored.Surv()] for more information.
#'
#' @export
uncensored <- function(obj){
  UseMethod("uncensored")
}
