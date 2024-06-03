##### Methods Surv S3 class #####

#' Get observation time
#'
#' Find the observation time from [survival::Surv()] object.
#'
#' @param obj [survival::Surv()] object
#'
#' @return numeric time value in object
#' @export
times.Surv <- function(obj){
  #Return the time information in the survival object
  return(unname(obj[,1]))
}

#.............................................................................#
#.............................................................................#

#' Get censoring information
#'
#' Find censoring information from [survival::Surv()] object.
#'
#' @param obj [survival::Surv()] object
#'
#' @return boolean `TRUE` if observation is not censored
#' @export
uncensored.Surv <- function(obj){
  #Return the censoring information in the survival object
  # 1 - uncensored, 0 - censored
  return(unname(as.logical(obj[,2])))
}

#.............................................................................#
#.............................................................................#
