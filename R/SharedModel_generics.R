##### Generic functions for SharedModel not defined in base R #####

#' Asymptotic covariance matrix
#'
#' This is an S3 generic. See [asymptvar.SharedModel()] for more information.
#'
#' @export
asymptvar <- function(obj, ...){
  UseMethod("asymptvar")
}

#' Draw from the asymptotic distribution of the estimate
#'
#' This is an S3 generic. See [asympt_draw.SharedModel()] for more information.
#'
#' @export
asympt_draw <- function(obj, ...){
  UseMethod("asympt_draw")
}

#' Baseline hazard functions
#'
#' This is an S3 generic. See [base_hazard.SharedModel()] for more information.
#'
#' @export
base_hazard <- function(obj, ...){
  UseMethod("base_hazard")
}

#' Baseline survival functions
#'
#' This is an S3 generic. See [base_survival.SharedModel()] for more information.
#'
#' @export
base_survival <- function(obj, ...){
  UseMethod("base_survival")
}

#'  Berndt–Hall–Hall–Hausman Hessian
#'
#' This is an S3 generic. See [BHHH_hessian.SharedModel()] for more information.
#'
#' @export
BHHH_hessian <- function(obj, ...){
  UseMethod("BHHH_hessian")
}

#' Check whether there is sample for numerical integration
#'
#' This is an S3 generic. See [check_intsample.SharedModel()] for more information.
#'
#' @export
check_intsample <- function(obj) {
  UseMethod("check_intsample")
}

#' Test relevance of recurrent event temporal distribution
#'
#' This is an S3 generic. See [dist_relevance_test.SharedModel()] for more information.
#'
#' @export
dist_relevance_test <- function(obj, ...){
  UseMethod("dist_relevance_test")
}

#'  Number of recurrent events for each individual
#'
#' This is an S3 generic. See [get_rec_number.SharedModel()] for more information.
#'
#' @export
get_rec_number <- function(obj){
  UseMethod("get_rec_number")
}

#' Hazard given history
#'
#' This is an S3 generic. See [hazard_given_history.SharedModel()] for more information.
#'
#' @export
hazard_given_history <- function(obj, ...){
  UseMethod("hazard_given_history")
}

#' Largest time observation in the sample
#'
#' This is an S3 generic. See [max_timespan.SharedModel()] for more information.
#'
#' @export
max_timespan <- function(obj) {
  UseMethod("max_timespan")
}

#' Optimization details
#'
#' This is an S3 generic. See [optim_details.SharedModel()] for more information.
#'
#' @export
optim_details <- function(obj, ...) {
  UseMethod("optim_details")
}

#' Coefficients of covariates
#'
#' This is an S3 generic. See [param_coef.SharedModel()] for more information.
#'
#' @export
param_coef <-function(obj){
  UseMethod("param_coef")
}

#' Parameters relative to the frailty term
#'
#' This is an S3 generic. See [param_frailty.SharedModel()] for more information.
#'
#' @export
param_frailty <-function(obj){
  UseMethod("param_frailty")
}

#' Parameters of the baseline hazards
#'
#' This is an S3 generic. See [param_hazard.SharedModel()] for more information.
#'
#' @export
param_hazard <-function(obj){
  UseMethod("param_hazard")
}

#' Parameter tables for printing and latex output
#'
#' This is an S3 generic. See [param_tables.SharedModel()] for more information.
#'
#' @export
param_tables <- function(obj, ...){
  UseMethod("param_tables")
}

#' Plot with terminal event risk predictions
#'
#' This is an S3 generic. See [predict_plot.SharedModel()] for more information.
#'
#' @export
predict_plot <- function(obj, ...){
  UseMethod("predict_plot")
}

#' Boxplot with number of recurrent events
#'
#' This is an S3 generic. See [rec_boxplot.SharedModel()] for more information.
#'
#' @export
rec_boxplot <- function(obj, ...) {
  UseMethod("rec_boxplot")
}

#' Scores of the model
#'
#' This is an S3 generic. See [scores.SharedModel()] for more information.
#'
#' @export
scores <- function(obj){
  UseMethod("scores")
}

#' Survival function given history
#'
#' This is an S3 generic. See [surv_given_history.SharedModel()] for more information.
#'
#' @export
surv_given_history <- function(obj, ...){
  UseMethod("surv_given_history")
}

#' Wald tests for the parameters of the model
#'
#' This is an S3 generic. See [Wald_test.SharedModel()] for more information.
#'
#' @export
Wald_test <- function(obj, ...){
  UseMethod("Wald_test")
}
