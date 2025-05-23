% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SharedModel_methods.R
\name{toLatex.SharedModel}
\alias{toLatex.SharedModel}
\title{Latex tables}
\usage{
\method{toLatex}{SharedModel}(
  obj,
  filename = "tab.tex",
  cols = c("est", "hr", "pvals"),
  col_names = NULL,
  merged = TRUE,
  digits = 3,
  horizontal = TRUE,
  append = FALSE,
  BHHH = TRUE
)
}
\arguments{
\item{obj}{\code{\link[=SharedModel]{SharedModel()}} object}

\item{filename}{Path where tables are saved}

\item{cols}{character vector with values to display. A subset of
\code{c('est', 'est_SE', 'est_CI', 'hr', 'hr_SE', 'hr_CI', 'pvals')}.}

\item{col_names}{Names of the columns.}

\item{merged}{\code{TRUE} (default) to build a single table. \code{FALSE} creates
separated tables for fit results.}

\item{digits}{Digits to round results.}

\item{horizontal}{\code{TRUE} (default) to order the results regarging the terminal
and recurrent events horizontally.}

\item{append}{\code{TRUE} to append results to \code{filename}. Default is \code{FALSE}.}

\item{BHHH}{\code{TRUE} (default) to approximate the Information Matrix by the
sum of squared scores. See \code{\link[=asymptvar.SharedModel]{asymptvar.SharedModel()}}.}
}
\value{
Silent. Creates files with latex tables.
}
\description{
Creates tables with fit information and exports them in latex format.
}
