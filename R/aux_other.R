##### Other auxiliary functions #####

#list w/ body lines from textable
body_lines<-function(textable){

  #Split the table in lines
  tab_split<-stringr::str_split(textable, '\n')

  #The body of the table is between the 2nd and 3rd hline
  ind<-which(stringr::str_detect(tab_split, 'hline'))
  bseq<-seq(ind[2]+1, ind[3]-1)

  #return the body
  return(tab_split[bseq])
}

#.............................................................................#
#.............................................................................#

# Checks if certain suggested packages are installed
check_suggested_packages <- function(packages){

  for (package in packages){
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(paste('Package ', package,
                 ' must be installed to use this functionality.', sep = ''),
        call. = FALSE)
    }
  }

  # Return TRUE if all where installed
  return(TRUE)
}

#.............................................................................#
#.............................................................................#

# Find confidence intervals
conf_intervals <- function(par, se, level = 0.95, HR=FALSE){

  #Get the z value
  z<-qnorm(1-(1-level)/2)

  #Get the confidence interval
  ci <- par + c(-1, 1) * z * se

  #If HR, return the HR
  if (HR){
    ci<-exp(ci)
  }

  # Return as string
  return(paste0('(', paste0(round(ci, 3), collapse = ', '), ')'))
}

#.............................................................................#
#.............................................................................#

# Get variable names from textable line
get_line_variable <- function(line){
  varname <- strsplit(line, '&')[[1]][1]

  return(varname)
}

#.............................................................................#
#.............................................................................#

#get recurrent event times from survival object
get_rec_times <- function(hist){

  #Get uncensored times (i.e., all but the last observation, which is death/cens)
  uncens<-uncensored(hist)

  #Get times in survival object
  survtimes<-times(hist)

  #return
  return(survtimes[uncens])
}

#.............................................................................#
#.............................................................................#

#From a formula object, where the dependent variable is a Survival object,
#get the times variable and the events variable
get_time_event_variables <- function(formula){

  #Get terms
  terms<-terms(formula)

  #Get the response as string
  resp_str<-as.character(attr(terms, "variables")[2] )

  #Split the character
  #Elements where to split: "(", ", " or ")"
  #The split gives [1]="Surv"
  #                [2]=Name of the times variable
  #                [3]=Name of the event variable
  resp_split<-strsplit(resp_str, '\\(|, |\\)')[[1]]

  return(list(times=resp_split[2], event=resp_split[3]))
}

#.............................................................................#
#.............................................................................#

join_body_lines <- function(list1, list2){

  # Get list variable names
  vars1 <- sapply(list1, get_line_variable)
  vars2 <- sapply(list2, get_line_variable)

  # Get the correspondent line in list2 of each element in list1
  corresp = match(vars1, vars2)

  # Variables only in the second list
  only2 <- setdiff(seq_along(vars2), corresp)

  # Take out end of line characters '\\'
  list1 <- lapply(list1, gsub, pattern = '\\\\', replacement = '', fixed = TRUE)
  list2 <- lapply(list2, gsub, pattern = '\\\\', replacement = '', fixed = TRUE)

  # Remove First cell from list2 for the ones that are matched
  list2 <- lapply(list2,
                  function(l) paste(strsplit(l, split = '&')[[1]][-1],
                                    collapse = '&'))

  # Line for variables that do not match
  NAline <- paste('- ',
                  paste(rep(' & - ', stringr::str_count(list2[[1]], '&')),
                        collapse = ''),
                  sep = '')

  # Fill the list with lines (first the ones in list1)
  outlist <- vector(mode = 'list', length = length(list1) + length(only2))
  for (l in seq_along(list1)){

    # If it is not matched
    if (is.na(corresp[l])){
      matched_line <- NAline
    } else {
      matched_line <- list2[[corresp[l]]]
    }

    outlist[[l]] <- paste(list1[[l]],
                          '&',
                          matched_line,
                          '\\\\',
                          sep = '')

  }

  # Fill the ones in list2 that are not matched
  for (l in seq_along(only2)){
    outlist[[length(list1) + l]] <- paste(vars2[only2[l]],
                                          '&',
                                          NAline,
                                          '&',
                                          list2[[only2[l]]],
                                          '\\\\',
                                          sep = '')
  }

  return(outlist)

}

#.............................................................................#
#.............................................................................#

#Write the dataset in the coresponding time scale (gap or calendar)
# Used to fit defaults in hazard_<...>_defaults
modify_timescale <- function(OBSL, timescale){

  #If it is Poisson, just leave it unchanged
  if (timescale == "Poisson"){

    rec_times_mod<-lapply(OBSL, function(obs) obs$rec_times)

    #If renewal, take the difference and turn it into a survival object
  } else if (timescale %in% c("renewal", "piecewise-renewal")){

    rec_times_mod<-lapply(OBSL, function(obs)
      survival::Surv(diff(c(0, times(obs$rec_times))),
                     uncensored(obs$rec_times)) )
  }

  #Return the results as a vector
  return(do.call(c, rec_times_mod))
}

#.............................................................................#
#.............................................................................#

# gradient for renewal model
# (used at surv_given_history.SharedModel, rec_piece_Surv, rec_piece_weights)
renewal_grad <- function(Sj, Sj_grad, weights = 1){

  # Multiply all the value of S for gap times
  S <- prod(Sj^weights)

  # For each gap, compute the derivative to value ratio
  dv_ratio <- weights * Sj_grad/replicate(n=ncol(Sj_grad), Sj)

  # Compute the sum across rows (ie., sum in gap times)
  dv_ratio_sum <- colSums(dv_ratio)

  # Multiply by the product across gap times
  grad <- S * dv_ratio_sum
}

#.............................................................................#
#.............................................................................#

# Index and exponent for each recurrent time
# (used at rec_piece_Surv)
rec_piece_exponent <- function(ts, model){

  # Remove first time in ts, which corresponds to zero
  ts <- ts[-1]

  # Get piecewise times and constants
  taus <- c(0, model$funtion_shapes$rec_piecewise_ts)
  pw_C <- c(1, param_hazard(model)$piecewise_Cs)

  ind_C <- sapply(ts, function(t) sum(taus < t))

  # Remove the first
  return(list(ind_C = ind_C, exp_C = pw_C[ind_C]))
}

#.............................................................................#
#.............................................................................#

# Compute Poisson (piecewise exponential) gap-time survivals
# (used at surv_given_history.SharedModel)
rec_piece_Surv <- function(ts, model, gradient = FALSE){

  # Find the gaps between times
  gaps <- diff(ts)

  # Survival function at the gaps
  Sb <- base_survival(gaps, obj=model, process = 'recurrent', gradient = gradient)

  # Find Poisson part weights for each gap
  weight <- rec_piece_weights(ts, model=model, gradient=gradient)

  # Get exponents for each gap
  expo <- rec_piece_exponent(ts, model=model)

  # Find conditional survival for each t
  Sjs <- weight * (Sb ^ expo$exp_C)

  if (gradient){

    # Get gradients
    Sb_grad <- attr(Sb, 'gradient')
    w_grad <- attr(weight, 'gradient')
    w_grad_C <- attr(weight, 'gradient_C')

    # Initialize
    Sjs_grad <- matrix(0, nrow = length(Sjs), ncol = length(model$par_pos$a_r))
    Sjs_grad_C <- matrix(0, nrow = length(Sjs), ncol = length(model$par_pos$pw_C))

    for (j in seq_along(Sjs)){

      # Gradient of product w.r.t. hazard parameters
      Sjs_grad[j,] <- renewal_grad(c(weight[j], Sb[j]),
                                   rbind(w_grad[j,], Sb_grad[j,]),
                                   c(1, expo$exp_C[j])
                                   )

      # Gradient of product w.r.t. piecewise constants
      w_deriv <- (Sb[j] ^ expo$exp_C[j]) * w_grad_C[j,]
      Sb_deriv_ind <- seq(from = 2, to = ncol(Sjs_grad_C) + 1) == expo$ind_C[j]
      Sjs_grad_C[j,] <-  w_deriv + weight[j] * log(Sb[j]) * Sb_deriv_ind
    }

    # Store
    attr(Sjs, 'gradient') <- Sjs_grad
    attr(Sjs, 'gradient_C') <- Sjs_grad_C
  }

  return(Sjs)
}

#.............................................................................#
#.............................................................................#

# Compute Poisson (piecewise exponential) weights for a given t
# (used at rec_piece_Surv)
rec_piece_weights <- function(ts, model, gradient = FALSE){

  # Get piecewise times and constants
  taus <- c(0, model$funtion_shapes$rec_piecewise_ts)
  pw_C <- c(1, param_hazard(model)$piecewise_Cs)

  ws <- numeric(length(ts) - 1)

  if(gradient){
    ws_grad <- matrix(0, nrow = length(ts) - 1, ncol = length(model$par_pos$a_r))
    ws_grad_C <- matrix(0, nrow = length(ts) - 1, ncol = length(model$par_pos$pw_C))
  }

  # For each subinterval
  for (j in seq(length(ts) - 1)){

    # Indexes for taus that fall in the subinterval
    # here index 1 never ocurs since taus[1] = 0
    tau_ind <- which(ts[j] < taus & taus < ts[j+1])

    if (length(tau_ind) == 0){
      ws[j] <- 1
      next
    }

    # Find the gaps between taus and time
    tau_gaps <- taus[tau_ind] - ts[j]

    # Find survival at gaps
    S_tau_gaps <- base_survival(model, t=tau_gaps, process = 'recurrent',
                                gradient = gradient)

    # Find the corrsponding constants
    c_expon <-  pw_C[tau_ind] - pw_C[tau_ind - 1]

    # Exponentiate each s_gap and multiply them
    ws[j] <- prod(S_tau_gaps^c_expon)

    if (gradient){

      # Get gradient wrt hazard of baseline survival
      S_tau_grad <- attr(S_tau_gaps, 'gradient')

      # Gradient wrt hazard parameters
      ws_grad[j,] <- renewal_grad(S_tau_gaps, S_tau_grad, c_expon)

      # Gracient wrt piecewise constants
      log_diff <- diff(c(0, log(S_tau_gaps)))
      ws_grad_C[j, tau_ind - 1] <- ws[j] * log_diff
    }
  }

  # Add gradients as attributed
  if (gradient){
    attr(ws, 'gradient') <- ws_grad
    attr(ws, 'gradient_C') <- ws_grad_C
  }

  return(ws)
}

#.............................................................................#
#.............................................................................#

textable_vertical <- function(Tabs, TexTabs, namesD, namesR){

  #Write preamble of table and assign Latex class
  textable<-paste('% latex table generated in ', version$version.string,
                  ' with the help of xtable package.\n',
                  '% ', Sys.time(), '\n', sep = '')
  class(textable)<-'Latex'

  #Get the info in one of the xtables to adapt the table
  tab_ex<-xtable::xtable(Tabs[[1]])

  #Alignmment and number of columns
  col_align <- attr(tab_ex, "align")
  ncol<-length(col_align)+1

  #Column names
  col_names<-attr(tab_ex, "names")

  ##Building the table
  #Table preamble
  textable<-paste(textable, '\\begin{table}[ht]\n', '\\centering\n',
                  '\\begin{tabular}{r', paste(col_align, collapse = ''), '}\n',
                  '\\hline\n', sep='')

  #First row w/ column names
  textable<-paste(textable, ' & & ', paste(xtable::sanitize(col_names),
                                           collapse = ' & '),
                  ' \\\\', '\n', '\\hline\n', sep='')

  #Terminal event
  tab_body_beta<-body_lines(TexTabs$beta_d)
  tab_body_a<-body_lines(TexTabs$a_d)
  total_lines<-length(tab_body_beta)+length(tab_body_a)

  #name of the event
  textable<-paste(textable, '\\multirow{', total_lines ,'}{*}{',
                  namesD$event, '}', sep='')

  #Write lines: regression coefficients
  for (l in seq_along(tab_body_beta)){
    textable<-paste(textable, ' & ', tab_body_beta[[l]], '\n', sep = '')
  }

  #line to separate the hazard coefficients
  textable<-paste(textable, '\\cline{2-', ncol ,'} \n', sep='')

  #Write lines: hazard coefficients
  for (l in seq_along(tab_body_a)){
    textable<-paste(textable, ' & ', tab_body_a[[l]], '\n', sep = '')
  }

  #hline
  textable<-paste(textable, '\\hline\n', sep='')

  #Recurrent event
  tab_body_beta<-body_lines(TexTabs$beta_r)
  tab_body_a<-body_lines(TexTabs$a_r)
  total_lines<-length(tab_body_beta)+length(tab_body_a)

  #name of the event
  textable<-paste(textable, '\\multirow{', total_lines ,'}{*}{',
                  namesR$event, '}', sep='')

  #Write lines: regression coefficients
  for (l in seq_along(tab_body_beta)){
    textable<-paste(textable, ' & ', tab_body_beta[[l]], '\n', sep = '')
  }

  #line to separate the hazard coefficients
  textable<-paste(textable, '\\cline{2-', ncol ,'} \n', sep='')

  #Write lines: hazard coefficients
  for (l in seq_along(tab_body_a)){
    textable<-paste(textable, ' & ', tab_body_a[[l]], '\n', sep = '')
  }

  #hline
  textable<-paste(textable, '\\hline\n', sep='')

  #Frailty
  tab_body_f<-body_lines(TexTabs$frail)

  #name of the event
  textable<-paste(textable, '\\multirow{', length(tab_body_f) ,'}{*}{',
                  'Frailty', '}', sep='')

  #Write lines: regression coefficients
  for (l in seq_along(tab_body_f)){
    textable<-paste(textable, ' & ', tab_body_f[[l]], '\n', sep = '')
  }

  #End table
  textable<-paste(textable, '\\hline\n\\end{tabular}\n\\end{table}', sep='')

  return(textable)
}

#.............................................................................#
#.............................................................................#

textable_horizontal <- function(Tabs, TexTabs, namesD, namesR){

  #Write preamble of table and assign Latex class
  textable<-paste('% latex table generated in ', version$version.string,
                  ' with the help of xtable package.\n',
                  '% ', Sys.time(), '\n', sep = '')
  class(textable)<-'Latex'

  #Get the info in one of the xtables to adapt the table
  tab_ex<-xtable::xtable(Tabs[[1]])

  #Alignmment and number of columns
  col_align <- attr(tab_ex, "align")[-1]
  eachcol <- length(col_align)
  ncol <- 2 * eachcol + 1

  #Column names
  col_names<-attr(tab_ex, "names")

  ##Building the table
  #Table preamble
  textable<-paste(textable, '\\begin{table}[ht]\n', '\\centering\n',
                  '\\begin{tabular}{r',
                  paste(rep(col_align, 2), collapse = ''),
                  '}\n',
                  sep='')

  # First row w/ event names
  textable <- paste(textable,
                    ' & ',
                    '\\multicolumn{', eachcol, '}{c}{', namesD$event, '} &',
                    '\\multicolumn{', eachcol, '}{c}{', namesR$event, '} \\\\',
                    '\n', '\\cmidrule(rl){2-', 1 + eachcol, '}',
                    '\\cmidrule(rl){', 2 + eachcol, '-', ncol, '}',
                    '\n',
                    sep = '')

  # Second row w/ column names
  textable<-paste(textable,
                  ' & ',
                  paste(xtable::sanitize(rep(col_names, 2)), collapse = ' & '),
                  ' \\\\', '\n', '\\hline\n', sep='')

  # Terminal and recurrent event coefficients (joint)
  body_beta_d <- body_lines(TexTabs$beta_d)
  body_beta_r <- body_lines(TexTabs$beta_r)
  body_beta <- join_body_lines(body_beta_d, body_beta_r)

  for (line in body_beta){
    textable <- paste(textable, line, '\n', sep = '')
  }

  # Hazards (joint)
  body_a_d <- body_lines(TexTabs$a_d)
  body_a_r <- body_lines(TexTabs$a_r)
  body_a <- join_body_lines(body_a_d, body_a_r)

  # Preamble
  textable <- paste(textable,
                    '\\hline\n',
                    '\\emph{Hazards:} & ',
                    paste(rep('&', ncol - 2), collapse = ' '),
                    '\\\\', '\n',
                    sep = '')

  # Parameters
  for (line in body_a){
    textable <- paste(textable, line, '\n', sep = '')
  }

  # Frailty
  body_frail <- body_lines(TexTabs$frail)

  # Preamble
  textable <- paste(textable,
                    '\\hline\n',
                    '\\emph{Frailty:} & ',
                    paste(rep('&', ncol - 2), collapse = ' '),
                    '\\\\', '\n',
                    sep = '')

  #Write lines: regression coefficients
  for (line in body_frail){

    # Remove line break
    line <- gsub(pattern = '\\\\', replacement = '', x=line, fixed = TRUE)
    line <- paste(line,
                  paste(rep('&', eachcol), collapse = ' '),
                  '\\\\',
                  sep = '')

    textable <- paste(textable, line, '\n', sep = '')
  }

  #End table
  textable<-paste(textable, '\\hline\n\\end{tabular}\n\\end{table}', sep='')

  return(textable)
}

#.............................................................................#
#.............................................................................#

#small implementation of z test w/ estimated parameters
z_test <- function(par, se, par0 = 0, one.sided=FALSE){

  #Fing Z-score
  z<-abs(par-par0)/se

  pval<- (2-as.integer(one.sided)) * stats::pnorm(-z, lower.tail = TRUE)

  return(list(z=z, pval=pval))
}

#.............................................................................#
#.............................................................................#
