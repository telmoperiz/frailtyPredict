# IN THIS README

- Brief description of the scripts
- How to implement a new S3 method for SharedModel
- How to implement a new frailty distribution
- How to implement a new hazard specification
- How to implement a new recurrent event model

## BRIEF DESCRIPTION OF THE SCRIPTS

- aux_fit.R

    Auxiliary functions to fit the model.

- aux_int.R

    Auxiliary functions for numerical integration.

- aux_other.R

    Other auxiliary	functions: gradient computation, table generation, etc.

- aux_pred.R

    Auxiliary functions for prediction.

- fun_haz_Weibull.R

    Functions required to model the baseline hazards as Weibull.

- fun_pdf_gamma.R

    Functions required to model the frailty distribution as Gamma.

- shared_frailty_fit.R

    Main function to fit the model.

- shared_frailty_pred.R

    Main function to predict the risk of terminal event.

- SharedModel_generics.R

    S3 generic function declaration for SharedModel object methods.

- SharedModel_methods.R

    efinition of S3 methods for SharedModel object.

- Surv_generics.R

    S3 generic function declaration for Surv object methods.

- Surv_methods.R

    Definition of S3 methods for Surv object.

## HOW TO IMPLEMENT A NEW S3 METHOD FOR SHAREDMODEL

1) Implement the function (method) logic in `SharedModel_methods.R`.

2) If the method is a SharedModel implementation of an already existing method (e.g., `print()`), DONE.

3) If it is not, code the generic S3 function in `SharedModel_generics.R`.

## HOW TO IMPLEMENT A NEW FRAILTY DISTRIBUTION

1) Create R file: `fun_pdf_<name>.R`, where <name> stands for the name of the new distribution.

2) In that file, implement the functions:

    a) `pdf_<name>()` returns pdf at a point.

    b) `pdf_<name>_defaults()` returns initial conditions for the parameters.

    c) `pdf_<name>_GQ()` returns Gaussian Quadrature nodes and weights.

    d) `pdf_<name>_inv()` returns quantile at point.

    e) `pdf_<name>_posit()` returns which parameters must be strictly positive.

    f) `pdf_<name>_random()` random generator of draws.

    g) `pdf_<name>_score()` returns the score (derivative of the log-likelihood).

3) In `shared_frailty_fit.R`, add <name> to the supported distributions (add to variable `frailty_supp`). Add also to function documentation.

4) Add brief script description to `README.md`.

The implementation of the gamma distribution in `fun_pdf_gamma.R` can be used as a model.

## HOW TO IMPLEMENT A NEW HAZARD SPECIFICATION

1) Create R file: `fun_haz_<name>.R`, where <name> stands for the name of the new model.

2) In that file, implement the functions:

    a) `hazard_<name>()` returns value of the hazard function.

    b) `hazard_<name>_defaults()` returns initial conditions for the parameters.

    c) `hazard_<name>_gradient()` returns the derivative of the hazard function.

    d) `hazard_<name>_names()` returns the names of the parameters.

    e) `hazard_<name>_posit()` returns which parameters must be strictly positive.

    f) `surv_<name>()` returns the value of the survival function.

    g) `surv_<name>_gradient()` returns the derivative of the survival function.

3) In `ShareModel_methods.R`, add case to the conditional in the function `dist_relevance_test.SharedModel()`.

4) In `shared_frailty_fit.R`, add <name> to the supported models (add to variable `hazard_supp`). Add also to function documentation.

5) Add brief script description to `README.md`.

The implementation of the Weibull model in `fun_haz_Weibull.R` can be used as a model.

## HOW TO IMPLEMENT A NEW RECURRENT EVENT MODEL

Need to modify the following functions:

1) If needed, add additional parameters to `shared_frailty_fit()`.

2) `hazard_given_history.SharedModel()` add the new case to the conditional and implement logic.

3) `surv_given_history.SharedModel()` add the new case to the conditional and implement logic.

4) If additional parameters, add to gradient computation in `individual_score()` (`aux_fit.R`).

5) `dist_relevance_test.SharedModel()` check whether the test applies for the model If it does, implement logic.

6) Add case to the supported models (add to variable 
`timescale_supp` in `shared_frailty_fit`).

7) Add to function documentation.