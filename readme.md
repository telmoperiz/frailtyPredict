# frailtyPredict package

This package provides tools to estimate joint/shared frailty models for terminal and recurrent events. Prediction of the risk of terminal event given the recurrent event history (and covariates) is implemented.

> [!WARNING] 
> This package is in development stage. Some functionality may be unstable and it may contain bugs. Bug reports are more than welcome.

## Estimation
The estimation interface is provided by the `shared_frailty_fit function`. It allows to fit different specifications
* Recurrent event model: Poisson or renewal.
* Frailty distribution: Gamma.
* Baseline hazards: Weibull.

> [!CAUTION] 
> Estimation with analytic Hessian (Berndt–Hall–Hall–Hausman algorith) is unstable and may lead to wrong results.

## Prediction
The prediction interface is provided by the `predict` and `predict_plot` methods of `SharedModel` objects. 
* `predict` computes the prediction for a given individual and time window.
* `predict_plot` computes prediction for a moving time window and a collection of individuals. It then uses `matplotlib` to plot the results.

> [!CAUTION] 
> Package `reticulate` has issues with subprocesses on Windows and RStudio (see [here](https://github.com/rstudio/reticulate/issues/518)). Posible solution: use `predict_plot` with `render_plot=FALSE` to generate a temporal JSON file with prediction information and then call `shared_frailty_plot` from `inst/shared_frailty_plot.py`.

## TO DO
* Analytic gradient for `rec_timescale = 'piecewise-renewal' may not be correctly computed.
* Need to implement print/toLatex methods for `rec_timescale = 'piecewise-renewal'.
* Estimation with `rec_timescale = 'piecewise-renewal' is NOT documented. 