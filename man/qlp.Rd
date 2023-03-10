% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_metrics.R
\name{qlp}
\alias{qlp}
\title{p-Quantile Loss (QLp)}
\usage{
qlp(Y, X, p = 0.5)
}
\arguments{
\item{Y}{target vector (response variable) [N x 1]}

\item{X}{predictions matrix [N x D]}

\item{p}{quantile [0 < p < 1], Default: 0.5}
}
\value{
p-quantile loss [scalar]
}
\description{
This function calculates the p-quantile loss (QLp) or p-risk given a set of predictions
'X', the target 'Y', and a quantile 'p' (0, 1).
}
\references{
Seeger, M.W., Salinas, D., Flunkert, V., 2016. Bayesian Intermittent Demand Forecasting
for Large Inventories, in: Lee, D.D., Sugiyama, M., Luxburg, U. V, Guyon, I., Garnett, R.
(Eds.), Advances in Neural Information Processing Systems 29. Curran Associates, Inc.,
pp. 4646-4654.

#' Li, S., Jin, X., Xuan, Y., Zhou, X., Chen, W., Wang, Y.-X., Yan, X., 2019. Enhancing
the Locality and Breaking the Memory Bottleneck of Transformer on Time Series Forecasting.
arXiv.

D. Salinas, V. Flunkert, J. Gasthaus et al., 2019, DeepAR: Probabilistic forecasting
with autoregressive recurrent networks. International Journal of Forecasting,
https://doi.org/10.1016/j.ijforecast.2019.07.001.
}
