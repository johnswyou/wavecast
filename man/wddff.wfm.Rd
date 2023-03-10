% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wddff.R
\name{wddff.wfm}
\alias{wddff.wfm}
\title{Model Agnostic WDDFF}
\usage{
wddff.wfm(ts, W_ts, wfm, lags, maxlag, nbc)
}
\arguments{
\item{ts}{time series/matrix/vector [N x 1]}

\item{W_ts}{wavelet decomposed time series [N x J+1]}

\item{wfm}{wavelet-based forecasting method [string; 'single', 'within', 'across','single_hybrid','within_hybrid','across_hybrid']}

\item{lags}{lag length to include as predictors [integer >= 0]}

\item{maxlag}{maximum lag length [1 < integer <= lags]}

\item{nbc}{number of boundary-effected coefficients [integer >= 0]}
}
\value{
list of input-output dataset(s) according to wfm [N x D x J+1]
}
\description{
This function creates input-output (IO) datasets for a given wavelet-based
forecasting method, the number of time series lags to consider as predictors,
and the number of boundary-effected coefficients to remove.
}
\references{
Quilty, J., and J. Adamowski (2018), Addressing the incorrect usage of wavelet-
based hydrological and water resources forecasting models for real-world
applications with best practices and a new forecasting framework, J. Hydrol.,
doi:10.1016/j.jhydrol.2018.05.003.
}
