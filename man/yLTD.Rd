% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{yLTD}
\alias{yLTD}
\title{Lag Target Variables}
\usage{
yLTD(x, leadtime = 1, timedelay = 1)
}
\arguments{
\item{x}{time series [N x 1]}

\item{leadtime}{forecast lead time [scalar >= 0], Default: 1}

\item{timedelay}{number of time delays [scalar >= 0], Default: 1}
}
\value{
forecast lead time and time delay matrix [N x timedelay+1]
}
\description{
This function time lags the target variable based on the lead time and
required time delay.
}
\details{
Based on the internal 'yLeadTimeDelay' function (used only for forecasting);See source code. This
function can be used for both modeling and forecasting.
}
