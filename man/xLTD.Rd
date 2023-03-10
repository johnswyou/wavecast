% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{xLTD}
\alias{xLTD}
\title{Lag Input Variables}
\usage{
xLTD(x, leadtime = 1, timedelays = NULL)
}
\arguments{
\item{x}{input time series [N x D]}

\item{leadtime}{forecast lead time [scalar >= 0], Default: 1}

\item{timedelays}{number of time delays [vector (of scalar >= 0)], Default: NULL}
}
\value{
input matrix based forecast lead time and time delays [N x D + sum(timedelays)]
}
\description{
This function time lags the input variables based on the lead time and
required time delay.
}
\details{
Based on the internal 'xLeadTimeDelay' function (used only for forecasting); See source code. This
function can be used for both modeling and forecasting.
}
