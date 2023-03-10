% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_metrics.R
\name{mpi}
\alias{mpi}
\title{Mean Prediction Interval (MPI)}
\usage{
mpi(X, alpha = 0.5)
}
\arguments{
\item{X}{predictions matrix [N x D]}

\item{alpha}{specifies 1 - alpha confidence level [0 < alpha < 1], Default: 0.5}
}
\value{
mean prediction interval [scalar]
}
\description{
This function calculates the mean prediction interval (MPI) given a
set of predictions 'X', the target 'Y', and a confidence level of
1 - alpha.
}
\references{
Shrestha, D.L., N. Kayastha, D. Solomatine, and R. Price (2014),
Encapsulation of parametric uncertainty statistics by various
predictive machine leanring models: MLUE method, J. Hydroinformatics,
16, 95-113.
}
