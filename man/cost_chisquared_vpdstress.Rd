% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cost_functions.R
\name{cost_chisquared_vpdstress}
\alias{cost_chisquared_vpdstress}
\title{Chi-squared cost function}
\usage{
cost_chisquared_vpdstress(par, obs, drivers, inverse = FALSE)
}
\arguments{
\item{par}{parameters}

\item{obs}{observed values}

\item{drivers}{drivers}

\item{inverse}{invert the function}
}
\value{
The Chi-squared value on VPD
}
\description{
Cost function using the chi-squared statistic,
after (Keenan et al., 2012 GCB)
}
