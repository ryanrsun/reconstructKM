% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonparam_rmst.R
\name{nonparam_rmst}
\alias{nonparam_rmst}
\title{Calculate RMST for each arm as well as contrast}
\usage{
nonparam_rmst(dat, tau, alpha = 0.05)
}
\arguments{
\item{dat}{Data frame of time-to-event data which MUST have the columns
'time', 'arm', and 'status}

\item{tau}{How long of a follow-up to consider, i.e. we integrate the survival
functions from 0 to tau}

\item{alpha}{Confidence interval is given for (alpha/2, 1-alpha/2) percentiles}
}
\value{
A list including data.frame of results in each arm (RMST, RMTL, SE, pvalue, CI)
as well as data.frame of results for Arm1 - Arm0 RMST.
}
\description{
Non-parametric RMST function that allows for
the tau (follow-up time) to be arbitrarily large. Uno package
restricts it to be min(last observed event in either arm).
Provides estimate, SE, CI for each arm. Provides same for
difference in arms (and also p-value).
}
\examples{
time <- rnorm(100)
status <- rbinom(n=100, size=1, prob=0.5)
arm <- c( rep(1, 50), rep(0, 50))
dat <- data.frame(time=time, status=status, arm=arm)
nonparam_rmst(dat=dat, tau=1, alpha=0.05)

}
