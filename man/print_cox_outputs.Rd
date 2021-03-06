% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_cox_outputs.R
\name{print_cox_outputs}
\alias{print_cox_outputs}
\title{Print outputs from Cox regression}
\usage{
print_cox_outputs(cox_fit, print_output = TRUE)
}
\arguments{
\item{cox_fit}{A model fitted with coxph()}

\item{print_output}{Print summary to screen if TRUE}
}
\value{
A list including beta, HR, SE, and CI
}
\description{
Just a wrapper to get quantities out of a call to coxph()
}
\examples{
time <- rnorm(100)
status <- rbinom(n=100, prob=0.5, size=1)
arm <- c(rep(1,50), rep(0,50))
temp_cox <- survival::coxph(survival::Surv(time, status) ~ arm)
print_cox_outputs(temp_cox)

}
