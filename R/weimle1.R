#' Fit Weibull distribution parameters using MLE
#'
#' Fit the shape and scale parameters for a Weibull distribution
#' to the time-to-event data using MLE.
#'
#' @param time A vector of event times
#' @param status A vector of 0-1 censoring status, 0 for censored, 1 for observed
#'
#' @return A list including out (the return from mle()), shape, and scale
#'
#' @export
#'
#' @examples
#' time <- rexp(100)
#' status <- rbinom(n=100, size=1, prob=0.5)
#' weimle1(time=time, status=status)
#'
weimle1 <- function(time, status){
  #--- exp(a): shape
  #--- exp(b): scale
  #--- minuslog likelihood function ---#
  ll<-function(a, b) {
    -sum(status*( log(exp(a)) + (exp(a)-1)*log(time) - exp(a)*log(exp(b)) ) -(time/exp(b))^exp(a))}
  est<-stats4::mle(minuslog=ll, start=list(a=0, b=0))
  shape<-as.numeric(exp(attributes(est)$coef[1]))
  scale<-as.numeric(exp(attributes(est)$coef[2]))
  return(list(out=est, shape=shape, scale=scale))
}
