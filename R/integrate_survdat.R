#' Integrate area under curve for single arm
#'
#' Calculate nonparametric RMST for a single arm up to tau for data.frame with time and status
#'
#' @param dat Data frame of time-to-event data which MUST have the columns
#' 'time' and 'status' exactly
#' @param tau The cutoff time, a scalar
#' @param alpha Level for confidence interval
#'
#' @return data.frame with rows for RMST and RMTL and columnns for estimate, std err, pvalue, and CI
#' @import dplyr magrittr survival survminer stats4
#' @importFrom stats integrate pchisq pweibull qnorm quantile runif time
#' @importFrom graphics hist
#' @importFrom rlang .data
#' @export
#' @examples
#'
#' time <- rnorm(100)
#' status <- rbinom(n=100, size=1, prob=0.5)
#' dat <- data.frame(time=time, status=status)
#' integrate_survdat(dat=dat, tau=2)
#'
integrate_survdat <- function(dat, tau, alpha=0.05) {
  #make sure data has necessary arms
  if ( length(which(c('time', 'status') %in% colnames(dat))) != 2) {
    stop('Column names must include time and status exactly.')
  }

  # get the KM table using prebuilt functions
  # it only tabulates the unique times
  # disregard those after \tau, then add a row for \tau
  KMfit <- survival::survfit(Surv(time, status) ~ 1, data=dat)
  n <- nrow(dat)
  nUnique <- length(KMfit$time)
  KMtab <- data.frame(time = KMfit$time, surv = KMfit$surv, nRisk = KMfit$n.risk, nEvent = KMfit$n.event) %>%
    dplyr::filter(.data$time <= tau) %>%
    dplyr::add_row(time=tau, surv=KMfit$surv[nUnique]) %>%
    # use lead() to add the "next time" column
    dplyr::mutate(nextTime = lead(.data$time)) %>%
    dplyr::mutate(diffTime = .data$nextTime - .data$time) %>%
    # multiply the survival by the difference in times to get each auc chunk
    dplyr::mutate(aucChunk = .data$diffTime * .data$surv) %>%
    # filter out the last NA row now that it's served its purpose
    dplyr::filter(!is.na(.data$aucChunk)) %>%
    # need this part for the variance calculation
    dplyr::mutate(cumAUC = cumsum(.data$aucChunk)) %>%
    dplyr::mutate(cumAUC = .data$cumAUC + .data$time[1]) %>%
    dplyr::mutate(intervalAUC = .data$cumAUC[length(.data$cumAUC)] - .data$cumAUC + .data$aucChunk)

  # have to add the first chunk to the auc
  aucTot <- KMtab$cumAUC[nrow(KMtab)]

  # variance estimate
  multiplier <- ifelse(KMtab$nRisk - KMtab$nEvent == 0, 0, KMtab$nEven / (KMtab$nRisk * (KMtab$nRisk - KMtab$nEvent)))
  varHat <- sum( KMtab$intervalAUC^2 * multiplier )

  # return data.frame
  returnDF <- data.frame(Stat = c("RMST", "RMTL"), Est = c(aucTot, tau - aucTot),
                         se = c(sqrt(varHat), sqrt(varHat)),
                         pval = c(1 - pchisq(aucTot^2 / varHat, df=1), 1 - pchisq((tau - aucTot)^2 / varHat, df=1)),
                         CIlower = c(aucTot - qnorm(1 - alpha/2) * sqrt(varHat), tau - aucTot - qnorm(1 - alpha/2) * sqrt(varHat)),
                         CIupper = c(aucTot + qnorm(1 - alpha/2) * sqrt(varHat), tau - aucTot + qnorm(1 - alpha/2) * sqrt(varHat)))

  return(returnDF)
}
