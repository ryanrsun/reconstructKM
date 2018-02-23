#' Fit the shape and scale parameters for a Weibull distribution
#' to the time-to-event data using MLE
#'
#' @param time A vector of event times
#' @param status A vector of 0-1 censoring status, 0 for censored, 1 for observed
#'
#' @return A list including out (the return from mle()), shape, and scale
#'
#' @export
#'
#' @examples
#' time <- rnorm(100)
#' status <- rbinom(n=100, size=0.5)
#' weimle1(time=time, status=status)
#'
weimle1 <- function(time, status){
    #--- exp(a): shape
    #--- exp(b): scale
    #--- minuslog likelihood function ---#
    ll<-function(a, b) {
        -sum(status*( log(exp(a)) + (exp(a)-1)*log(time) - exp(a)*log(exp(b)) ) -(time/exp(b))^exp(a))}
    est<-mle(minuslog=ll, start=list(a=0, b=0))
    shape<-as.numeric(exp(attributes(est)$coef[1]))
    scale<-as.numeric(exp(attributes(est)$coef[2]))
    return(list(out=est, shape=shape, scale=scale))
}


#' RMST for time-to-event data under parametric Weibull fit for data in each
#' arm separately. Also can provide CI for RMST estimate and difference in RMST.
#'
#' @param num_boots Number of bootstrap iterations
#' @param dat Data frame of time-to-event data which MUST have the columns
#' 'time', 'arm', and 'status
#' @param tau How long of a follow-up to consider, i.e. we integrate the survival
#' functions from 0 to tau
#' @param alpha Confidence interval is given for (alpha/2, 1-alpha/2) percentiles
#' @param seed For reproducibility
#'
#' @return A list including out_tab (estimate and CI in both arms), trt_rmst,
#' pbo_rmst, diff_rmst, trt_CI, pbo_CI, diff_CI. Assumes trt coded as arm 1 and
#' placebo coded as arm 0.
#'
#' @export
#'
#' @examples
#' time <- rnorm(100)
#' status <- rbinom(n=100, size=0.5)
#' arm <- c( rep(1, 50), rep(0, 50))
#' dat <- data.frame(time=time, status=status, arm=arm)
#' weibull_rmst(dat=dat, tau=1, alpha=0.05)
#'
weibull_rmst <- function(num_boots=500, dat, tau, alpha, seed=NULL) {
    if(is.numeric(seed)) {set.seed(seed)}

    # survival function of weibull, to be integrated from 0 to \tau
    weib_surv <- function(x, shape, scale) {
        1 - pweibull(q=x, shape=shape, scale=scale)
        #pweibull(q=x, shape=shape, scale=scale, lower.tail=FALSE)
    }

    if ( length(which(c('time', 'status', 'arm') %in% colnames(dat))) != 3) {
        error('Column names must include time, status, arm')
    }

    # make one dataset for each arm
    trt_dat <- dat[which(dat$arm == 1), ]
    pbo_dat <- dat[which(dat$arm == 0), ]
    boot_trt <- rep(NA, num_boots)
    boot_pbo <- rep(NA, num_boots)
    boot_diff <- rep(NA, num_boots)

    # resample
    for (i in 1:num_boots) {
        trt_samp <- dplyr::sample_n(trt_dat, nrow(trt_dat), replace=TRUE)
        pbo_samp <- dplyr::sample_n(pbo_dat, nrow(pbo_dat), replace=TRUE)

        # fit weibull parameters
        temp_trt_est <- weimle1(time=trt_samp$time, status=trt_samp$status)
        temp_pbo_est <-  weimle1(time=pbo_samp$time, status=pbo_samp$status)

        # calculate AUC
        boot_trt[i] <- integrate(f=weib_surv, shape=temp_trt_est$shape,
                                 scale=temp_trt_est$scale, lower=0, upper=tau)$value
        boot_pbo[i] <- integrate(f=weib_surv, shape=temp_pbo_est$shape,
                                 scale=temp_pbo_est$scale, lower=0, upper=tau)$value
        boot_diff[i] <- boot_trt[i] - boot_pbo[i]
    }

    # fit weibull and calculate auc for unperturbed data
    trt_est <- weimle1(time=trt_dat$time, status=trt_dat$status)
    pbo_est <- weimle1(time=pbo_dat$time, status=pbo_dat$status)
    trt_rmst <- integrate(f=weib_surv, shape=trt_est$shape, scale=trt_est$scale, lower=0, upper=tau)$value
    pbo_rmst <- integrate(f=weib_surv, shape=pbo_est$shape, scale=pbo_est$scale, lower=0, upper=tau)$value
    diff_rmst <- trt_rmst - pbo_rmst

    # also find weibull mean time-to-event
    trt_time_event_mean <- trt_est$scale * gamma(1 + 1/trt_est$shape)
    pbo_time_event_mean <- pbo_est$scale * gamma(1 + 1/pbo_est$shape)

    # find CI
    trt_CI <- quantile(boot_trt, probs=c(alpha/2, 1-alpha/2))
    pbo_CI <- quantile(boot_pbo, probs=c(alpha/2, 1-alpha/2))
    diff_CI <- quantile(boot_diff, probs=c(alpha/2, 1-alpha/2))

    # format output into df
    out_df <- data.frame( cbind(c(trt_rmst, pbo_rmst, diff_rmst),
                                rbind(trt_CI, pbo_CI, diff_CI)) )
    colnames(out_df) <- c('Est', alpha/2, 1-alpha/2)
    rownames(out_df) <- c('Arm1', 'Arm0', 'Arm1-Arm0')
    return( list(tab=out_df, trt_rmst=trt_rmst, pbo_rmst=pbo_rmst, diff_rmst=diff_rmst,
                 trt_CI=trt_CI, pbo_CI=pbo_CI, diff_CI=diff_CI,
                 trt_time_event_mean=trt_time_event_mean, pbo_time_event_mean=pbo_time_event_mean,
                 trt_param=list(shape=trt_est$shape, scale=trt_est$scale),
                 pbo_param=list(shape=pbo_est$shape, scale=pbo_est$scale)) )
}
