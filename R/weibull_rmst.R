#' RMST using Weibull fit
#'
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
#' @param find_pval Boolean, if TRUE then does bootstrap under the null to find p-value of mean difference and RMST difference
#'
#' @return A list including out_tab (estimate and CI in both arms), trt_rmst,
#' pbo_rmst, diff_rmst, trt_CI, pbo_CI, diff_CI. Assumes trt coded as arm 1 and
#' placebo coded as arm 0.
#'
#' @export
#'
#' @examples
#' time <- rexp(100)
#' status <- rbinom(n=100, prob=0.5, size=1)
#' arm <- c( rep(1, 50), rep(0, 50))
#' dat <- data.frame(time=time, status=status, arm=arm)
#' weibull_rmst(dat=dat, tau=1, alpha=0.05, num_boots=200)
#'
weibull_rmst <- function(num_boots=1000, dat, tau, alpha, find_pval=FALSE, seed=NULL) {
    if(is.numeric(seed)) {set.seed(seed)}

    # survival function of weibull, to be integrated from 0 to \tau
    weib_surv <- function(x, shape, scale) {
        1 - pweibull(q=x, shape=shape, scale=scale)
        #pweibull(q=x, shape=shape, scale=scale, lower.tail=FALSE)
    }

    if ( length(which(c('time', 'status', 'arm') %in% colnames(dat))) != 3) {
        stop('Column names must include time, status, arm')
    }

    # make one dataset for each arm
    trt_dat <- dat[which(dat$arm == 1), ]
    pbo_dat <- dat[which(dat$arm == 0), ]

    boot_df <- data.frame(rmst_trt=rep(NA, num_boots), rmst_pbo=NA,
                          rmst_diff=NA, mean_trt=NA, mean_pbo=NA,
                          mean_diff=NA)
    # resample
    for (i in 1:num_boots) {
        trt_samp <- dplyr::sample_n(trt_dat, nrow(trt_dat), replace=TRUE)
        pbo_samp <- dplyr::sample_n(pbo_dat, nrow(pbo_dat), replace=TRUE)

        # fit weibull parameters
        temp_trt_est <- weimle1(time=trt_samp$time, status=trt_samp$status)
        temp_pbo_est <-  weimle1(time=pbo_samp$time, status=pbo_samp$status)

        # calculate AUC
        boot_df$rmst_trt[i] <- integrate(f=weib_surv, shape=temp_trt_est$shape,
                                 scale=temp_trt_est$scale, lower=0, upper=tau)$value
        boot_df$rmst_pbo[i] <- integrate(f=weib_surv, shape=temp_pbo_est$shape,
                                 scale=temp_pbo_est$scale, lower=0, upper=tau)$value
        boot_df$rmst_diff[i] <- boot_df$rmst_trt[i] - boot_df$rmst_pbo[i]

        # bootstrapped means
        boot_df$mean_trt[i] <- temp_trt_est$scale * gamma(1 + 1/temp_trt_est$shape)
        boot_df$mean_pbo[i] <-  temp_pbo_est$scale * gamma(1 + 1/temp_pbo_est$shape)
        boot_df$mean_diff[i] <- boot_df$mean_trt[i] - boot_df$mean_pbo[i]
    }

    # fit weibull and calculate auc for unperturbed data
    weib_trt_params <- weimle1(time=trt_dat$time, status=trt_dat$status)
    weib_pbo_params <- weimle1(time=pbo_dat$time, status=pbo_dat$status)
    rmst_trt <- integrate(f=weib_surv, shape=weib_trt_params$shape, scale=weib_trt_params$scale, lower=0, upper=tau)$value
    rmst_pbo <- integrate(f=weib_surv, shape=weib_pbo_params$shape, scale=weib_pbo_params$scale, lower=0, upper=tau)$value
    rmst_diff <- rmst_trt - rmst_pbo

    # find CI for RMST
    rmst_trt_CI <- quantile(boot_df$rmst_trt, probs=c(alpha/2, 1-alpha/2))
    rmst_pbo_CI <- quantile(boot_df$rmst_pbo, probs=c(alpha/2, 1-alpha/2))
    rmst_diff_CI <- quantile(boot_df$rmst_diff, probs=c(alpha/2, 1-alpha/2))

    # format RMST output into df
    rmst_df <- data.frame( cbind(c(rmst_trt, rmst_pbo, rmst_diff),
                                 rbind(rmst_trt_CI, rmst_pbo_CI, rmst_diff_CI)) )
    colnames(rmst_df) <- c('Est', alpha/2, 1-alpha/2)
    rownames(rmst_df) <- c('Arm1', 'Arm0', 'Arm1-Arm0')

    # find weibull mean time-to-event
    mean_trt <- weib_trt_params$scale * gamma(1 + 1/weib_trt_params$shape)
    mean_pbo <- weib_pbo_params$scale * gamma(1 + 1/weib_pbo_params$shape)
    mean_diff <- mean_trt - mean_pbo

    # find CI for mean
    mean_trt_CI <- quantile(boot_df$mean_trt, probs=c(alpha/2, 1-alpha/2))
    mean_pbo_CI <- quantile(boot_df$mean_pbo, probs=c(alpha/2, 1-alpha/2))
    mean_diff_CI <- quantile(boot_df$mean_diff, probs=c(alpha/2, 1-alpha/2))

    # format mean output into df
    mean_df <- data.frame( cbind(c(mean_trt, mean_pbo, mean_diff),
                                rbind(mean_trt_CI, mean_pbo_CI, mean_diff_CI)) )
    colnames(mean_df) <- c('Est', alpha/2, 1-alpha/2)
    rownames(mean_df) <- c('Arm1', 'Arm0', 'Arm1-Arm0')

    # only if we want p-values
    if (find_pval) {
        null_boot_df <- data.frame(rmst_trt=rep(NA, num_boots), rmst_pbo=NA,
                              rmst_diff=NA, mean_trt=NA, mean_pbo=NA,
                              mean_diff=NA)
        # resample under null for pvalue
        for (i in 1:num_boots) {
            temp_trt_assign <- sample(x=1:nrow(dat), size=nrow(trt_dat), replace=FALSE)
            trt_samp <- dat[temp_trt_assign, ] %>%
                mutate(arm = 1)
            pbo_samp <- dat[-temp_trt_assign, ] %>%
                mutate(arm = 0)

            # fit weibull parameters
            temp_trt_est <- weimle1(time=trt_samp$time, status=trt_samp$status)
            temp_pbo_est <-  weimle1(time=pbo_samp$time, status=pbo_samp$status)

            # calculate AUC
            null_boot_df$rmst_trt[i] <- integrate(f=weib_surv, shape=temp_trt_est$shape,
                                             scale=temp_trt_est$scale, lower=0, upper=tau)$value
            null_boot_df$rmst_pbo[i] <- integrate(f=weib_surv, shape=temp_pbo_est$shape,
                                             scale=temp_pbo_est$scale, lower=0, upper=tau)$value
            null_boot_df$rmst_diff[i] <- null_boot_df$rmst_trt[i] - null_boot_df$rmst_pbo[i]

            # bootstrapped means
            null_boot_df$mean_trt[i] <- temp_trt_est$scale * gamma(1 + 1/temp_trt_est$shape)
            null_boot_df$mean_pbo[i] <-  temp_pbo_est$scale * gamma(1 + 1/temp_pbo_est$shape)
            null_boot_df$mean_diff[i] <- null_boot_df$mean_trt[i] - null_boot_df$mean_pbo[i]
        }
        pval_df <- data.frame(rmst_diff_pside = length(which(null_boot_df$rmst_diff > rmst_diff)) / num_boots,
                              rmst_diff_nside = length(which(null_boot_df$rmst_diff < rmst_diff)) / num_boots,
                              rmst_diff_2side = length(which(abs(null_boot_df$rmst_diff) > abs(rmst_diff))) / num_boots,
                              mean_diff_pside = length(which(null_boot_df$mean_diff > mean_diff)) / num_boots,
                              mean_diff_nside = length(which(null_boot_df$mean_diff < mean_diff)) / num_boots,
                              mean_diff_2side = length(which(abs(null_boot_df$mean_diff) > abs(mean_diff))) / num_boots)
    } else {
        pval_df <- NULL
    }

    return( list(rmst_df=rmst_df, mean_df=mean_df, pval_df=pval_df, rmst_trt=rmst_trt, rmst_pbo=rmst_pbo,
                 rmst_diff=rmst_diff,
                 rmst_trt_CI=rmst_trt_CI, rmst_pbo_CI=rmst_pbo_CI, rmst_diff_CI=rmst_diff_CI,
                 mean_trt=mean_trt, mean_pbo=mean_pbo,
                 trt_param=list(shape=weib_trt_params$shape, scale=weib_trt_params$scale),
                 pbo_param=list(shape=weib_pbo_params$shape, scale=weib_pbo_params$scale)) )
}
