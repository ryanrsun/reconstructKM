#' Calculate nonparametric RMST for data.frame with time, status, arm

#' @param dat Data frame of time-to-event data which MUST have the columns
#' 'time', 'arm', and 'status
#'
#' @return list with rmst_pbo (arm=0) and rmst_trt (arm=1)
#' @export
#' @examples
#' time <- rnorm(100)
#' status <- rbinom(n=100, size=0.5)
#' arm <- c( rep(1, 50), rep(0, 50))
#' dat <- data.frame(time=time, status=status, arm=arm)
#' integrate_survdat(dat=dat)
#'
integrate_survdat <- function(dat, tau) {
    #make sure data has necessary arms
    if ( length(which(c('time', 'status', 'arm') %in% colnames(dat))) != 3) {
        error('Column names must include time, status, arm')
    }

    # get the KM table using prebuilt functions
    KM_fit <- survival::survfit(Surv(time, status) ~ arm, data=dat)

    # make sure the arms are coded as 0/1 to get the correct KM table
    num_zero <- KM_fit$strata[1]
    num_one <- KM_fit$strata[2]
    KM_pbo <- data.frame(time = KM_fit$time[1:num_zero],
                         surv = KM_fit$surv[1:num_zero])
    KM_trt <- data.frame(time = KM_fit$time[(num_zero+1):(num_zero+num_one)],
                         surv = KM_fit$surv[(num_zero+1):(num_zero+num_one)])

    # depending on \tau, we need to add or remove rows
    if (tau %in% KM_trt$time) {
        KM_trt <- KM_trt %>% filter(time <= tau)
    } else if (tau > KM_trt$time[num_one]) {
        KM_trt <- rbind(KM_trt, c(tau, KM_trt$surv[num_one]))
    } else if (tau < KM_trt$time[num_one]) {
        KM_trt <- KM_trt %>% filter(time < tau)
        KM_trt <- rbind(KM_trt, c(tau, KM_trt$surv[nrow(KM_trt)]))
    }

    if (tau %in% KM_pbo$time) {
        KM_pbo <- KM_pbo %>% filter(time <= tau)
    } else if (tau > KM_pbo$time[num_zero]) {
        KM_pbo <- rbind(KM_pbo, c(tau, KM_pbo$surv[num_zero]))
    } else if (tau < KM_pbo$time[num_zero]) {
        KM_pbo <- KM_pbo %>% filter(time < tau)
        KM_pbo <- rbind(KM_pbo, c(tau, KM_pbo$surv[nrow(KM_pbo)]))
    }

    # sum intervals in dplyr
    KM_pbo <- KM_pbo %>%
        mutate(next_time = lead(time)) %>%
        mutate(diff_time = next_time - time) %>%
        mutate(surv_chunk = diff_time * surv)
    # this table doesn't start at 0, have to add that first interval
    rmst_pbo <- sum(KM_pbo$surv_chunk[1:(nrow(KM_pbo) - 1)]) +
        KM_pbo$time[1]

    KM_trt <- KM_trt %>%
        mutate(next_time = lead(time)) %>%
        mutate(diff_time = next_time - time) %>%
        mutate(surv_chunk = diff_time * surv)
    # this table doesn't start at 0, have to add that first interval
    rmst_trt <- sum(KM_trt$surv_chunk[1:(nrow(KM_trt)-1)]) +
        KM_trt$time[1]

    return(list(rmst_pbo=rmst_pbo, rmst_trt=rmst_trt))
}

#' My personal non-parametric RMST function that allows for
#' the tau (follow-up time) to be arbitrarily large. Uno package
#' restricts it to be min(last observed event in either arm).
#' Provides estimate, SE, CI for each arm. Provides same for
#' difference in arms (and also p-value).
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
#' nonparam_rmst(dat=dat, tau=1, alpha=0.05)
#'
nonparam_rmst <- function(num_boots=1000, dat, tau, alpha, find_pval=FALSE, seed=NULL) {

     if(is.numeric(seed)) {set.seed(seed)}

    # make sure data has necessary arms
    if ( length(which(c('time', 'status', 'arm') %in% colnames(dat))) != 3) {
        error('Column names must include time, status, arm')
    }

    # make one dataset for each arm
    trt_dat <- dat[which(dat$arm == 1), ]
    pbo_dat <- dat[which(dat$arm == 0), ]

    boot_df <- data.frame(rmst_trt=rep(NA, num_boots), rmst_pbo=NA,
                          rmst_diff=NA)
    # resample under the distribution of the data (not the null)
    for (i in 1:num_boots) {
        trt_samp <- dplyr::sample_n(trt_dat, nrow(trt_dat), replace=TRUE)
        pbo_samp <- dplyr::sample_n(pbo_dat, nrow(pbo_dat), replace=TRUE)

        # put it back together, calculate AUC
        samp_dat <- rbind(trt_samp, pbo_samp)
        AUC_output <- integrate_survdat(dat=samp_dat, tau=tau)

        # calculate AUC
        boot_df$rmst_trt[i] <- AUC_output$rmst_trt
        boot_df$rmst_pbo[i] <- AUC_output$rmst_pbo
        boot_df$rmst_diff[i] <- AUC_output$rmst_trt - AUC_output$rmst_pbo
    }

    # find CI for RMST
    rmst_trt_CI <- quantile(boot_df$rmst_trt, probs=c(alpha/2, 1-alpha/2))
    rmst_pbo_CI <- quantile(boot_df$rmst_pbo, probs=c(alpha/2, 1-alpha/2))
    rmst_diff_CI <- quantile(boot_df$rmst_diff, probs=c(alpha/2, 1-alpha/2))

    # RMST for observed data
    AUC_output <- integrate_survdat(dat=dat, tau=tau)
    rmst_trt <- AUC_output$rmst_trt
    rmst_pbo <- AUC_output$rmst_pbo
    rmst_diff <- rmst_trt - rmst_pbo

    # format RMST output into df
    rmst_df <- data.frame( cbind(c(rmst_trt, rmst_pbo, rmst_diff),
                                 rbind(rmst_trt_CI, rmst_pbo_CI, rmst_diff_CI)) )
    colnames(rmst_df) <- c('Est', alpha/2, 1-alpha/2)
    rownames(rmst_df) <- c('Arm1', 'Arm0', 'Arm1-Arm0')

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

            # put it back together, calculate AUC
            samp_dat <- rbind(trt_samp, pbo_samp)
            AUC_output <- integrate_survdat(dat=samp_dat, tau=tau)

            # calculate AUC
            null_boot_df$rmst_trt[i] <- AUC_output$rmst_trt
            null_boot_df$rmst_pbo[i] <- AUC_output$rmst_pbo
            null_boot_df$rmst_diff[i] <- AUC_output$rmst_trt - AUC_output$rmst_pbo
        }
        pval_df <- data.frame(rmst_diff_pside = length(which(null_boot_df$rmst_diff > rmst_diff)) / num_boots,
                              rmst_diff_nside = length(which(null_boot_df$rmst_diff < rmst_diff)) / num_boots,
                              rmst_diff_2side = length(which(abs(null_boot_df$rmst_diff) > abs(rmst_diff))) / num_boots)
    } else {
        pval_df <- NULL
    }

    return( list(rmst_df=rmst_df, pval_df=pval_df))
}
