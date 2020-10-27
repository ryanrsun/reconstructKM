#' Calculate nonparametric RMST for a single arm up to tau for data.frame with time and status

#' @param dat Data frame of time-to-event data which MUST have the columns
#' 'time' and 'status
#' @param tau The cutoff time, a scalar
#' @param alpha Level for confidence interval
#'
#' @return data.frame with rows for RMST and RMTL and columnns for estimate, std err, pvalue, and CI
#' @export
#' @examples
#' time <- rnorm(100)
#' status <- rbinom(n=100, size=0.5)
#' dat <- data.frame(time=time, status=status, tau=2)
#' integrate_survdat(dat=dat)
#'
integrate_survdat <- function(dat, tau, alpha=0.05) {
    #make sure data has necessary arms
    if ( length(which(c('time', 'status') %in% colnames(dat))) != 2) {
        error('Column names must include time and status exactly.')
    }

    # get the KM table using prebuilt functions
    # it only tabulates the unique times
    # disregard those after \tau, then add a row for \tau
    KMfit <- survival::survfit(Surv(time, status) ~ 1, data=dat)
    n <- nrow(dat)
    nUnique <- length(KM_fit$time)
    KMtab <- data.frame(time = KMfit$time, surv = KMfit$surv, nRisk = KMfit$n.risk, nEvent = KMfit$n.event) %>%
        filter(time <= tau) %>%
        add_row(time=tau, surv=KM_fit$surv[nUnique]) %>%
        # use lead() to add the "next time" column
        mutate(nextTime = lead(time)) %>%
        mutate(diffTime = nextTime - time) %>%
        # multiply the survival by the difference in times to get each auc chunk
        mutate(aucChunk = diffTime * surv) %>%
        # filter out the last NA row now that it's served its purpose
        filter(!is.na(aucChunk)) %>%
        # need this part for the variance calculation
        mutate(cumAUC = cumsum(aucChunk)) %>%
        mutate(cumAUC = cumAUC + time[1]) %>%
        mutate(intervalAUC = cumAUC[nrow(.)] - cumAUC + aucChunk)

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

#' My personal non-parametric RMST function that allows for
#' the tau (follow-up time) to be arbitrarily large. Uno package
#' restricts it to be min(last observed event in either arm).
#' Provides estimate, SE, CI for each arm. Provides same for
#' difference in arms (and also p-value).
#'
#' @param dat Data frame of time-to-event data which MUST have the columns
#' 'time', 'arm', and 'status
#' @param tau How long of a follow-up to consider, i.e. we integrate the survival
#' functions from 0 to tau
#' @param alpha Confidence interval is given for (alpha/2, 1-alpha/2) percentiles
#'
#' @return A list including outTab (estimate and CI in both arms), trt_rmst,
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
nonparam_rmst <- function(dat, tau, alpha = 0.05) {

    # make sure data has necessary arms
    if ( length(which(c('time', 'status', 'arm') %in% colnames(dat))) != 3) {
        error('Column names must include time, status, arm exactly.')
    }
    uniqueArms <- sort(unique(dat$arm))
    if ( length(uniqueArms) > 2) {
        print("More than two arms, contrasts will only be between first two numerically sorted arms.")
    }

    # loop through arms
    outList <- list()
    for (arm_it in 1:length(uniqueArms)) {
        tempDat <- dat %>% filter(arm == uniqueArms[arm_it])
        tempIntegrate <- integrate_survdat(dat = tempDat, tau = tau, alpha=alpha)
        outList[[arm_it]] <- tempIntegrate
    }
    # name the list
    names(outList) <- paste0("Arm", uniqueArms)

    # contrasts
    # rmst difference
    rmstDiff10 <- outList[[2]]$Est[1] - outList[[1]]$Est[1]
    rmstDiff10se <- sqrt(outList[[2]]$se[1]^2 + outList[[1]]$se[1]^2)
    rmstDiff10low <- rmstDiff10 - qnorm(1 - alpha/2) * rmstDiff10se
    rmstDiff10up <- rmstDiff10 + qnorm(1 - alpha/2) * rmstDiff10se
    rmstDiff10pval <- 1 - pchisq( (rmstDiff10 / rmstDiff10se)^2, df=1 )
    contrastDF <- data.frame(Contrast = "Arm1 - Arm0", Est = rmstDiff10, se = rmstDiff10,
                             Lower = rmstDiff10low, Upper = rmstDiff10up, pval = rmstDiff10pval)

    return( list(oneArmList = outList, contrastDF = contrastDF))
}
