#' Calculate RMST for each arm as well as contrast
#'
#' Non-parametric RMST function that allows for
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
#' @return A list including data.frame of results in each arm (RMST, RMTL, SE, pvalue, CI)
#' as well as data.frame of results for Arm1 - Arm0 RMST.
#'
#' @export
#'
#' @examples
#' time <- rnorm(100)
#' status <- rbinom(n=100, size=1, prob=0.5)
#' arm <- c( rep(1, 50), rep(0, 50))
#' dat <- data.frame(time=time, status=status, arm=arm)
#' nonparam_rmst(dat=dat, tau=1, alpha=0.05)
#'
nonparam_rmst <- function(dat, tau, alpha = 0.05) {

    # make sure data has necessary arms
    if ( length(which(c('time', 'status', 'arm') %in% colnames(dat))) != 3) {
        stop('Column names must include time, status, arm exactly.')
    }
    uniqueArms <- sort(unique(dat$arm))
    if ( length(uniqueArms) > 2) {
        stop("More than two arms, please run this function with only two arms at a time.")
    }

    # loop through arms
    outList <- list()
    for (arm_it in 1:length(uniqueArms)) {
        tempDat <- dat %>% dplyr::filter(.data$arm == uniqueArms[arm_it])
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
    contrastDF <- data.frame(Contrast = "Arm1 - Arm0", Est = rmstDiff10, se = rmstDiff10se,
                             Lower = rmstDiff10low, Upper = rmstDiff10up, pval = rmstDiff10pval)

    return( list(oneArmList = outList, contrastDF = contrastDF))
}
