#' Reconstruct digitized Kaplan-Meier curves and generate invididual patient data
#'
#' Reconstruct individual-level data from augmented survival table and
#' augmented NAR table, with augmentation performed by format_raw_tabs().
#'
#' @param aug_NAR A data frame processed through format_raw_tabs().
#' @param aug_surv A data frame processed through format_raw_tabs().
#'
#' @return A list including IPD_time, IPD_event, n_hat=n_hat,
#' KM_hat, n_cen, n_event, int_censor
#'
#' @export
#'
#' @examples
#' data(pembro_NAR)
#' data(pembro_clicks)
#' augTabs <- format_raw_tabs(raw_NAR=pembro_NAR, raw_surv=pembro_clicks)
#' KM_reconstruct(aug_NAR=augTabs$aug_NAR, aug_surv=augTabs$aug_surv)
#'
KM_reconstruct <- function(aug_NAR, aug_surv) {

    # info from NAR table
    TAR <- aug_NAR$time
    NAR <- aug_NAR$NAR
    lower <- aug_NAR$lower
    upper <- aug_NAR$upper

    # make sure the time/survival is nonincreasing
    t_surv <- aug_surv$time
    surv <- aug_surv$surv
    if ( is.unsorted(t_surv) | is.unsorted(rev(surv)) ) {stop('aug_surv unsorted')}

    # number of intervals
    total_ints <- length(NAR)
    # number of event times
    total_e_times <- upper[total_ints]

    # number censored on each interval (except last)
    int_censor <- rep(0, total_ints-1)
    # last value of t where we had an event
    last_event <- rep(1, total_ints)

    # estimated subjects remaining at each k
    n_hat <- rep(NAR[1]+1, total_e_times)
    # number censored at each k
    n_cen <- rep(0, total_e_times)
    # number of events at each k
    n_event <- rep(0, total_e_times)
    # S(t) at each k
    KM_hat <- rep(1, total_e_times)

    # loop through intervals
    for (int_idx in 1:(total_ints-1)) {

        # it's possible that surv[lower[int_idx]] = 0 if the KMC goes to 0
        if (surv[lower[int_idx]] == 0) {
            int_censor[int_idx] <- 0
        } else {
            # first approximation of no. censored on interval int_idx
            int_censor[int_idx] <- round(NAR[int_idx] * surv[lower[int_idx+1]] /
                                             surv[lower[int_idx]] - NAR[int_idx+1])
        }

        # adjust int_censor[int_idx] until n_hat = NAR at the start of the next interval
        # if we have too many events, then just add more censoring
        # if too few events and no. censored > 0, then remove censoring
        # if too few events and no.censored <=0, then stuck, just move on
        while ( n_hat[lower[int_idx+1]] > NAR[int_idx+1] |
                (n_hat[lower[int_idx+1]] < NAR[int_idx+1]&&int_censor[int_idx]>0) ) {

            # can't have negative censoring
            if (int_censor[int_idx] <= 0) {
                n_cen[lower[int_idx]:upper[int_idx]] <- 0
                int_censor[int_idx] <- 0
            } else {

                # evenly distribute censoring times
                cen_times <- t_surv[lower[int_idx]] + (1:int_censor[int_idx])*
                    (t_surv[lower[int_idx+1]] - t_surv[lower[int_idx]]) / (int_censor[int_idx]+1)
                n_cen[lower[int_idx]:upper[int_idx]] <- hist(cen_times,
                    breaks=t_surv[lower[int_idx]:lower[int_idx+1]], plot=F)$counts
            }

            # now account for all events in the interval
            n_hat[lower[int_idx]] <- NAR[int_idx]
            last <- last_event[int_idx]
            for (click_idx in lower[int_idx]:upper[int_idx]) {
                # initial row
                if (click_idx == 1) {
                    n_event[click_idx] <- 0
                    KM_hat[click_idx] <- 1
                } else {
                    # have to check if our KMC goes to zero
                    if (KM_hat[last] == 0) {
                        n_event[click_idx] <- 0
                        KM_hat[click_idx] <- 0
                    } else {
                        # KM_hat and S are ideally the same, but since we are rounding/estimating,
                        # there will be small differences
                        n_event[click_idx] <- round(n_hat[click_idx] * (1-(surv[click_idx] / KM_hat[last])))
                        KM_hat[click_idx] <- KM_hat[last] * (1-(n_event[click_idx] / n_hat[click_idx]))
                    }
                }

                # fill in next n_hat
                n_hat[click_idx+1] <- n_hat[click_idx] - n_event[click_idx] - n_cen[click_idx]
                # update last
                if (n_event[click_idx] != 0) {last <- click_idx}
            }

            # update amount of censoring we need
            int_censor[int_idx] <- int_censor[int_idx] + (n_hat[lower[int_idx+1]] - NAR[int_idx+1])
        } # end while loop through one interval

        # if ended the interval with fewer estimated at risk than published number,
        # that means there int_censor[int_idx] was not positive, so nobody else to redistribute,
        # so we need to change the number at risk from the published number to our estimated
        # number before continuing the estimation
        if (n_hat[lower[int_idx+1]] < NAR[int_idx+1]) {NAR[int_idx+1] <- n_hat[lower[int_idx+1]]}

        last_event[int_idx+1] <- last
    } # end looping through intervals

    # record events and event times
    IPD_event <- rep(0, NAR[1])
    IPD_event[1:sum(n_event)] <- 1
    IPD_time <- c()
    for (click_idx in 1:total_e_times) {
        IPD_time <- c(IPD_time, rep(t_surv[click_idx], n_event[click_idx]))
    }

    # record censoring times
    for (click_idx in 1:(total_e_times-1)) {
        IPD_time <- c(IPD_time, rep((t_surv[click_idx]+t_surv[click_idx+1])/2, n_cen[click_idx]))
    }

    # fill rest of events as censored at max(t_surv)
    ppl_remain <- length(IPD_event) - length(IPD_time)
    if (ppl_remain < 0) {
        stop('Algorithm failed, ended up with too many people')
    } else {
        IPD_time <- c(IPD_time, rep(max(t_surv), ppl_remain))
    }

    # return separated results
    return( list(IPD_time=IPD_time, IPD_event=IPD_event, n_hat=n_hat,
           KM_hat=KM_hat, n_cen=n_cen, n_event=n_event, int_censor=int_censor) )
}
