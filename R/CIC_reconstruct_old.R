#' Using reconstruct individual-level data from KM_reconstruct and
#' digitized inputs from cumulative incidence curves, assign a cause
#' for each failure in competing risks situations. Only consider competing
#' risks, not semi-competing risks. One or two cumulative incidence curves.
#' One arm at a time.
#'
#' @param dat A data frame with arm, time, event.
#' @param aug_CIC A data frame processed through format_CIC().
#'
#' @return An augmented version of dat that additionally gives the cause
#' of the event (cause 1 or cause 2).
#'
#' @export
#'
#' @examples
#' data(TTfields_pfs_trt_NAR)
#' data(TTfields_pfs_trt_clicks)
#' augmented_NAR <- format_NAR_tab(rawNAR=TTfields_pfs_trt_NAR, rawSurv=TTfields_pfs_trt_clicks)
#' KM_reconstruct(aug_NAR=augmented_NAR$aug_NAR, aug_surv=augmented_NAR$aug_surv)
#'
CIC_reconstruct <- function(IPD, clicks1, clicks2=NULL) {

    # get standard KM outputs from reconstructed IPD.
    # only times where an event occured, corresponding to jumps in
    # cumulative incidence.
    num_sub <- nrow(IPD)
    km_dat <- IPD %>% group_by(time) %>%
        summarise(deaths = sum(status), timeRisk=n()) %>%
        mutate(cens = timeRisk - deaths) %>%
        mutate(atRisk = num_sub - cumsum(timeRisk) + timeRisk) %>%
        mutate(haz = deaths / atRisk) %>%
        mutate(cumhaz = cumsum(haz)) %>%
        mutate(surv = cumprod(1-haz)) %>%
        filter(deaths > 0)

    # clicks1 and clicks2 should have time and est_cum_inc.
    # we need to turn clicks1 and clicks2 into aug_CIC.
    # aug_CIC should have the same time pts as IPD,
    # so if the clicks have more time points, need to shrink them down.
    # if the clicks have less time points, insert some time points.
    ntimes_km <- nrow(km_dat)

    # check if time points line up across nsplits intervals,
    # better to do local matching of times.
    # we'd like to have at least five data points in each split
    nsplits <- max(round(ntimes_km / 7), 5)
    seg_length <- ceiling(ntimes_km / nsplits)
    remainder <- ntimes_km - seg_length * (nsplits - 2)
    km_dat_split <- km_dat %>%
        mutate(split = c(rep(1:(nsplits-2), each=ceiling(ntimes_km / nsplits)),
                         rep((nsplits-1), ceiling(remainder/2)),
                         rep(nsplits, floor(remainder/2))) )

    # time split endpoints
    min_times <- km_dat_split %>% group_by(split) %>%
        summarise(min_time = min(time)) %>% select(min_time) %>% unlist(.) %>% as.numeric(.)
    max_times <- km_dat_split %>% group_by(split) %>%
        summarise(max_time = max(time)) %>% select(max_time) %>% unlist(.) %>% as.numeric(.)
    split_times <- c(0, max_times[1:(length(max_times)-1)] +
                         (min_times[2:length(min_times)] - max_times[1:(length(max_times)-1)]) / 2,
                     max(max_times) * 2)

    # match clicks1
    matched_clicks1 <- c()
    prev_times <- c(0, max_times[1:(length(max_times)-1)])
    prev_cuminc <- 0
    for (split_it in 1:nsplits) {
        # split
        temp_km <- km_dat_split %>% filter(split == split_it)
        temp_clicks <- clicks1 %>% filter(time > split_times[split_it] &
                                              time <= split_times[split_it + 1])
        # match and append
        new_temp_clicks <- match_clicks(clicks_dat=temp_clicks, km_dat=temp_km,
                                        prev_row = c(prev_times[split_it], prev_cuminc))
        matched_clicks1 <- c(matched_clicks1, new_temp_clicks)
        prev_cuminc <- max(new_temp_clicks$est_cum_inc)
    }

    # match clicks2
    if (!is.null(clicks2)) {
        ntimes2 <- nrow(clicks2)
        split_clicks2 <- split(clicks2, rep(1:ceiling(ntimes2 / nsplits),
                                            each=nsplits, length.out=ntimes1))
        matched_clicks2 <- c()
        for (split_it in 1:nsplits) {
            # split
            temp_clicks <- split_clicks2[[split_it]]
            temp_km <- split_km[[split_it]]
            # match and append
            new_temp_clicks <- match_clicks(click_dat=temp_clicks, km_dat=temp_km)
            matched_clicks2 <- c(matched_clicks2, new_temp_clicks)
        }
    }

    # if both curves scanned
    if (!is.null(clicks2)) {
        clicks_CIC <- data.frame(time = matched_clicks1$time,
                                 est_cum_inc1 = matched_clicks1$est_cum_inc,
                                 est_cum_inc2 = matched_clicks2$est_cum_inc)

        # clicks_CIC should have time, est_cum_inc1, est_cum_inc2
        aug_CIC <- clicks_CIC %>%
            mutate(surv_prev = c(1, km_dat$surv[(nrow(km_dat)-1)])) %>%
            mutate(t1_inc = c(est_cum_inc1[1], diff(est_cum_inc1))) %>%
            mutate(t1_haz = t1_inc / surv_prev) %>%
            mutate(t2_inc = c(est_cum_inc2[1], diff(est_cum_inc2))) %>%
            mutate(t2_haz = t2_inc / surv_prev)
    } else {    # just scanned one curve
        clicks_CIC <- matched_clicks1

        # get the hazard as total hazard - t1 hazard
        aug_CIC <- clicks_CIC %>%
            mutate(surv_prev = c(1, km_dat$surv[(nrow(km_dat)-1)])) %>%
            mutate(t1_inc = c(est_cum_inc1[1], diff(est_cum_inc1))) %>%
            mutate(t1_haz = t1_inc / surv_prev) %>%
            mutate(t2_haz = max(km_dat$haz - t1_haz, 0)) %>%
            mutate(t2_inc = t2_haz * surv_prev)
    }

    # aug_CIC should have time, t1_haz, t2_haz, deaths
    # assume the time points are matched with IPD
    imputed <- aug_CIC %>%
        mutate(perc1 <- t1_haz / (t1_haz + t2_haz)) %>%
        mutate(type1 <- round(perc1 * deaths)) %>%
        mutate(type2 <- deaths - type1)

    # loop through times and put in type of event
    IPD_aug <- IPD %>% mutate(type = 0)
    imputed_events <- imputed %>% filter(deaths > 0)
    for (time_it in 1:length(unique(imputed$time))) {
        temp_time <- imputed$time[time_it]
        temp_type <- c(rep(1, imputed$type1[time_it]), rep(2, imputed$type2[time_it]))
        IPD_aug <- IPD_aug %>% mutate(type = ifelse(time == temp_time & event == 1,
                                     temp_type, type))
    }

    return(list(IPD_aug=IPD_aug, imputed=imputed))
}


#' Match clicks from cumulative incidence curve to clicks from KM curve.
#' We need the times from both sets of clicks to match exactly to guess
#' the number of events due to each cause.
#'
#' @param clicks_dat A data frame with time and est_cum_inc
#' @param km_dat A data frame with time (and usually surv but that's not used here)
#' @param prev_row The last row of the previous split (c(0,0) if not splitting and matching)
#'
#' @return new_temp_clicks, the time-matched cum incidence clicks
#'
#' @export
#'
#' @examples
#' data(TTfields_pfs_trt_NAR)
#' data(TTfields_pfs_trt_clicks)
#' augmented_NAR <- format_NAR_tab(rawNAR=TTfields_pfs_trt_NAR, rawSurv=TTfields_pfs_trt_clicks)
#' KM_reconstruct(aug_NAR=augmented_NAR$aug_NAR, aug_surv=augmented_NAR$aug_surv)

match_clicks <- function(km_dat, clicks_dat, prev_row=c(0,0)) {
    # size difference between splits
    temp_clicks_rows <- nrow(clicks_dat)
    temp_km_rows <- nrow(km_dat)
    diff_times <- temp_km_rows - temp_clicks_rows

    # new data
    new_temp_clicks <- data.frame(time = km_dat$time,
                                  est_cum_inc=NA)

    # same number of clicks, perfect just match them up
    if (diff_times == 0) {
        new_temp_clicks$est_cum_inc <- clicks_dat$est_cum_inc
    } else if (diff_times > 0) {           # too few clicks

        clicks_dat <- rbind(prev_row, clicks_dat) %>% set_colnames(c("time", "est_cum_inc"))
        temp_clicks_rows <- nrow(clicks_dat)
        temp_click_idx <- 2
        for (time_it in 1:temp_km_rows) {
            # generally number of km left is more
            temp_km_left <- temp_km_rows - time_it + 1
            temp_clicks_left <- temp_clicks_rows - temp_click_idx + 1

            # forced, same number left
            if (temp_km_left <= temp_clicks_left) {
                new_temp_clicks$est_cum_inc[time_it] <-
                    clicks_dat$est_cum_inc[temp_click_idx]
                temp_click_idx <- temp_click_idx + 1
            } else if (temp_clicks_left <= 0) {
                # no more clicks, fill remaining with the previous cum inc
                new_temp_clicks$est_cum_inc[time_it] <-
                    clicks_dat$est_cum_inc[temp_click_idx - 1]
            } else {
                # wiggle room, can match closer ones
                # is the current km time closer to the prev cum inc time or
                # the next cum inc time
                curr_diff <- abs(new_temp_clicks$time[time_it] -
                                         clicks_dat$time[temp_click_idx])
                next_diff <- abs(new_temp_clicks$time[time_it + 1] -
                                     clicks_dat$time[temp_click_idx])
                if (next_diff < curr_diff) {
                    new_temp_clicks$est_cum_inc[time_it] <-
                        clicks_dat$est_cum_inc[temp_click_idx - 1]
                } else {
                    new_temp_clicks$est_cum_inc[time_it] <-
                        clicks_dat$est_cum_inc[temp_click_idx]
                    temp_click_idx <- temp_click_idx + 1
                }
            }   # done with this time_it
        }   # done with all time_its
    }   else if (diff_times < 0) {        # too many clicks, shrink, shouldn't really happen

        temp_km_idx <- 1
        for (time_it in 1:temp_clicks_rows) {
            # generally number of clicks left is more
            temp_km_left <- temp_km_rows - temp_km_idx + 1
            temp_clicks_left <- temp_clicks_rows - time_it + 1

            # forced, same number left
            if (temp_clicks_left <= temp_km_left) {
                new_temp_clicks$est_cum_inc[temp_km_idx] <-
                    clicks_dat$est_cum_inc[time_it]
                temp_km_idx <- temp_km_idx + 1
            } else if (temp_km_left <= 0) {
                # no more km spots left, break out of loop
                break
            } else {
                # wiggle room, can match closer ones
                # is the current km time closer to the current cum_inc time or the
                # next cum_inc time
                curr_diff <- new_temp_clicks$time[temp_km_idx] - clicks_dat$time[time_it]
                next_diff <- new_temp_clicks$time[temp_km_idx] - clicks_dat$time[time_it+1]
                if (prev_diff < next_diff) {
                    new_temp_clicks$est_cum_inc[temp_km_idx] <-
                        clicks_dat$est_cum_inc[time_it]
                    temp_km_idx <- temp_km_idx + 1
                } else {
                    next
                }
            }   # done with this time_it
        }   # done with all time_its
    }   # done with if (diff_times < 0)

    return(new_temp_clicks)
}


