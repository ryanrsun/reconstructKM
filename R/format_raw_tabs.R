#' Format raw survival and NAR tables so they are ready for reconstruction algorithm
#'
#' Augment a raw number at risk table with the necessary information to run
#' the reconstruction algorithm.
#'
#' @param raw_NAR A data frame with the columns 'time' and NAR' at least.
#' @param raw_surv A data frame with the columns 'time' and 'survival' at least.
#' @param tau End of follow-up time, defaults to last time in NAR table.
#'
#' @return A list with aug_NAR and aug_surv, properly cleaned tables that can be used as input in KM_reconstruct().
#'
#' @export
#'
#' @examples
#' data(pembro_clicks)
#' data(pembro_NAR)
#' augTabs <- format_raw_tabs(raw_NAR=pembro_NAR, raw_surv=pembro_clicks)
#'
format_raw_tabs <- function(raw_NAR, raw_surv, tau=NULL) {

    # check clicks tab has correct columns
    has_col <- length(which(colnames(raw_surv) %in% c('time', 'survival')))
    if (has_col != 2) { stop('raw_surv must have columns named time and survival exactly') }

    # subset and order clicks
    raw_surv <- dplyr::select(raw_surv, .data$time, .data$survival) %>%
        dplyr::arrange(time)

    # make sure survival is non-increasing, starts with (0,1), ends with an event
    if (is.unsorted(rev(raw_surv$survival))) {
        stop('survival must be non-increasing in time') }
    if (raw_surv[1, 1] != 0 | raw_surv[1,2] != 1) {
        stop('Your raw_clicks_tab did not have a t=0,S=1 row, check your work') }
    last_click_row <- nrow(raw_surv)
    last_click_t <- raw_surv$time[last_click_row]
    last_surv <- raw_surv$survival[last_click_row]
    if (last_click_t <= raw_surv$time[last_click_row-1] |
        last_surv >= raw_surv$survival[last_click_row-1]) {
        stop('Your last click should have been at the end of a vertical (not
             horizontal) segment')
    }

    # check NAR tab has correct columns
    has_col <- length(which(colnames(raw_NAR) %in% c('time', 'NAR')))
    if (has_col != 2) { stop('raw_NAR must have columns named time and NAR exactly') }

    # subset and order NAR
    raw_NAR <- dplyr::select(raw_NAR, .data$time, .data$NAR) %>%
        dplyr::arrange(time)

    # make sure NAR is non-increasing
    if (is.unsorted(rev(raw_NAR$NAR))) {
        stop('NAR must be non-increasing in time') }

    # follow-up end is the last NAR time unless otherwise specified (e.g. surv goes to 0)
    if (is.null(tau)) {tau=max(raw_NAR$time)}

    # match NAR intervals with raw_clicks rows - remember last row done manually
    ints <- data.frame(lower=rep(NA, nrow(raw_NAR)-1), upper=NA)
    for (int_idx in 1:nrow(ints)) {
        temp_rows <- which(raw_surv$time >= raw_NAR$time[int_idx] &
                               raw_surv$time < raw_NAR$time[int_idx+1])
        if (length(temp_rows) == 0) {
            next
        } else {
            ints$lower[int_idx] <- min(temp_rows)
            ints$upper[int_idx] <- max(temp_rows)
        }
    }

    # augment NAR, remove NA rows
    aug_NAR <- dplyr::bind_cols(raw_NAR[-nrow(raw_NAR), ], ints) %>%
        dplyr::filter(!is.na(.data$lower))

    # manually add last row to NAR and clicks tables
    last_NAR_row <- data.frame(time=max(raw_NAR$time),
                           NAR=min(raw_NAR$NAR), lower=aug_NAR$upper[nrow(aug_NAR)]+1,
                           upper=aug_NAR$upper[nrow(aug_NAR)]+1)
    last_surv_row <- data.frame(time=tau, survival=last_surv)

    aug_NAR <- dplyr::bind_rows(aug_NAR, last_NAR_row)
    aug_surv <- dplyr::bind_rows(raw_surv, last_surv_row)

    return(list(aug_NAR=aug_NAR, aug_surv=aug_surv))

}
