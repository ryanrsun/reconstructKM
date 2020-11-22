#' Reconstruct cumulative incidence curves
#'
#' In competing risks situations, papers may provide one overall KM plot for
#' the composite outcome of event 1 or event 2 as well as cumulative incidence
#' plots for the each event separately. We can use these three plots to reconstruct
#' individual level data with event-specific labels (censored, event 1, or event 2).
#' Can also handle the case when the CIC for event 2 is not given.
#' Run this separately for each arm.
#'
#' @param overallIPD The individual patient data from the overall (composite outcome) plot
#' that has already been processed through reconstructKM. Should have three columns: time,
#' status, and arm.
#' @param clicks1 A data.frame with "time" and "cuminc" columns that are output from the
#' digitizing software, similar to what you would input for reconstructKM except it's a
#' cumulative incidence function for a specific event, not a survival function (make sure first click is (0,0)).
#' @param arm The arm corresponding to clicks1 and possibly clicks2.
#' @param clicks2 Same as clicks1 but for the second event if it's provided. Default is null.
#'
#' @return An augmented version of overallIPD that additionally gives the cause
#' of the event (cause 1 or cause 2) as a fourth "event" column.
#'
#' @export
#'
#' @examples
#' data(pembro_clicks)
#' data(pembro_NAR)
#' augTabs <- format_raw_tabs(raw_NAR=pembro_NAR, raw_surv=pembro_clicks)
#' reconstruct <- KM_reconstruct(aug_NAR=augTabs$aug_NAR, aug_surv=augTabs$aug_surv)
#' IPD <- data.frame(arm=1, time=reconstruct$IPD_time, status=reconstruct$IPD_event)
#' clicks1 <- dplyr::mutate(pembro_clicks, cuminc=1-survival)
#' CIC_reconstruct(overallIPD = IPD, clicks1 = clicks1, arm=1, clicks2=NULL)
CIC_reconstruct <- function(overallIPD, clicks1, arm, clicks2=NULL) {

  # make sure the time/survival is nondecreasng
  if (is.unsorted(clicks1$time) | is.unsorted(rev(clicks2$surv)) ) {stop('clicks1 is not nondecreasing.')}
  if (!is.null(clicks2)) {
    if (is.unsorted(clicks2$time) | is.unsorted(rev(clicks2$surv)) ) {stop('clicks2 is not nondecreasing.')}
  }

  # get standard KM outputs from reconstructed IPD.
  # only times where an event occured, corresponding to jumps in
  # cumulative incidence.
  overallIPD <- overallIPD %>% filter(arm == arm)
  num_sub <- nrow(overallIPD)
  kmDat <- overallIPD %>% group_by(time) %>%
    summarise(events = sum(.data$status), eventsPlusCens=n()) %>%
    # number censored at each time point
    mutate(cens = .data$eventsPlusCens - .data$events) %>%
    # number at risk at each time point
    mutate(atRisk = num_sub - cumsum(.data$eventsPlusCens) + .data$eventsPlusCens) %>%
    # estimated hazard at each time point
    mutate(haz = .data$events / .data$atRisk) %>%
    # estimated cumulative hazard and survival
    mutate(cumhaz = cumsum(.data$haz)) %>%
    mutate(surv = cumprod(1-.data$haz)) %>%
    filter(.data$events > 0)

  # we need to make sure that both clicks1 and clicks2 have a number of clicks equal
  # to the number of rows in km_dat plus 1 (km_dat won't have the initial t=0 point).
  ntimesKM <- nrow(kmDat)
  nclicks1 <- nrow(clicks1)
  # if we need to add clicks - usually it's this one because overall includes all clicks from both events,
  # so it'll have more clicks than either event alone.
  if (nclicks1 < ntimesKM + 1) {
	  nAdd <- ntimesKM + 1 - nclicks1
	  # add clcks here
	  clicks1 <- add_clicks(clicksDF = clicks1, targetTimes = unlist(kmDat$time), nAdd = nAdd)
  } else if (nclicks1 > ntimesKM + 1) {
	  nRemove <- nclicks1 - ntimesKM - 1
    # or remove clicks here
    clicks1 <- remove_clicks(clicksDF = clicks1, targetTimes = unlist(kmDat$time), nRemove = nRemove)
  }

  # same procedure with clicks2
  if (!is.null(clicks2)) {
	  nclicks2 <- nrow(clicks2)
    if (nclicks2 < ntimesKM + 1) {
      nAdd <- ntimesKM + 1 - nclicks2
      # add to clicks2
      clicks2 <- add_clicks(clicksDF = clicks2, targetTimes = unlist(kmDat$time), nAdd = nAdd)
	  } else if (nclicks2 > ntimesKM + 1) {
      nRemove <- nclicks2 - ntimesKM - 1
      # remove from clicks2
      clicks2 <- remove_clicks(clicksDF = clicks2, targetTimes = unlist(kmDat$time), nRemove = nRemove)
	  }
  }

  # merge event clicks times with overall click times, figure out how
  # many of each event at each time point.
  clicks1$time <- c(0, kmDat$time)
  if (is.null(clicks2)) {
	  clicksAug <- clicks1 %>% mutate(time = c(0, kmDat$time))
	    # note we just assume the two sets of clicks both correspond to the same times, specifically
	    # the times from the overall clicks. that is, the 10th click from clicks1 corresponds to the same
	    # time as the 10th click from overall reconstructed data. this might not hold depending on where
	    # add_clicks() or remove_clicks() operated.
	  clicksAug <- merge(clicksAug, kmDat, by="time", all.x=TRUE) %>%
	    arrange(time) %>%
	    mutate(prevSurv = c(1, 1, .data$surv[2:(length(.data$surv) - 1)])) %>%
	    # difference in cuminc from clicks tells us how many events are of type 1
	    mutate(f1Hat = c(0, diff(.data$cuminc))) %>%
	    mutate(lambda1Hat = .data$f1Hat / .data$prevSurv) %>%
	    mutate(d1Hat = round(.data$lambda1Hat * .data$atRisk)) %>%
	    mutate(d1Hat = ifelse(.data$d1Hat > .data$events, .data$events, .data$d1Hat)) %>%
	    # all other events are of type 2
	    mutate(d2Hat = .data$events - .data$d1Hat)
  } else {

    # if two clicks provided.
    # again just assume the overall reconstructed times are the correct ones
    clicksAug <- data.frame(time=c(0, kmDat$time), cuminc1 = clicks1$cuminc, cuminc2 = clicks2$cuminc)
    clicksAug <- merge(clicksAug, kmDat, by="time", all.x=TRUE) %>%
      arrange(time) %>%
      mutate(prevSurv = c(1, 1, .data$surv[2:(length(.data$surv) - 1)])) %>%
      mutate(f1Hat = c(0, diff(.data$cuminc1))) %>%
	    mutate(f2Hat = c(0, diff(.data$cuminc2))) %>%
      mutate(coinFlip = runif(n=length(.data$f2Hat))) %>%
      # flip a coin if no movement in either cuminc but there are events
      # (could happen if needed to add a row at this time in both clicks)
      mutate(f1Hat = ifelse(.data$f1Hat == 0 & .data$f2Hat == 0 & .data$coinFlip < 0.5, .data$events, .data$f1Hat)) %>%
      mutate(f2Hat = ifelse(.data$f1Hat == 0 & .data$f2Hat == 0 & .data$coinFlip >= 0.5, .data$events, .data$f2Hat)) %>%
      # look at relative differences in cuminc jumps to decide where events should be allocated - no randomness here.
	    mutate(percent1 = .data$f1Hat / (.data$f1Hat + .data$f2Hat)) %>%
	    mutate(d1Hat = round(.data$events * .data$percent1)) %>%
	    mutate(d2Hat = .data$events - .data$d1Hat)
  }

  # now go insert the event type into overallIPD
  eventVec <- rep(0, nrow(overallIPD))
  for (time_it in 2:nrow(clicksAug)) {
	  tempTime <- clicksAug$time[time_it]
	  nReplace1 <- clicksAug$d1Hat[time_it]
	  nReplace2 <- clicksAug$d2Hat[time_it]
	  replaceIdx <- which(overallIPD$time == tempTime & overallIPD$status == 1)
	  # add type 1 or 2
	  if (nReplace1 > 0) {
	    replace1idx <- replaceIdx[1:nReplace1]
	    eventVec[replace1idx] <- 1
	  }
	  if (nReplace2 > 0) {
	    replace2idx <- replaceIdx[(nReplace1+1):(nReplace1+nReplace2)]
	    eventVec[replace2idx] <- 2
	  }
  }

  # return
  overallIPD <- overallIPD %>% mutate(event = eventVec)
  return(overallIPD)
}


