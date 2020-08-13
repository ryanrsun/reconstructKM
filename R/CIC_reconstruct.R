#' In competing risks situations, papers may provide one overall KM plot for 
#' the composite outcome of event 1 or event 2 as well as cumulative incidence
#' plots for the each event separately. We can use these three plots to reconstruct
#' individual level data with event-specific labels (censored, event 1, or event 2).
#' Can also handle the case when the CIC for event 2 is not given.
#'
#' @param overallIPD The individual patient data from the overall (composite outcome) plot
#' that has already been processed through reconstructKM. Should have three columns: time, 
#' status, and arm.
#' @param clicks1 A data.frame with "time" and "cuminc" columns that are output from the
#' digitizing software, similar to what you would use for reconstructKM except it's a 
#' cumulative incidence function, not a survival function (make sure first click is (0,0)).
#' @param arm The arm corresponding to clicks1 and possibly clicks2.
#' @param clicks2 Same as clicks1 but for the second event. Default is null.
#'
#' @return An augmented version of overallIPD that additionally gives the cause
#' of the event (cause 1 or cause 2) as a fourth "event" column.
#'
#' @export
#'
#' @examples
#' data(TTfields_pfs_trt_NAR)
#' data(TTfields_pfs_trt_clicks)
#' augmented_NAR <- format_NAR_tab(rawNAR=TTfields_pfs_trt_NAR, rawSurv=TTfields_pfs_trt_clicks)
#' reconstruct <- KM_reconstruct(aug_NAR=augmented_NAR$aug_NAR, aug_surv=augmented_NAR$aug_surv)
#' IPD <- data.frame(arm=1, time=reconstruct$IPD_time, status=reconstruct$IPD_event)
#' CIC_reconstruct(overallIPD = IPD, clicks1 = TTfields_pfs_trt_clicks, arm=1, clicks2=NULL)
CIC_reconstruct <- function(overallIPD, clicks1, arm, clicks2=NULL) {

  # make sure the time/survival is nonincreasing
  if (is.unsorted(clicks1$time) | is.unsorted(rev(clicks2$surv)) ) {stop('clicks1 unsorted')}
  if (!is.null(clicks2)) {
    if (is.unsorted(clicks2$time) | is.unsorted(rev(clicks2$surv)) ) {stop('clicks2 unsorted')}
  }
  
  # get standard KM outputs from reconstructed IPD.
  # only times where an event occured, corresponding to jumps in
  # cumulative incidence.
  overallIPD <- overallIPD %>% filter(arm == arm)
  num_sub <- nrow(overallIPD)
  kmDat <- overallIPD %>% group_by(time) %>%
    summarise(events = sum(status), eventsPlusCens=n()) %>%
    mutate(cens = eventsPlusCens - events) %>%
    mutate(atRisk = num_sub - cumsum(eventsPlusCens) + eventsPlusCens) %>%
    mutate(haz = events / atRisk) %>%
    mutate(cumhaz = cumsum(haz)) %>%
    mutate(surv = cumprod(1-haz)) %>%
    filter(events > 0)

  # we need to make sure that clicks1 and clicks2 have a number of clicks equal
  # to the number of rows in km_dat plus 1 (km_dat won't have the initial t=0 point).
  ntimesKM <- nrow(kmDat)
  nclicks1 <- nrow(clicks1)
  if (nclicks1 < ntimesKM + 1) {
	  nAdd <- ntimesKM + 1 - nclicks1
	  cat("Need to add ", nAdd, " time points to clicks1")
	  clicks1 <- add_clicks(clicksDF = clicks1, targetTimes = unlist(kmDat$time), nAdd = nAdd)
  } else if (nclicks1 > ntimesKM + 1) {
	  nRemove <- nclicks1 - ntimesKM - 1
    cat("Need to remove ", nRemove, " time points to clicks1")
    clicks1 <- remove_clicks(clicksDF = clicks1, targetTimes = unlist(kmDat$time), nRemove = nRemove)
  }

  # same with clicks2
  if (!is.null(clicks2)) {
	  nclicks2 <- nrow(clicks2)
    if (nclicks2 < ntimesKM + 1) {
      nAdd <- ntimesKM + 1 - nclicks2
      cat("Need to add ", nAdd, " time points to clicks2")
      clicks2 <- add_clicks(clicksDF = clicks2, targetTimes = unlist(kmDat$time), nAdd = nAdd)
	  } else if (nclicks2 > ntimesKM + 1) {
      nRemove <- nclicks2 - ntimesKM - 1
      cat("Need to remove ", nRemove, " time points to clicks2")
      clicks2 <- remove_clicks(clicksDF = clicks2, targetTimes = unlist(kmDat$time), nRemove = nRemove)      
	  }
  }

  # sync up the times
  clicks1$time <- c(0, kmDat$time)
  if (is.null(clicks2)) {
	  clicksAug <- clicks1 %>% mutate(time = c(0, kmDat$time)) %>%
	    merge(., kmDat, by="time", all.x=TRUE) %>%
	    arrange(time) %>%
	    mutate(prevSurv = c(1, 1, surv[2:(nrow(.) - 1)])) %>%
	    mutate(f1Hat = c(0, diff(cuminc))) %>%
	    mutate(lambda1Hat = f1Hat / prevSurv) %>%
	    mutate(d1Hat = lambda1Hat * atRisk) %>%
	    mutate(d1Hat = ifelse(d1Hat > events, events, d1Hat)) %>%
	    mutate(d2Hat = events - d1Hat)    
  } else {

    # if two clicks provided
    clicksAug <- data.frame(time=c(0, kmDat$time), cuminc1 = clicks1$cuminc, cuminc2 = clicks2$cuminc) %>%
	    merge(., kmDat, by="time", all.x=TRUE) %>%
      arrange(time) %>%
      mutate(prevSurv = c(1, 1, surv[2:(nrow(.) - 1)])) %>%
      mutate(f1Hat = c(0, diff(cuminc1))) %>%
	    mutate(f2Hat = c(0, diff(cuminc2))) %>%
      mutate(coinFlip = runif(n=nrow(.))) %>%
      # flip a coin if no movement in either cuminc but there are events 
      # (could happen if needed to add a row at this time in both clicks)
      mutate(f1Hat = ifelse(f1Hat == 0 & f2Hat == 0 & coinFlip < 0.5, events, f1Hat)) %>%
      mutate(f2Hat = ifelse(f1Hat == 0 & f2Hat == 0 & coinFlip >= 0.5, events, f2Hat)) %>%
	    mutate(percent1 = f1Hat / (f1Hat + f2Hat)) %>%
	    mutate(d1Hat = round(events * percent1)) %>%
	    mutate(d2Hat = events - d1Hat)
  }
    
  # now go insert the event type into overallIPD
  eventVec <- rep(0, nrow(overallIPD)) 
  for (time_it in 2:nrow(clicksAug)) {
	  tempTime <- clicksAug$time[time_it]
	  nReplace1 <- clicksAug$d1Hat[time_it]
	  nReplace2 <- clicksAug$d2Hat[time_it]
	  replaceIdx <- which(overallIPD$time == tempTime & overallIPD$status == 1)
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


