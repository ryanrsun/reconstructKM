#' When there are more clicks in the composite (overall) outcome curve,
#' we need to add them to the subdistribution curves. Find the time points
#' in the composite data that are furthest away from the times in clicksDF,
#' add these times to clicksDF with 0 jumps in cuminc.
#' @param clicksDF A data frame with the two columns, time and cuminc. 
#' @param targetTimes A vector of times from the composite KM plot.
#' @param nAdd Number of times to add to clicksDF.
#'
#' @return An augmented clicksDF with extra rows (no cuminc jumps in those extra times).
#'
#' @export
#'
#' @examples
#' clicksDF <- data.frame(time=0:10, cuminc=seq(from=0, to=1, by=0.1))
#' add_clicks(clicksDF, targetTimes = runif(n=14, min=0, max=10), nAdd=5)

add_clicks <- function(clicksDF, targetTimes, nAdd) {
  
  # add the targetTime that is furthest away from existing times in clicksDF
  for (add_it in 1:nAdd) {
    # minimum distance to clicks times
    minDists <- sapply(X=targetTimes, FUN=function(x, searchVec){min(abs(x - searchVec))},
                       searchVec = unlist(clicksDF$time))
    # maximum min distance to clicks times
    furthestTime <- targetTimes[which.max(minDists)]
    # add this time in
    prevRow <- max(which(clicksDF$time < furthestTime))
    clicksDF <- clicksDF %>% add_row(time=furthestTime, cuminc=clicksDF$cuminc[prevRow]) %>%
      arrange(time)
  }
  
  return(clicksDF)
}



#' When there are fewer clicks in the composite (overall) outcome curve,
#' we need to remove them from the subdistribution curves. Find the time points
#' in the subdistribution data that are furthest away from the composite curve times,
#' remove those times.
#' @param clicksDF A data frame with the two columns time and cuminc. 
#' @param targetTimes A vector of times from the composite KM plot.
#' @param nRemove Number of times to remove from clicksDF.
#'
#' @return A clicksDF with fewer rows.
#'
#' @export
#'
#' @examples
#' clicksDF <- data.frame(time=0:10, cuminc=seq(from=0, to=1, by=0.1))
#' remove_clicks(clicksDF, targetTimes = runif(n=7, min=0, max=10), nRemove=3)

remove_clicks <- function(clicksDF, targetTimes, nRemove) {
  
  # remove the clicksDF time that is closest to another clicksDF time
  for (remove_it in 1:nRemove) {
    # distance to next time
    distanceToNext <- diff(clicksDF$time)
    removeRow <- order(distanceToNext)[1]
    if (removeRow == 1) {removeRow <- 2}
    
    # remove
    clicksDF <- rbind(clicksDF[1:(removeRow-1), ], clicksDF[(removeRow+1):nrow(clicksDF), ] )
  }
  
  return(clicksDF)
}