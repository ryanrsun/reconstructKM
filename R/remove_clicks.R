#' Remove clicks from subdistribution curves for reconstructing CIC
#'
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
