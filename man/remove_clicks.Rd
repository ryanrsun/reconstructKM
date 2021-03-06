% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_clicks.R
\name{remove_clicks}
\alias{remove_clicks}
\title{Remove clicks from subdistribution curves for reconstructing CIC}
\usage{
remove_clicks(clicksDF, targetTimes, nRemove)
}
\arguments{
\item{clicksDF}{A data frame with the two columns time and cuminc.}

\item{targetTimes}{A vector of times from the composite KM plot.}

\item{nRemove}{Number of times to remove from clicksDF.}
}
\value{
A clicksDF with fewer rows.
}
\description{
When there are fewer clicks in the composite (overall) outcome curve,
we need to remove them from the subdistribution curves. Find the time points
in the subdistribution data that are furthest away from the composite curve times,
remove those times.
}
\examples{
clicksDF <- data.frame(time=0:10, cuminc=seq(from=0, to=1, by=0.1))
remove_clicks(clicksDF, targetTimes = runif(n=7, min=0, max=10), nRemove=3)
}
