#' Pembrolizumab example OS KM reconstruction clicks - pembrolizumab arm
#'
#' A dataset containing the clicks used to reconstruct the pembrolizumab OS KM curve.
#'
#' @format A data frame with 97 rows and 2 variables:
#' \describe{
#'   \item{time}{event time, months}
#'   \item{survival}{probability of progression-free survival}
#' }
#' @source \url{NEJM. 2018;378(22):2078-2092}
"pembro_clicks"


#' Pembrolizumab example OS NAR table - pembrolizumab arm
#'
#' A dataset containing the number at risk information for the pembrolizumab OS KM curve.
#'
#' @format A data frame with 8 rows and 2 variables:
#' \describe{
#'   \item{time}{event time, months}
#'   \item{NAR}{number of subjects at risk}
#'   ...
#' }
#' @source \url{NEJM. 2018;378(22):2078-2092}
"pembro_NAR"


#' Pembrolizumab example OS KM reconstruction clicks - placebo arm
#'
#' A dataset containing the clicks used to reconstruct the placebo OS KM curve.
#'
#' @format A data frame with 96 rows and 2 variables:
#' \describe{
#'   \item{time}{event time, months}
#'   \item{survival}{probability of progression-free survival}
#' }
#' @source \url{NEJM. 2018;378(22):2078-2092}
"pbo_clicks"



#' Pembrolizumab example OS NAR table - placebo arm
#'
#' A dataset containing the number at risk information for the placebo OS KM curve.
#'
#' @format A data frame with 8 rows and 2 variables:
#' \describe{
#'   \item{time}{event time, months}
#'   \item{NAR}{number of subjects at risk}
#'   ...
#' }
#' @source \url{NEJM. 2018;378(22):2078-2092}
"pbo_NAR"
