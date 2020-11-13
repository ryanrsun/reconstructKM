#' Tumor-Treating Fields example PFS KM reconstruction clicks - Temozolomide alone arm
#'
#' A dataset containing the clicks used to reconstruct the Temozolomide alone KM curve.
#'
#' @format A data frame with 114 rows and 2 variables:
#' \describe{
#'   \item{time}{event time, months}
#'   \item{survival}{probability of progression-free survival}
#' }
#' @source \url{JAMA. 2018;318(23):2306-2316}
"PFS_noTTF_clicks"


#' Tumor-Treating Fields example PFS NAR table - Temozolomide alone arm
#'
#' A dataset containing the number at risk information for the Temozolomide alone KM curve.
#'
#' @format A data frame with 6 rows and 2 variables:
#' \describe{
#'   \item{time}{event time, months}
#'   \item{NAR}{number of subjects at risk}
#'   ...
#' }
#' @source \url{JAMA. 2018;318(23):2306-2316}
"PFS_noTTF_NAR"


#' Tumor-Treating Fields example PFS KM reconstruction clicks - TTFields + temozolomide arm
#'
#' A dataset containing the clicks used to reconstruct the TTFields + temozolomide KM curve.
#'
#' @format A data frame with 223 rows and 2 variables:
#' \describe{
#'   \item{time}{event time, months}
#'   \item{survival}{probability of progression-free survival}
#' }
#' @source \url{JAMA. 2018;318(23):2306-2316}
"PFS_TTF_clicks"



#' Tumor-Treating Fields example PFS NAR table - TTFields + temozolomide arm
#'
#' A dataset containing the number at risk information for the TTFields + temozolomide KM curve.
#'
#' @format A data frame with 6 rows and 2 variables:
#' \describe{
#'   \item{time}{event time, months}
#'   \item{NAR}{number of subjects at risk}
#'   ...
#' }
#' @source \url{JAMA. 2018;318(23):2306-2316}
"PFS_TTF_NAR"
