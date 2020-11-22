#' Print outputs from Cox regression
#'
#' Just a wrapper to get quantities out of a call to coxph()
#'
#' @param cox_fit A model fitted with coxph()
#' @param print_output Print summary to screen if TRUE
#'
#' @return A list including beta, HR, SE, and CI
#'
#' @export
#'
#' @examples
#' time <- rnorm(100)
#' status <- rbinom(n=100, prob=0.5, size=1)
#' arm <- c(rep(1,50), rep(0,50))
#' temp_cox <- survival::coxph(survival::Surv(time, status) ~ arm)
#' print_cox_outputs(temp_cox)
#'
print_cox_outputs <- function(cox_fit, print_output=TRUE) {
    cox_beta <- summary(cox_fit)$coefficients[1,1]
    cox_HR <- summary(cox_fit)$coefficients[1,2]
    cox_beta_SE <- summary(cox_fit)$coefficients[1,3]
    cox_pvalue <- summary(cox_fit)$coefficients[1,5]

    return( list(beta=cox_beta, HR=cox_HR, se=cox_beta_SE, pvalue=cox_pvalue,
                 CI=c(exp(cox_beta - 1.96*cox_beta_SE), exp(cox_beta + 1.96*cox_beta_SE))) )
}
