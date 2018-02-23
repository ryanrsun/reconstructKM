#' Just a wrapper to get quantities out of a call to coxph()
#'
#' @param cox_fit A model fitted with coxph()
#'
#' @return A list including IPD_time, IPD_event, n_hat=n_hat,
#' KM_hat, n_cen, n_event, int_censor
#'
#' @export
#'
#' @examples
#' time <- rnorm(100)
#' status <- rbinom(n=100, size=0.5)
#' arm <- c(rep(1,50), rep(0,50))
#' temp_cox <- coxph(Surv(time, status) ~ arm)
#' print_cox_outputs(temp_cox)
#'
print_cox_outputs <- function(cox_fit, print_output=TRUE) {
    cox_beta <- summary(cox_fit)$coefficients[1,1]
    cox_HR <- summary(cox_fit)$coefficients[1,2]
    cox_beta_SE <- summary(cox_fit)$coefficients[1,3]
    cat('Hazard Ratio: ', cox_HR, '\n')
    cat('95% CI: ',
        c(exp(cox_beta - 1.96*cox_beta_SE), exp(cox_beta + 1.96*cox_beta_SE)), '\n')
    cat('P-value: ',
        summary(cox_fit)$coefficients[1,5], '\n')

    return( list(beta=cox_beta, HR=cox_HR, se=cox_beta_SE,
                 CI=c(exp(cox_beta - 1.96*cox_beta_SE), exp(cox_beta + 1.96*cox_beta_SE))))
}
