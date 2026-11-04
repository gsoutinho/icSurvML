#' Fit Survival Curves on Imputed Data
#' 
#' @description Fit Survival Curves on Imputed Data
#' @param times Numeric vector containing the survival times.
#' @param status Numeric vector containing the event indicator
#' (1 = event, 0 = censored)
#' @return A survfit object.
#' 
#' @examples
#' \dontrun{
#' data(breastCancer)
#' head(breastCancer)
#' dim(breastCancer)
#' 
#' breastCancer[1:6, c("rectime", "censrec", "survtime", "censdead", "left_rec", 
#'                     "right_rec", "left", "right", "age", "menopause")]
#' 
#' imputed_mid <- impute_ic(Surv(left, right) ~ 1, 
#'                data = breastCancer, 
#'                method = "midpoint")
#' head(imputed_mid[, c("left", "right", "time_imp", "status_imp")])
#' 
#' fit_pm<- fit_ic(imputed_mid$time_imp, imputed_mid$status_imp)
#' summary(fit_pm, times = c(365, 730, 1095, 1460))
#' 
#' imputed_lnorm<- impute_ic(Surv(left, right) ~ age + menopause + 
#'                 hormone + size + grade + nodes + prog_recp + 
#'                 estrg_recp, data = breastCancer, method = "lognormal", 
#'                 n_rep = 1000)
#' head(imputed_lnorm[, c("left", "right", "time_imp", "status_imp")])
#' fit_lnorm<- fit_ic(imputed_lnorm$time_imp, imputed_lnorm$status_imp)
#' summary(fit_lnorm, times = c(365, 730, 1095, 1460))
#' }

fit_ic <- function(times, status) {
  
  formula<-Surv(times, status) ~ 1 
  
  fit <- survival::survfit(formula) 
  
  return(fit)
}