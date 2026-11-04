#' Log-Rank-Type Tests for Interval-Censored Data
#'
#' @description
#' Compares survival distributions between two or more groups
#' using statistical tests for interval-censored data.
#'
#' @param left Numeric vector containing the left endpoints of the
#'   observation intervals.
#' @param right Numeric vector containing the right endpoints of
#'   the observation intervals. Use \code{Inf} or \code{NA} for
#'   right-censored observations.
#' @param group Numeric, character, or factor vector identifying
#'   the groups to be compared. This argument is required.
#' @param method Character string specifying the type of test.
#'   Supported options are \code{"logrank"} (default) and
#'   \code{"Wilcoxon"}.
#'
#' @return An object of class \code{ictest} containing the results
#'   of the statistical test, including the test statistic and
#'   p-value.
#'
#' @details
#' This function wraps \code{interval::ictest} to facilitate
#' comparisons of survival distributions for interval-censored
#' data within the package framework.
#'
#' @examples
#' 
#' data(brca)
#' head(brca)
#' dim(brca)
#' 
#' set.seed(123)
#' 
#' times <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
#' 
#' ints2 <- get_intervals(brca$left, brca$right2, group = brca$treatment, 
#' nboot = 5)
#' 
#' tb_surv2 <- survival_ic(ints2, times = times, conf = TRUE)
#' 
#' logTest <- test_ic(brca$left, brca$right2, 
#'                   group = brca$treatment,
#'                   method = "logrank")
#'                   
#' logTest$p.value
#' 
#' logTest2 <- test_ic(brca$left, brca$right2, 
#'                   group = brca$treatment,
#'                   method = "wilcoxon")
#'                   
#' logTest2$p.value
#' 
#' @importFrom interval ictest
#' @importFrom survival Surv

test_ic<- function(left, right, group, method = "logrank1") {
  
  
  method<-ifelse(method == "logrank", "logrank1", "wmw")
  
  if (is.null(group)) {
    stop("A 'group' vector must be provided to perform the comparison.")
  }
  
  left[is.infinite(left)] <- NA
  right[is.infinite(right)] <- NA
  
  test_result <- interval::ictest(
    survival::Surv(left, right, type = "interval2") ~ group, 
    scores = method
  )
  
  return(test_result)
}
