#' Compute Survival Estimates and Confidence Intervals
#'
#' @description
#' Estimates survival probabilities at specified time points using
#' previously computed Turnbull intervals. Supports stratified
#' analyses and percentile-based bootstrap confidence intervals.
#'
#' @param data Object returned by \code{get_intervals}, containing
#'   the Turnbull estimates and, when available, bootstrap
#'   replications.
#' @param times Numeric vector specifying the time points at which
#'   survival probabilities are estimated.
#' @param conf Logical. If \code{TRUE}, computes confidence intervals
#'   using the bootstrap replications stored in \code{data}.
#'   Defaults to \code{FALSE}.
#' @param conf.level Numeric value between 0 and 1 specifying the
#'   confidence level. Defaults to \code{0.95}.
#'
#' @return An object of class \code{flexIC}, implemented as a data
#'   frame with the following columns:
#'   \itemize{
#'     \item \code{time}: Time points at which survival probabilities
#'     are estimated.
#'     \item \code{survival}: Estimated survival probabilities.
#'     \item \code{group}: Group identifier associated with each
#'     estimate.
#'     \item \code{ci_lower}: Lower confidence limit.
#'     \item \code{ci_upper}: Upper confidence limit.
#'   }
#' @examples
#' 
#' data(brca)
#' head(brca)
#' dim(brca)
#' 
#' set.seed(123)
#' ints <- get_intervals(brca$left, brca$right2, group = NULL, nboot = 5)
#' 
#' times <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
#' tb_surv <- survival_ic(ints, times = times, conf = TRUE)
#' tb_surv 
#' 
#' @importFrom graphics abline axis legend lines matplot par plot polygon
#'
#' @importFrom stats approxfun as.formula binomial delete.response glm
#' @importFrom stats model.frame model.matrix model.offset model.response
#' @importFrom stats na.omit pchisq pnorm predict quantile rnorm
#' @importFrom stats rpois sd terms weighted.mean
#'
#' @importFrom survival Surv coxph strata survfit untangle.specials
#'
#' @importFrom utils capture.output flush.console head
#'
#' @importFrom mstate LMAJ cutLMms events msdata2etm msfit
#' @importFrom mstate msprep probtrans to.trans2 trans2Q
#' @importFrom mstate transMat xsect
#'
#' @importFrom stats median runif
#' @importFrom survidm tprob
#' 
#' @importFrom xgboost xgb.train xgb.DMatrix setinfo
#' 
#' @importFrom randomForestSRC rfsrc
#' 
#' @importFrom mboost glmboost gamboost bbs bols CoxPH boost_control survFit
#' 
#' @importFrom survivalsvm survivalsvm
#'  
#' @export get_intervals
#' @export survival_ic
#' @export test_ic
#' @export plot_ic
#' @export impute_ic
#' @export fit_ic
#' @export cond_surv_ic
#' @export cond_aZb_surv_ic
#' @export tp_ic
#' 
survival_ic<- function(data, times, conf = FALSE, conf.level = 0.95) {
  
  group_names <- sapply(data$original, function(x) x$group)
  
  all_results <- lapply(seq_along(data$original), function(i) {
    grp_name <- data$original[[i]]$group
    grp_data <- data$original[[i]]$data
    
    surv_est <- .single.flexIC_survival(grp_data, times)
    
    res <- data.frame(time = times, survival = surv_est, group = grp_name)
    
    if (conf) {
      if (is.null(data$bootstrap)) stop("Bootstrap data not found.")
      
      boot_samples <- lapply(data$bootstrap, function(b) b[[i]]$data)
      
      boot_matrix <- matrix(NA, nrow = length(times), ncol = length(boot_samples))
      for (b in seq_along(boot_samples)) {
        boot_matrix[, b] <- .single.flexIC_survival(boot_samples[[b]], times)
      }
      
      alpha <- (1 - conf.level) / 2
      res$ci_lower <- apply(boot_matrix, 1, quantile, probs = alpha, na.rm = TRUE)
      res$ci_upper <- apply(boot_matrix, 1, quantile, probs = 1 - alpha, na.rm = TRUE)
    }
    return(res)
  })
  
  final_df <- do.call(rbind, all_results)
  class(final_df) <- c("data.frame", "flexIC")
  return(final_df)
}


.single.flexIC_survival <- function(data, times) {
  data <- data[order(data$left), ]
  
  min_t <- min(data$left, na.rm = TRUE)
  data_extended <- rbind(data.frame(left = -Inf, right = min_t, weight = 0, survival = 1), data)
  
  survival_estimates <- numeric(length(times))
  
  for (i in seq_along(times)) {
    t <- times[i]
    
    if (t <= min_t) {
      survival_estimates[i] <- 1
    } else if (t >= max(data$right, na.rm = TRUE)) {
      survival_estimates[i] <- data$survival[nrow(data)]
    } else {
      
      idx <- which(data$left <= t & data$right >= t)
      
      if (length(idx) > 0) {
        idx <- idx[1] 
        s_prev <- if(idx == 1) 1 else data$survival[idx-1]
        t_start <- data$left[idx]
        t_end <- data$right[idx]
        
        if (t_end == t_start || is.infinite(t_end)) {
          survival_estimates[i] <- data$survival[idx]
        } else {
          survival_estimates[i] <- s_prev + (data$survival[idx] - s_prev) * (t - t_start) / (t_end - t_start)
        }
      } else {
        last_idx <- max(which(data$left <= t))
        survival_estimates[i] <- data$survival[last_idx]
      }
    }
  }
  return(survival_estimates)
}