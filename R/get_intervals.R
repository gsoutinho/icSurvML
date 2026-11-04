#' Constructs Turnbull intervals from the observed censoring intervals
#' 
#' @description
#' Constructs Turnbull intervals from the observed censoring intervals
#' and estimates the corresponding event probabilities using the
#' nonparametric maximum likelihood estimator (NPMLE). The estimation
#' is performed iteratively using a self-consistency algorithm, which
#' updates the probability assigned to each Turnbull interval based on
#' the observed interval-censoring information. The resulting estimates
#' define a stepwise survival function. The function also supports
#' stratified analyses and bootstrap estimation of uncertainty.
#'
#' @param left Numeric vector containing the left endpoints of the
#'   observation intervals.
#' @param right Numeric vector containing the right endpoints of the
#'   observation intervals. Use \code{Inf} or \code{NA} for
#'   right-censored observations.
#' @param group Optional numeric, character, or factor vector defining
#'   the groups for stratified analysis. If provided, Turnbull
#'   intervals are computed separately for each group. Defaults to
#'   \code{NULL}.
#' @param nboot Integer specifying the number of bootstrap
#'   replications. Defaults to \code{1}, which corresponds to
#'   estimation without bootstrapping.
#'
#' @return An object consisting of a list with
#'   the following components:
#'   \itemize{
#'     \item \code{original}: A list of data frames, one per group,
#'     containing the estimated Turnbull intervals, probability
#'     weights, and survival estimates for the original dataset.
#'     \item \code{bootstrap}: A list of bootstrap replications if
#'     \code{nboot > 1}; otherwise, \code{NULL}.
#'   }
#'
#' @examples
#' 
#' data(brca)
#' head(brca)
#' dim(brca)
#' 
#' set.seed(123)
#' ints <- get_intervals(brca$left, brca$right2, group = NULL, nboot = 5)
#' ints$original
#' ints$bootstrapp
#' 
#' @importFrom interval icfit
#' @importFrom survival Surv

get_intervals<- function(left, right, group = NULL, nboot = 1) {
  
  if (length(left) != length(right)) stop("Length of 'left' and 'right' must be
                                          different.")
  
  
  if (is.null(group)) {
    group <- rep("All", length(left))
  }
  
  unique_groups <- unique(group)
  
  process_all_groups <- function(l, r, g) {
    lapply(unique_groups, function(grp) {
      idx <- which(g == grp)
      res <- .process_single_run2(l[idx], r[idx])
      list(group = grp, data = res)
    })
  }
  
  original_result <- process_all_groups(left, right, group)
  
  bootstrap_results <- NULL
  if (nboot > 1) {
    bootstrap_results <- replicate(nboot, {
      indices <- sample(seq_along(left), replace = TRUE)
      process_all_groups(left[indices], right[indices], group[indices])
    }, simplify = FALSE)
  }
  
  result <- list(original = original_result, bootstrap = bootstrap_results)
  class(result) <- c("list", "TB")
  return(invisible(result))
}

.process_single_run2 <- function(left, right) {
  
  left[is.infinite(left)] <- NA
  right[is.infinite(right)] <- NA
  
  
  res <- interval::icfit(survival::Surv(left, right, type = "interval2") ~ 1)
  
  sink(tempfile())
  summary_res <- summary(res)
  sink()
  
  ob <- data.frame(
    left = as.numeric(gsub("\\(([^,]+),.*", "\\1", summary_res$Interval)),
    right = suppressWarnings(as.numeric(gsub(".*,(.*)\\]", "\\1", 
                                             summary_res$Interval))),
    weight = as.numeric(summary_res$Probability),
    survival = 1 - cumsum(as.numeric(summary_res$Probability))
  )
  ob[is.na(ob)] <- Inf
  
  
  if (sum(ob$weight) > 1) {
    ex <- sum(ob$weight) - 1
    len <- nrow(ob)
    ob$weight[len] <- ob$weight[len] - ex
    ob$survival[len] <- ob$survival[len] + ex
  }
  
  return(ob)
}