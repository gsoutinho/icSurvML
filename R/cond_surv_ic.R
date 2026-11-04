#' Estimate a Conditional Survival Function
#'
#' @description
#' Estimates the conditional survival function at a specified covariate
#' value using either the conditional Turnbull estimator or the
#' midpoint-based Beran estimator.
#'
#' @param L Numeric vector of left endpoints of the observed censoring
#' intervals.
#' @param R Numeric vector of right endpoints of the observed censoring
#' intervals.
#' @param Z Numeric vector containing the continuous covariate used to
#' estimate the conditional survival function.
#' @param z0 Numeric value specifying the covariate value at which the
#' conditional survival function is estimated.
#' @param times Numeric vector of time points at which the conditional
#' survival function is evaluated.
#' @param method Character string specifying the estimation method:
#' \code{"Turnbull"} for the conditional Turnbull estimator or
#' \code{"Beran"} for the midpoint-based Beran estimator.
#' Defaults to \code{"Turnbull"}.
#' @param h Optional numeric value specifying the bandwidth used for
#' kernel weighting by the Beran estimator. Defaults to \code{NULL}.
#' @param kernel Kernel function. Defaults to
#' \code{epanechnikov_kernel}. Its use depends on the implementation
#' of the selected estimation method.
#' @param tol Numeric convergence tolerance. Defaults to \code{1e-7}.
#' @param max_iter Integer specifying the maximum number of iterations.
#' Defaults to \code{10000}.
#' @param init Initialization method for the estimation algorithm.
#' Defaults to \code{"uniform"}.
#' @param verbose Logical value indicating whether progress or diagnostic
#' information should be displayed. Defaults to \code{FALSE}.
#'
#' @return For \code{method = "Turnbull"}, a list containing
#' \code{survival}, the estimated conditional survival probabilities,
#' and \code{fit}, the fitted model returned by the conditional
#' Turnbull estimator. For \code{method = "Beran"}, the result returned
#' by \code{beran_midpoint_estimator()}.

  cond_surv_ic<- function(L, R, Z, z0, times, method='Turnbull', h = NULL,
                          kernel = epanechnikov_kernel,
                          tol = 1e-7,
                          max_iter = 10000,
                          init = "uniform",
                          verbose = FALSE) {
    
    if(method=='Turnbull'){
      
      results<- conditional_turnbull_survival_maximal(
        L = L,
        R = R,
        Z = Z,
        z0 = z0,
        times = times)
      
      output <- list(survival = results$survival, fit=results$fit)
      
      return(output)
  
    }
    
    if(method=='Beran'){
      beran_midpoint_estimator(
             L=L, R=R, X=Z, x0=z0, times=times, h = h)  
  
      }
  }

#-------------------------------------------------------------------------------
#conditional_turnbull_survival_maximal:
#-------------------------------------------------------------------------------

make_turnbull_intervals_maximal <- function(L, R) {
  validate_interval_data(L, R)
  finite_endpoints <- sort(unique(c(L, R[is.finite(R)])))
  if (length(finite_endpoints) < 2L) {
    stop("Not enough finite endpoints to construct Turnbull intervals.")
  }
  
  elementary <- data.frame(
    a = finite_endpoints[-length(finite_endpoints)],
    b = finite_endpoints[-1L]
  )
  if (any(is.infinite(R))) {
    elementary <- rbind(
      elementary,
      data.frame(a = max(finite_endpoints), b = Inf)
    )
  }
  
  incidence <- make_alpha_matrix(L, R, elementary)
  occupied <- colSums(incidence) > 0
  elementary <- elementary[occupied, , drop = FALSE]
  incidence <- incidence[, occupied, drop = FALSE]
  if (!nrow(elementary)) {
    stop("No compatible elementary intervals were constructed.")
  }
  
  # Merge connected elementary cells with exactly the same incidence vector.
  merged <- vector("list", nrow(elementary))
  merged_incidence <- vector("list", nrow(elementary))
  group <- 1L
  merged[[group]] <- elementary[1L, , drop = FALSE]
  merged_incidence[[group]] <- incidence[, 1L]
  
  if (nrow(elementary) > 1L) {
    for (j in 2:nrow(elementary)) {
      adjacent <- identical(merged[[group]]$b, elementary$a[j])
      same_incidence <- identical(merged_incidence[[group]], incidence[, j])
      if (adjacent && same_incidence) {
        merged[[group]]$b <- elementary$b[j]
      } else {
        group <- group + 1L
        merged[[group]] <- elementary[j, , drop = FALSE]
        merged_incidence[[group]] <- incidence[, j]
      }
    }
  }
  
  merged <- do.call(rbind, merged[seq_len(group)])
  incidence <- do.call(cbind, merged_incidence[seq_len(group)])
  if (group == 1L) {
    incidence <- matrix(incidence, ncol = 1L)
  }
  
  dominated <- logical(group)
  if (group > 1L) {
    for (j in seq_len(group)) {
      for (k in seq_len(group)) {
        if (j != k && all(incidence[, j] <= incidence[, k]) &&
            any(incidence[, j] < incidence[, k])) {
          dominated[j] <- TRUE
          break
        }
      }
    }
  }
  
  intervals <- merged[!dominated, , drop = FALSE]
  rownames(intervals) <- NULL
  intervals
}


conditional_turnbull_survival_maximal <- function(L, R, Z, z0, times,
                                                  h = NULL,
                                                  kernel = epanechnikov_kernel,
                                                  tol = 1e-7,
                                                  max_iter = 10000,
                                                  init = "uniform",
                                                  verbose = FALSE) {
  if (is.null(h)) {
    h <- rule_of_thumb_bandwidth(Z)
  }
  fit <- weighted_turnbull_em_core(
    L, R, Z, z0, h,
    support_function = make_turnbull_intervals_maximal,
    kernel = kernel, tol = tol, max_iter = max_iter,
    init = init, verbose = verbose
  )
  list(fit = fit, survival = predict_conditional_survival(fit, times))
}


#---------------buscar funcoes:

## Shared utilities for conditional Turnbull estimators.

epanechnikov_kernel <- function(x) {
  out <- 0.75 * (1 - x^2)
  out[abs(x) > 1] <- 0
  out
}


make_alpha_matrix <- function(L, R, intervals) {
  alpha <- matrix(0, nrow = length(L), ncol = nrow(intervals))
  for (i in seq_along(L)) {
    for (j in seq_len(nrow(intervals))) {
      alpha[i, j] <- as.numeric(
        L[i] <= intervals$a[j] && intervals$b[j] <= R[i]
      )
    }
  }
  alpha
}


validate_interval_data <- function(L, R) {
  if (length(L) != length(R)) {
    stop("L and R must have the same length.")
  }
  if (length(L) == 0L || anyNA(L) || anyNA(R)) {
    stop("L and R must be non-empty and must not contain missing values.")
  }
  if (any(!is.finite(L))) {
    stop("All left endpoints L must be finite.")
  }
  if (any(R <= L)) {
    stop("All intervals must satisfy L < R.")
  }
  invisible(NULL)
}


predict_conditional_survival <- function(fit, times) {
  survival <- vapply(times, function(t) {
    if (is.infinite(t)) {
      return(0)
    }
    1 - sum(fit$mass[fit$intervals$b <= t])
  }, numeric(1))
  
  data.frame(time = times, survival = survival)
}


rule_of_thumb_bandwidth <- function(x) {
  scale_x <- min(
    stats::sd(x, na.rm = TRUE),
    stats::IQR(x, na.rm = TRUE) / 1.349
  )
  h <- 1.06 * scale_x * length(x)^(-1 / 5)
  if (!is.finite(h) || h <= 0) {
    stop("The automatic bandwidth is not positive; supply h explicitly.")
  }
  h
}


#-----------buscar funcoes:


## Conditional Turnbull estimator using the original endpoint-cell support.
## Source R/turnbull_utils.R before this file.

make_turnbull_intervals_original <- function(L, R) {
  validate_interval_data(L, R)
  finite_endpoints <- sort(unique(c(L, R[is.finite(R)])))
  if (length(finite_endpoints) < 2L) {
    stop("Not enough finite endpoints to construct Turnbull intervals.")
  }
  
  intervals <- data.frame(
    a = finite_endpoints[-length(finite_endpoints)],
    b = finite_endpoints[-1L]
  )
  compatible <- vapply(seq_len(nrow(intervals)), function(j) {
    any(L <= intervals$a[j] & intervals$b[j] <= R)
  }, logical(1))
  intervals <- intervals[compatible, , drop = FALSE]
  
  if (any(is.infinite(R))) {
    last_left <- max(finite_endpoints)
    if (any(L <= last_left & is.infinite(R))) {
      intervals <- rbind(intervals, data.frame(a = last_left, b = Inf))
    }
  }
  rownames(intervals) <- NULL
  intervals
}


weighted_turnbull_em_core <- function(L, R, Z, z0, h, support_function,
                                      kernel = epanechnikov_kernel,
                                      tol = 1e-8, max_iter = 5000,
                                      init = c("uniform", "weighted"),
                                      verbose = FALSE) {
  init <- match.arg(init)
  if (length(L) != length(R) || length(L) != length(Z)) {
    stop("L, R and Z must have the same length.")
  }
  if (!is.finite(h) || h <= 0) {
    stop("Bandwidth h must be finite and positive.")
  }
  
  intervals <- support_function(L, R)
  alpha <- make_alpha_matrix(L, R, intervals)
  if (any(rowSums(alpha) == 0)) {
    warning("Some subjects have no compatible support interval.")
  }
  
  w <- kernel((Z - z0) / h)
  if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
    stop("Kernel weights must be finite, non-negative, and have positive sum.")
  }
  
  m <- nrow(intervals)
  if (init == "uniform") {
    mass <- rep(1 / m, m)
  } else {
    mass <- colSums(alpha * w)
    mass <- if (sum(mass) > 0) mass / sum(mass) else rep(1 / m, m)
  }
  
  loglik_trace <- numeric(max_iter)
  converged <- FALSE
  for (iteration in seq_len(max_iter)) {
    denominator <- as.vector(alpha %*% mass)
    denominator[denominator <= 0] <- NA_real_
    posterior <- alpha * matrix(mass, nrow = length(L), ncol = m, byrow = TRUE)
    posterior <- posterior / denominator
    posterior[is.na(posterior)] <- 0
    
    mass_new <- colSums(posterior * w) / sum(w)
    if (sum(mass_new) > 0) {
      mass_new <- mass_new / sum(mass_new)
    }
    
    fitted_probability <- as.vector(alpha %*% mass_new)
    positive <- fitted_probability > 0
    loglik_trace[iteration] <- sum(
      w[positive] * log(fitted_probability[positive])
    )
    difference <- max(abs(mass_new - mass))
    mass <- mass_new
    
    if (verbose && iteration %% 100L == 0L) {
      cat("Iteration:", iteration, "difference:", difference, "\n")
    }
    if (difference < tol) {
      converged <- TRUE
      break
    }
  }
  
  if (!converged) {
    warning("Maximum number of iterations reached before convergence.")
  }
  
  structure(
    list(
      intervals = intervals,
      mass = mass,
      weights = w,
      alpha = alpha,
      z0 = z0,
      h = h,
      iterations = iteration,
      converged = converged,
      final_convergence_criterion = difference,
      max_iteration_reached = iteration >= max_iter,
      loglik = loglik_trace[seq_len(iteration)]
    ),
    class = "conditional_turnbull"
  )
}


conditional_turnbull_survival_original <- function(L, R, Z, z0, times,
                                                   h = NULL,
                                                   kernel = epanechnikov_kernel,
                                                   tol = 1e-8,
                                                   max_iter = 5000,
                                                   init = "uniform",
                                                   verbose = FALSE) {
  if (is.null(h)) {
    h <- rule_of_thumb_bandwidth(Z)
  }
  fit <- weighted_turnbull_em_core(
    L, R, Z, z0, h,
    support_function = make_turnbull_intervals_original,
    kernel = kernel, tol = tol, max_iter = max_iter,
    init = init, verbose = verbose
  )
  list(fit = fit, survival = predict_conditional_survival(fit, times))
}

#-------------------------------------------------------------------------------
#beran estimator
#-------------------------------------------------------------------------------


beran_midpoint_estimator <- function(L, R, X, x0, times, h,
                                     kernel = epanechnikov_kernel) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for the midpoint Beran estimator.")
  }
  if (length(L) != length(R) || length(L) != length(X)) {
    stop("L, R and X must have the same length.")
  }
  if (!is.finite(h) || h <= 0) {
    stop("Bandwidth h must be finite and positive.")
  }
  
  pseudo_time <- ifelse(is.finite(R), (L + R) / 2, L)
  event <- as.numeric(is.finite(R))
  weights <- kernel((X - x0) / h)
  if (sum(weights) <= 0) {
    stop("All kernel weights are zero.")
  }
  
  fit <- survival::survfit(
    survival::Surv(pseudo_time, event) ~ 1,
    weights = weights
  )
  as.numeric(summary(fit, times = times, extend = TRUE)$surv)
}


