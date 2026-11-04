#' Estimate a Conditional Survival Function Given an Interval Constraint
#'
#' @description
#' Estimates the conditional survival function of a terminal event
#' given that the intermediate event time lies within the interval
#' \eqn{(a,b]}. The estimation method can be selected between the
#' reconstruction-based Turnbull estimator and the midpoint-based
#' Beran estimator.
#'
#' @param L Numeric vector of left endpoints of the observed censoring
#' intervals for the intermediate event time.
#' @param R Numeric vector of right endpoints of the observed censoring
#' intervals for the intermediate event time.
#' @param T Numeric vector of observed terminal event times, which may
#' be right-censored.
#' @param Delta Numeric or logical vector indicating whether the
#' terminal event was observed (1 or \code{TRUE}) or right-censored
#' (0 or \code{FALSE}).
#' @param a Numeric value specifying the lower bound of the conditioning
#' interval for the intermediate event time.
#' @param b Numeric value specifying the upper bound of the conditioning
#' interval for the intermediate event time.
#' @param times Numeric vector of time points at which the conditional
#' survival function is evaluated.
#' @param method Character string specifying the estimation method:
#' \code{"Turnbull"} for the reconstruction-based estimator or
#' \code{"Beran"} for the midpoint-based estimator. Defaults to
#' \code{"Turnbull"}.
#'
#' @return For \code{method = "Turnbull"}, a list containing
#' \code{survival}, the estimated conditional survival probabilities.
#' For \code{method = "Beran"}, returns the result of
#' \code{midpoint_conditional_survival()}.

cond_aZb_surv_ic<- function(L, R, T, Delta, a, b, times, method='Turnbull'){
  
  if(method=='Turnbull'){
    
      results<- conditional_survival_aZb_corrected(  
      LZ = L,
      RZ = R,
      Tobs = T,
      Delta = Delta,
      a = a,
      b = b,
      times = times,
      tol = 1e-7,
      max_iter = 10000,
      init = "uniform",
      verbose = FALSE
    )
      
      output <- list(survival = results$survival)
      
      return(output)
    
  }
  
  
  if(method=='Beran'){
    
    midpoint_conditional_survival(
      LZ = L,
      RZ = R,
      Tobs = T,
      Delta = Delta,
      a = a,
      b = b,
      times = times
    )
    
  }
  
}


turnbull_marginal_em_maximal <- function(L, R, tol = 5e-7,
                                         max_iter = 10000,
                                         init = c("uniform", "compatible"),
                                         verbose = FALSE) {
  init <- match.arg(init)
  validate_interval_data(L, R)
  if (!is.finite(tol) || tol <= 0 || length(tol) != 1L) {
    stop("tol must be a finite positive scalar.")
  }
  if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1 ||
      max_iter != as.integer(max_iter)) {
    stop("max_iter must be a positive integer.")
  }
  
  intervals <- make_turnbull_intervals_maximal(L, R)
  alpha <- make_alpha_matrix(L, R, intervals)
  if (any(rowSums(alpha) == 0)) {
    stop("At least one observation has no compatible maximal intersection.")
  }
  
  n <- length(L)
  m <- nrow(intervals)
  mass <- if (init == "uniform") {
    rep(1 / m, m)
  } else {
    initial <- colSums(alpha)
    initial / sum(initial)
  }
  
  loglik_trace <- numeric(max_iter)
  converged <- FALSE
  final_criterion <- Inf
  
  for (iteration in seq_len(max_iter)) {
    observed_probability <- as.vector(alpha %*% mass)
    if (any(observed_probability <= 0)) {
      stop("A non-positive observed-interval probability occurred in the EM algorithm.")
    }
    
    posterior <- alpha * matrix(mass, nrow = n, ncol = m, byrow = TRUE)
    posterior <- posterior / observed_probability
    mass_new <- colSums(posterior) / n
    mass_new <- mass_new / sum(mass_new)
    
    fitted_probability <- as.vector(alpha %*% mass_new)
    loglik_trace[iteration] <- sum(log(fitted_probability))
    final_criterion <- max(abs(mass_new - mass))
    mass <- mass_new
    
    if (verbose && iteration %% 100L == 0L) {
      cat("Iteration:", iteration,
          "criterion:", final_criterion,
          "log-likelihood:", loglik_trace[iteration], "\n")
    }
    if (final_criterion < tol) {
      converged <- TRUE
      break
    }
  }
  
  max_iteration_reached <- iteration >= max_iter
  if (!converged) {
    warning("Maximum number of iterations reached before convergence.")
  }
  
  structure(
    list(
      intervals = intervals,
      support_points = intervals$b,
      mass = mass,
      alpha = alpha,
      iterations = iteration,
      converged = converged,
      final_convergence_criterion = final_criterion,
      max_iteration_reached = max_iteration_reached,
      tolerance = tol,
      max_iterations = as.integer(max_iter),
      loglik = loglik_trace[seq_len(iteration)]
    ),
    class = "turnbull_marginal_maximal"
  )
}


turnbull_probability_right_endpoint <- function(fit, lower, upper) {
  if (length(lower) != 1L || length(upper) != 1L || is.na(lower) ||
      is.na(upper) || lower >= upper) {
    return(0)
  }
  support <- fit$support_points
  sum(fit$mass[support > lower & support <= upper])
}


compute_pi_aZb_right_endpoint <- function(fit_Z, LZ, RZ, a, b) {
  support <- fit_Z$support_points
  mass <- fit_Z$mass
  
  vapply(seq_along(LZ), function(i) {
    denominator <- sum(mass[support > LZ[i] & support <= RZ[i]])
    if (!is.finite(denominator) || denominator <= 0) {
      return(0)
    }
    lower <- max(a, LZ[i])
    upper <- min(b, RZ[i])
    numerator <- if (lower < upper) {
      sum(mass[support > lower & support <= upper])
    } else {
      0
    }
    value <- numerator / denominator
    min(max(value, 0), 1)
  }, numeric(1))
}


reverse_km_G_minus <- function(Tobs, Delta) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }
  n <- length(Tobs)
  if (length(Delta) != n || n == 0L || anyNA(Tobs) || anyNA(Delta)) {
    stop("Tobs and Delta must be non-empty, equal-length, and complete.")
  }
  if (any(!is.finite(Tobs)) || any(Tobs < 0)) {
    stop("Tobs must contain finite non-negative times.")
  }
  if (any(!Delta %in% c(0, 1))) {
    stop("Delta must contain only 0 and 1.")
  }
  if (all(Delta == 1)) {
    return(rep(1, n))
  }
  
  fit <- survival::survfit(survival::Surv(Tobs, 1 - Delta) ~ 1)
  vapply(Tobs, function(t) {
    before <- which(fit$time < t)
    if (!length(before)) 1 else fit$surv[max(before)]
  }, numeric(1))
}


validate_conditional_survival_aZb_inputs <- function(LZ, RZ, Tobs, Delta,
                                                     a, b, times) {
  validate_interval_data(LZ, RZ)
  n <- length(LZ)
  if (length(Tobs) != n || length(Delta) != n) {
    stop("LZ, RZ, Tobs and Delta must have the same length.")
  }
  if (length(a) != 1L || length(b) != 1L || !is.finite(a) ||
      !is.finite(b) || a >= b) {
    stop("a and b must be finite scalars satisfying a < b.")
  }
  if (!length(times) || anyNA(times) || any(!is.finite(times)) ||
      any(times < 0)) {
    stop("times must be finite, non-negative, and non-empty.")
  }
  if (anyNA(Tobs) || any(!is.finite(Tobs)) || any(Tobs < 0)) {
    stop("Tobs must contain finite non-negative times.")
  }
  if (anyNA(Delta) || any(!Delta %in% c(0, 1))) {
    stop("Delta must contain only 0 and 1.")
  }
  invisible(NULL)
}


conditional_survival_aZb_corrected <- function(LZ, RZ, Tobs, Delta,
                                               a, b, times,
                                               tol = 5e-7,
                                               max_iter = 10000,
                                               init = "uniform",
                                               verbose = FALSE) {
  validate_conditional_survival_aZb_inputs(
    LZ, RZ, Tobs, Delta, a, b, times
  )
  
  fit_Z <- turnbull_marginal_em_maximal(
    LZ, RZ, tol = tol, max_iter = max_iter,
    init = init, verbose = verbose
  )
  denominator <- turnbull_probability_right_endpoint(fit_Z, a, b)
  if (!is.finite(denominator) || denominator <= 0) {
    stop("The estimated probability P(a < Z <= b) is not positive.")
  }
  
  pi_hat <- compute_pi_aZb_right_endpoint(fit_Z, LZ, RZ, a, b)
  G_minus <- reverse_km_G_minus(Tobs, Delta)
  invalid_G <- Delta == 1 & G_minus <= 0
  if (any(invalid_G)) {
    stop("G(Tobs-) is zero for at least one observed terminal event.")
  }
  failure_weight <- numeric(length(Tobs))
  observed_event <- Delta == 1
  failure_weight[observed_event] <-
    1 / (length(Tobs) * G_minus[observed_event])
  
  weighted_cdf <- vapply(times, function(t) {
    sum(failure_weight * as.numeric(Tobs <= t) * pi_hat)
  }, numeric(1))
  survival_raw <- 1 - weighted_cdf / denominator
  survival <- pmin(pmax(survival_raw, 0), 1)
  
  structure(
    list(
      survival = data.frame(
        time = times,
        survival = survival,
        survival_raw = survival_raw,
        weighted_failure_probability = weighted_cdf,
        denominator = denominator
      ),
      fit_Z = fit_Z,
      pi_hat = pi_hat,
      G_minus = G_minus,
      failure_weights = failure_weight,
      convergence = list(
        converged = fit_Z$converged,
        iterations = fit_Z$iterations,
        final_convergence_criterion = fit_Z$final_convergence_criterion,
        max_iteration_reached = fit_Z$max_iteration_reached
      ),
      a = a,
      b = b,
      probability_convention = "right_endpoint_discrete_mass"
    ),
    class = "conditional_survival_aZb_corrected"
  )
}


## Manuscript numerator retained only for controlled methodological comparison.
## It deliberately reuses the same maximal Turnbull fit, D, pi, and G as the
## corrected estimator so that only the survival representation differs.
conditional_survival_aZb_manuscript_from_corrected_fit <- function(fit,
                                                                   Tobs,
                                                                   Delta,
                                                                   times) {
  numerator <- vapply(times, function(t) {
    sum(fit$failure_weights * as.numeric(Tobs > t) * fit$pi_hat)
  }, numeric(1))
  survival_raw <- numerator / fit$survival$denominator[1]
  data.frame(
    time = times,
    survival = pmin(pmax(survival_raw, 0), 1),
    survival_raw = survival_raw,
    numerator = numerator,
    denominator = fit$survival$denominator[1]
  )
}


midpoint_conditional_survival <- function(LZ, RZ, Tobs, Delta,
                                          a, b, times) {  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }
  midpoint <- ifelse(is.finite(RZ), (LZ + RZ) / 2, NA_real_)
  selected <- !is.na(midpoint) & midpoint > a & midpoint <= b
  if (sum(selected) < 5L) {
    return(rep(NA_real_, length(times)))
  }
  fit <- survival::survfit(
    survival::Surv(Tobs[selected], Delta[selected]) ~ 1
  )
  as.numeric(summary(fit, times = times, extend = TRUE)$surv)
}

