#' Plot Survival Curves
#'
#' @description
#' Plots estimated survival curves. Confidence intervals can be displayed when
#' available.
#'
#' @param x An object containing estimated
#'   survival probabilities and, when available, confidence intervals.
#' @param main Character string specifying the plot title.
#'   Defaults to \code{""}.
#' @param xlab Character string specifying the x-axis label.
#'   Defaults to \code{"Time"}.
#' @param ylab Character string specifying the y-axis label.
#'   Defaults to \code{"Survival"}.
#' @param ylim Numeric vector of length two specifying the y-axis
#'   limits. Defaults to \code{c(0, 1)}.
#' @param conf Logical. If \code{TRUE}, displays confidence intervals
#'   when available. Defaults to \code{TRUE}.
#' @param showlegend Logical. If \code{TRUE}, displays the group
#'   legend. Defaults to \code{TRUE}.
#' @param ... Additional graphical parameters passed to the plotting
#'   function.
#'
#' @examples
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
#' plot_ic(tb_surv, conf = TRUE)
#' 
#' ints2 <- get_intervals(brca$left, brca$right2, group = brca$treatment, 
#' nboot = 5)
#' 
#' tb_surv2 <- survival_ic(ints2, times = times, conf = TRUE)
#' 
#' plot_ic(tb_surv2, ylim=c(0, 1), conf = TRUE, showlegend = TRUE)
#'
#' @return A survival plot displayed in the active graphics device.
#'
#' @importFrom graphics plot lines polygon grid legend

plot_ic<- function(x, main = "",
                         xlab = "Time", ylab = "Survival",  ylim = c(0, 1),
                         conf = TRUE, showlegend = TRUE, ...) {
  
  if (!inherits(x, "flexIC")) {
    stop("The object must be of class 'flexIC'.")
  }
  
  return(.plot_surv_static(x = x, main = main, xlab = xlab, ylab = ylab, 
                           ylim = ylim, showlegend = showlegend, 
                           conf = conf, ...))
}


.plot_surv_static <- function(x, main, xlab, ylab, showlegend, ylim, 
                              conf = TRUE, ...) {
  groups <- unique(x$group)
  cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
  group_cols <- cols[1:length(groups)]
  names(group_cols) <- groups
  
  plot(NULL, xlim = range(x$time), ylim = ylim, 
       main = main, xlab = xlab, ylab = ylab, ...)
  grid()
  
  for (g in groups) {
    df_g <- x[x$group == g, ]
    color <- group_cols[g]
    
    if (conf && "ci_lower" %in% colnames(df_g)) {
      polygon(c(df_g$time, rev(df_g$time)), 
              c(df_g$ci_lower, rev(df_g$ci_upper)),
              col = paste0(color, "33"), border = NA)
    }
    
    lines(df_g$time, df_g$survival, col = color, lwd = 2)
  }
  
  if (showlegend) {
    legend("topright", legend = groups, col = group_cols, lty = 1, lwd = 2, 
           bty = "n")
  }
}
