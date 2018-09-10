#' Plot two-way analysis cost-effectiveness plane based on bootstrap data
#'
#' This function plots the cost-efectiveness plane of a two-way analysis of a single study.
#' @param x Object of class study2way.
#' @param y Should be omitted.
#' @param xlab x axis label. Defaults to 'x'.
#' @param ylab y axis label. Defaults to 'y'.
#' @param xlim numeric vector length 2 with limits for the x axis. Defaults to range of x$x.bt.
#' @param ylim numeric vector length 2 with limits for the y axis. Defaults to range of x$y.bt.
#' @param col color for individual studies. Defaults to 1.
#' @param pch integer point character for bootstrap points and confidence area. Defaults to 2.
#' @param lwd integer line width for individual studies confidence area. Defaults to 1.
#' @param points logical whether bootstrap points should be plotted. Defaults to TRUE.
#' @param max.points integer the maximum number of bootstrap points to be plotted per study and meta-analysis. Bootstrap samples 1 to max.points will be plotted if available. Use NA or NULL to plot all points. Defaults to 2000.
#' @param quadrant.proportions logical whether the proportions of bootstrap points in each quadrant should be printed in the plot corners. Defaults to FALSE.
#' @param add logical whether to add the plot to existing plot. Defaults to FALSE.
#' @keywords cost-effectiveness; bootstrapping
#' @export
#' @return
#' Returns NULL invisibly.
#' @examples
#' s1 <- study2way(treatments$Abrahams2018.x, treatments$Abrahams2018.y, 'Abrahams 2018')
#' plot(s1, xlab='Delta effects', ylab='Delta costs')

plot.study2way <- function(x, y, xlab="x", ylab="y", xlim=range(x$x.bt), ylim=range(x$y.bt), col=1, pch=3,
                           center.pch=15, lwd=1, lty=1,
                           points=TRUE, max.points=2000, quadrant.proportions=FALSE, add=FALSE){

  # function to draw ellipse (simplified function inspired by car package)
  dataEllipse <- function (x, y, center.pch = 19, center.cex = 2, col = 1, lwd = 2, lty=1, level=.95, segments=51){
    dfn <- 2
    dfd <- length(x) - 1
    v <- cov.wt(cbind(x, y))
    shape <- v$cov
    center <- v$center
    radius <- sqrt(dfn * qf(level, dfn, dfd))
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    Q <- chol(shape, pivot = TRUE)
    order <- order(attr(Q, "pivot"))
    ellipse <- t(center + radius * t(unit.circle %*% Q[, order]))
    colnames(ellipse) <- c("x", "y")
    lines(ellipse, col = col, lwd = lwd, lty=lty)
    if ((center.pch != FALSE) && (!is.null(center.pch)))
      points(center[1], center[2], pch = center.pch, cex = center.cex, col = col)
  }

  if(!add){
    # create empty plot
    plot(NA, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
    abline(h=0, lty=2)
    abline(v=0, lty=2)
  }

  # determine which rows to add
  r <- 1:x$n.bt
  if(!is.null(max.points)){
    if(!is.na(max.points)){
      if(x$n.bt > max.points) r <- 1:max.points
    }
  }

  # add bootstrap points
  if(points){
    points(x$x.bt[r], x$y.bt[r], col=col, pch=pch, cex=.3)
  }

  # add ellipse
  dataEllipse(x$x.bt, x$y.bt, col=col, lwd=lwd, lty=lty, center.pch=center.pch)

}
