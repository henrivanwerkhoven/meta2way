#' Plot two-way meta-analysis cost-effectiveness plane based on bootstrap data
#'
#' This function plots the cost-efectiveness plane of a two-way meta-analysis.
#' @param x Object of class meta2way.
#' @param y Should be omitted.
#' @param do.fixed whether to plot the fixed effects meta-analysis. Defaults to x$do.fixed.
#' @param do.random whether to plot the random effects meta-analysis. Defaults to x$do.random.
#' @param xlab x axis label. Defaults to 'x'.
#' @param ylab y axis label. Defaults to 'y'.
#' @param col color for individual studies. Defaults to 2:(x$n.studies+1).
#' @param pch integer point character for individual studies bootstrap points and confidence area. Length should equal number of included studies. Defaults to 2:(x$n.studies+1).
#' @param lwd integer line width for individual studies confidence area. Defaults to 1.
#' @param points logical whether bootstrap points of individual studies should be plotted. Defaults to TRUE.
#' @param meta.col color for meta-analysis bootstrap points and confidence area. Vector length one or two; if two, the first will be used for the fixed effects and the second for the random effects meta-analysis if performed. Defaults to 1.
#' @param meta.pch integer point character for meta-analysis points: vector length one or two (like meta.col). Defaults to c(15,17) because I like these characters most.
#' @param meta.lwd integer line width for meta-analysis confidence area. Defaults to 2.
#' @param meta.lty integer line type for meta-analysis confidence area. Defaults to c(2,3).
#' @param meta.points logical whether bootstrap points of meta-analyses should be plotted. Defaults to FALSE.
#' @param max.points integer the maximum number of bootstrap points to be plotted per study and meta-analysis. Bootstrap samples 1 to max.points will be plotted if available. Use NA or NULL to plot all points. Defaults to 2000.
#' @param quadrant.proportions logical whether the proportions of meta-analysis bootstrap points in each quadrant should be printed in the plot corners. Defaults to FALSE.
#' @param quadrant.format format for the quadrant proportions. If do.fixed and do.random are both true, two numbers are required in the format. Default depends on do.fixed and do.random.
#' @param legend logical whether a legend of the studies should be printed. The legend will be printed below the plot. Defaults to TRUE.
#' @param legend.ncol integer number of columns of the legend. Defaults to 2.
#' @param legend.cex relative font size of the legend. Defaults to 1.
#' @param legend.I2 logical whether to print the I2 statistics in the legend. Defaults to TRUE.
#' @param legend.tau2 logical whether to print the tau2 statistics in the legend. Defaults to TRUE.
#' @param legend.ci logical whether to print the 95% confidence interval of the I2 or tau2 statistics. Only in effect if legend.I2==TRUE | legend.tau2==TRUE. Defaults to TRUE.
#' @keywords meta-analysis; cost-effectiveness; bootstrapping
#' @export
#' @return
#' Returns NULL invisibly.
#' @examples
#' m <- meta2way(study2way(treatments$Abrahams2018.x, treatments$Abrahams2018.y, 'Abrahams 2018'),
#'               study2way(treatments$Baruch2018.x, treatments$Baruch2018.y, 'Baruch 2018'))
#' plot(m, xlab='delta effects', ylab='delta costs')
plot.meta2way <- function(x, y, do.fixed=x$do.fixed, do.random=x$do.random,
                          xlab="x", ylab="y", col=NULL, pch=NULL, lwd=1, lty=1, points=TRUE,
                          meta.col=1, meta.pch=c(17,15), meta.lwd=2, meta.lty=c(2,3), meta.points=FALSE,
                          label.fixed.estimate='Fixed estimate', label.fixed.conf.area='Fixed conf area',
                          label.random.estimate='Random estimate', label.random.conf.area='Random conf area',
                          max.points=2000, quadrant.proportions=FALSE, quadrant.format=if(do.fixed & do.random) '%.1f%% (FE model)\n%.1f%% (RE model)' else '%.1f%% (meta-analysis)',
                          legend=TRUE, legend.ncol=2, legend.cex=1, legend.I2=TRUE, legend.tau2=TRUE, legend.ci=TRUE){

  # check arguments
  if(!missing(y)) stop("Argument y should not be used.")

  if(is.null(col)) col <- 2:(x$n.studies+1)
  if(is.null(pch)) pch <- 2:(x$n.studies+1)

  # these should be named vectors with fixed and random or else forse it
  if(length(meta.lty)==2) names(meta.lty) <- c('fixed','random')
  else if(length(meta.lty) == 1) meta.lty <- c(fixed=meta.lty, random=meta.lty)
  else stop('Argument meta.lty should be a vector length 1 or 2')

  if(length(meta.col)==2) names(meta.col) <- c('fixed','random')
  else if(length(meta.col) == 1) meta.col <- c(fixed=meta.col, random=meta.col)
  else stop('Argument meta.col should be a vector length 1 or 2')

  if(length(meta.pch)==2) names(meta.pch) <- c('fixed','random')
  else if(length(meta.pch) == 1) meta.pch <- c(fixed=meta.pch, random=meta.pch)
  else stop('Argument meta.pch should be a vector length 1 or 2')

  if(do.fixed & !x$do.fixed){
    warning("Argument do.fixed should only be set to TRUE if x$do.fixed is also TRUE")
    do.fixed <- FALSE
  }
  if(do.random & !x$do.random){
    warning("Argument do.random should only be set to TRUE if x$do.random is also TRUE")
    do.random <- FALSE
  }

  oldpar <- par(mar=c(6,3,.5,.5), mgp=c(2,1,0))

  plot(NA, type="n", xlim=range(x$x$bt), ylim=range(x$y$bt), xlab=xlab, ylab=ylab)
  abline(h=0, lty=2)
  abline(v=0, lty=2)

  # add each study
  for(i in 1:x$n.studies){
    plot(study2way(x$x$bt[[i]], x$y$bt[[i]], x$study.names[i]), col=col[i], pch=pch[i], center.pch=NULL, lwd=lwd, lty=lty, points=points, max.points=max.points, add=TRUE)
  }

  # add meta-analysis
  if(do.fixed)
    plot(study2way(x$x$meta.fixed, x$y$meta.fixed, 'FE'), col=meta.col['fixed'], pch=meta.pch['fixed'], center.pch=meta.pch['fixed'], lwd=meta.lwd, lty=meta.lty['fixed'], points=meta.points, max.points=max.points, add=TRUE)
  if(do.random)
    plot(study2way(x$x$meta.random, x$y$meta.random, 'RE'), col=meta.col['random'], pch=meta.pch['random'], center.pch=meta.pch['random'], lwd=meta.lwd, lty=meta.lty['random'], points=meta.points, max.points=max.points, add=TRUE)

  if(quadrant.proportions){
    quadrants <- c('bottomleft','bottomright','topleft','topright')
    arg <- list(quadrant.format)
    if(do.fixed){
      # same order as quadrants
      arg <- c(arg, list(c(
        mean(x$x$meta.fixed < 0 & x$y$meta.fixed < 0),
        mean(x$x$meta.fixed > 0 & x$y$meta.fixed < 0),
        mean(x$x$meta.fixed < 0 & x$y$meta.fixed > 0),
        mean(x$x$meta.fixed > 0 & x$y$meta.fixed > 0)
      )*100))
    }
    if(do.random){
      arg <- c(arg, list(c(
        mean(x$x$meta.random < 0 & x$y$meta.random < 0),
        mean(x$x$meta.random > 0 & x$y$meta.random < 0),
        mean(x$x$meta.random < 0 & x$y$meta.random > 0),
        mean(x$x$meta.random > 0 & x$y$meta.random > 0)
      )*100))
    }
    q.label <- do.call(sprintf, arg)
    names(q.label) <- quadrants
    sapply(quadrants, function(q){
      legend(q, legend = q.label[q], bty='n', cex=legend.cex)
    })
  }

  par(mar=c(.1,.1,.1,.1), new=TRUE)
  plot(NA, type="n", xlim=0:1, ylim=0:1, yaxt='n', xaxt='n', xlab='', ylab='', bty ="n")

  if(legend.I2 | legend.tau2){
    l <- function(d){
      t <- ''
      if(legend.tau2){
        t <- paste0(t, "tau2 = ",formatC(d$tau2, format='g', digits=2))
        if(legend.ci) t <- paste0(t, sprintf(" (%.2g-%.2g)", d$tau2.ci[1], d$tau2.ci[2]))
      }
      if(legend.tau2 & legend.I2) t <- paste0(t, "; ")
      if(legend.I2){
        t <- paste0(t, "I2 = ",sprintf("%.1f%%",d$I2*100))
        if(legend.ci) t <- paste0(t, sprintf(" (%.1f%%-%.1f%%)", d$I2.ci[1]*100, d$I2.ci[2]*100))
      }
      t
    }
    text(0,0,labels=paste0(xlab, ": ", l(x$x), "\n", ylab, ": ", l(x$y)), cex=legend.cex, adj=0)
  }

  if(legend){
    l.txt <- x$study.names
    l.lty <- rep(1,x$n.studies)
    l.pch <- rep(NA, x$n.studies)
    l.col <- c(col, rep(meta.col['fixed'], do.fixed*2), rep(meta.col['random'], do.random*2))
    if(do.fixed){
      l.txt <- c(l.txt, c(label.fixed.estimate, label.fixed.conf.area))
      l.lty <- c(l.lty, c(NA, meta.lty['fixed']))
      l.pch <- c(l.pch, c(meta.pch['fixed'], NA))
    }
    if(do.random){
      l.txt <- c(l.txt, c(label.random.estimate, label.random.conf.area))
      l.lty <- c(l.lty, c(NA, meta.lty['random']))
      l.pch <- c(l.pch, c(meta.pch['random'], NA))
    }
    legend('bottomright', lwd=2, bty='n', cex=legend.cex, pt.cex=1.5, ncol=legend.ncol,
           legend = l.txt, lty=l.lty, col = l.col, pch = l.pch)
  }

  par(oldpar)
  invisible(NULL)
}
