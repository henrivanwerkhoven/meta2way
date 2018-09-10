#' Plot two-way meta-analysis cost-effectiveness plane based on bootstrap data
#'
#' This function plots the cost-efectiveness plane of a two-way meta-analysis.
#' @param x Object of class meta2way.
#' @param y Should be omitted.
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
#' @param legend logical whether a legend of the studies should be printed. The legend will be printed below the plot. Defaults to TRUE.
#' @param legend.ncol integer number of columns of the legend. Defaults to 2.
#' @param legend.cex relative font size of the legend. Defaults to 1.
#' @keywords meta-analysis; cost-effectiveness; bootstrapping
#' @export
#' @return
#' Returns NULL invisibly.
#' @examples
#' m <- meta2way(study2way(treatments$Abrahams2018.x, treatments$Abrahams2018.y, 'Abrahams 2018'),
#'               study2way(treatments$Baruch2018.x, treatments$Baruch2018.y, 'Baruch 2018'))
#' plot(m, xlab='delta effects', ylab='delta costs')
plot.meta2way <- function(x, y, xlab="x", ylab="y", col=NULL, pch=NULL, lwd=1, lty=1, points=TRUE,
                          meta.col=1, meta.pch=c(17,15), meta.lwd=2, meta.lty=c(2,3), meta.points=FALSE,
                          max.points=2000, quadrant.proportions=FALSE,
                          legend=TRUE, legend.ncol=2, legend.cex=1){

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

  oldpar <- par(mar=c(6,3,.5,.5), mgp=c(2,1,0))

  plot(NA, type="n", xlim=range(x$x$bt), ylim=range(x$y$bt), xlab=xlab, ylab=ylab)
  abline(h=0, lty=2)
  abline(v=0, lty=2)

  # add each study
  for(i in 1:x$n.studies){
    plot(study2way(x$x$bt[[i]], x$y$bt[[i]], x$study.names[i]), col=col[i], pch=pch[i], center.pch=NULL, lwd=lwd, lty=lty, points=points, max.points=max.points, add=TRUE)
  }

  # add meta-analysis
  if(x$do.fixed)
    plot(study2way(x$x$meta.fixed, x$y$meta.fixed, 'FE'), col=meta.col['fixed'], pch=meta.pch['fixed'], center.pch=meta.pch['fixed'], lwd=meta.lwd, lty=meta.lty['fixed'], points=meta.points, max.points=max.points, add=TRUE)
  if(x$do.random)
    plot(study2way(x$x$meta.random, x$y$meta.random, 'RE'), col=meta.col['random'], pch=meta.pch['random'], center.pch=meta.pch['random'], lwd=meta.lwd, lty=meta.lty['random'], points=meta.points, max.points=max.points, add=TRUE)

  if(quadrant.proportions){
    q.fix <- q.ran <- c(bottomleft=NA,bottomright=NA,topleft=NA,topright=NA)
    if(x$do.fixed){
      q.fix <- c(
        bottomleft = mean(x$x$meta.fixed < 0 & x$y$meta.fixed < 0),
        bottomright = mean(x$x$meta.fixed > 0 & x$y$meta.fixed < 0),
        topleft = mean(x$x$meta.fixed < 0 & x$y$meta.fixed > 0),
        topright = mean(x$x$meta.fixed > 0 & x$y$meta.fixed > 0)
      )*100
    }
    if(x$do.random){
      q.ran <- c(
        bottomleft = mean(x$x$meta.random < 0 & x$y$meta.random < 0),
        bottomright = mean(x$x$meta.random > 0 & x$y$meta.random < 0),
        topleft = mean(x$x$meta.random < 0 & x$y$meta.random > 0),
        topright = mean(x$x$meta.random > 0 & x$y$meta.random > 0)
      )*100
    }
    sapply(names(q.fix), function(pos){
      txt <- c(
        sprintf("%.1f%% (FE model)", q.fix[pos]),
        sprintf("%.1f%% (RE model)", q.ran[pos])
      )[c(x$do.fixed, x$do.random)]
      legend(pos,legend = txt, bty='n', cex=legend.cex)
    })
  }

  par(mar=c(.1,.1,.1,.1), new=TRUE)
  plot(NA, type="n", xlim=0:1, ylim=0:1, yaxt='n', xaxt='n', xlab='', ylab='', bty ="n")

  text(0,0,labels=paste0(xlab, ": tau2 = ",formatC(x$x$tau2, format='e', digits=2),"; I2 = ",sprintf("%.1f%%",x$x$I2*100),"\n",
                         ylab, ": tau2 = ",formatC(x$y$tau2, format='e', digits=2),"; I2 = ",sprintf("%.1f%%",x$y$I2*100)), cex=legend.cex, adj=0)

  if(legend){
    l.txt <- x$study.names
    l.lty <- rep(1,x$n.studies)
    l.pch <- rep(NA, x$n.studies)
    l.col <- c(col, rep(meta.col['fixed'], x$do.fixed*2), rep(meta.col['random'], x$do.random*2))
    if(x$do.fixed){
      l.txt <- c(l.txt, c('FE estimate','FE conf area'))
      l.lty <- c(l.lty, c(NA, meta.lty['fixed']))
      l.pch <- c(l.pch, c(meta.pch['fixed'], NA))
    }
    if(x$do.random){
      l.txt <- c(l.txt, c('RE estimate','RE conf area'))
      l.lty <- c(l.lty, c(NA, meta.lty['random']))
      l.pch <- c(l.pch, c(meta.pch['random'], NA))
    }
    legend('bottomright', lwd=2, bty='n', cex=legend.cex, pt.cex=1.5, ncol=legend.ncol,
           legend = l.txt, lty=l.lty, col = l.col, pch = l.pch)
  }

  par(oldpar)
  invisible(NULL)
}
