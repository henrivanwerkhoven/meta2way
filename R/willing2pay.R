#' Willingness to pay plot
#'
#' This function allows you to plot a willingness to pay plot based on a two-way meta analysis object.
#' @param object object of class meta2way.
#' @param type character length 1, either 'fixed' or 'random'. If object contains only one type of meta-analysis, this argument can be omitted.
#' @param xlim numeric length 2 with the x-axis limits. Defaults to NA which means the argument xlim.quantile is used to determine xlim.
#' @param xlim.quantile numeric length 1: which upper extreme of costs per effect will be plotted. Defaults to .999. Setting the value to 1 means the most extreme value will be plotted.
#' @param ylim numeric length 2 with the y-axis limits. Defaults to c(0,1)
#' @param xlab character x-axis label. Defaults to 'Costs per effect'
#' @param ylab character y-axis label. Defaults to 'Probability of treshold'
#' @param main character plot title. Defaults to 'Willingness to Pay'
#' @param col line color. Defaults to 1.
#' @param lty line type. Defaults to 1.
#' @param add logical whether to add the line to the previous plot. Defaults to FALSE.
#' @param res integer how many data points should be calculated between xlim[1] and xlim[2]. Defaults to 1000.
#' @keywords Willingness to pay; cost-effectiveness; meta-analysis
#' @export
#' @return
#' Returns NULL invisibly.
#' @examples
#' m <- meta2way(study2way(treatments$Abrahams2018.x, treatments$Abrahams2018.y, 'Abrahams 2018'),
#'               study2way(treatments$Baruch2018.x, treatments$Baruch2018.y, 'Baruch 2018'))
#' willing2pay(m, 'fixed')
willing2pay <- function(object, type, xlim=NA, xlim.quantile=.999, ylim=0:1,
                        xlab="Costs per effect", ylab="Probability of acceptance", main="Willingness to Pay",
                        col=1, lty=1, add=FALSE, res=1e3){
  # check arguments and extract data
  if(class(object) == 'study2way'){
    x <- object$x.bt
    y <- object$y.bt
  }else if(class(object) == "meta2way"){
    if(missing(type)){
      if(object$do.fixed & !object$do.random) type <- 'fixed'
      else if(object$do.random & !object$do.fixed) type <- 'random'
      else type <- 'not provided'
    }
    if(type=='fixed'){ x <- object$x$meta.fixed; y <- object$y$meta.fixed }
    else if(type=='random'){ x <- object$x$meta.random; y <- object$y$meta.random }
    else stop("Argument type should be either of 'fixed' or 'random' or can be omitted if object contains only one of these meta-analysis types.")
  }else{
    stop("Argument object should be of class meta2way or study2way")
  }
  in.dominated <- x <= 0 & y >= 0 & !(x == 0 & y == 0)
  in.costsaving <- x >= 0 & y <= 0 & !(x == 0 & y == 0)
  in.topright <- x > 0 & y > 0
  in.bottomleft <- x < 0 & y < 0
  in.cross <- x==0 & y==0
  if(is.na(xlim[1])) xlim <- c(0, quantile(c(rep(-Inf, sum(in.dominated | in.costsaving | in.cross)), y[in.topright | in.bottomleft] / x[in.topright | in.bottomleft]), xlim.quantile))
  else if(length(xlim) != 2 | !is.numeric(xlim) | xlim[2] <= xlim[1]) stop("Argument xlim should be numeric with length 2 with second number > first number.")
  tresholds <- seq(xlim[1], xlim[2], length.out = res)
  probability <- sapply(tresholds, function(treshold){
    mean(in.costsaving | (in.topright & (y / x) <= treshold) | (in.bottomleft & (y / x) > treshold) )
  })
  if(add){
    lines(tresholds, probability, col=col, lty=lty)
  }
  else{
    plot(tresholds, probability, type="l", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=col, lty=lty, main=main)
  }
  invisible(NULL)
}
