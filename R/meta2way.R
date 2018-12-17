#' Two-way meta-analysis of bootstrapped effect and cost data
#'
#' This function allows you to perform a two-way meta-analysis of bootstrapped effect and cost data
#' @param x.bt list of each studies x axis bootstrap values.
#' @param y.bt list of each studies y axis bootstrap values.
#' @param type character which type of meta-analysis to perform. Acceptes either of 'fixed', 'random' or 'both'. Defaults to 'both'.
#' @param study.names character with readable names of the studies. Defaults to NULL, in which case the names of x.bt will be used (with a warning if names of y.bt are different or not available). If x.bt has no names, 'Study A','Study B', etc will be used.
#' @keywords Two-way meta-analysis; Cost-effectiveness
#' @export
#' @return
#' meta2way returns an object of class "meta2way" which represents a two-way meta-analysis of studies, e.g. cost-effectiveness studies.
#' The function plot can be used to plot the two dimensions of the meta-analysis, such as a cost-effectiveness plane. The function print will display the summary results of each study and the meta-analysis for the two dimentions. The function willing2pay can be used to plot a willingness to pay curve.
#' An object of class "meta2way" is a named list containing the following components:
#' \describe{
#'   \item{x}{named list with values for x-variable:
#'     \describe{
#'       \item{bt}{copy of x.bt; if study.names was provided, these will be used as names of bt}
#'       \item{meta.fixed}{fixed effects meta analysis result for x}
#'       \item{meta.random}{random effects meta-analysis result for x}
#'       \item{Q}{value of Q for x}
#'       \item{tau2}{value of tau2 for x}
#'       \item{I2}{value of I2 for x}
#'     }}
#'   \item{y}{named list with values for y-variable, same as x}
#'   \item{study.names}{character vector of study names}
#'   \item{do.fixed}{logical whether or not fixed effects meta-analysis is performed}
#'   \item{do.random}{logical whether or not random effects meta-analysis is performed}
#'   \item{n.studies}{integer length 1 containing number of studies}
#'   \item{n.bt}{integer length 1 containing number of bootstrap samples per study}
#' }
#' @examples
#' # perform meta-analysis using list of x-values and list of y-values
#' m <- meta2way(list(Abrahams=treatments$Abrahams2018.x, Baruch=treatments$Baruch2018.x),
#'               list(Abrahams=treatments$Abrahams2018.y, Baruch=treatments$Baruch2018.y))
#' m
#' # this is a trick to do meta-analysis of the inverse comparison
#' m <- meta2way(list(Abrahams=-treatments$Abrahams2018.x, Baruch=-treatments$Baruch2018.x),
#'               list(Abrahams=-treatments$Abrahams2018.y, Baruch=-treatments$Baruch2018.y))
#' m
#'
#' # meta-analysis can also be performed using study2way objects
#' s1 <- study2way(treatments$Abrahams2018.x, treatments$Abrahams2018.y, 'Abrahams 2018')
#' s2 <- study2way(treatments$Baruch2018.x, treatments$Baruch2018.y, 'Baruch 2018')
#' m <- meta2way(s1, s2)
#' m

meta2way <- function(...){
  UseMethod('meta2way')
}

#' Two-way meta-analysis of bootstrapped effect and cost data
#'
#' For the study2way objects
#' @export
meta2way.study2way <- function(..., type='both'){
  studies <- lapply(match.call(expand.dots=FALSE)$`...`, function(x) if(is.name(x)){eval(x)}else{NULL})
  if(!all(sapply(studies, class) == 'study2way'))
    stop("All arguments except argument type have to be of class study2way.")
  study.names <- sapply(studies, function(x) x$study.name)
  x.bt <- lapply(studies, function(x) x$x.bt)
  y.bt <- lapply(studies, function(x) x$y.bt)
  meta2way(x.bt=x.bt, y.bt=y.bt, type=type, study.names=study.names)
}

#' Two-way meta-analysis of bootstrapped effect and cost data
#'
#' For list objects
#' @export
meta2way.list <- function(x.bt, y.bt, type='both', study.names=NULL){
  # Q calculation Function
  Q <- function(est, se, tau2=0){
    # w.fix <- 1/se^2
    # sum(w.fix * est^2) - sum(w.fix * est)^2 / sum(w.fix)
    w <- 1 / (tau2 + se^2)
    u <- sum(w * est) / sum(w)
    sum((est - u) ^ 2 / (tau2 + se^2))
  }
  Qp <- function(est, se){
    Q <- Q(est, se)
    df <- length(est)-1
    pchisq(Q, df = df, lower.tail = FALSE)
  }
  # I2 calculation Function
  I2 <- function(est, se){
    df <- length(est)-1
    Q <- Q(est, se)
    #C <- sum(1 / se^2) - sum((1 / se^2)^2) / sum(1 / se^2)
    #tau2 <- ifelse(Q <= df, 0, (Q - df) / C)
    I2 <- max(0,(Q - df) / Q)
    I2
  }
  # get I2 based on Tau2
  I2tau2 <- function(tau2, se){
    df <- length(se) - 1
    1 - df / (df + tau2 * (sum(1/se^2) - sum(1/se^2^2) / sum(1/se^2)))
  }
  # Tau2 calculation Function
  tau2 <- function(est, se){
    w.fix <- 1/se^2
    df <- length(est) - 1
    C <- sum(w.fix) - sum(w.fix^2) / sum(w.fix)
    Q <- Q(est, se)
    tau2 <- ifelse(Q <= df, 0, (Q - df) / C)
    tau2
  }
  # get confidence interval boundary of Tau2
  tau2CIbound <- function(est, se, p){
    # if there is no variance in the estimates, tau2 CI is always 0
    if(var(est) == 0) return(0)
    # chi2 statistic of the boundary
    chi2 <- qchisq(p=1-p, df=length(est)-1)
    # check if the boundary is <0 (if so, return 0)
    Qnull <- Q(est, se, 0)
    if(Qnull < chi2) return(0)
    # find positive boundary by using the optimizer function
    opt <- optimize(function(tau2){
      abs(Q(est, se, tau2) - chi2)
    }, interval=c(0,1e9*mean(se^2)), tol=1e-9)
    opt$minimum
  }

  # for confidence intervals of tau2 try method of
  # Viechtbauer Stat Med 2007 (PMID 16463355)
  # I2 can also be calculated as tau2 / (vi.avg + tau2), where vi.avg is ??? see confint function in metafor package
  # then I2 upper bound can be calculated as tau2.upp / (vi.avg + tau2.upp), etc.
  # zie metafor rma en confint functies (in files gezet)

  # check arguments
  if(length(type) != 1 | !(type %in% c('both','fixed','random'))) stop("Type should be one of 'fixed', 'random' or 'both'.")
  if(!missing(x.bt) & !missing(y.bt)){
    if(!is.list(x.bt) | !is.list(y.bt)) stop("x.bt and y.bt should be a list containing each studies bootstrap results")
    if(length(x.bt) != length(y.bt)) stop("x.bt and y.bt should be of equal length")
    if(var(c(sapply(x.bt,length),sapply(y.bt,length))) != 0) stop("Number of bootstrap samples should be the same for each study and for both dimensions")
  }

  # prepare variables
  n.studies <- length(x.bt)
  if(is.null(study.names)){
    if(!is.null(names(x.bt))){
      if(!identical(names(x.bt), names(y.bt))) warning("Names of x.bt and y.bt are not identical; names of x.bt are used as study names.")
      study.names <- names(x.bt)
    }else{
      study.names <- paste('study', LETTERS[1:n.studies])
    }
  }
  names(x.bt) <- names(y.bt) <- study.names
  n.bt <- length(x.bt[[1]])
  do.fixed <- type %in% c('both','fixed')
  do.random <- type %in% c('both','random')

  x.est <- sapply(x.bt, median)
  x.se <- sapply(x.bt, sd)
  x.weight.fixed <- 1 / x.se^2
  x.Q <- Q(x.est, x.se)
  x.Qp <- Qp(x.est, x.se)
  x.tau2 <- tau2(x.est, x.se)
  x.tau2.ci <- c(tau2CIbound(x.est, x.se, .025), tau2CIbound(x.est, x.se, .975))
  x.I2 <- I2(x.est, x.se)
  x.I2.ci <- c(I2tau2(x.tau2.ci[1], x.se), I2tau2(x.tau2.ci[2], x.se))

  y.est <- sapply(y.bt, median)
  y.se <- sapply(y.bt, sd)
  y.weight.fixed <- 1 / y.se^2
  y.Q <- Q(y.est, y.se)
  y.Qp <- Qp(y.est, y.se)
  y.tau2 <- tau2(y.est, y.se)
  y.tau2.ci <- c(tau2CIbound(y.est, y.se, .025), tau2CIbound(y.est, y.se, .975))
  y.I2 <- I2(y.est, y.se)
  y.I2.ci <- c(I2tau2(y.tau2.ci[1], y.se), I2tau2(y.tau2.ci[2], y.se))

  if(do.fixed){
    x.meta.fixed <- rowSums(sapply(1:n.studies, function(i) x.bt[[i]] * x.weight.fixed[i])) / sum(x.weight.fixed)
    y.meta.fixed <- rowSums(sapply(1:n.studies, function(i) y.bt[[i]] * y.weight.fixed[i])) / sum(y.weight.fixed)
  }else{
    x.meta.fixed <- NA
    y.meta.fixed <- NA
  }
  if(do.random){
    # calculate the weighting and the required inflation due to the random effect
    x.weight.random <- 1 / (x.se^2 + x.tau2)
    x.se.inflation <- sqrt(1 / sum(x.weight.random)) / sqrt(1 / sum(x.weight.fixed))
    y.weight.random <- 1 / (y.se^2 + y.tau2)
    y.se.inflation <- sqrt(1 / sum(y.weight.random)) / sqrt(1 / sum(y.weight.fixed))
    # calculate the weighted sum of bootstrap samples
    x.meta.prep <- rowSums(sapply(1:n.studies, function(i) x.bt[[i]] * x.weight.random[i])) / sum(x.weight.random)
    y.meta.prep <- rowSums(sapply(1:n.studies, function(i) y.bt[[i]] * y.weight.random[i])) / sum(y.weight.random)
    # apply the inflation
    x.meta.random <- (1 - x.se.inflation) * median(x.meta.prep) + x.se.inflation * x.meta.prep
    y.meta.random <- (1 - y.se.inflation) * median(y.meta.prep) + y.se.inflation * y.meta.prep
  }else{
    x.meta.random <- NA
    y.meta.random <- NA
  }
  object <- list(
      x = list(
        bt=x.bt,
        meta.fixed=x.meta.fixed,
        meta.random=x.meta.random,
        Q=x.Q,
        Qp=x.Qp,
        tau2=x.tau2,
        tau2.ci=x.tau2.ci,
        I2=x.I2,
        I2.ci=x.I2.ci
      ),
      y = list(
        bt=y.bt,
        meta.fixed=y.meta.fixed,
        meta.random=y.meta.random,
        Q=y.Q,
        Qp=y.Qp,
        tau2=y.tau2,
        tau2.ci=y.tau2.ci,
        I2=y.I2,
        I2.ci=y.I2.ci
      ),
      study.names=study.names,
      do.fixed=do.fixed,
      do.random=do.random,
      n.studies=n.studies,
      n.bt=n.bt
  )
  class(object) <- "meta2way"
  object
}
