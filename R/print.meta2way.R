#' Print the two-way meta-analysis results based on bootstrap data
#'
#' This function allows you to print the two-way meta-analysis results based on bootstrap data. The DerSimonian and Laird estimator is used for calculation of the tau2 and I2 confidence intervals.
#' @param object Object of class meta2way.
#' @keywords meta-analysis; cost-effectiveness; bootstrapping
#' @export
#' @return
#' Returns NULL invisibly.
#' @examples
#' m <- meta2way(study2way(treatments$Abrahams2018.x, treatments$Abrahams2018.y, 'Abrahams 2018'),
#'               study2way(treatments$Baruch2018.x, treatments$Baruch2018.y, 'Baruch 2018'))
#' m
#' print(m)
print.meta2way <- function(object){
  # function to summarize
  summarizedim <- function(xy){
    q <- t(sapply(xy$bt, function(bt) formatC(quantile(bt, c(.5,.025,.975)), format='g', digits=4)))
    q <- rbind(q, `-----`=rep('-----',ncol(q)))
    if(object$do.fixed){
      q <- rbind(q, `Fixed effects model`= formatC(quantile(xy$meta.fixed, c(.5,.025,.975)), format='g', digits=4))
    }
    if(object$do.random){
      q <- rbind(q, `Random effects model`= formatC(quantile(xy$meta.random, c(.5,.025,.975)), format='g', digits=4))
    }
    colnames(q) <- c('Estimate', 'Lower 95% CI limit', 'Upper 95% CI limit')
    as.data.frame(q)
  }
  psci <- function(p,digits=4){
    if(p < 10^-digits) return(paste0("<",10^-digits))
    if(p > 1-10^-digits) return(paste0(">", 1-10^-digits))
    return(formatC(p, digits, format="f"))
  }

  cat("Two-way meta-analysis of", object$n.studies, "studies.\n")
  cat("Each study contains", object$n.bt, "bootstrap samples.\n")

  cat("\n==========\nX variable:\n")
  print(summarizedim(object$x))
  cat("\n")
  cat("Q: ", formatC(object$x$Q, format='f', digits=4), " (df=", object$n.studies-1, "; p-value=", psci(object$x$Qp),")\n", sep="")
  cat("I2: ", sprintf("%.2f%% (95%% CI %.2f%% - %.2f%%)", object$x$I2 * 100, object$x$I2.ci[1] * 100, object$x$I2.ci[2] * 100), "\n", sep="")
  cat("Tau2: ", sprintf("%.3g (95%% CI %.3g - %.3g)", object$x$tau2, object$x$tau2.ci[1], object$x$tau2.ci[2]), "\n", sep="")
  #cat("Tau2: ", formatC(object$x$tau2, format='e', digits=2), "\n", sep="")

  cat("\n==========\nY variable:\n")
  print(summarizedim(object$y))
  cat("\n")
  cat("Q: ", formatC(object$y$Q, format='f', digits=4), " (df=", object$n.studies-1, "; p-value=", psci(object$y$Qp),")\n", sep="")
  cat("I2: ", sprintf("%.2f%% (95%% CI %.2f%% - %.2f%%)", object$y$I2 * 100, object$y$I2.ci[1] * 100, object$y$I2.ci[2] * 100), "\n", sep="")
  cat("Tau2: ", sprintf("%.3g (95%% CI %.3g - %.3g)", object$y$tau2, object$y$tau2.ci[1], object$y$tau2.ci[2]), "\n", sep="")
  #cat("Tau2: ", formatC(object$y$tau2, format='e', digits=2), "\n", sep="")

  invisible(NULL)
}
