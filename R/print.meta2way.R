#' Print the two-way meta-analysis results based on bootstrap data
#'
#' This function allows you to print the two-way meta-analysis results based on bootstrap data.
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
    q <- t(sapply(xy$bt, function(bt) formatC(quantile(bt, c(.5,.025,.975)), format='e', digits=2)))
    q <- rbind(q, `-----`=rep('-----',ncol(q)))
    if(object$do.fixed){
      q <- rbind(q, `Fixed effects model`= formatC(quantile(xy$meta.fixed, c(.5,.025,.975)), format='e', digits=2))
    }
    if(object$do.random){
      q <- rbind(q, `Random effects model`= formatC(quantile(xy$meta.random, c(.5,.025,.975)), format='e', digits=2))
    }
    colnames(q) <- c('Estimate', 'Lower 95% CI limit', 'Upper 95% CI limit')
    as.data.frame(q)
  }

  cat("Two-way meta-analysis of", object$n.studies, "studies.\n")
  cat("Each study contains", object$n.bt, "bootstrap samples.\n")

  cat("\n==========\nX variable:\n")
  print(summarizedim(object$x))
  cat("\n")
  cat("Q: ", formatC(object$x$Q, format='e', digits=2), "\n")
  cat("I2: ", sprintf("%.1f%%", object$x$I2 * 100), "\n")
  cat("Tau2: ", formatC(object$x$tau2, format='e', digits=2), "\n")

  cat("\n==========\nY variable:\n")
  print(summarizedim(object$y))
  cat("\n")
  cat("Q: ", formatC(object$y$Q, format='e', digits=2), "\n")
  cat("I2: ", sprintf("%.1f%%", object$y$I2 * 100), "\n")
  cat("Tau2: ", formatC(object$y$tau2, format='e', digits=2), "\n")

  invisible(NULL)
}
