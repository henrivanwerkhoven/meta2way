#' Two-way analysis of bootstrapped effect and cost data
#'
#' This function allows you to perform a two-way analysis of bootstrapped effect and cost data of a single study
#' @param x.bt vector of the x axis bootstrap values.
#' @param y.bt vector of the y axis bootstrap values.
#' @param study.name character with readable name of the studies.
#' @keywords Two-way meta-analysis; Cost-effectiveness
#' @export
#' @return
#' study2way returns an object of class "study2way" which represents a two-way analysis of studies, e.g. cost-effectiveness studies.
#' The function plot can be used to plot the two dimensions of the analysis, such as a cost-effectiveness plane. The function print will display the results of the study for the two dimentions. The function willing2pay can be used to plot a willingness to pay curve.
#' An object of class "study2way" is a named list containing the following components:
#' \describe{
#'   \item{x.bt}{copy of x.bt}
#'   \item{y.bt}{copy of y.bt}
#'   \item{study.name}{character vector of study name}
#'   \item{n.bt}{integer length 1 containing number of bootstrap samples per study}
#' }
#' @examples
#' # plot two-way analysis of study 1
#' s1 <- study2way(treatments$Abrahams2018.x, treatments$Abrahams2018.y, 'Abrahams 2018')
#' plot(s1)
#'
#' # Add study 2 and perform meta-analysis
#' s2 <- study2way(treatments$Baruch2018.x, treatments$Baruch2018.y, 'Baruch 2018')
#' m <- meta2way(s1, s2)
#' plot(m)

study2way <- function(x.bt, y.bt, study.name){
  # check variables
  if(!is.numeric(x.bt) | !is.numeric(y.bt))
    stop("Arguments x.bt and y.bt should be numeric")
  if((n.bt <- length(x.bt)) != length(y.bt))
    stop("Arguments x.bt and y.bt should have equal length.")
  if(length(study.name) != 1 | !is.character(study.name))
    stop("Argument study.name should be character vector length 1.")
  object <- list(
    x.bt = x.bt,
    y.bt = y.bt,
    study.name = study.name,
    n.bt = n.bt
  )
  class(object) <- 'study2way'
  object
}
