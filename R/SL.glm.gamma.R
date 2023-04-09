#' Generalized linear model with gamma family and user-specified link function
#'
#' @param Y A numeric outcome variable
#' @param X A \code{data.frame} of covariates constituting the training sample
#' @param newX A \code{data.frame} with the same column names and format as \code{X} constituting
#' the validation sample.
#' @param family Gaussian only
#' @param obsWeights Observation-level weights (not currently used)
#' @param link Link function
#' @param ... Other arguments (not currently used)
#'
#' @export
#' @importFrom stats glm
#' @return
#' \describe{
#'  \item{\code{pred}}{Predicted outcomes based on predictors in \code{newX}}
#'  \item{\code{fit}}{A list with named entries \code{object} (the fitted \code{glm} regression object)}
#' }
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit gamma glm with log link model
#' fit_glm.gamma <- SL.glm.gamma(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_glm.gamma <- predict(fit_glm.gamma$fit, newdata = cost_data[,c("female", "race")])
#'
SL.glm.gamma <- function(Y, X, newX, family, obsWeights, link = 'log', ...){
  if(family$family=="gaussian"){
    fit.glm <- stats::glm(Y ~ ., data = X, family = Gamma(link = link), weights = obsWeights,
                          control=list(maxit = 100))
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm" # can use predict.SL.glm
    out <- list(pred = pred, fit = fit)
    return(out)
  }else{
    stop("SL.glm.gammalog not written for binomial family")
  }
}



#' Generalized linear model with gamma family and identity link function
#'
#' See \code{\link[twostageSL]{SL.glm.gamma}}.
#'
#' @param ... Other arguments passed to \code{SL.glm.gamma}.
#' @param link Link function
#'
#' @export
#' @importFrom stats glm
#' @return
#' \describe{
#'  \item{\code{pred}}{Predicted outcomes based on predictors in \code{newX}}
#'  \item{\code{fit}}{A list with named entries \code{object} (the fitted \code{glm} regression object)}
#' }
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit gamma glm with identity link
#' fit_glm.gammaid <- SL.glm.gammaid(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_glm.gammaid <- predict(fit_glm.gammaid$fit, newdata = cost_data[,c("female", "race")])
#'
SL.glm.gammaid <- function(..., link = 'identity'){
  SL.glm.gamma(..., link = link)
}
