#' Generalized linear model with Gaussian family and log link
#'
#' The \code{glm.fit} algorithm is notoriously unstable for this fitting. If
#' convergence issues arise, try different starting values.
#' @param Y A numeric outcome variable
#' @param X A \code{data.frame} of covariates constituting the training sample
#' @param newX A \code{data.frame} with the same column names and format as \code{X} constituting
#' the validation sample.
#' @param family Gaussian only
#' @param obsWeights Observation-level weights (not currently used)
#' @param link Link function
#' @param ... Other arguments (not currently used)
#' @return
#' \describe{
#'  \item{\code{pred}}{Predicted outcomes based on predictors in \code{newX}}
#'  \item{\code{fit}}{A list with named entries \code{object} (the fitted \code{glm} regression object)}
#' }
#' @export
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit gaussian glm with log link
#' fit_glm.gaussianlog <- SL.glm.gaussianlog(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_glm.gaussianlog <- predict(fit_glm.gaussianlog$fit, newdata = cost_data[,c("female", "race")])

SL.glm.gaussianlog <- function(Y, X, newX, family, obsWeights,
                               link = "log", ...){
  if(family$family=="binomial"){
    stop("SL.glm.gaussianlog not written for binomial family")
  }else{
    myFamily <- gaussian(link="log")
    startValues <- c(log(mean(Y)), rep(0, ncol(X)))
    fit.glm <- glm(Y ~ ., data = X, family = gaussian(link = link),
                   weights = obsWeights, start = startValues,
                   control = list(maxit=100))
    pred <- predict(fit.glm, newdata = newX, type="response")

    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
  }
}
