#' Generalized linear model with Gaussian family and log link - Smearing Estimate
#'
#' Implements the "smearing" estimator of Duan (1983), which
#' fits a linear model on the log transformed outcome and transforms to an estimate of the
#' conditional mean on the original scale.
#'
#' @param Y A numeric outcome variable
#' @param X A \code{data.frame} of covariates constituting the training sample
#' @param newX A \code{data.frame} with the same column names and format as \code{X} constituting
#' the validation sample.
#' @param family Gaussian only
#' @param obsWeights Observation-level weights
#' @param ... Other arguments (not currently used)
#' @return
#' \describe{
#'  \item{\code{pred}}{Predicted outcomes based on predictors in \code{newX}}
#'  \item{\code{fit}}{A list with named entries \code{object} (the fitted \code{glm} regression object)
#'                    and \code{smear_factor} (the estimated smearing factor)}
#' }
#' @references Duan, N. (1983) Smearing estimate: a nonparametric retransformation method. Journal of the American Statistical Association, 78, 605-610.
#' @export
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit smear estimator
#' fit_smearglm <- SL.smearglm(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_smearglm <- predict(fit_smearglm$fit, newdata = cost_data[,c("female", "race")])

SL.smearglm <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    logY <- log(Y)
    fit.logGLM <- glm(logY ~ ., data = X, family = family, weights = obsWeights)

    mu <- predict(fit.logGLM, type = "response", newdata = X)
    resid <- logY - mu

    pred <- exp(predict(fit.logGLM, type = "response", newdata = newX)) * mean(exp(resid))
    fit <- list(object = fit.logGLM, smear_factor = mean(exp(resid)))
    class(fit) <- "SL.smearglm"
  }else{
    stop("SL.logGLM.smear not written for binomial family")
  }

  out <- list(fit=fit, pred=pred)
  return(out)
}


#' Prediction method for \code{SL.smearglm}
#'
#' @export
#'
#' @param object An object of class \code{SL.smearglm}
#' @param newdata A \code{data.frame} to generate predictions for
#' @param ... Other arguments (unused)
#' @return A numeric vector of predictions
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit smear estimator
#' fit_smearglm <- SL.smearglm(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_smearglm <- predict(fit_smearglm$fit, newdata = cost_data[,c("female", "race")])

predict.SL.smearglm <- function(object, newdata, ...){
  mu <- predict(object$object, newdata=newdata, type="response")
  return(exp(mu) * object$smear_factor)
}
