#' Cox proportional hazard regression models back-transformed to the conditional mean
#'
#' This function implements Cox proportional hazard regression models (David 1972),
#' and back-transforms to the conditional mean scale (see also Basu 2005).
#'
#' @param Y A numeric outcome variable
#' @param X A \code{data.frame} of covariates constituting the training sample
#' @param newX A \code{data.frame} with the same column names and format as \code{X} constituting
#' the validation sample.
#' @param family Gaussian only
#' @param obsWeights Observation-level weights (not currently used)
#' @param ... Other arguments (unused)
#' @return
#' \describe{
#'  \item{\code{pred}}{Predicted outcomes based on predictors in \code{newX}}
#'  \item{\code{fit}}{A list with named entries \code{object} (the fitted \code{coxph} model object),
#'                    \code{randomUpper}.}
#' }
#' @export
#' @references Basu A and Rathouz PJ. Estimating marginal and incremental effects on health outcomes using flexible link and variance function models. Biostatistics 2005; 6(1): 93–109.
#' @importFrom survival Surv coxph
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit cox ph model
#' fit_coxph <- SL.coxph(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_coxph <- predict(fit_coxph$fit, newdata = cost_data[,c("female", "race")])

SL.coxph  <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family == "gaussian"){
    library(survival)
    fit.coxph <- survival::coxph(survival::Surv(Y, rep(1,length(Y))) ~ ., data = X)
    fit <- list(object = fit.coxph)
    class(fit) <- "SL.coxph"
    pred <- predict(object = fit, newdata = newX)
  }else{
    stop("SL.coxph not implemented for binominal family")
  }
  return(list(fit = fit, pred = pred))
}

#' Prediction method for cox proportional hazard regression
#' @param object An object of class \code{SL.coxph}
#' @param newdata A \code{data.frame} to generate predictions for
#' @param ... Other arguments (unused)
#' @export
#' @return A numeric vector of predictions
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit cox ph model
#' fit_coxph <- SL.coxph(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_coxph <- predict(fit_coxph$fit, newdata = cost_data[,c("female", "race")])
#' @importFrom survival survfit
predict.SL.coxph <- function(object, newdata, ...){
  # use surv.fit to get survival estimate and because by default it uses
  # nelson-aalen hazard, easy to convert back to an estimate of the mean
  surv.fit <- survival::survfit(object$object, newdata = newdata)
  pred <- colSums(
    diff(c(0,surv.fit$time))*rbind(
      rep(1,dim(surv.fit$surv)[2]),
      surv.fit$surv[1:(dim(surv.fit$surv)[1]-1),]
    )
  )
  return(pred)
}
