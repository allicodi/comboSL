#' Accelerated failure time models
#'
#' This function implements accelerated failure time regression models to estimate
#' the conditional survival, which is back-transformed to estimate the conditional mean.
#' It leverages the \code{flexsurv} package to fit accelerated failure time models and
#' uses numerical integration to back-transform. The numerical integration step can
#' be unstable and we have tried to build in checks to guard against this. In particular, we
#' first try to integrate with upper limit = Inf, but if that fails move to 1e8,
#'  which sometimes is able to provide a sane answer when upper limit = Inf fails.
#' The function keeps trying smaller and smaller values, but will not go smaller than 1e6.
#' In that case it returns a random number between 0 and \code{randomUpper} (default is \code{maxY}).
#'
#' @param Y A numeric outcome variable
#' @param X A \code{data.frame} of covariates constituting the training sample
#' @param newX A \code{data.frame} with the same column names and format as \code{X} constituting
#' the validation sample.
#' @param family Gaussian only
#' @param obsWeights Observation-level weights (not currently used)
#' @param dist Distribution for accelerated failure time model (defaults to generalized Gamma)
#' @param randomUpper If numeric integration fails, for the purposes of stability while fitting
#' the \code{SuperLearner}, a random number is selected uniformly between 0 and \code{randomUpper}
#' (defaults to \code{max(Y)}).
#' @return
#' \describe{
#'  \item{\code{pred}}{Predicted outcomes based on predictors in \code{newX}}
#'  \item{\code{fit}}{A list with named entries \code{object} (the fitted \code{flexsurv} model object),
#'                    \code{randomUpper}.}
#' }
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit AFT model
#' fit_flexsurv <- SL.flexsurv(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_flexsurv <- predict(fit_flexsurv$fit, newdata = cost_data[,c("female", "race")])
#' @export
#' @importFrom flexsurv flexsurvreg
SL.flexsurv <- function(Y, X, newX, family, obsWeights,
                        dist = "gengamma",
                        randomUpper = max(Y), ...){
  .slcost.require("flexsurv")

  if(family$family == "gaussian"){
    fit.flexSurv <- flexsurv::flexsurvreg(
      as.formula(paste0("Surv(Y, rep(1, length(Y))) ~", paste(colnames(X), collapse = " + "))),
      data = X, dist = dist
    )

    pred <- predict.SL.flexsurv(object = list(object = fit.flexSurv, randomUpper = randomUpper),
                                newdata = newX,
                                randomUpper = randomUpper)

    fit <- list(object = fit.flexSurv, randomUpper = randomUpper)
    class(fit) <- "SL.flexsurv"
    out <- list(fit = fit, pred = pred)
  }else{
    stop("SL.flexsurv not implemented for binominal family")
  }

  return(out)
}

#' Accelerated failure time model with generalized gamma family.
#'
#' @param ... Other arguments passed to \code{\link[twostageSL]{SL.flexsurv}}.
#' @param link Link function
#' @references Manning WG, Basu A and Mullahy J. Generalized modeling approaches to risk adjustment of skewed outcomes data. Journal of Health Economics 2005; 24(3): 465–488.
#' @export

SL.flexsurv.gengamma <- function(..., dist="gengamma"){
  SL.flexsurv(..., dist=dist)
}

#' Accelerated failure time model with Weibull family.
#'
#' @param ... Other arguments passed to \code{\link[twostageSL]{SL.flexsurv}}.
#' @param link Link function
#' @references Manning WG, Basu A and Mullahy J. Generalized modeling approaches to risk adjustment of skewed outcomes data. Journal of Health Economics 2005; 24(3): 465–488.
#' @export
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit AFT model with weibull family
#' fit_flexsurv.weibull <- SL.flexsurv.weibull(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_flexsurv.weibull <- predict(fit_flexsurv.weibull$fit, newdata = cost_data[,c("female", "race")])

SL.flexsurv.weibull <- function(..., dist="weibull"){
  SL.flexsurv(..., dist=dist)
}

#' Accelerated failure time model with log-Normal family.
#'
#' See Manning (2005) and \code{\link[twostageSL]{SL.flexsurv}}.
#'
#' @param ... Other arguments passed to \code{SL.flexsurv}.
#' @param link Link function
#' @references Manning WG, Basu A and Mullahy J. Generalized modeling approaches to risk adjustment of skewed outcomes data. Journal of Health Economics 2005; 24(3): 465–488.
#' @export
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit AFT model with weibull family
#' fit_flexsurv.lognormal <- SL.flexsurv.lognormal(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_flexsurv.lognormal <- predict(fit_flexsurv.lognormal$fit, newdata = cost_data[,c("female", "race")])

SL.flexsurv.lognormal <- function(..., dist = "lognormal"){
  SL.flexsurv(..., dist = dist)
}

#' Prediction method for accelerated failure time models
#'
#' The numerical integration step can
#' be unstable and we have tried to build in checks to guard against this. In particular, we
#' first try to integrate with upper limit = Inf, but if that fails move to 1e8,
#' which sometimes is able to provide a sane answer when upper limit = Inf fails.
#' The function keeps trying smaller and smaller values, but will not go smaller than 1e6.
#' In that case it returns a random number between 0 and \code{object$randomUpper}.
#'
#' @param object An object of class \code{SL.smearglm}
#' @param newdata A \code{data.frame} to generate predictions for
#' @param ... Other arguments (unused)
#' @export
#' @return A vector of numeric predictions
#' @importFrom stats integrate
predict.SL.flexsurv <- function(object, newdata, ...){
  # function to return survival probability based on flexsurv object
  .getSurv <- function(x, fit, thisnewdata){
    summary(fit, t = x, B = 0, newdata = thisnewdata)[[1]][,2]
  }

  pred <- as.numeric(apply(matrix(1:nrow(newdata)), 1, function(i){
    upper <- Inf
    out <- NA; class(out) <- "try-error"
    while(class(out)=="try-error" & upper > 1e6){
      out <- try(stats::integrate(.getSurv, fit=object$object,
                                  thisnewdata = newdata[i,],
                                  lower = 0,
                                  upper = upper)$value, silent=TRUE)
      if(upper==Inf){
        upper <- 1e8
      }else{
        upper <- upper/2
      }
    }
    if(class(out)=="try-error"){
      warning(paste0("Unable to integrate survival function. Returning random number between 0 and ", object$randomUpper))
      out <- runif(1, 0, object$randomUpper)
    }
    out
  }))
  return(pred)
}


#' @export
.slcost.require <- function (package, message = paste("loading required package (",
                                                      package, ") failed", sep = ""))
{
  if (!require(package, character.only = TRUE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}
