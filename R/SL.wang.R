#' Quantile regression method of Wang and Zhou (2009)
#'
#' This function implements the estimator of Wang and Zhou (2010).
#' The estimator estimates the conditional mean costs by modeling the conditional quantiles of a transformed cost.
#' Using the equivariance property of quantiles to monotone transformations, the quantile estimators may be
#' transformed back to the conditional mean cost on the original scale.
#'
#' @param Y A numeric outcome variable
#' @param X A \code{data.frame} of covariates constituting the training sample
#' @param newX A \code{data.frame} with the same column names and format as \code{X} constituting
#' the validation sample.
#' @param family Gaussian only
#' @param obsWeights Observation-level weights (not currently used)
#' @param g Transformation to apply to \code{Y} before quantile regression is used. Choices
#' are \code{"log"} or \code{"sqrt"}
#' @param m Number of quantiles to compute
#' @param c Constant used to determine truncation level for transforming quantiles to conditional
#' mean
#' @param b Constant used to determine truncation level for transforming quantiles to conditional
#' mean
#' @param ... Other arguments (not currently used)
#'
#' @references Wang HJ, Zhou X (2010). “Estimation of the retransformed conditional mean in health care cost studies.” Biometrika, 97(1), 147–158.
#' @export
#' @return
#' \describe{
#'  \item{\code{pred}}{Predicted outcomes based on predictors in \code{newX}}
#'  \item{\code{fit}}{A list with named entries \code{object} (the fitted \code{rq} regression object),
#'                    \code{alpha} (the controlled level of trimming based), and \code{g_inv}
#'                    (the inverse function of the inputted \code{g})}
#' }
#' @importFrom Rdpack reprompt
#' @importFrom quantreg rq rearrange
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit wang model
#' fit_wang <- SL.wang(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                     newX = cost_data[, c("female", "race")])
SL.wang <- function(Y, X, newX,
                    family = gaussian(),
                    obsWeights = rep(1, length(Y)),
                    g = "log", # transformation of Y
                    m = length(Y), # number of quantiles
                    c = 0.2, # for calculating truncated mean
                    b = 0.05,# for calculating truncated mean
                    ...){
  .slcost.require("quantreg")
  if(family$family=="gaussian"){
    n <- length(Y)
    # calculate alpha_n for calculating truncated mean
    alpha <- c*n^(-1/(1+4*b))
    tau <- seq(alpha, 1-alpha, length = m)

    # transform Y
    if(g == "log"){
      gY <- log(Y)
      g_inv <- function(x){ exp(x) }
    }else if(g == "sqrt"){
      gY <- sqrt(Y)
      g_inv <- function(x){ x^2 }
    }else{
      stop("SL.wang only implemented for g = 'log' or 'sqrt'")
    }

    # get quantile regressions
    suppressWarnings(
      fm <- quantreg::rq(formula = as.formula("gY ~ ."), tau = tau,
                         weights = obsWeights, data = data.frame(gY, X))
    )
    # get predictions
    QhatList <- predict(fm, newdata = newX, stepfun = TRUE, type = "Qhat")

    QhatRearrange <- lapply(QhatList, quantreg::rearrange)
    # transform to means
    pred <- unlist(lapply(QhatRearrange, FUN = function(Q){
      Qw <- g_inv(environment(Q)$y[-which(duplicated(environment(Q)$x))])
      1/(1-2*alpha) * sum(Qw * diff(c(0,tau)))
    }))
  }else{
    stop("SL.wang not written for binomial family")
  }
  fit <- list(object = fm, alpha = alpha, g_inv = g_inv)
  class(fit) <- "SL.wang"
  out <- list(pred = pred, fit = fit)
  return(out)
}


#' Prediction method for \code{SL.wang}
#'
#' @export
#'
#' @param object An object of class \code{SL.wang}
#' @param newdata A \code{data.frame} to generate predictions for
#' @param ... Other arguments (unused)
#' @return A numeric vector of predictions
#' @export
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit wang model
#' fit_wang <- SL.wang(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_wang <- predict(fit_wang$fit, newdata = cost_data[,c("female", "race")])
predict.SL.wang <- function(object, newdata, ...){
  .slcost.require(quantreg)

  QhatList <- predict(object$object, newdata = newdata, stepfun = TRUE, type = "Qhat")

  QhatRearrange <- lapply(QhatList, quantreg::rearrange)

  pred <- mapply(Q = QhatRearrange, dt = diff(c(0,object$object$tau)), function(Q,dt){
    Qw <- do.call(object$g_inv, args = list(x = environment(Q)$y[-which(duplicated(environment(Q)$x))]))
    1/(1-2*object$alpha) * sum(Qw * dt)
  })

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
