#' Adaptive generalized linear model of Manning (2001)
#'
#' This function implements the estimator of Manning and Mullahy (2001).
#' The estimator adaptively selects a link function and family based on the skewness and
#' kurtosis of the data. Some of the family/link functions recommended by Manning are
#' numerically unstable and so this function returns Duan's smearing estimate if the GLM
#' fitting step breaks down.
#'
#' @param Y A numeric outcome variable
#' @param X A \code{data.frame} of covariates constituting the training sample
#' @param newX A \code{data.frame} with the same column names and format as \code{X} constituting
#' the validation sample.
#' @param family Gaussian only
#' @param obsWeights Observation-level weights (not currently used)
#' @param kCut Cut point for kurtosis
#' @param lambdaCut Cut points for skewness
#' @param startNLS Starting values for the non-linear least squares
#' @param ... Other arguments (not currently used)
#'
#' @export
#' @importFrom moments kurtosis
#' @importFrom stats nls
#' @references Manning WG, Mullahy J (2001). “Estimating log models: to transform or not to transform?” Journal of Health Economics, 20(4), 461–494.
#' @return
#' \describe{
#'  \item{\code{pred}}{Predicted outcomes based on predictors in \code{newX}}
#'  \item{\code{fit}}{A list with named entries \code{object} (the fitted regression model object)}
#' }
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit manning model
#' fit_glm.manning <- SL.glm.manning(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_glm.manning <- predict(fit_glm.manning$fit, newdata = cost_data[,c("female", "race")])

SL.glm.manning <- function(Y, X, newX, family = gaussian(),
                           obsWeights = rep(1, length(Y)),
                           kCut = 3, # kurtosis cutpoint
                           lambdaCut = c(0.5, 1.5, 2.5), # skew cutpoint
                           startNLS = 0, # starting values for NLS
                           ...){
  if(family$family=="gaussian"){
    .slcost.require("moments")
    # first do ols on log scale
    logY <- log(Y)
    fit.logGLM <- glm(logY ~ ., data = X, family = family, weights = obsWeights)

    mu <- predict(fit.logGLM, type = "response", newdata = X)
    resid <- logY - mu

    # check kurtosis of residuals
    k <- moments::kurtosis(resid)

    # by default use these methods
    # some of the other GLMs are unstable and if they fail, this
    # algorithm returns log OLS + smearing estimate
    pred <- exp(predict(fit.logGLM, type = "response", newdata = newX)) * mean(exp(resid))
    fit <- list(object = fit.logGLM, mean(exp(resid)))
    class(fit) <- "SL.smearglm"

    try({
      if(k < kCut){
        # perform park's test
        fit.initGLM <- stats::glm(Y ~ ., data = X, weights = obsWeights,
                                  family = gaussian(link="log"))
        muPark <- predict(fit.initGLM, type = "response", newdata = X)
        resid2Park <- (Y - muPark)^2
        fit.parkGLM <- stats::glm(resid2Park ~ muPark, family = gaussian(link="log"))
        lambda1 <- fit.parkGLM$coef[2]
        # use NLS
        if(lambda1 < lambdaCut[1]){
          xNames <- colnames(X)
          d <- length(xNames)
          bNames <- paste0("b",1:d)
          form <- apply(matrix(1:d), 1, function(i){
            paste(c(bNames[i],xNames[i]), collapse = "*")
          })
          formula <- paste(form, collapse = " + ")
          try({
            fit.nls <- stats::nls(as.formula(paste0("Y ~ exp(b0 +",formula,")")),
                                  data = data.frame(Y, X),
                                  start=eval(parse(text=paste0(
                                    "list(b0=0.5,",paste(paste0(bNames, "=", startNLS),collapse=","),")"
                                  )))
            )
          })
          pred <- predict(fit.nls, newdata = newX)
          fit <- list(object = fit.nls)
          class(fit) <- "SL.glm.manning"
        }else if(lambda1 < lambdaCut[2] & lambda1 >= lambdaCut[1]){
          # use poisson glm
          fit.poisGLM <- suppressWarnings(
            stats::glm(Y ~ ., data = X, weights = obsWeights, family = "poisson", control = list(maxit=100))
          )
          pred <- predict(fit.poisGLM, newdata = newX, type = "response")
          fit <- list(object = fit.poisGLM)
          class(fit) <- "SL.glm.manning"
        }else if(lambda1 < lambdaCut[3] & lambda1 >= lambdaCut[2]){
          # use gamma glm
          fit.gammaGLM <- stats::glm(Y ~ ., data = X, weights = obsWeights, family = Gamma(link = 'log'),
                                     control = list(maxit=100))
          pred <- predict(fit.gammaGLM, newdata = newX, type = "response")
          fit <- list(object = fit.gammaGLM)
          class(fit) <- "SL.glm.manning"
        }else if(lambda1 > lambdaCut[3]){
          # use inverse gaussian glm -- not very stable
          fit.invGaussianGLM <- glm(Y ~ ., data = X, weights = obsWeights,
                                    family = inverse.gaussian(link = "log"),
                                    control = list(maxit=100))
          pred <- predict(fit.invGaussianGLM, newdata = newX, type = "response")
          fit <- list(object = fit.invGaussianGLM)
          class(fit) <- "SL.glm.manning"
        }
      }
    }, silent=TRUE)
  }else{
    stop("SL.glm.manning doesn't work with binomial family.")
  }
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Prediction method for \code{SL.glm.manning}
#'
#' @export
#'
#' @param object An object of class \code{SL.glm.manning}
#' @param newdata A \code{data.frame} to generate predictions for
#' @param ... Other arguments (unused)
#' @export
#' @return A numeric vector of predictions
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit manning model
#' fit_glm.manning <- SL.glm.manning(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_glm.manning <- predict(fit_glm.manning$fit, newdata = cost_data[,c("female", "race")])

predict.SL.glm.manning <- function(object, newdata,...){
  pred <- predict(object = object$object, newdata = newdata, type = "response")
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
