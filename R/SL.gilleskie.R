#' Adaptive Hazard method of Gilleskie and Mroz
#'
#' This function implements the estimator of Gilleskie and Mroz (2004),
#' which estimates the conditional mean costs by adaptively modeling the conditional hazard
#' function and back-transforming into an estimate of the conditional mean. For more description,
#' we refer users to the original paper.
#'
#' @param Y A numeric outcome variable
#' @param X A \code{data.frame} of covariates constituting the training sample
#' @param newX A \code{data.frame} with the same column names and format as \code{X} constituting
#' the validation sample.
#' @param family Gaussian only
#' @param obsWeights Observation-level weights (not currently used)
#' @param kValues Number of intervals to bin the variable into in order to estimate the discrete
#' hazard function
#' @param maxPoly The largest degree polynomial to be used in fitting the discrete hazard
#' function
#' @param pValThresh The threshold for p-values used in the covariate selection process
#' @param ... Other arguments (not currently used)
#'
#' @references Gilleskie DB, Mroz TA (2004). “A flexible approach for estimating the effects of covariates on health expenditures.” Journal of Health Economics, 23(2), 391–418.
#'
#' @export
#'
#' @importFrom Rdpack reprompt
#' @importFrom MASS ginv
#' @importFrom stats quantile
#' @importFrom sandwich vcovHC
#'
#' @return
#' \describe{
#'  \item{\code{pred}}{Predicted outcomes based on predictors in \code{newX}}
#'  \item{\code{fit}}{A list with named entries \code{object} (the fitted hazard regression object),
#'                    \code{maxK} (the selected number of partitions), and \code{hK}
#'                    (the mean outcome in each partition)}
#' }
#' @examples
#' # load cost data
#' data(cost_data)
#'
#' fit_gilleskie <- SL.gilleskie(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])

SL.gilleskie <- function(Y, X, newX, family = gaussian(),
                         obsWeights = rep(1, length(Y)),
                         kValues = c(5, 15, 25), # number of intervals
                         maxPoly = 2, # maximum polynomial in hazard regressions
                         pValThresh = 0.05,
                         ...){
  if(family$family == "binomial"){
    stop("SL.gilleskie written only for gaussian family.")
  }
  # need the sandwich package for covariate selection algorithm
  .slcost.require("sandwich")

  # maxPoly describes the polynomial used for the partition variable
  # in the hazard regression. Choosing the number of partitions to be
  # less than this value leads to rank deficiency in glm()
  if(any(kValues < maxPoly)){
    warning("kValue specified that is less than maxPoly. These kValues will be ignored")
    kValues <- kValues[kValues > maxPoly]
  }

  #====================================================
  # get hazard fit over different partitions of data
  #====================================================
  outList <- lapply(split(kValues, 1:length(kValues)), FUN = function(K){
    # break up Y into K+1 partitions
    Ytilde <- cut(Y, breaks = stats::quantile(Y, p = seq(0, 1, length = K+1)),
                  labels = FALSE, include.lowest = TRUE)

    # make a long version data set
    longDat <- data.frame(Ytilde, X, id = 1:length(Y))[rep(1:length(Y),Ytilde),]

    # assign parition number variable
    row.names(longDat)[row.names(longDat) %in% paste(row.names(data.frame(Ytilde,X)))] <-
      paste(row.names(data.frame(Ytilde,X)),".0",sep="")

    # string split the row names of longDat to get bin indicator
    longDat$k <- as.numeric(paste(unlist(strsplit(row.names(longDat),".",fixed = TRUE))[seq(2,nrow(longDat)*2,2)]))+1

    # indicator of falling in a particular partition
    longDat$indY <- as.numeric(longDat$k == longDat$Ytilde)

    # loop to do covariate selection
    pVal <- Inf
    d <- maxPoly
    while(pVal > pValThresh & d >= 1){
      # generate the regression equation
      rhs <- NULL
      for(i in 1:(ncol(X)-1)){
        rhs <- c(rhs, paste0("poly(",colnames(X)[i],",",
                             ifelse(length(unique(X[,i]))>d, d, length(unique(X[,i]))-1),
                             ")*poly(k,",d,")*",colnames(X)[(i+1):(ncol(X))],collapse="+"))
      }
      rhs <- c(rhs, paste0("poly(",colnames(X)[ncol(X)],",",
                           ifelse(length(unique(X[,i]))>d, d, length(unique(X[,i]))-1),
                           ")*poly(k,",d,")"))

      # fit the hazard regression
      suppressWarnings(
        fm <- glm(as.formula(paste0("indY ~ ",paste0(rhs, collapse = "+"))),
                  data = longDat, family = "binomial")
      )
      # get coefficients of degree d
      dropNum <- NULL
      for(cn in colnames(X)){
        dropNum <- c(dropNum, grep(paste0(cn,", ",d,")",d), names(fm$coef[!is.na(fm$coef)])))
      }
      if(length(dropNum) > 0){
        dropCoef <- fm$coef[!is.na(fm$coef)][dropNum]

        # get sandwich covariance matrix for all those coefficients
        fullCov <- sandwich::vcovHC(fm, type = "HC0")
        dropCov <- fullCov[dropNum, dropNum]
        # test significance of those coefficients
        SigInv <- tryCatch(solve(dropCov), error = function(e){
          MASS::ginv(dropCov)
        })
        chi2Stat <- tryCatch({
          t(dropCoef) %*% SigInv %*% dropCoef
        },error = function(e){
          return(0)
        })
        pVal <- pchisq(chi2Stat, lower.tail = FALSE, df = length(dropCoef))
      }
      d <- d-1
    }
    # after done dropping polynomial terms, get hazard predictions
    suppressWarnings(
      longDat$haz <- predict(fm, newdata = longDat, type = "response")
    )

    # calculate likelihood
    tmp <- by(longDat, factor(longDat$id), FUN = function(x){
      prod((c(1,cumprod(1 - x$haz[x$k < x$Ytilde])) * x$haz)^x$indY)
    })
    LKR <- sum(log(as.numeric(tmp))) + length(Y)*log(K)
    # return the likelihood ratio and estimate hazard regression
    return(list(LKR = LKR, fm = fm))
  })

  # figure out which one had highest likelihood ratio
  LKRs <- unlist(lapply(outList, "[[", 1))
  maxLKR <- which.max(LKRs)
  maxK <- kValues[maxLKR]
  thisOut <- outList[[maxLKR]]

  # get mean in each partition for transforming back to mean-scale
  Ytilde <- cut(Y, breaks = quantile(Y, p = seq(0, 1, length = maxK+1)), labels = FALSE,
                include.lowest = TRUE)
  hK <- apply(matrix(1:maxK), 1, function(x){
    mean(Y[Ytilde==x])
  })

  # calculate mean by calculating density of each partition
  pred <- apply(matrix(1:length(newX[,1])), 1, FUN=function(x){
    suppressWarnings(
      haz <- predict(thisOut$fm,
                     newdata=data.frame(newX[x,], k = 1:maxK),
                     type = "response")
    )
    dens <- c(1, cumprod(1 - haz[1:(maxK-1)])) * haz
    sum(hK*dens)
  })

  fit <- list(object = thisOut$fm, maxK = maxK, hK = hK)
  class(fit) <- "SL.gilleskie"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Prediction method for \code{SL.gilleskie}
#'
#' @export
#'
#' @param object An object of class \code{SL.gilleskie}
#' @param newdata A \code{data.frame} for which predictions are desired
#' @param ... Other arguments (unused)
#'
#' @examples
#' # load cost data
#' data(cost_data)
#' # fit gilleskie model
#' fit_gilleskie <- SL.gilleskie(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
#'                               newX = cost_data[, c("female", "race")])
#' # get back predictions
#' pred_gilleskie <- predict(fit_gilleskie$fit, newdata = cost_data[,c("female", "race")])
predict.SL.gilleskie <- function(object, newdata,...){
  pred <- apply(matrix(1:length(newdata[,1])), 1, FUN=function(x){
    suppressWarnings(
      haz <- predict(object$object,newdata = data.frame(newdata[rep(x,object$maxK),],k = 1:object$maxK),
                     type = "response")
    )
    dens <- c(1, cumprod(1 - haz[1:(object$maxK-1)])) * haz
    sum(object$hK*dens)
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
