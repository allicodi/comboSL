#' Method for generating weights (coefficients) under Huber Loss: scaled quadratic programming
#'
#' @param Z A \code{data.frame} with predictions on test data for each learner
#' @param Y A numeric outcome variable
#' @param libraryNames A list of learners to use for prediction
#' @param verbose True if printing fit (not currently used)
#' @param lambda A numeric robustification parameter to regulate influence of outliers
#' @param obsWeights Observation-level weights
#' @param errorsInLibrary (not currently used)
#' @param ... Other arguments (not currently used)
#'
#' @export
#' @importFrom stats glm
#' @return
#' \describe{
#'  \item{\code{computeCoef}}{List with three elements: (i) coef: weights for each algorithm, (ii) cvRisk: the V-fold CV risk for each algorithm, (iii) optimizer: result object from weight optimization}
#'  \item{\code{computePred}}{Super learner predicted values}
#' }
#'
method.CC_HUBER <- function() {
  computeCoef = function(Z, Y, libraryNames, verbose,lambda,
                         obsWeights=rep(1, length(Y)),
                         errorsInLibrary = NULL, ...) {



    # look at Y minus Yhat residuals, figure out lambda based on residuals   
    # nlambda = input? specifying number of lambdas to search over
    # 
    # lambda_grid <- 
    #    # find CV errors (|y - yhat|)
    #    apply(Z, 2, function(x){ abs(x - Y) }) %>% 
    #    # turn into vector
    #    c() %>% 
    #    # look at quantiles of errors
    #    quantile(p = seq(0.01, 0.99, length = nlambda))                 
    # compute cvRisk
    cvRisk <- apply(Z, 2, function(x) 
      mean(ifelse((abs(Y-x) > 10000*lambda),
                  10000*lambda*(obsWeights*abs(Y-x) - 0.5*10000*lambda),
                  0.5*(obsWeights*(x-Y)^2))))
    names(cvRisk) <- libraryNames
    
    modZ <- Z
    # check for columns of all zeros. assume these correspond
    # to errors that SuperLearner sets equal to 0. not a robust
    # solution, since in theory an algorithm could predict 0 for
    # all observations (e.g., SL.mean when all Y in training = 0)
    naCols <- which(apply(Z, 2, function(z){ all(z == 0 ) }))
    anyNACols <- length(naCols) > 0
    if(anyNACols){
      # if present, throw warning identifying learners
      warning(paste0(paste0(libraryNames[naCols],collapse = ", "), " have NAs.",
                     "Removing from super learner."))
    }
    # check for duplicated columns
    # set a tolerance level to avoid numerical instability
    tol <- 8
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    if(anyDupCols){
      # if present, throw warning identifying learners
      warning(paste0(paste0(libraryNames[dupCols],collapse = ", "),
                     " are duplicates of previous learners.",
                     " Removing from super learner."))
    }
    # remove from Z if present
    if(anyDupCols | anyNACols){
      rmCols <- unique(c(naCols,dupCols))
      modZ <- Z[,-rmCols]
    }
    
    # compute coefficients on remaining columns -- Use CVXR 
    # Variables minimized over
    beta <- Variable(ncol(modZ))
    # Problem definition
    objective <- Minimize(sum(huber((Y - modZ%*%beta)/10000, lambda)))
    constraints <- list(beta >= 0, sum(beta)==1)
    prob <- Problem(objective,constraints)
    # Problem solution
    result <- solve(prob,solver="ECOS", MAXIT=as.integer(2000))
    coef <- as.vector(result$getValue(beta))
    
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] = 0
    }
    # add in coefficients with 0 weights for algorithms with NAs
    if(anyDupCols | anyNACols){
      ind <- c(seq_along(coef), rmCols - 0.5)
      coef <- c(coef, rep(0, length(rmCols)))
      coef <- coef[order(ind)]
    }
    # Set very small coefficients to 0 and renormalize.
    coef[coef < 1.0e-4] <- 0
    coef <- coef / sum(coef)
    if(!sum(coef) > 0) warning("All algorithms have zero weight", call. = FALSE)
    list(cvRisk = cvRisk, coef = coef, optimizer = result)
  }
  
  computePred = function(predY, coef, ...) {
    predY %*% matrix(coef)
  }
  out <- list(require = "CVXR",
              computeCoef = computeCoef,
              computePred = computePred)
  invisible(out)
}
