

method.CC_HUBERLAMBDA <- function() {
    computeCoef = function(Z, Y, libraryNames, verbose, nlambda = 100,
                         obsWeights=rep(1, length(Y)),
                         errorsInLibrary = NULL, ...) {

        # (1) compute cvRisk
            
        cvRisk <- apply(Z, 2, function(x){ mean((Y-x)^2) })

        #CVRisk with prespecified lambda (Jeremey defined)
        #cvRisk <- apply(Z, 2, function(x) 
        #mean(ifelse((abs(Y-x) > 10000*lambda),
        #            10000*lambda*(obsWeights*abs(Y-x) - 0.5*10000*lambda),
        #            0.5*(obsWeights*(x-Y)^2))))

        names(cvRisk) <- libraryNames

        modZ <- Z

        # (2) compute coefficients

        #create grid of lambdas
        cv_err <- c(apply(Z, 2, function(x){ abs(Y-x) }))
        lambda_grid <- quantile(cv_err, probs=seq(0.01, 1, length = nlambda))

        #keep track coefficients for each choice of lambda
        coef_grid <- NULL

        for(lambda in lambda_grid){

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

            #matrix to track coefficients for each unique lambda
            coef_grid <- cbind(coef_grid, coef)

            if(!sum(coef) > 0) warning("All algorithms have zero weight", call. = FALSE)

        }

        colnames(coef_grid) <- lambda_grid

        #reweight predictions by coefficients for each unique lambda
        #produces matrix with k cross val rows, nlambda columns
        Ypred_grid <- Z %*% coef_grid 

        #recalculate MSE and find lambda that minimizes MSE
        mse <- apply(Ypred_grid, 2, function(x){ mean((Y-x)^2) })
        min_mse <- which(mse == min(mse))

        #handles case if multiple lambdas with identical MSE 
        if(length(min_mse) > 1){
            min_mse <- min_mse[1]
        }

        #return coefficients corresponding to smallest mse
        optimal_coef <- coef_grid[,min_mse]

        # do we also want to return the lambda that minimized cvRisk?? are we allowed to return more things here without breaking it?
        list(cvRisk = cvRisk, coef = optimal_coef, optimizer = result)
    }

    computePred = function(predY, coef, ...) {
        predY %*% matrix(coef)
    }

    out <- list(require = "CVXR",
                computeCoef = computeCoef,
                computePred = computePred)
    invisible(out)

}