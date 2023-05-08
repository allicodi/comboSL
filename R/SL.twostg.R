
#' @title Create two stage functions
#' 
#' @description Used to create two stage functions for SuperLearner
#' 
#' @param stg1_fn list of stage 1 function(s)
#' @param stg2_fn list of stage 2 function(s)
#' @param environment to create functions within, default to global environment
#' 
#' @return list of newly created library names in format "SL.stg1.STG1NAME__stg2.STG2NAME"
#' 
#' @export 
#' 
SL.twostg <- function(stg1_fn = c("SL.glm", "SL.ranger"), stg2_fn = c("SL.glm"), twostg_env = .GlobalEnv) {
    
    #get template function
    template_fn <- make_template()

    #create libraries
    SL.library_names <- NULL
    for (s1 in stg1_fn){
        for (s2 in stg2_fn){   

            #gsub names into the template function and assign name to function
            wrapper_fn <- gsub("STG1_FN", s1, template_fn)
            wrapper_fn <- gsub("STG2_FN", s2, wrapper_fn)
            wrapper_fn_eval <- eval(parse(text = wrapper_fn))
            wrapper_fn_name <- paste0("SL.stg1.", s1, "__stg2.", s2)
            
            #assign function to the twostage environment
            assign(wrapper_fn_name, wrapper_fn_eval, env = twostg_env)
            SL.library_names <- c(SL.library_names, wrapper_fn_name)

        }
    }

    #return names of libraries
    return(SL.library_names)
}

#' @title Helper function for creating 2 stage libraries
#' 
#' @description SL.twostg uses the string returned by this function as a template for creating two stage learner functions.
#' Gsub is used to substitute the names of the desired learners in for "STG1_FN" and "STG2_FN".
#' 
#' @return character string to be used as function
#' 
#' @noexport 
#'
make_template <- function(){
    template_fn <- "
    new_fn <- function(Y, X, newX, family, obsWeights, ...){
        stg1_output <- STG1_FN(Y = as.numeric(Y > 0), X = X, newX = newX, family = binomial(), obsWeights = obsWeights, ...)
        stg2_output <- STG2_FN(Y = Y[Y > 0], X = X[Y > 0, , drop = FALSE], newX = newX, family = gaussian(), obsWeights = obsWeights[Y > 0],  ...)
        # E[Y | X] = P(Y > 0 | X) * E[Y | Y > 0, X]
        pred <- stg1_output$pred * stg2_output$pred
        fit <- list(stg1_fit = stg1_output$fit,
            stg2_fit = stg2_output$fit)
        class(fit) <- 'SL.twostg'             
        return(list(pred = pred, fit = fit))
    }
    "
    return(template_fn)
}


#' 
#' @title Prediction for SL.twostg
#' @description Prediction for SL.twostg
#' 
#' @param object SL.twostg object
#' @param newdata Dataframe to generate predictions
#' @param ... unused additional arguments
#' 
#' @export
predict.SL.twostg <- function(object, newdata, ...) {
    # predict from object$stg1_fit
    stg1_predict <- predict(
        object$stg1_fit, newdata = newdata, family = binomial()
    )
    # predict from object$stg2_fit
    stg2_predict <- predict(
        object$stg2_fit, newdata = newdata, family = gaussian()
    )
    # multiply return
    pred <- stg1_predict * stg2_predict
    return(pred)
}
