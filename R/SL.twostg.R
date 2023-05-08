

SL.makeTwoStg <- function(stg1_fn = c("SL.glm", "SL.ranger"), stg2_fn = c("SL.glm"), template_fn) {
    #create new environment for two stage functions
    twostg_env <- new.env() #maybe move
    assign("twostg_env", twostg_env, .GlobalEnv)

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

    #attach environment so custom learners can be accessed
    attach(twostg_env)
    return(SL.library_names)

}

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
