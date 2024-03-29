% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method.huberlambda.R
\name{method.CC_HUBERLAMBDA}
\alias{method.CC_HUBERLAMBDA}
\title{Method for generating weights (coefficients) under Huber Loss: scaled quadratic programming}
\usage{
method.CC_HUBERLAMBDA()
}
\arguments{
\item{Z}{A \code{data.frame} with predictions on test data for each learner}

\item{Y}{A numeric outcome variable}

\item{libraryNames}{A list of learners to use for prediction}

\item{verbose}{True if printing fit (not currently used)}

\item{lambda}{A numeric robustification parameter to regulate influence of outliers}

\item{nlambda}{If lambda not specified, number of lambdas test for optimality}

\item{obsWeights}{Observation-level weights}

\item{errorsInLibrary}{(not currently used)}

\item{...}{Other arguments (not currently used)}
}
\value{
\describe{
 \item{\code{computeCoef}}{List with three elements: (i) coef: weights for each algorithm, (ii) cvRisk: the V-fold CV risk for each algorithm, (iii) optimizer: result object from weight optimization}
 \item{\code{computePred}}{Super learner predicted values}
}
}
\description{
Generates weights under pre-specified lambda or finds optimal lambda from quantiles of residuals
}
