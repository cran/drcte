\name{predict.drcte}
\alias{predict.drcte}
\alias{predict.list}

\title{Prediction}

\description{
  Predicting probabilities for models of class 'drcte'.
}

\usage{
  \method{predict}{drcte}(object, newdata, se.fit = FALSE,
  interval = FALSE,
  level = 0.95, na.action = na.pass,
  npmle.type = c("interpolation", "left", "right", "midpoint"),
  robust = FALSE, units = NULL, B = 200, ...)
  \method{predict}{list}(object, newdata, coefs, vcov. = NULL, ...)
}

\arguments{
  \item{object}{an object of class 'drcte' or a list representing a time-to-event model (e.g.: LL.3(), LN.3()... )}
  \item{newdata}{A data frame in which to look for variables with which to predict. The first variable, must always be the time at which predictions are sought; if necessary, other covariates can be added in succeeding columns, in the same order as they are specified by the fitted model. For models where time is the only predictor, 'newdata' can also be provided as a vector.}
  \item{se.fit}{logical. If TRUE standard errors are provided.}
  \item{interval}{logical. If TRUE confidence intervals are provided.}
  \item{level}{Confidence level.}
  \item{na.action}{function determining what should be done with missing values in 'newdata'. The default is to predict 'NA'.}
  \item{npmle.type}{character string relating to the type of prediction for NPMLE. For interval censored observations there is not a unique MLE for a specific time, and, therefore, predictions can be obtained by three methods: "interpolation" (default; predictions are sought on the line connecting the two points bounding the non-unique MLE interval), "left" (taking the left side of the non-unique MLE interval) or "right" (taking the right side of the non-unique MLE interval. The argument is neglected for parametric or kernel based time-to-event models.}
  \item{robust}{a logical value. If TRUE, robust confidence interval are provided (see below for detail)}
  \item{units}{a vector coding for randomisation units or other clustering structures. It causes a cluster robust standard error to be provided.}
  \item{B}{Number of bootstrap resamples for bootstrap based inference}
  \item{coefs}{The coefficients of a parametric fit}
  \item{vcov.}{Variance-covariance matrix}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{With parametric models, ‘naive’ (asymptotic) standard errors are obtained by using the 'delta method, in association to the second derivative of the likelihood function. For seed germination and emergence assays, the observational units are the individual seeds and, most often, they are clustered within randomisation units (e.g., Petri dishes or other types of containers), to which the experimental treatments are allocated. For these and other cases of clustered data, cluster-robust SEs can be obtained by using the so-called ‘sandwich estimator’ (Carroll et al. 1998), which has proven reliable for clustered survival data (Yu and Peng 2008). Robust SEs are obtained by exploiting the facilities provided in the ‘sandwich’ package (Zeileis et al. 2020).}

\value{
A vector (when 'newdata' is missing and 'se.fit' = 'interval' = FALSE) or a data.frame with as many rows as in 'newdata' (or as in the original dataset, when 'newdata' is missing) and a set of columns containing predictors, predictions, standard errors, lower and upper limits of confidence intervals.
}

\references{
Carroll RJ, Wang S, Simpson DG, Stromberg AJ, Ruppert D (1998) The sandwich (robust covariance matrix) estimator. Technical report, Department of statistics. A&M University. Texas (USA)

Onofri, A., Mesgaran, M., & Ritz, C. (2022). A unified framework for the analysis of germination, emergence, and other time-to-event data in weed science. Weed Science, 1-13. doi:10.1017/wsc.2022.8

Yu B, Peng Y (2008) Mixture cure models for multivariate survival data. Comput Stat Data Anal 52:1524–1532

Zeileis A, Koell S, Graham N (2020) Various Versatile Variances: An Object-Oriented Implementation of Clustered Covariances in R. J Stat Softw 95

}

\author{Andrea Onofri, borrowing code from Christian Ritz}

\examples{
library(drcte)
data(chickweed)
modb <- drmte(count ~ start + end,
               data=chickweed,
               fct = LL.3())
predict(modb, se.fit = TRUE, interval = TRUE,
newdata = data.frame(time = c(1, 10, 100, 250)))
}

\keyword{seed germination}
\keyword{time-to-event model}
\concept{prediction}
