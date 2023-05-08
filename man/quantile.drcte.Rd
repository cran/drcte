\name{quantile.drcte}
\alias{quantile.drcte}

\title{Estimating quantiles for time-to-event models}

\description{
This function estimates the quantiles for time-to-events data, either for the whole sample (restricted = FALSE) or only for the fraction of individuals with events (restricted = TRUE). The quantiles for rates (the inverse of times) are also available, by using the 'rate = TRUE' option. The quantiles represent the time/rate to reach a given fraction of individuals with event}

\usage{
  \method{quantile}{drcte}(x, probs, restricted = FALSE, rate = FALSE,
                           interval = FALSE,
                           level = ifelse(!(interval == "none"), 0.95, NULL),
                           robust = FALSE, B = 999, units = NULL,
                           display = TRUE, ...)
}

\arguments{
  \item{x}{an object of class 'drcte'.}
  \item{probs}{a vector with the sought probabilities; the values must be included between 0 and 1.}
  \item{restricted}{a logical value: if FALSE, the quantiles are calculated for the whole sample (the default); if TRUE, only the fraction of individuals with events is considered.}
    \item{rate}{a logical value: if FALSE quantiles for time-to-events are reported (default), if TRUE, rates (inverse of times) are reported.}

  \item{interval}{a logical value: if TRUE, confidence intervals are supplied. The default is FALSE. See Details below for more explanation.}
  \item{level}{numeric. The level for the confidence intervals. The default is 0.95.}
  \item{robust}{a logical value: if TRUE, robust standard error and confidence intervals are reported. The default is FALSE.}
  \item{B}{Number of bootstrap resamples for bootstrap based inference}
  \item{units}{a vector coding for randomisation units or other clustering structures. It causes a cluster robust standard error to be provided.}
  \item{display}{a logical value: if TRUE (default), the results are displayed at the end of calculations.}
  \item{\dots}{further arguments passed to or from other methods. In particular, it is used to pass environmental covariates, wherever necessary}
}

\details{
Quantiles for parametric models are calculated by using the inverse of the cumulative probability function. For NPMLE fits, it is assumed that events are evenly spread within each Turnbull interval and the quantiles are the values corresponding to the abscissas of the points where the horizontal line starting from the selected percentile crosses the time-to-event line. For KDE fits, quantiles are calculated from the time-to-event curve by the bisection method, which can be rather slow. If requested, confidence intervals are calculated, according to the type of time-to-event model. For parametric models, standard errors are calculated by using the delta method with Wald-type confidence intervals (normal distribution). If robust = TRUE, a sandwich estimator of the variance-covariance matrix is used, which is in its cluster-robust form, if a 'units' variable is provided, containing the coding fro Petri dishes, other grouping structures. For NPMLE fits, bootrapped confidence intervals are provided, based on B resamples (default to 1000); if the 'units' argument is provided, the cluster-robust version of the bootstrap is used. For KDE fits, standard errors and confidence intervals are not yet provided.}

\value{
  An invisible matrix containing the shown matrix with two or more columns, containing the estimates and the corresponding estimated standard errors and possibly lower and upper confidence limits.}

\seealso{
The 'ED()' function.
}

\author{Andrea Onofri}

\examples{
start <- c(0, 6, 7, 10, 13, 16, 22, 23)
end <- c(6, 7, 10, 13, 16, 22, 23, Inf)
count <- c(3, 9, 3, 8, 1, 1, 5, 5)

# Fitting a non-parametric time-to-event model
mod <- drmte(count ~ start + end, fct = NPMLE())

# Getting the times required to reach three different fractions
# of individuals with events (i.e.: 0.1, 0.3 and 0.5), either including
# the whole lot of individuals, or only the individuals with events
quantile(mod, probs = c(0.10, 0.30, 0.50))
quantile(mod, probs = c(0.10, 0.30, 0.50), restricted = TRUE)
}

\keyword{models}
\keyword{time-to-event}
