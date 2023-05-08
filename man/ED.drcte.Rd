\name{ED.drcte}
\alias{ED.drcte}

\title{Estimating the quantiles for time-to-event models}

\description{ The function \code{ED} estimates the expected time to obtain a certain fraction (percentage) of individuals with events, for one or more time-to-event curves. With 'type = "relative"' the above fraction is regarded as relative to the maximum observed fraction. This function is included for compatibility with the 'drc' package, but the use of the 'quantile()' function is recommended, instead.
}

\usage{
  \method{ED}{drcte}(object, respLev, interval = c("none", "delta", "boot"),
  clevel = NULL, level = ifelse(!(interval == "none"), 0.95, NULL),
  reference = c("control", "upper"), type = c("relative", "absolute"),
  lref, uref, bound = TRUE, od = FALSE, vcov. = vcov, display = TRUE,
  pool = TRUE, logBase = NULL, multcomp = FALSE,
  intType = "confidence", rate = FALSE, B = 200,
  seed = 1234, units = NULL, ...)
}

\arguments{
  \item{object}{an object of class 'drcte'.}
  \item{respLev}{a numeric vector containing the response levels. For 'type = "absolute"' 'respLev' must be included between 0 and 1, while for 'type = "relative"', 'respLev' expresses the percentage of maximum response and must be included between 0 and 100\%. }
  \item{interval}{character string specifying the type of confidence intervals to be supplied. The default is "none".
  See Details below for more explanation.}
  \item{clevel}{character string specifying the curve id in case on estimates for a specific curve or compound is requested. By default estimates
  are shown for all curves.}
  \item{level}{numeric. The level for the confidence intervals. The default is 0.95.}
  \item{reference}{character string. Is the upper limit or the control level the reference?}
  \item{type}{character string. Whether the specified response levels are absolute or relative (default).}
  \item{lref}{numeric value specifying the lower limit to serve as reference.}
  \item{uref}{numeric value specifying the upper limit to serve as reference (e.g., 100\%).}
  \item{bound}{logical. If TRUE only ED values between 0 and 100\% are allowed. FALSE is useful for hormesis models.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{vcov.}{function providing the variance-covariance matrix. \code{\link{vcov}} is the default,
  but \code{sandwich} is also an option (for obtaining robust standard errors).}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}
  \item{pool}{logical. If TRUE curves are pooled. Otherwise they are not. This argument only works for models with independently fitted curves as specified in \code{\link{drm}}.}
  \item{logBase}{numeric. The base of the logarithm in case logarithm transformed dose values are used.}
  \item{multcomp}{logical to switch on output for use with the package multcomp (which needs to be activated first). Default is FALSE (corresponding to the original output).}
  \item{intType}{string specifying the type of interval to use with the predict method in case the type of confidence
  interval chosen with the argument "type" is "inverse regression."}
  \item{rate}{Logical: if FALSE quantiles for time-to-events are reported (default), if TRUE, rates (inverse of times) are reported.}
  \item{B}{Number of bootstrap resamples for bootstrap based inference}
  \item{seed}{seed for sampling methods}
  \item{units}{The name of a variable, coding for the experimental units, within which the observational units are clustered (e.g., Pretri dishes in germination assays)}
  \item{...}{see the details section below.}
}

\details{
There are several options for calculating confidence intervals through the argument \code{interval}. The option "none" (default) does not provide confidence intervals. The option "delta" results in asymptotical Wald-type confidence intervals (using the delta method and the normal distribution) and it is only available for parametric time-to-event models. The option "boot" provides bootraspped confidence intervals and it is only available for NPMLE or KDE curves; it may be rather slow. For (\code{\link{KDE}} models, the additional arguments \code{lower} and \code{upper} may be supplied. These arguments specify the lower and upper limits of the bisection method used to find the ED values.  The lower and upper limits need to be smaller/larger than the EDx level to be calculated. The default limits are 0 and 1}

\value{ An invisible matrix containing the shown matrix with two or more columns, containing the estimates and the corresponding estimated standard errors and possibly lower and upper confidence limits. Or, alternatively, a list with elements that may be plugged directly into \code{parm} in the package \emph{multcomp} (in case the argument \code{multcomp} is TRUE).
}

\seealso{
The 'quantile()' function
}

\author{Andrea Onofri}

\examples{
start <- c(0, 6, 7, 10, 13, 16, 22, 23)
end <- c(6, 7, 10, 13, 16, 22, 23, Inf)
count <- c(3, 9, 3, 8, 1, 1, 5, 5)

# Fitting a non-parametric time-to-event model
mod <- drmte(count ~ start + end, fct = NPMLE())

# Getting the times required to reach three different fractions
# of individuals with events (i.e.: 0.1, 0.3 and 0.5)
ED(mod, respLev = c(0.10, 0.30, 0.50), type = "absolute")

# Getting the times required to reach three 10, 30 and 50%
# of the maximum fraction of individuals with events
ED(mod, respLev = c(10, 30, 50), type = "relative")
}

\keyword{models}
\keyword{time-to-event}
