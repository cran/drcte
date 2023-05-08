\name{NPMLE}
\alias{NPMLE}
\alias{NPMLE.fun}
\title{
NPMLE for time-to-event data
}
\description{
This function provides the non-parametric maximum likelihood estimator for time-to-event data for the distribution function for interval‐grouped data.
}
\usage{
NPMLE()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  This function has no arguments
}
\details{
The function 'KDE()' provides the NPMLE for the cumulative distribution function of time-to-event data. It uses the algorithms proposed in the 'interval' package (Fay and Shaw, 2010)
}
\value{
The 'KDE()' function returns a list containing the nonlinear function and other facilities, which are internally used by the 'drmte()' function.
}
\references{
Michael P. Fay, Pamela A. Shaw (2010). Exact and Asymptotic Weighted Logrank
  Tests for Interval Censored Data: The interval R Package. Journal of Statistical
  Software, 36(2), 1-34. URL https://www.jstatsoft.org/v36/i02/.
}
\author{
Andrea Onofri. Codes for the estimation have been largely taken from the 'interval' package (Fay and Shaw, 2010).
}
\examples{
data(chickweed)
mod <- drmte(count ~ start + end, data=chickweed,
              fct = NPMLE())
summary(mod)
plot(mod, ylim = c(0, 0.25), log = "")
}
\keyword{ time-to-event models }
\keyword{ kernel density estimator }
\keyword{ seed germination }
