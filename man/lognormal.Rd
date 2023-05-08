\name{lognormal}
\alias{lognormal}
\alias{lognormalSurv}
\title{
Log-normal distribution of times to an event
}
\description{
This function provides the truncated log-normal cumulative distribution of event times, to be used to fit time-to-event methods to the time-course of seed germination for a seed lot.
}
\usage{
lognormal(fixed = c(NA, NA, NA), names = c("b", "d", "e"))
lognormalSurv(fixed = c(NA, NA, NA), names = c("b", "d", "e"))
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed. It is often used to fix the 'd' parameter to 1, which provides the usual log-logistic distribution, with no ungerminated fraction.
}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.
}
}
\details{
The log-normal distribution of germination times is parameterised as:

\deqn{ P(t) = d * pnorm(b*(log(t + 0.000001) - log(d)))}


where 't' is the time and 'P(t)' is the proportion of germinated seeds at time = t. LogNormalSurv is parameterised as usual in survival analysis:

\deqn{ P(t) = d * pnorm(b*(log(t + 0.000001) - d))}

i.e. with 'e' on a log-scale and a log-link on 'b'. With no cured fraction (d = 1), logNormalSurv gives same results as 'survereg' in the 'survival package', with 'dist = "lognormal"'.
}
\value{
The value returned is a list containing the nonlinear function, the self starter function and the parameter names.
}
\references{
Kleinbaum, D.G., Klein, M., 2005. Survival analysis. Springer Science, USA.
}
\author{
Andrea Onofri
}
\examples{
library(drcte)
data(chickweed)
modb <- drmte(count ~ start + end,
               data=chickweed,
               fct = lognormal())
}
