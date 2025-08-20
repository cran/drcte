\name{exponential}
\alias{exponential}
\alias{exponentialSurv}
\title{
(Shifted) - Exponential distribution of germination times within a seed lot
}
\description{
This function provides the truncated and shifted exponential cumulative distribution of germination times, to be used to fit time-to-event methods to the time-course of seed germination for a seed lot.
}
\usage{
exponential(fixed = c(NA, NA, NA), names = c("b", "d", "shift"))
exponentialSurv(fixed = c(NA, NA, NA), names = c("b", "d", "shift"))
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed. It is used to fix the shifting parameter ('shift') to 0, so that this becomes the usual exponential cumulative probability function. It is also used to fix the 'd' parameter to 1, which provides the usual (non-truncated) exponential distribution, with no final fraction of individuals without the event.
}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.
}
}
\details{
The exponential distribution of event times is parameterised as:

\deqn{ P(t) = d\left[ 1 - \exp(-b (t - shift)) \right]}


where 't' is the time and 'P(t)' is the proportion of germinated seeds at time = t.
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
}
