\name{KDE}
\alias{KDE}
\alias{KDE.fun}
\title{
Kernel estimator for the cumulative distribution function
}
\description{
This function provides the kernel estimator for the distribution function for interval‐grouped data.
}
\usage{
KDE(bw = c("AMISE", "boot"))
KDE.fun(x, start, end, count, h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{for time-to-event data, the predictor variable is time}
  \item{start}{for each assessment interval, this is the start time}
  \item{end}{for each assessment interval, this is the end time}
  \item{count}{for each assessment interval, this is the number of events}
  \item{h}{bandwidth value}
  \item{bw}{Method for bandwidth selection}
}
\details{
The function 'KDE.fun()' provides the kernel estimator for the distribution function, as modified for interval-grouped data. The equation is:

\deqn{F(x) = \sum {w_i * pnorm( \frac{x - t_i}{h}) }}

where 'w_i' is the observed proportion of individuals getting the event in the i-th interval, 't_i' is the mid-point of the i-th time interval, 'x' is the time and 'h' is the bandwidth. The function KDE() is meant to be used with the 'drmte()' function, to estimate the optimal bandwidth, according to a bootstrap selector (Barreiro-Ures et al., 2019).
}
\value{
The 'KDE.fun()' functions returns the CDF, for given values of time. The KDE() function returns a list containing the nonlinear function and other facilities, which are internally used by the 'drmte()' function.
}
\references{
Barreiro‐Ures, D., Francisco‐Fernández, M., Cao, R., Fraguela, B.B., Doallo, R., González‐Andújar, J.L., Reyes, M., 2019. Analysis of interval‐grouped data in weed science: The binnednp Rcpp package. Ecol Evol 9, 10903–10915.

Cao, R., Francisco-Fernández, M., Anand, A., Bastida, F., González-Andújar, J.L., 2013. Modeling Bromus diandrus Seedling Emergence Using Nonparametric Estimation. Journal of Agricultural, Biological, and Environmental Statistics 18, 64–86.

Cao, R., Francisco-Fernández, M., Anand, A., Bastida, F., González-Andújar, J.L., 2011. Computing statistical indices for hydrothermal times using weed emergence data. The Journal of Agricultural Science 149, 701–712.

Gonzalez-Andujar, J.L., Francisco-Fernandez, M., Cao, R., Reyes, M., Urbano, J.M., Forcella, F., Bastida, F., 2016. A comparative study between nonlinear regression and nonparametric approaches for modelling Phalaris paradoxa seedling emergence. Weed Research 56, 367–376.
}
\author{
Andrea Onofri. Codes for the estimation of optimal bandwidth have been largely taken from the 'binnednp' package (Barreiro-Ures et al., 2019).
}
\examples{
data(chickweed)
mod <- drmte(count ~ start + end, data=chickweed,
              fct = KDE())
summary(mod)
plot(mod, ylim = c(0, 0.25), log = "")
}
\keyword{ time-to-event models }
\keyword{ kernel density estimator }
\keyword{ seed germination }
