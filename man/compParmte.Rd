\name{compParmte}
\alias{compParmte}

\title{Comparison of parameters for time-to-event models}

\description{
Compare parameters from different assays, either by means of ratios or differences. It is only usable with parametric models
}

\usage{
compParmte(object, strVal, operator = "/", robust = FALSE,
         units = NULL, display = TRUE)
}

\arguments{
  \item{object}{an object of class 'drcte'}
  \item{strVal}{a name of parameter to compare.}
  \item{operator}{a character. If equal to "/" (default) parameter ratios are compared. If equal to "-" parameter differences are compared.}
  \item{robust}{Logical: if TRUE, robust sandwich standard errors are printed instead of the asymptotic formula. Defaults to TRUE if 'units' or 'type' is given.}
  \item{units}{Optional vector that identifies groups of subjects, used in computing a cluster robust standard error. Like model variables, this is searched for in the dataset pointed to by the data argument.}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}
}

\value{
  A matrix with columns containing the estimates, estimated standard errors, values of t-statistics and p-values for the null hypothesis that the ratio equals 1 or that the difference equals 0 (depending on the \code{operator} argument).
}

\details{
  The function compares actual parameter estimates, and therefore the results depend on the parameterisation used.
}

\author{Andrea Onofri and Christian Ritz}

\examples{
data(lotus)
mod <- drmte(count ~ timeBef + timeAf, curveid = Stage,
                fct = loglogistic(), data = lotus,
             upperl = c(NA, 1, NA))
compParmte(mod, "d", robust = TRUE)
}
\keyword{models}
\keyword{time-to-event}
