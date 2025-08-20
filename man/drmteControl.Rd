\name{drmteControl}

\alias{drmteControl}

\title{Sets control arguments}

\description{
  Set control arguments in the control argument in the function 'drmte()'.
}

\usage{
  drmteControl(constr = FALSE, errorm = TRUE, maxIt = 500, method="BFGS",
  noMessage = FALSE, relTol = 1e-07, rmNA=FALSE, useD = FALSE,
  trace = FALSE, otrace = FALSE, warnVal = -1, dscaleThres = 1e-15, rscaleThres = 1e-15,
  conCheck = TRUE)
}

\arguments{
  \item{constr}{logical. If TRUE optimisation is constrained, only yielding non-negative parameters.}
  \item{errorm}{logical specifying whether failed convergence in \code{\link{drmte}} should result
    in an error or only a warning.}
  \item{maxIt}{numeric. The maximum number of iterations in the optimisation procedure.}
  \item{method}{character string. The method used in the optimisation procedure.
    See \code{\link{optim}} for available methods.}
  \item{noMessage}{logical, specifying whether or not messages should be displayed.}
  \item{relTol}{numeric. The relative tolerance in the optimisation procedure.}
  \item{rmNA}{logical. Should NAs be removed from sum of squares used for estimation?
    Default is FALSE (not removed).}
  \item{useD}{logical. If TRUE derivatives are used for estimation (if available).}
  \item{trace}{logical. If TRUE the trace from \code{\link{optim}} is displayed.}
  \item{otrace}{logical. If TRUE the output from \code{\link{optim}} is displayed.}
  \item{warnVal}{numeric. If equal to 0 then the warnings are stored and displayed at the end.
    See under 'warn' in \code{\link{options}}. The default results in suppression of warnings.}
  \item{dscaleThres}{numeric value specifying the threshold for dose scaling.}
  \item{rscaleThres}{numeric value specifying the threshold for response scaling.}
  \item{conCheck}{logical, switching on/off handling of control measurements.}
}

\value{
  A list with 8 components, one for each of the above arguments.
}

\author{Christian Ritz}

\note{none, yet.}

\examples{
### Displaying the default settings
drmteControl()
}
\keyword{models}
\keyword{nonlinear}
