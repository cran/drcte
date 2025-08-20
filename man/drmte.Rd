\name{drmte}
\alias{drmte}
\title{Fitting time-to-event models for seed science}

\description{
A general platform for the analysis of seed germination/emergence and other time-to-event data in agriculture.
}

\usage{
  drmte(formula, curveid, pmodels, data = NULL, subset, fct,
  start, na.action = na.omit, logDose = NULL, type = "event",
  control = drmteControl(), lowerl = NULL, upperl = NULL, separate = FALSE,
  pshifts = NULL, varcov = NULL)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. I must be in the form 'count \eqn{~} start + end', where count is the number
  of events observed in the interval between 'start' and 'end'. Other variables can be specified after the 'start' and 'end', as needed to fit,
  e.g., hydro-time or thermal-time models.}
  \item{fct}{a list with three or more elements specifying the distribution of event times, the accompanying self starter function, the names
  of the parameter in the non-linear function and, optionally, the first and second derivatives as well as information used for the calculation
  of quantiles. Currently available distributions are the exponential Cumulative Distribution Function \code{\link{exponential}}, the log-normal
  CDF \code{\link{lognormal}} and the log-logistic CDF \code{\link{loglogistic}}. All these CDF permit the inclusion of a fraction without
  the event  at the end of the assay ('cured' fraction). Several nonlinear functions in the 'drc' package can also be used for the distribution
  of event-times, such as the three-parameter log-logistic model \code{\link[drc]{LL.3}}, the three-parameter Weibull type-1 model
  \code{\link[drc]{W1.3}} and the Weibull type-2 model \code{\link[drc]{W2.3}}.}
  \item{curveid}{a numeric vector or factor containing the grouping of the data.}
  \item{pmodels}{a data frame with a many columns as there are parameters in the non-linear function.
  Or a list containing a formula for each parameter in the nonlinear function.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{type}{a character string specifying the distribution of the data (parameter estimation will depend on the assumed distribution as different log likelihood functions will be used). The default is "event", implying a multinomial distribution.}
  \item{start}{an optional numeric vector containing starting values for all mean parameters in the model.
  Overrules any self starter function.}

\item{na.action}{a function for treating mising values ('NA's). Default is \code{\link{na.omit}}.}

\item{logDose}{a numeric value or NULL. If log doses value are provided the base of the logarithm should be specified (exp(1) for the natural logarithm
  and 10 for 10-logarithm).}

\item{control}{a list of arguments controlling constrained optimisation (zero as boundary),
  maximum number of iteration in the optimisation,
  relative tolerance in the optimisation, warnings issued during the optimisation.}

  \item{lowerl}{a numeric vector of lower limits for all parameters in the model
  (the default corresponds to minus infinity for all parameters).}

  \item{upperl}{a numeric vector of upper limits for all parameters in the model
  (the default corresponds to plus infinity for all parameters).}

  \item{separate}{logical value indicating whether curves should be fit separately (independent of each other).}

  \item{pshifts}{a matrix of constants to be added to the matrix of parameters. Default is no shift for all parameters.}

   \item{varcov}{an optional user-defined known variance-covariance matrix for the responses. Default is the identity matrix (NULL), corresponding to independent response values with a common standard deviation, which will be estimated from the data.}
}

\details{
  This function relies on the general optimiser function \code{\link{optim}} for the minimisation of negative log likelihood function.
The control arguments are specified using the function \code{\link[drc]{drmc}}.

Setting \code{lowerl} and/or \code{upperl} automatically invokes constrained optimisation.

  The columns of a data frame argument to \code{pmodels} are automatically converted into factors.
  This does not happen if a list is specified.
}

\value{
  An object of class 'drcte' and 'drc'.
}

%\references{ ~put references to the literature/web site here ~ }

\author{Andrea Onofri}

% \note{}

\examples{
data(verbascum)
modVerb <- drmte(nSeeds ~ timeBef + timeAf, curveid = Species,
             fct = NPMLE(), data = verbascum)
}

\keyword{seed science}
\keyword{time-to-event regression}
