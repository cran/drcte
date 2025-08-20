\name{plot.drcte}
\alias{plot.drcte}

\title{Plotting fitted time-to-event models}

\description{
  \code{plot} displays fitted curves and observations in the same plot window, distinguishing between curves by different plot symbols and line types.}

\usage{

  \method{plot}{drcte}(x, ...,
  npmle.type = c("interpolation", "midpoint", "right", "left", "none"),
  npmle.points = FALSE, kde.points = TRUE,
  shading = TRUE,
  bp)
}

\arguments{
  \item{x}{an object of class 'drcte'.}
  \item{...}{additional graphical arguments. For instance, use \code{lwd=2} or \code{lwd=3} to increase the width of plot symbols.}
  \item{shading}{For NPMLE, it prints shading on the graphs for the points where the likelihood value is not unique. Defaults to TRUE}
\item{npmle.type}{the NPMLE of the cumulative density function is only specified at the end of each inspection interval, while it is not unique within each interval. This argument specifies how the CDF increases within each interval: possible values are "interpolation" (it is assumed that the CDF increases progressively), "left" (the CDF increases at the beginning of each interval), "right" (the CDF increases at the end of each interval; it is very common in survival analysis) and "midpoint" (the CDF increases in the middle of each interval; it is very common in survival analysis). This argument is neglected with parametric and KDE fits.}
  \item{npmle.points}{If F, with NPMLE fits, plotting the NPMLEs at the end of each interval as symbols is suppressed. Not implemented, yet.}
  \item{kde.points}{If F, with KDE fits, plotting the NPMLE at the end of each interval is suppressed}
  \item{bp}{numeric value specifying the break point below which the dose is zero (the amount of stretching on
  the dose axis above zero in order to create the visual illusion of a logarithmic scale \emph{including} 0).
  The default is the base-10 value corresponding to the rounded value of the minimum of the log10 values of
  all positive dose values. This argument is only working for logarithmic dose axes.}
}

\details{The plot method for 'drcte' objects inherits from the plot method for 'drc' objects and adds functionalities for nonparametric time-to-event fits. For parametric time-to-event models, the fitted curve is presented, together with symbols, corresponding to the NPMLE estimator of the cumulative distribution function at the end of each observation interval. For NPMLE fits, the cumulative CDF at the end of each interval is presented and the way that the density increases during the interval is specified through the argument 'npmle.type'. For KDE fits, the cumulative density is presented as a line, while the observed NPMLE of the CDF at the end of each interval can be either presented as symbols or suppressed. For all other information, please consult ?plot.drc. This method does not work for models with environmental covariates.
}

\value{
An invisible data frame with the values used for plotting the fitted curves. The first column contains the dose values, and the following columns (one for each curve) contain the fitted response values.
}

\author{Andrea Onofri}

\references{
  Weimer, M., Jiang, X., Ponta, O., Stanzel, S., Freyberger, A., Kopp-Schneider, A. (2012)
  The impact of data transformations on concentration-response modeling.
  \emph{Toxicology Letters}, \bold{213}, 292--298.
}

%\note{ ~~further notes~~ }

%\seealso{\code{\link{plotraw}} plots the observations only.}

\examples{
library(drcte)
data(verbascum)
mod <- drmte(nSeeds ~ timeBef + timeAf, fct = NPMLE(),
                 data = verbascum, curveid = Species)
plot(mod)
}
\keyword{plot}

