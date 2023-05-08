\name{plotData}
\alias{plotData}

\title{Data for plotting fitted time-to-event models}

\description{
  \code{plotData} returns the data for plotting a regression model with other plotting packages, such as ggplot2}

\usage{
  plotData(x, xlim, confidence.level = 0.95, gridsize = 100,
  type = c("average", "all", "bars", "none", "obs", "confidence"),
  npmle.type = c("interpolation", "midpoint", "right", "left","none"))
  }

\arguments{
  \item{x}{an object of class 'drcte' or 'drc'.}
  \item{xlim}{as in the basic plot method: the x-range for plotting}
  \item{confidence.level}{For confidence intervals}
  \item{gridsize}{the gridsize for predictions}
  \item{type}{it has been left for analogy with the plot method for 'drc' objects, but it is neglected in the case of 'drcte' objects.}
  \item{npmle.type}{the NPMLE of the cumulative density function is only specified at the end of each inspection interval, while it is not unique within each interval. This argument specifies how the CDF increases within each interval: possible values are "interpolation" (it is assumed that the CDF increases progressively), "left" (the CDF increases at the beginning of each interval), "right" (the CDF increases at the end of each interval; it is very common in survival analysis) and "midpoint" (the CDF increases in the middle of each interval; it is very common in survival analysis). This argument is neglected with parametric and KDE fits.}
 }

\details{The plotData method for 'drcte' objects is used to get the data for plotting with ggplot2.}

\value{
A list with two elements. The first element (plotPoints) contains the data for plotting the observed values, the second element (plotFits) contains the data for plotting curves.
}

\author{Andrea Onofri}

% \references{ }

%\note{ ~~further notes~~ }

%\seealso{\code{\link{plotraw}} plots the observations only.}

\examples{
library(drcte)
data(verbascum)
mod <- drmte(nSeeds ~ timeBef + timeAf, fct = LL.3(),
                 data = verbascum, curveid = Species)
plotData(mod)
}
\keyword{plot}

