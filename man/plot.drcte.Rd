\name{plot.drcte}
\alias{plot.drcte}

\title{Plotting fitted time-to-event models}

\description{
  \code{plot} displays fitted curves and observations in the same plot window, distinguishing between curves by different plot symbols and line types.}

\usage{

  \method{plot}{drcte}(x, ..., add = FALSE, level = NULL, shading = TRUE,
  type = "all",
  npmle.type = c("interpolation", "midpoint", "right", "left", "none"),
  npmle.points = FALSE, kde.points = TRUE,
  broken = FALSE, bp, bcontrol = NULL, conName = NULL, axes = TRUE,
  gridsize = 100, log = "", xtsty, xttrim = TRUE,
  xt = NULL, xtlab = NULL, xlab, xlim,
  yt = NULL, ytlab = NULL, ylab, ylim,
  cex, cex.axis = 1, col = FALSE, lty, pch,
  legend, legendText, legendPos, cex.legend = 1,
  normal = FALSE, normRef = 1, confidence.level = 0.95)
}

\arguments{
  \item{x}{an object of class 'drcte'.}
  \item{...}{additional graphical arguments. For instance, use \code{lwd=2} or \code{lwd=3} to increase the width of plot symbols.}
  \item{add}{logical. If TRUE then add to already existing plot.}
  \item{level}{vector of character strings. To plot only the curves specified by their identification number.}
  \item{shading}{For NPMLE, it prints shading on the graphs for the points where the likelihood value is not unique. Defaults to TRUE}
  \item{type}{it has been left for analogy with the plot method for 'drc' objects, but it is neglected in the case of 'drcte' objects.}
  \item{npmle.type}{the NPMLE of the cumulative density function is only specified at the end of each inspection interval, while it is not unique within each interval. This argument specifies how the CDF increases within each interval: possible values are "interpolation" (it is assumed that the CDF increases progressively), "left" (the CDF increases at the beginning of each interval), "right" (the CDF increases at the end of each interval; it is very common in survival analysis) and "midpoint" (the CDF increases in the middle of each interval; it is very common in survival analysis). This argument is neglected with parametric and KDE fits.}
  \item{npmle.points}{If F, with NPMLE fits, plotting the NPMLEs at the end of each interval as symbols is suppressed. Not implemented, yet.}
  \item{kde.points}{If F, with KDE fits, plotting the NPMLE at the end of each interval is suppressed}
  \item{broken}{logical. If TRUE the x axis is broken provided this axis is logarithmic
  (using functionality in the CRAN package 'plotrix').}
  \item{bp}{numeric value specifying the break point below which the dose is zero (the amount of stretching on
  the dose axis above zero in order to create the visual illusion of a logarithmic scale \emph{including} 0).
  The default is the base-10 value corresponding to the rounded value of the minimum of the log10 values of
  all positive dose values. This argument is only working for logarithmic dose axes.}
  \item{bcontrol}{a list with components \code{factor}, \code{style} and \code{width}.
  Controlling the appearance of the break (in case \code{broken} is \code{TRUE}).
  The component \code{factor} is the distance from the control to the break as a
  multiple of the value of \code{bp} (default is 2).
  The component \code{style} can take the values: \code{gap}, \code{slash} or \code{zigzag}.
  The component \code{width} is the width of the break symbol (default is 0.02).}
  \item{conName}{character string. Name on x axis for dose zero. Default is '"0"'.}
  \item{axes}{logical indicating whether both axes should be drawn on the plot.}
  \item{gridsize}{numeric. Number of points in the grid used for plotting the fitted curves.}
  \item{log}{a character string which contains '"x"' if the x axis is to be logarithmic, '"y"' if the y axis is to be logarithmic and '"xy"' or
  '"yx"' if both axes are to be logarithmic. The default is "x". The empty string "" yields the original axes.}
  \item{xtsty}{a character string specifying the dose axis style for arrangement of tick marks. By default ("base10")
  For a logarithmic axis by default only base 10 tick marks are shown ("base10"). Otherwise sensible
  equidistantly located tick marks are shown ("standard"), relying on \code{\link{axTicks}}.}
  \item{xttrim}{logical specifying if the number of tick marks should be trimmed in case too many tick marks
  are initially determined.}
  \item{xt}{a numeric vector containing the positions of the tick marks on the x axis.}
  \item{xtlab}{a vector containing the tick marks on the x axis.}
  \item{xlab}{an optional label for the x axis.}
  \item{xlim}{a numeric vector of length two, containing the lower and upper limit for the x axis.}
  \item{yt}{a numeric vector, containing the positions of the tick marks on the y axis.}
  \item{ytlab}{a vector containing the tick marks on the y axis.}
  \item{ylab}{an optional label for the y axis.}
  \item{ylim}{a numeric vector of length two, containing the lower and upper limit for the y axis.}
  \item{cex}{numeric or numeric vector specifying the size of plotting symbols and text
  (see \code{\link{par}} for details).}
  \item{cex.axis}{numeric value specifying the magnification to be used for axis annotation
  relative to the current setting of cex.}
  \item{col}{either logical or a vector of colours. If TRUE default colours are used.
  If FALSE (default) no colours are used.}
  \item{legend}{logical. If TRUE a legend is displayed.}
  \item{legendText}{a character string or vector of character strings specifying the legend text
  (the position of the upper right corner of the legend box).}
  \item{legendPos}{numeric vector of length 2 giving the position of the legend.}
  \item{cex.legend}{numeric specifying the legend text size.}
  \item{lty}{a numeric vector specifying the line types.}
  \item{pch}{a vector of plotting characters or symbols (see \code{\link{points}}).}
  \item{normal}{logical. If TRUE the plot of the normalized data and fitted curves are shown
  (for details see Weimer et al. (2012) for details).}
  \item{normRef}{numeric specifying the reference for the normalization (default is 1).}
  \item{confidence.level}{confidence level for error bars. Defaults to 0.95.}
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

