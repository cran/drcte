\name{compCDF}
\alias{compCDF}
\title{
Compare time-to-event curves
}
\description{
This function is used to compare time-to-event curves for different experimental treatments. The apprach depends on the type of time-to-event curves, i.e., parametric, NPMLE and KDE (see details)
}
\usage{
compCDF(obj, scores = c("wmw", "logrank1", "logrank2"),
B = 199, type = c("naive", "permutation"), units = NULL, upperl,
lowerl, display = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{obj}{
A time-to-event model fitted with 'drmte' (drcte object)
}
\item{scores}{
This argument is only meaningful for NPMLEs of time-to-event curves and it specifies the type of log-rank statistic to compare the curves. Possible values are 'wmw' (Wilcoxon-Mann-Whitney) , 'logrank1' (Sun's scores), 'logrank2' (Finkelstein's scores) and it defaults to 'wmw'. See details in Fay and Shaw, 2010.
}
\item{B}{
Number of permutations for permutation tests
}
\item{type}{
For parametric models, comparisons are based on a Likelihood Ratio Test (LRT). The argument type specifies whether P-levels should be calculated by using a 'naive' chi-square approximation, or a 'permutation' approach. It defaults to 'naive' and it is neglected for NPMLE and KDE fits, which always use the permutation approach.
}
\item{units}{
A vector for the clustering units (e.g. Petri dishes or other containers for seed germination/emergence assays). When type = 'permutation' and 'units' is given, the permutation approach is performed at the units level and not at the individuals level.
}
\item{upperl}{Only for parametric models: a numeric vector of upper limits for all parameters in the model (the default corresponds to infinity for all parameters).
}
\item{lowerl}{Only for parametric models: a numeric vector of lower limits for all parameters in the model (the default corresponds to minus infinity for all parameters).
}
\item{display}{
logical. If TRUE, results are displayed, otherwise, they are not (useful in simulations).
}
}
\details{
This function is used to compare time-to-event curves fitted with the function 'drmte()' in the 'drcte' package, where treatment factors have been included by using the 'curveid' argument. The type of test statistic depends on the fitted model, whether it is a parametric time-to-event model or a NPMLE or a KDE. With parametric time-to-event models a Likelihood Ratio test is used, with NPMLE, weighted log-rank statistics are used, while for KDEs, a Cramér‐von Mises type distance among curves is used. In all cases, this function is used to determine the p-level for the null (i.e. the curves are not significantly different from each other). WARNING: a permutation based approach may not necessarily work with parametric and KDE models; indeed, the permuted sample can be such that a model cannot be fitted into it. We are working a solution for these cases, but, in the meantime, we recommend that the user resorts to the naive inference approach whenever a permutation based approach does not work.
}
\value{
\item{method}{The test statistic used for the comparison}
\item{scores}{The individual scores for the statistic (when applicable)}
\item{val0}{The observed value for the test statistic}
\item{vali}{The observed values of the test statistic for the permutational samples}
\item{pval}{The P-value for the LRT statistic, by using the chi square approximation}
\item{pvalb}{The P-value for the test statistic, by using permutation methods}
\item{U}{The score values for each group (only for Weighted Log-rank methods)}
\item{N}{Sample size for each group (only for Weighted Log-rank methods)}
}
\references{
Bolker, BM (2008) Ecological models and data in R. Princeton University Press

Fay MP, Shaw PA (2010) Exact and Asymptotic Weighted Logrank Tests for Interval Censored Data: The interval R Package. J Stat Softw 36:1–34

Barreiro-Ures D, Francisco-Fernández M, Cao R, Fraguela BB, Doallo R, González-Andújar JL, Reyes M (2019) Analysis of interval-grouped data in weed science: The binnednp Rcpp package. Ecol Evol 9:10903–10915
}
\author{
Andrea Onofri
}
\note{
Some code is taken from the 'interval' package and from the 'binnednp' package (citations below).
}
\examples{
data(verbascum)
\donttest{
modVerb <- drmte(nSeeds ~ timeBef + timeAf, curveid = Species,
              fct = NPMLE(), data = verbascum)
compCDF(modVerb, units = verbascum$Dish)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ time-to-event methods }
\keyword{ seed germination }
\keyword{ emergence }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
