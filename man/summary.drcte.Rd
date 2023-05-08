\name{summary.drcte}
\alias{summary.drcte}

\title{Summarising time-to-event fits}

\description{
  'summary' returns a comprehensive summary for model parameters from a parametric time-to-event fit object of class 'drcte'. For non-parametric models (NPLME), this function reports the set of Turnbull's interval and related masses, while for KDEs it reports the AMISE or bootstrap bandwidth value.
}

\usage{
  \method{summary}{drcte}(object, robust = FALSE, units = NULL,
  type = c("sandwich", "bootstrap", "jackknife"), ...)
}

\arguments{
  \item{object}{an object of class 'drcte'.}
  \item{robust}{Logical: if TRUE, robust sandwich standard errors are printed instead of the asymptotic formula. Defaults to TRUE if 'units' or 'type' is given. It is neglected for non-parametric time-to-event models}
  \item{units}{Optional vector that identifies groups of subjects, used in computing a cluster robust standard error. Like model variables, this is searched for in the dataset pointed to by the data argument. It is neglected for non-parametric time-to-event models.}
  \item{type}{Type of robust standard error. It defaults to 'sandwich' that is the only available method at moment. Bootstrap and jackknife standard errors are under study. It is neglected for non-parametric time-to-event models}
  \item{...}{additional arguments.}
}

\value{
  A list of summary statistics that includes parameter estimates and estimated standard errors.
}
\details{
With parametric models, ‘naive’ (asymptotic) standard errors are obtained from the second derivative of the likelihood function. For seed germination and emergence assays, the observational units are the individual seeds and, most often, they are not independent from one another, as they are clustered within randomisation units (e.g., Petri dishes or other types of containers), to which the experimental treatments are allocated. For these and other cases of clustered data, cluster-robust SEs can be obtained by using the so-called ‘sandwich estimator’ (Carroll et al. 1998), which has proven reliable for clustered survival data (Yu and Peng 2008). Robust SEs are obtained by exploiting the facilities provided in the ‘sandwich’ package (Zeileis et al. 2020).
}
\references{
Carroll RJ, Wang S, Simpson DG, Stromberg AJ, Ruppert D (1998) The sandwich (robust covariance matrix) estimator. Technical report, Department of statistics. A&M University. Texas (USA)

Onofri, A., Mesgaran, M., & Ritz, C. (2022). A unified framework for the analysis of germination, emergence, and other time-to-event data in weed science. Weed Science, 1-13. doi:10.1017/wsc.2022.8

Yu B, Peng Y (2008) Mixture cure models for multivariate survival data. Comput Stat Data Anal 52:1524–1532

Zeileis A, Koell S, Graham N (2020) Various Versatile Variances: An Object-Oriented Implementation of Clustered Covariances in R. J Stat Softw 95
}
\author{Andrea Onofri, Christian Ritz}

%\examples{}

\keyword{models}
\keyword{time-to-event}
