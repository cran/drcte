\name{reshape_te}
\alias{melt_te}
\alias{decumulate_te}
\alias{group_te}
\alias{ungroup_te}

\title{
Reshaping time-to-event datasets}

\description{
This is a set of functions that can be used to transform time-to-event datasets into a form that is amenable for data analyses with the 'drc' and 'drcte' packages}

\usage{
melt_te(data = NULL, count_cols, treat_cols, monitimes,
n.subjects = NULL, grouped = TRUE)
decumulate_te(data = NULL, resp, treat_cols, monitimes,
units, n.subjects, type = c("count", "proportion"))
group_te(data)
ungroup_te(data, counts)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame to be reshaped. In this data frame, each row represents a container (e.g., a Petri dish) and each column represents a piece of information about that container, including the counts of germinated seeds at each inspection time.}
  \item{count_cols}{a numeric vector, specifying the positions (not names) of the columns containing the counts in 'data'}
  \item{treat_cols}{a vector, specifying the positions/names of the columns containing the treatment levels in 'data'}
  \item{monitimes}{a numeric vector of monitoring times. Must be of same length as 'count_cols'}
  \item{n.subjects}{a numeric vector listing the number of viable seeds for each container (e.g., Petri dish). It may be either of length one (if it is the same for all randomisation units) or it must have same length as the number of rows in 'data' and 'treat_cols'. With the function 'melt_te()', it is also possible to pass a reference to the name of a variable in 'data', containing the number of individuals. If missing, it is assumed that all individuals experienced the event and there is no 'cured' fraction.}
  \item{grouped}{logical: specifying whether the output should be in LONG GROUPED or LONG UNGROUPED format}
  \item{counts}{the name of a numeric vector listing the counts of events in each interval. It should be available in data}
  \item{resp}{For 'decumulate_te(): the response variable, i.e. a vector of cumulative counts/proportions}
  \item{units}{For 'decumulate_te(): a label for Petri dishes or other containers}
  \item{type}{For 'decumulate_te(): character string, it is "count" or "proportion", depending on what it is given in the argument 'resp'}
}
\details{
Details are given in \url{https://www.statforbiology.com/2021/stat_drcte_3-reshapingdata/}
}
\value{
Returns a data frame
}
\references{
See \url{https://www.statforbiology.com/2021/stat_drcte_3-reshapingdata/}
}
\author{
Andrea Onofri
}

\examples{
# Transform a dataset from WIDE to LONG GROUPED
library(drcte)
data(lotusOr)
datasetG <- melt_te(lotusOr, count_cols = 3:17, treat_cols = Stage,
            monitimes = 1:15, n.subjects = rep(25,12))
head(datasetG, 16)

# Transform a dataset from WIDE to LONG UNGROUPED
datasetU <- melt_te(lotusOr, count_cols = 3:17, treat_cols = 1,
                monitimes = 1:15, n.subjects = rep(25,12),
                grouped = FALSE)
head(datasetU, 16)

# From LONG GROUPED to LONG UNGROUPED
datasetU2 <- ungroup_te(datasetG, count)[,-c(5, 6)]
head(datasetU2, 16)

# From LONG UNGROUPED to LONG GROUPED
datasetG2 <- group_te(datasetU)
head(datasetG2, 16)

# Chunk 7
# Decumulate a dataset with cumulative counts
data(lotusCum)
dataset_sd <- decumulate_te(lotusCum,
                            resp = nCum,
                            treat_cols = Stage,
                            monitimes = Time,
                            units = Dish,
                            n.subjects = rep(25, 12),
                            type = "count")

dataset_sd <- decumulate_te(lotusCum,
                            resp = Prop,
                            treat_cols = "Stage",
                            monitimes = Time,
                            units = Dish,
                            n.subjects = rep(25, 12),
                            type = "proportion")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ germination }
\keyword{ seed }
