resamp <- function(data = NULL, cluster = NULL, replace) {
  # New functions for bootstrap resampling (02/01/2022)
  # Code was taken and modified from:
  # https://biostat.app.vumc.org/wiki/Main/HowToBootstrapCorrelatedData
  # This function resamples the data with replacement
  # data is a data.frame, having a field for each level of hierarchy
  # cluster is a character vector that identifies the hierarchy in order from top to bottom
  # replace is a logical vector that indicates whether sampling should be with
  # or without replacement at the corresponding level of hierarchy
  # if replace is NA, resampling is not carried out at a certain level
  if(is.null(data)) stop("No data has been given")
  if(!is.data.frame(data)) stop("Data is not a data frame")
  noCall <- try(is.null(cluster), silent = T)
  if(!is.null(cluster)){
    clusterPos <- cluster
    inDF <- try(dplyr::select(data, {{ cluster }}), silent = T)
    if(is(inDF, "try-error")) stop("Group variables not in data")
    cluster <- inDF
  }

  # Start resampling/permuting
  nStrata <- length(cluster)
  if(is.null(cluster)){
    if(!is.na(replace[1])) {
      sel <- sample(1:length(data[,1]), replace = replace[1])
      ret <- data[sel, ]
    } else {
      ret <- data
    }
  } else {
    if(!is.na(replace[1])) {
      cls <- sample(unique(cluster[,1]), replace = replace[1])
    } else {
      cls <- unique(cluster[,1])
    }
    # subset on the sampled clustering factors
    sub <- lapply(cls, function(b) subset(data, cluster[,1] == b))
    if(nStrata > 1) {
      if(!is.na(replace[2])) sub <- lapply(sub, function(j) resamp(j, cluster = clusterPos[-1], replace = replace[-1]))
    }
    if(nStrata == 1) {
      if(!is.na(replace[2])) sub <- lapply(sub, function(j) resamp(j, cluster = NULL, replace = replace[2]))
    }

    # join and return samples
    ret <- do.call(rbind, sub)
  }
  if(is.vector(ret)){
    ret <- data.frame(ret)
  } else {
    if(!is.vector(data)) row.names(ret) <- 1:length(ret[,1]) else row.names(ret) <- 1:length(ret[,1])
  }
  # if(is.vector(data)) return(invisible(as.numeric(unlist(ret)))) else return(invisible(ret))
  if(is.vector(data)) return(as.numeric(unlist(ret))) else return(ret)
}


resample.cens <- function(dat, cluster = NULL, replace) {
  # Old function for bootstrap resampling (23/7/2021)
  # Now superseeded by newer resamp.R
  # This function resamples the data with replacement
  # dat is a data.frame, containing the survival data, in the form left/right
  # One row per individual
  # cluster is a vector (external or internal) codes for groups
  # Only one group is permitted
  # replace is a list, in the form unts/groups
  # the reverse with respect to resamp.R

  if(is.null(cluster)){
    sel <- sample(1:length(dat[,1]), replace = replace[[1]])
    ret <- dat[sel, ]
  } else {
    cls <- sample(unique(cluster), replace = replace[[2]])
    # subset on the sampled clustering factors
    sub <- lapply(cls, function(b) subset(dat, cluster == b))
    sub
    # sample lower levels of hierarchy (if any)
    sub <- lapply(sub, resample.cens, replace = replace[[1]])
    # join and return samples
    ret <- do.call(rbind, sub)
  }
  if(is.vector(ret)){
    ret <- data.frame(ret)
    names(ret) <- names(dat)
    row.names(ret) <- 1:length(dat[,1])
  } else {
    # print(ret)
    row.names(ret) <- 1:length(dat[,1])
  }
  ret
}

simulateTE <- function(start, end, count, B = 1, groups = NULL){
  # This function produces a new time-to-event dataset, by bootstrap
  # resampling the given dataset. It returns a list of datasets
  # set.seed(seed)
  # Superseded !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # But still used by quantile(), with KDE fits
  res <- vector("list", length = B)
  for(i in 1:B){
    # bootSam <- sample(1:sum(count), replace = T)
    # startSam <- rep(start, count)[bootSam]
    # endSam <- rep(end, count)[bootSam]
    # bootSam <- sample(1:sum(count), replace = T)
    L <- rep(start, count)
    R <- rep(end, count)
    # print(R); print(L); print(count); print(groups)
    # stop()
    if(!is.null(groups)) groups <- rep(groups, count)
    df <- data.frame(L, R)

    if(!is.null(groups)){
      cat("\r Cluster Resampling:", i)
      newSam <- resample.cens(df, cluster = groups, replace = list(TRUE, TRUE))
      # print(newSam)
      stop()
    } else {
      cat("\r Resampling:", i)
      newSam <- resample.cens(df, replace = T)
    }

    startSam <- newSam$L; endSam <- newSam$R
    obj <- getNPMLE(survival::Surv(startSam, endSam, type = "interval2") ~ 1)
    newStart <- obj$intmap[1,]
    newEnd <- obj$intmap[2,]
    n <- sum(count)
    newCount <- obj$pf*n
    df <- data.frame(startTime = newStart, endTime = newEnd,
                     count = newCount)
    res[[i]] <- df
  }
  return(res)
}

