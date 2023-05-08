melt_te <- function(data = NULL, count_cols, treat_cols, monitimes,
                    n.subjects = NULL, grouped = TRUE){
  # this function reshapes a common field-book for germination assays (WIDE format)
  # into one of two LONG formats (LONG GROUPED and LONG UNGROUPED)

  # data is the dataframe to be processed;
  # count_cols is a numeric vector, specifying the positions of the columns containing the counts
  # treat_cols is a vector, specifying the positions/names of the columns containing the treatment levels
  # monitimes is a numeric vector of monitoring times. Same length as number of columns in dataset
  # n.subjects is a numeric vector listing the number of viable seeds for each Petri dish (same length as
  # number of rows in dataset and treat). If missing, it is assumed that all individuals
  # experienced the event.
  # grouped is a logic element, speciefying whether the output should be LONG GROUPED or UNGROUPED

  # Check a dataset has been passed as argument
  if(is.null(data)) stop("Original dataset not found")

  # Check that moniTimes is of correct length
  if(length(count_cols) != length(monitimes))
  stop("The number of monitoring times must be equal to the number
       of columns with counts")


  # Load the counts as positions (not as a reference)
  # Edited: 4/7/22: tibble needs to be transformed as data.frame, to avoid errors
  counts <- as.data.frame(data[, count_cols])
  colnames(counts) <- monitimes

  # treat_cols: loads as position or name
  treat <- dplyr::select(data, {{ treat_cols }})
  treat <- as.data.frame(treat)
  # treat <- data[, treat]
  # if(is.vector(treat) | is.factor(treat)){
  #   treat <- data.frame(treat)
  #   colnames(treat) <- names(data)[treat_cols]
  # }

  # n.subjects: uses data column or external vector
  tmp1 <- try(is.vector(n.subjects), silent = T)
  if(!is(tmp1, "try-error")){
      if(tmp1){

        if(length(n.subjects) == 1) nViable <- rep(n.subjects, length(counts[,1])) else nViable <- n.subjects
       } else {
       nViable <- apply(counts, 1, sum)
       }}
  else {
    nam.n.subject <- deparse(substitute(n.subjects))
    tmp2 <- try(dplyr::select(data, {{ nam.n.subject }}), silent = T)
    if(is(tmp2, "try-error")) {
      stop("Variable for the number of subjects not found in data")
      } else {
      nViable <- tmp2[,1]
      }
    }

  df <- makeDrm.drcte(counts = counts, treat = treat, nViable = nViable,
                moniTimes = monitimes)

  if(grouped == F){
    df <- df[df$count > 0,]
    frequency <- df$count
    df_surv <- df[rep(seq_len(nrow(df)), frequency),]
    row.names(df_surv) <- 1:(length(df_surv[,1]))
    df <- df_surv
    df <- df[,1:(length(treat) + 3)]
  }
  row.names(df) <- 1:nrow(df)
  return(df)
}

decumulate_te <- function(data = NULL, resp, treat_cols, monitimes, units, n.subjects,
                          type = c("count", "proportion")){
  # This function transform a dataset as cumulative germinations in a dataset for
  # use with DRM type = "event". Dish is a factor that identifies seeds
  # Group is a data.frame

  if(is.null(data)) stop("Original dataset not found")
  type <- match.arg(type)

  cumulative <- T

  # Response variable
  # respName <- substitute(resp)
  # if(length(respName) > 1) stop("The response variable is not unique")
  # if(is.numeric(respName)){
  #   count <- data[,respName]
  # } else {
  #   if(is.symbol(respName)) respName <- deparse(respName)
  #   count <- data[ ,which(names(data) == respName)]
  # }
  tmp <- try(dplyr::select(data, {{ resp }}), silent = T)
  if(is(tmp, "try-error")){
    count <- resp
  } else {
    if(length(tmp[1,]) > 1) stop("The response variable is not unique")
    count <-  as.numeric(tmp[,1])
  }

  # Time variable
  # tmp.name <- substitute(monitimes)
  # if(length(tmp.name) > 1) stop("The time variable is not unique")
  # if(is.numeric(tmp.name)){
  #   moniTime <- data[,tmp.name]
  # } else {
  #   if(is.symbol(tmp.name)) tmp.name <- deparse(tmp.name)
  #   moniTime <- data[ ,which(names(data) == tmp.name)]
  # }
  tmp <- try(dplyr::select(data, {{ monitimes }}), silent = T)
  if(is(tmp, "try-error")){
    moniTime <- monitimes
  } else {
    if(length(tmp[1,]) > 1) stop("The time variable is not unique")
    moniTime <-  as.numeric(tmp[,1])
  }

  # Treatment variables
  tmpGroups <- try(dplyr::select(data, {{ treat_cols }}), silent = T)
  # print(class(tmpGroups))
  # print(treat_cols)
  if(is(tmpGroups, "try-error")){
    treatGroups <- treat_cols
  } else {
    treatGroups <-  tmpGroups
  }

  # Units could be an external variable.
  tmp <- try(dplyr::select(data, {{ units }}), silent = T)
  if(is(tmp, "try-error")){
    Dish <- units
  } else {
    Dish <-  as.factor(tmp[,1])
  }

  # Number of viable seeds per dish. It must be an external variable.
  # It must be either length one,
  # or it must be the same length as the number of levels in Dish
  tmp <- factor(Dish)
  nLev <- length(levels(tmp))
  tmp <- factor(tmp, levels = 1:nLev)
  if(length(n.subjects) == 1) n.subjects <- rep(n.subjects, nLev)
  if(length(n.subjects) != nLev) stop("The length of the vector for subjects is not equal to the number of units")
  nViable <- n.subjects[tmp]

  temp <- data.frame()
  result <- data.frame()
  # Dish <- factor(Dish)
  temp <- data.frame(moniTime, count, nViable, Dish, treatGroups)
  for(i in 1:length(levels(factor(temp$Dish)))){
    temp2 <- temp[temp$Dish == levels(factor(temp$Dish))[i],]
    timeBef <- c(0, temp2$moniTime)
    timeAf <- c(temp2$moniTime, Inf)
    nViable <- c(temp2$nViable, temp2$nViable[1])
    if(cumulative == T){
    #Decumulate if necessary
      if(type == "count") {
        last <- temp2$nViable[1] - tail(temp2$count, 1)
        nSeeds <- c(temp2$count[1], diff(temp2$count), last)
        nCum <- c(temp2$count, NA)
        propCum <- nCum / temp2$nViable[1]
      } else {
        last <- temp2$nViable[1] - tail(temp2$count, 1) * temp2$nViable[1]
        nSeeds <- c(temp2$count[1] * temp2$nViable[1], diff(temp2$count)*temp2$nViable[1], last)
        nCum <- c(temp2$count * temp2$nViable[1], NA)
        propCum <- c(temp2$count, NA)
        }
      }
      else{
      #Cumulate, if necessary NOT NECESSARY!!!!!!
      nSeeds <- c(temp2$count, temp2$nViable[1] - sum(temp2$count))
      nCum <- cumsum(nSeeds)
      nCum[length(nCum)] <- NA
      }

    Dish <- c(temp2$Dish, i)
    GroupT <- temp2[,5:(length(temp2[1,]))]
    GroupT <- data.frame(GroupT)
    GroupT <- rbind(GroupT, tail(GroupT, 1))
    colnames(GroupT) <- colnames(treatGroups)

    # if(type == "proportion") nSeed <- nSeeds * nViable
    dataset_t <- data.frame()
    dataset_t <- data.frame(GroupT, Dish, timeBef, timeAf, count = nSeeds, nCum, propCum)
    result <- rbind(result, dataset_t)
    }
  result
}

# ungroup_te <- function(data, counts) {
#   #This function organises a dataset to be submitted to survival analysis
#   #i.e. one row per each seed.
#
#   anName <- deparse(substitute(counts))
#   dfr <- data["anName" > 0,]
#   frequency <- subset(dfr, select = anName)[,1]
#   nr <- nrow(dfr)
#   df_surv <- dfr[rep(seq_len(nr), frequency), ]
#   row.names(df_surv) <- 1:(length(df_surv[,1]))
#   toRem <- which(names(df_surv) == anName)
#   df_surv[,-toRem]
# }

ungroup_te <- function(data, counts) {
  #This function organises a dataset to be submitted to survival analysis
  #i.e. one row per each seed.
  data <- as.data.frame(data)
  dfr <- data["anName" > 0,]
  tmp <- try(dplyr::select(dfr, {{ counts }}), silent = T)
  if(is(tmp, "try-error")){
    # Counts is not in data
    cond <- "external"
    frequency <- counts
  } else {
    # Counts is in data
    cond <- "internal"
    anName <- deparse(substitute(counts))
    frequency <- subset(dfr, select = anName)[,1]
  }
  nr <- nrow(dfr)
  df_surv <- dfr[rep(seq_len(nr), frequency), ]
  row.names(df_surv) <- 1:(length(df_surv[,1]))

  if(cond == "internal"){
    toRem <- which(names(df_surv) == anName)
    dfSurv <- df_surv[,-toRem]
  } else {
    dfSurv <- df_surv[,-length(df_surv[1,])]
  }
  dfSurv
}

group_te <- function(data) {
  tmp <- dplyr::group_by(data, dplyr::across())
  tmp <- dplyr::summarise(tmp, count = dplyr::n())
  as.data.frame(tmp)
}

makeDrm.drcte <- function(counts, treat, nViable, moniTimes) {
  # this function reshapes a common field-book for germination assays (WIDE format)
  # into the kind of dataset required by the function drm() in the drc package
  # (LONG GROUPED format)

  #Creating objects to store the results
  dati <- as.data.frame(counts)
  tempi <- c(moniTimes, Inf)

  dati[is.na(dati)] <- 0
  dati2 <- t(apply(dati, 1, cumsum))

  numPetri <- length(dati[,1]); numTimes <- length(tempi); numTesi<-length(treat[1,])
  numRecords <- numPetri*numTimes
  nSeeds <- numeric(numRecords); nCum <- numeric(numRecords); Prop <- numeric(numRecords); Dish <- numeric(numRecords)
  nGerm <- numeric(numRecords); nGermPetri <- numeric(numRecords)
  nGermPetri <- apply(dati, 1, sum)
  tempi2 <- c(0, tempi[1:length(tempi)-1])
  group <- data.frame()
  final.date <- max(tempi2)
  Petri <- rep(c(1:numPetri), each=numTimes)
  timeAf <- rep(tempi, numPetri)
  timeBef <- rep(tempi2, numPetri)

  cont <- 1
  for (j in 1:numPetri){ #For each dish
    for(i in 1:(numTimes-1)) { #For each sampling
      Dish[cont] <- j
      nSeeds[cont] <- dati[j,i]
      nGerm[cont] <- nGermPetri[j]
      nCum[cont] <- dati2[j,i]
      Prop[cont] <- dati2[j,i]/nViable[j]
      cont= cont+1
    }
    Dish[cont] <- j
    if(!is.na(nViable[j])){
      nSeeds[cont] <- nViable[j] - nGermPetri[j]
    } else {
      nSeeds[cont] <- 0
    }
    nGerm[cont] <- nGermPetri[j]
    nCum[cont] <- NA
    Prop[cont] <- NA
    cont=cont + 1
  }
  group <- treat[rep(row.names(treat), rep(numTimes, length(treat[,1]))), 1:length(treat[1,])]
  if(is.vector(group) | is.factor(group)){
    group <- data.frame(group)
    colnames(group) <- names(treat)
  }
  datiFin <- data.frame(group, Units = Dish, timeBef=timeBef, timeAf=timeAf, count=nSeeds, nCum=nCum, propCum=Prop)
  return(datiFin)
}


makeSurv.drcte <- function(counts, treat, nViable, moniTimes) {
  # this function reshapes a common field-book for germination assays (WIDE format)
  # into the kind of dataset required by the function Surv() in the drc package
  # (LONG UNGROUPED format, type = "interval2")
}
