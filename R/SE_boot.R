boot.drcte <- function(obj, B = 499, units = NULL){
  # This is working for a drcte time-to-event object
  numFac <- length(attr(terms(as.formula(obj$call$formula)), "term.labels"))
  counts <- obj$data[,numFac + 1]
  L <- rep(obj$data[,1], counts)
  R <- rep(obj$data[,2], counts)
  numTreat <- length(levels(factor(obj$data[,numFac + 2])))
  if(numTreat == 1) {
    df <- data.frame(L, R)
  } else if(numTreat > 1) {
    treats <- rep(obj$data[, numFac + 3], counts)
    df <- data.frame(L, R, treats)
  }
  res <- matrix(NA, nrow = B, ncol = length(coef(obj)))
  for(i in 1:B){
    message("\r Resampling:", i)
    if(numTreat == 1){
      newSample <- resample.cens(df, replace = T) # resampling cases
      newSample$newCounts <- rep(1, length(newSample[,1]))
      newMod <- try(update(obj, formula = newCounts ~ L + R,
                          data = newSample), silent = T)
    } else if(numTreat > 1){
      dfList <- by(df, treats, I)
      newList <- lapply(dfList, function(x) resample.cens(x, replace = T)) # Resampling cases
      newSample <- do.call(rbind, newList)
      newSample$newCounts <- rep(1, length(newSample[,1]))
      newMod <- try(update(obj, formula = newCounts ~ L + R,
                          data = newSample, curveid = treats), silent = T)
    }
    if(any(class(newMod) == "try-error")){
      res[,i] <- NA
    } else {
      #if(all(coef(newMod) > 0))
      res[i,] <- coef(newMod)
    }
  }
  vals <- apply(res, 2, mean, na.rm = T)
  ses <- apply(res, 2, sd, na.rm = T)
  dfres <- cbind(vals, ses)

  class(obj) <- "drc"
  sumObj <- summary(obj)
  sumObj$coefficients <- data.frame("Estimate" = sumObj$coefficients[,1],
             "Bootstrap Mean" = vals,
             "Bootstrap SE" = ses)
  sumObj$varMat <- NA
  sumObj$robust <- paste("Bootstrap resampling (B = ", B, ")", sep = "")
  sumObj$resVar <- NULL
  sumObj$resamples <- res
  class(sumObj) <- c("summary.drcte", "summary.drc")
  return(sumObj)
}

