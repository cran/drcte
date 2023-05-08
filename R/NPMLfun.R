# NPMLE for time-to-event data
NPMLE <- function(){

  fct <- function(x, parm) {
    start <- x[,1]; end <- x[,2]; count <- x[,3]
    tmp <- NPcdf(start, end, count)
    tmp$SurvObjSum
  }
  ssfct <- function(data){
     k <- 0
     return(k)
  }
  names <- NULL
  text <- "NPML estimator for time-to-event data"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

# #Skeleton for DRC function ######################
# KDE.fun <- Vectorize( function(x, start, end, count, h) {
#   mod <- NPcdf(start, end, count)
#   y <- mod$Type1plot$time
#   w <- mod$Type1plot$pdf[-1]
#   t <- y[-length(y)] + diff(y)/2
#   sum(w * stats::pnorm((x - t)/h))
# }, "x")
#
# KDEfun <- function(){
#
#   fct <- function(x, parm) {
#     pred <- x[,1]; start <- x[,2]; end <- x[,3]; count <- x[,4]
#     h <- parm[1]
#     KDE.fun(x, start, end, count, h)
#   }
#   ssfct <- function(data){
#      # Self-starting code here
#      # x1 <- data[, 1]
#      # x2 <- data[, 2]
#      # x3 <- data[ ,3]
#      # k <- Kest(timeBef, timeAf, count, gplugin, type="N")
#     h <- 1
#   }
#   names <- c("h")
#   text <- "Kernel density estimator"
#
#   ## Returning the function with self starter and names
#   returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
#   class(returnList) <- "drcMean"
#   invisible(returnList)
# }

