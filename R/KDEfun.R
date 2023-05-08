# Kernel estimator for the distribution function
# Modified on 27/5/2021
KDE.fun <- Vectorize( function(x, start, end, count, h) {
  mod <- NPcdf(start, end, count)
  y <- mod$Type1plot$time
  w <- mod$Type1plot$pdf[-1]
  t <- y[-length(y)] + diff(y)/2
  sum(w * stats::pnorm((x - t)/h))
}, "x")

KDE <- function(bw = c("AMISE", "boot")){
  fct <- function(x, parm) {
    pred <- x[,1]; start <- x[,2]; end <- x[,3]; count <- x[,4]
    h <- parm[1]
    KDE.fun(x, start, end, count, h)
  }
  ssfct <- function(data){
    h <- 1
  }
  names <- c("h")
  text <- "Kernel estimator for the distribution function"
  bw <- match.arg(bw)

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text,
                     bw = bw)
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
#
