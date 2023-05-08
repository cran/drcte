# Function to predict the CDF for given times. It works
# with an icfit object with only one stratum (as obtained
# from a 'drmte' fit) and returns three types of preditions
# interpolation, right and left. The function was derived by
# modifying the 'getsurv' function in the 'interval' package
predictCDF <- function(obj, times, method = "interpolation"){

    ntimes <- length(times)
    method <- match.arg(method, c("interpolation","left","right", "midpoint"))
    p <- obj$pf
    L <- obj$intmap[1, ]
    R <- obj$intmap[2, ]

    ## function to get survival given L, R, and p
    ## does not use LRin attributes, so values exactly on intmap values may represent
    ## the limit approaching the intmap value

    Sout <- rep(NA, ntimes)
    mle <- rep(TRUE, ntimes)
    S <- c(1 - cumsum(p))

    # Core function: taken from getsurv in the interval package
    nonUMLE.func <- switch(method,
        interpolation = function(i, Time){
            if (i == 1){
                1 + ((Time - L[i])/(R[i] - L[i]))*(S[i]-1)
            } else if(R[i] != Inf) {
                S[i-1] + ((Time - L[i])/(R[i]-L[i]))*(S[i] - S[i-1])
              } else {
                S[i-1]
              }},
        right = function(i, Time){ ifelse(i <= 1, 1, S[i-1]) },
        left = function(i, Time){ S[i]},
        midpoint = function(i, Time){
          if (i == 1){
                ifelse(Time < (L[i] + R[i])/2, 1, S[i])
            } else {
                ifelse(Time < (L[i] + R[i])/2, S[i-1], S[i])
              }},
        )

    k<-length(p)
    for (i in 1:ntimes){
        if (any(times[i]==R)){ Sout[i]<-S[times[i]==R]
        } else if (times[i]<=L[1]){ Sout[i]<-1
        } else if (times[i]>=R[k]){ Sout[i]<-0
        } else {
            if  (times[i]>L[k]){
                Sout[i]<-nonUMLE.func(k,times[i])
                mle[i]<-FALSE
            } else {
                ## iLup is the index of the smallest L endpoint
                ## larger than times[i]
                iLup<-min((1:k)[L>=times[i]])
                if (R[iLup-1]<=times[i] | L[iLup]==times[i]) Sout[i]<-S[iLup-1]
                else {
                    Sout[i]<- nonUMLE.func(iLup-1,times[i])
                    mle[i]<-FALSE
                }
            }
       }
    }
    out <- list(cdf = 1 - Sout, S = Sout, times = times,
                unique.mle = mle, method = method)

    return(out)

}
