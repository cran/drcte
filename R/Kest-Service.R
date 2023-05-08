# Service functions (translated from C++) ###################
# Taken from binnednp package, Barreiro-Ures et al., 2019
# Ecology and Evolution, 9, 10903-10915.
# ##########################################################
Fg_c <- function(x, w, t, g){
  arg <- (x - t) / g
  suma <- sum( w * pnorm(arg) )
  return(suma)
}

biasFh_c <- function(x, n, t, w, p, g, h){
  arg <- (x - t) / h
  suma <- sum( pnorm(arg) * p )
  suma <- suma - Fg_c(x, w, t, g)
  return(suma)
}

varFh_c <- function(x, n, t, p, h){
	k <- length(t)
	invn <- 1.0/n
	suma1 <- 0
	suma2 <- 0
	for(i in 1:k){
		arg <- (x - t[i]) / h
		pnormarg <- pnorm(arg)
		suma1 <- suma1 + pnormarg * pnormarg * p[i] * (1 - p[i])
		l <- i + 1
		if(l <= k){
		for(j in (i+1):k){
			arg2 <- (x - t[j]) / h
			pnormarg2 <- pnorm(arg2)
			suma2 <- suma2 + pnormarg * pnormarg2 * p[i] * p[j]
		} }
	}
	invn * suma1 - 2.0 * invn * suma2
}

mise_Fh_c <- function(h, n, t, w, p, g, lgrid, lim1, lim2){
	cte <- (lim2 - lim1) / lgrid
	biasFhlim1 <- biasFh_c(lim1, n, t, w, p, g, h)
	biasFhlim2 <- biasFh_c(lim2, n, t, w, p, g, h)
	varFhlim1 <- varFh_c(lim1, n, t, p, h)
	varFhlim2 = varFh_c(lim2, n, t, p, h)
	suma <- 0.5 * (biasFhlim1 * biasFhlim1 + varFhlim1 + biasFhlim2 * biasFhlim2 + varFhlim2)
	for(i in 1:(lgrid-2)){
		xi <- lim1 + i * cte
		biasFhxi <- biasFh_c(xi, n, t, w, p, g, h)
		varFhxi <- varFh_c(xi, n, t, p, h)
		suma <- suma + biasFhxi * biasFhxi + varFhxi
	}
	suma <- suma * cte
	return(suma)
}

boot_bw_dist_c <- function(nit, h0, h1, rho, n, t, w, p, g, lgrid, lim1, lim2){
   mises <- c()
   hseq <- c()
	 j <- 2
	 newh0 <- h0; newh1 <- h1; newrho <- rho
	 for(it in 1:nit){
	   
		for(i in 1:5){
		  # print(i)
			hseq[i] = newh0 * newrho^(i - 1)
			mises[i] = mise_Fh_c(hseq[i], n, t, w, p, g, lgrid, lim1, lim2)
		}
	  
		j <- which.min(mises)
		# j = min(mises) - mises[1]
    # print(j)
		if(j == 1) {
			newh0 = hseq[1] / newrho;
			newh1 = hseq[2];
			newrho = exp((log(newh1) - log(newh0)) / 4);
		} else {
			if(j == 5) {
				newh0 = hseq[4];
				newh1 = hseq[5] * newrho;
				exp((log(newh1) - log(newh0)) / 4);
			} else {
				newh0 = hseq[j - 1];
				newh1 = hseq[j + 1];
				newrho = exp((log(newh1) - log(newh0)) / 4);
			}
		}
	 }
	 # print(mises, digits = 20)
	 # print(hseq)
	 # print(j)
	return(hseq[j])
}

dnorm_cpp_c <- function(x){
	1.0 / sqrt(2.0 * pi) * exp(-0.5 * x * x)
}

pnorm_cpp_fun_c <- function(x){
  1 - pnorm( - x )
}

dnorm_d1_cpp_c <- function(x){
	-x * dnorm_cpp_c(x)
}


Fh_combt_c <- function(x, t, w, h){
  dif = (x - t) / h
  sum (w * pnorm_cpp_fun_c(dif))
}


fh_combt_d1_c <- function(x, t, w, h){
		dif = (x - t) / h
		suma <- sum( w * dnorm_d1_cpp_c(dif) )
	  suma <- suma/ (h * h)
	  return(suma)
}

slope_cpp_c <- function(t, w, h, lim1, lim2, lgrid){
	cte <- (lim2 - lim1) / lgrid
	fhtd1lim1 <- fh_combt_d1_c(lim1, t, w, h)
	fhtd1lim2 <- fh_combt_d1_c(lim2, t, w, h)
	suma <- 0.5 * (fhtd1lim1 * fhtd1lim1 + fhtd1lim2 * fhtd1lim2)
	for(i in 1:(lgrid-1))
	{
		xi <- lim1 + i * cte
		fhtd1xi <- fh_combt_d1_c(xi, t, w, h)
		suma <- suma + fhtd1xi * fhtd1xi
	}
	suma <- suma * cte
	return(suma)
}

dicoto_lambda_dist_c <- function(lambda, nith, h0, h1, rho, emp,
                                 comby, combt, combw, lim1, lim2){

  output <- rep(NA, 2)
  hseq <- rep(NA, 5)
  slope <- rep(NA, 5)
  objective <- rep(NA, 5)
  newh0 <- h0
  newh1 <- h1
  newrho <- rho
  kcomby <- length(comby)

  for(it in 1:nith){
    for(i in 1:5){
      auxhi <- newh0 * newrho^(i - 1)
      hseq[i] <- auxhi

      suma = 0
      dif = 0
      for(ii in 1:kcomby){
        auxcombyii <- comby[ii]
        dif <- emp[ii] - Fh_combt_c(auxcombyii, combt, combw, auxhi)
        suma = suma + dif * dif
      }

      auxslopei <- slope_cpp_c(combt, combw, auxhi, lim1, lim2, 100)
      slope[i] <- auxslopei
      objective[i] <- suma + lambda * auxslopei

    }

    minind <- which.min(objective)
    gboot <- hseq[minind]
    output[1] <- gboot
    output[2] <- slope[minind]

    if(minind == 1)
    {
      newh1 <- hseq[2]
      newh0 <- gboot / newrho
      newrho <- exp((log(newh1)-log(newh0))/4)
    } else {
      if(minind == 5)
      {
        newh0 <- hseq[4]
        newh1 <- gboot * newrho;
        newrho <- exp((log(newh1)-log(newh0))/4);
      } else {
        newh0 <- hseq[minind - 1];
        newh1 <- hseq[minind + 1];
        newrho <- exp((log(newh1)-log(newh0))/4);
      }
    }

  }
  return(output)
}

zeta_hist_p_dist_c <- function(emp, combt, comby, combw, Af1_mixt,
                               l0, l1, h0, h1, lrho, rho, nitlambda,
                               nith, lim1, lim2) {
  # print(list(emp, combt, comby, combw, Af1_mixt,
  #                              l0, l1, h0, h1, lrho, rho, nitlambda,
  #                              nith, lim1, lim2))
  lseq <- rep(NA, 5)
  dicoto_return <- rep(NA, 2)
  ldist <- rep(NA, 5)
  bw <- rep(NA, 5)
  
  newl0 <- l0
  newl1 <- l1
  newlrho <- lrho
  gboot <- 0

  for(itlambda in 1:nitlambda){
    for(i in 1:5) {
      lseq[i] <- newl0 * newlrho^(i-1)
    }
    for(li in 1:5)
    {
      lambda <- lseq[li]
      dicoto_return <- dicoto_lambda_dist_c(lambda, nith, h0, h1, rho, emp, comby, combt, combw, lim1, lim2)
      bw[li] <- dicoto_return[1]
      ldist[li] <- abs(dicoto_return[2] - Af1_mixt)
    }
    
    lminind <- which.min(ldist)
    gboot <- bw[lminind]
    
    if(lminind == 1)
    {
      newl1 <- lseq[2]
      newl0 <- lseq[1] / newlrho
      newlrho <- exp((log(newl1)-log(newl0))/4);
    } else {
      if(lminind == 5)
      {
        newl0 <- lseq[4]
        newl1 <- lseq[5] * newlrho
        newlrho <- exp((log(newl1)-log(newl0))/4)
      } else {
        newl0 <- lseq[lminind - 1];
        newl1 <- lseq[lminind + 1];
        newlrho <- exp((log(newl1)-log(newl0))/4);
      }
    }

  }
  return(gboot)
}

