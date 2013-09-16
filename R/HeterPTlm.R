##' Bayesian Quantile Regression with Polya Tree
##' for univariate outcome.
##'
##' Fit a quantile regression for univariate response.
##'
##' @title Quantile Regression with Polya Tree
##' for univariate outcome
##' @param y [n] response
##' @param x [n, p] covariates matrix
##' @param mcmc mcmc parameters
##' @param prior priors
##' @param quan quantile requested. Multiple quantiles allowed.
##' @param method \itemize{
##' \item [normal] defaul method. Use normal priors
##' \item [ss] spike-slab priors
##' }
##' @return an object of class \code{HeterPTlm}
##' @author Minzhao Liu, Mike Daniels
##' @export
HeterPTlm <- function(y, X, mcmc, prior = NULL, quan = 0.5,
                      method = "normal",
                      den = FALSE){

  ## DATA
  nrec <- length(y)
  p <- dim(X)[2]

  ## PT
  ## mdzero = 0 => median fix at 0
  if (is.null(prior$maxm)) {
    maxm <- floor(log(nrec)/log(2)) }  else maxm <- prior$maxm
  if (is.null(prior$mdzero)) {mdzero <- 0}  else  mdzero <- prior$mdzero

  ## PRIOR
  if (is.null(prior$betapm)){
    betapm <- solve(t(X)%*%X)%*%t(X)%*%y
    res <- y - X%*%betapm
    betapv <- sqrt(diag(solve(t(X)%*%X)*sum(res^2)/(nrec - p)))
    gammapm <- rep(0, p)
    gammapv <- rep(100, p)
    a <- b <- 1
  } else {
    betapm <- prior$betapm
    betapv <- prior$betapv
    gammapm <- prior$gammapm
    gammapv <- prior$gammapv
    a <- prior$a
    b <- prior$b
  }

  ## MCMC
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nburn <- mcmc$nburn
  ndisp <- mcmc$ndisp
  arate <- mcmc$arate
  mcmc <- c(nburn, nskip, nsave, ndisp)

  ## QUAN
  nquan <- length(quan)

  ## SAVE
  betasave <- gammasave <- matrix(0, nsave, p)
  sigmasave <- alphasave <- rep(0, nsave)
  quansave <- matrix(0, nsave, nquan)
  deltabetasave <- deltagammasave <- matrix(0, nsave, p)

  ## grid
  ngrid <- 200
  f <- rep(0, ngrid)

  ## INITIAL
  beta <- as.vector(solve(t(X)%*%X)%*%t(X)%*%y)
  gamma <- c(1, rep(0, p-1))
  sigma <- 1
  alpha <- 1
  v <- as.vector((y-X%*%beta)/(X%*%gamma))

  isave <- 0
  skipcount <- 0
  dispcount <- 0
  nscan <- nburn + nskip * nsave

  ## initial for spike-slab
  pibeta <- rbeta(p, 1, 1)
  pigamma <- rbeta(p, 1, 1)

  ## deltabeta, deltagamma = 1 : from slab prior
  deltabeta <- rep(1, p)
  deltagamma <- rep(1, p)

  ## new grid

  left <- min(v) - 0.5*sd(v)
  right <- max(v) + 0.5*sd(v)
  grid <- seq(left, right, length=ngrid)

  ## WORKING
  whicho <- whichn <- rep(0, nrec)

  ## TUNE
  tunegamma <- tunebeta <- rep(0.3, p)
  tunesigma <- 0.3
  tunealpha <- 0.3
  attgamma <- accgamma <- attbeta <- accbeta <- rep(0, p)
  attsigma <- attalpha <- accsigma <- accalpha <- 0

#################################################

  ## first
  loglikeo <- ll(beta, gamma, sigma, alpha, mdzero, maxm, y, X)

  start <- proc.time()[3]

  ## MCMC roll
  for (iscan in 1:nscan) {
    ## beta
    for (i in 1:p){
      attbeta[i] <- attbeta[i] + 1
      betac <- beta
      betac[i] <- rnorm(1, beta[i], tunebeta[i])

      if (method == 'normal') {
        logpriorc <- dnorm(betac[i], betapm[i], betapv[i], log = T)
        logprioro <- dnorm(beta[i], betapm[i], betapv[i], log = T)
      } else if (method == 'ss') {
        logpriorc <- ifelse(deltabeta[i] == 0, log((1 - pibeta[i])*dnorm(betac[i], 0, betapv[i]/1000)), log(pibeta[i]) + dnorm(betac[i], betapm[i], betapv[i], log = T))
        logprioro <- ifelse(deltabeta[i] == 0, log((1 - pibeta[i])*dnorm(beta[i], 0, betapv[i]/1000)), log(pibeta[i]) + dnorm(beta[i], betapm[i], betapv[i], log = T))
      }

      loglikec <- ll(betac, gamma, sigma, alpha, mdzero, maxm, y, X)

      ratio <- loglikec + logpriorc - loglikeo - logprioro

      if (log(runif(1)) <= ratio) {
        accbeta[i] <- accbeta[i] + 1
        loglikeo <- loglikec
        beta <- betac
      }
    }

    ## gamma
    for (i in 2:p){
      attgamma[i] = attgamma[i] + 1
      gammac <- gamma
      gammac[i] <- rnorm(1, gamma[i], tunegamma[i])
      while (any(X%*%gammac < 0)) {
        gammac[i] <- rnorm(1, gamma[i], tunegamma[i])
      }

      if (method == 'normal') {
        logpriorc <- dnorm(gammac[i], gammapm[i], gammapv[i], log = T)
        logprioro <- dnorm(gamma[i], gammapm[i], gammapv[i], log = T)
      } else if (method == 'ss') {
        logpriorc <- ifelse(deltagamma[i] == 0, log((1 - pigamma[i])*dnorm(gammac[i], 0, gammapv[i]/1000)), log(pigamma[i]) + dnorm(gammac[i], gammapm[i], gammapv[i], log = T))
        logprioro <- ifelse(deltagamma[i] == 0, log((1 - pigamma[i])*dnorm(gamma[i], 0, gammapv[i]/1000)), log(pigamma[i]) + dnorm(gamma[i], gammapm[i], gammapv[i], log = T))
      }

      loglikec <- ll(beta, gammac, sigma, alpha, mdzero, maxm, y, X)

      loglikeaddc <- -sum(log(X%*%gammac))
      loglikeaddo <- -sum(log(X%*%gamma))

      ratio <- loglikec + logpriorc - loglikeo - logprioro + loglikeaddc - loglikeaddo

      if (log(runif(1)) <= ratio) {
        accgamma[i] <- accgamma[i] + 1
        loglikeo <- loglikec
        gamma <- gammac
      }
    }

    ## sigma
    attsigma <- attsigma + 1
    theta <- log(sigma)
    thetac <- rnorm(1, theta, tunesigma)
    logcgkc <- -theta
    logcgko <- -thetac
    sigmac <- exp(thetac)

    loglikec <- ll(beta, gamma, sigmac, alpha, mdzero, maxm, y, X)

    logpriorc <- dgamma(sigmac, a/2, b/2, log = T)
    logprioro <- dgamma(sigma, a/2, b/2, log = T)

    ratio <- loglikec + logpriorc - loglikeo - logprioro + logcgkc - logcgko

    if (log(runif(1)) <= ratio) {
      accsigma <- accsigma + 1
      loglikeo <- loglikec
      sigma <- sigmac
    }

    ## alpha
    attalpha <- attalpha + 1
    theta <- log(alpha)
    thetac <- rnorm(1, theta, tunealpha)
    logcgkc <- -theta
    logcgko <- -thetac
    alphac <- exp(thetac)

    loglikec <- ll(beta, gamma, sigma, alphac, mdzero, maxm, y, X)

    logpriorc <- dgamma(alphac, a/2, b/2, log = T)
    logprioro <- dgamma(alpha, a/2, b/2, log = T)

    ratio <- loglikec + logpriorc - loglikeo - logprioro + logcgkc - logcgko

    if (log(runif(1)) <= ratio) {
      accalpha <- accalpha + 1
      loglikeo <- loglikec
      alpha <- alphac
    }

    ## deltabeta and deltagamma
    if (method == 'ss') {
      for (i in 1:p){
        deltabeta[i] <- rbinom(1, 1, pibeta[i]*dnorm(beta[i], betapm[i], betapv[i])/(pibeta[i]*dnorm(beta[i], betapm[i], betapv[i]) + (1-pibeta[i])*dnorm(beta[i], 0, betapv[i]/1000)))
      }
      for (i in 2:p){
        deltagamma[i] <- rbinom(1, 1, pigamma[i]*dnorm(gamma[i], gammapm[i], gammapv[i])/(pigamma[i]*dnorm(gamma[i], gammapm[i], gammapv[i]) + (1-pigamma[i])*dnorm(gamma[i], 0, gammapv[i]/1000)))
      }

      ## pibeta and pigamma
      pibeta <- rbeta(p, 1 + deltabeta, 1 + (1 - deltabeta))
      pigamma <- rbeta(p, 1 + deltagamma, 1 + 1 - deltagamma)
    }

    ## TUNE
    if (attbeta[1] >= 100 & iscan < nburn) {
      tunegamma <- tunegamma*ifelse(accgamma/attgamma > arate,
                                    2, 0.5)
      tunebeta <- tunebeta*ifelse(accbeta/attbeta > arate,
                                    2, 0.5)
      tunesigma <- tunesigma*ifelse(accsigma/attsigma > arate, 2, 0.5)
      tunealpha <- tunealpha*ifelse(accalpha/attalpha > arate, 2, 0.5)
      attgamma <- accgamma <- attbeta <- accbeta <- rep(0, p)
      attsigma <- accsigma <- attalpha <- accalpha <- 0
    }

    ## save
    if (iscan > nburn) {
      skipcount = skipcount + 1
      if (skipcount >= nskip) {
        isave <- isave + 1
        dispcount <- dispcount + 1
        betasave[isave, ] <- beta
        gammasave[isave, ] <- gamma
        sigmasave[isave] <- sigma
        alphasave[isave] <- alpha
        deltabetasave[isave, ] <- deltabeta
        deltagammasave[isave, ] <- deltagamma
        quansave[isave, ] <- PostQuantile(beta, gamma, sigma, alpha, mdzero, maxm, y, X, quan)
        if (den) {
          f <- f + PostDensity(grid, beta, gamma, sigma, alpha, mdzero, maxm, y, X)
        }
        skipcount <- 0
        if (dispcount >= ndisp) {
          dispcount <- 0
          cat(isave, 'out of', nsave, proc.time()[3] - start, '\n')
        }
      }
    }
  }

  ans <- list(betasave=betasave,
              gammasave=gammasave,
              sigmasave=sigmasave,
              alphasave=alphasave,
              deltabetasave = deltabetasave,
              deltagammasave = deltagammasave,
              quansave=quansave,
              dens=den,
              f = f/nsave,
              grid=grid,
              mcmc=mcmc,
              prior=prior,
              n=nrec,
              p=p,
              quan=quan,
              y=y,
              X=X,
              ## ratesave=ratesave,
              tune = list(beta = tunebeta, gamma = tunegamma,
                sigma = tunesigma, alpha = tunealpha),
              ## hetersave=hetersave,
              method = method
              )

  class(ans) <- "HeterPTlm"

  return(ans)
}

##' @rdname HeterPTlm
##' @method coef HeterPTlm
##' @S3method coef HeterPTlm
coef.HeterPTlm <- function(mod, ...){
  nquan <- length(mod$quan)
  betasave <- mod$betasave
  gammasave <- mod$gammasave
  quansave <- mod$quansave
  quan <- mod$quan
  nsave <- dim(betasave)[1]

  deltabetaprop <- 1 - apply(mod$deltabetasave, 2, sum)/nsave
  deltagammaprop <- 1 - apply(mod$deltagammasave, 2, sum)/nsave

  betaMedian <- apply(betasave, 2, median)
  gammaMedian <- apply(gammasave, 2, median)
  quanMedian <- apply(quansave,2,median)
  betatauMedian <- matrix(0, nquan, mod$p)
  betatauCIlbd <- matrix(0, nquan, mod$p)
  betatauCIubd <- matrix(0, nquan, mod$p)
  rownames(betatauMedian) <- mod$quan
  rownames(betatauCIlbd) <- mod$quan
  rownames(betatauCIubd) <- mod$quan
  for (i in 1:nquan) {
    tmp <- betasave + gammasave*as.numeric(quansave[,i])
    betatauMedian[i, ] <- apply(tmp, 2, median)
    betatauCIlbd[i, ] <- apply(tmp, 2, function(x) quantile(x, probs = c(0.025)))
    betatauCIubd[i, ] <- apply(tmp, 2, function(x) quantile(x, probs = c(0.975)))
  }

  betaMean <- apply(betasave, 2, mean)
  gammaMean <- apply(gammasave, 2, mean)
  quanMean <- apply(quansave,2,mean)
  betatauMean <- matrix(0, nquan, mod$p)
  rownames(betatauMean) <- mod$quan
  for (i in 1:nquan) {
    tmp <- betasave + gammasave*as.numeric(quansave[,i])
    betatauMean[i, ] <- apply(tmp, 2, mean)
  }

  return(list(betaMedian = betaMedian, gammaMedian = gammaMedian,
              quanMedian = quanMedian, betatauMedian = betatauMedian,
              betaMean = betaMean, gammaMean = gammaMean,
              quanMean = quanMean, betatauMean = betatauMean,
              betatauCIlbd = betatauCIlbd, betatauCIubd = betatauCIubd,
              deltabetaprop = deltabetaprop, deltagammaprop = deltagammaprop))
}

##' @rdname HeterPTlm
##' @method plot HeterPTlm
##' @S3method plot HeterPTlm
plot.HeterPTlm <- function(obj, ...){
  for (i in 1:obj$p){
    title1 <- paste("Trace of beta" , i-1, sep=" ")
    title2 <- paste("Density of beta", i-1, sep=" ")
    plot(obj$betasave[,i], type='l', main=title1, xlab="MCMC scan", ylab=" ")
    plot(density(obj$betasave[,i]), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')
  }

  for (i in 2:obj$p){
    title1 <- paste("Trace of gamma" , i-1, sep=" ")
    title2 <- paste("Density of gamma", i-1, sep=" ")
    plot(obj$gammasave[,i], type='l', main=title1, xlab="MCMC scan", ylab=" ")
    plot(density(obj$gammasave[,i]), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')
  }

  title1 <- "Trace of sigma"
  title2 <- "Density of sigma"
  plot(obj$sigmasave, typ='l', main=title1, xlab="MCMC scan", ylab=" ")
  plot(density(obj$sigmasave), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')

  title1 <- "Trace of alpha"
  title2 <- "Density of alpha"
  plot(obj$alphasave, typ='l', main=title1, xlab="MCMC scan", ylab=" ")
  plot(density(obj$alphasave), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')

  if (obj$den) {
    title1 <- "Predictive Error Density"
    plot(obj$grid, obj$f, ylab="density", main=title1, type='l', lwd=2, xlab="values")
  }
}

##' @rdname HeterPTlm
##' @method summary HeterPTlm
##' @S3method summary HeterPTlm
summary.HeterPTlm <- function(mod, ...){
  n <- mod$n
  quan <- mod$quan
  p <- mod$p
  method <- mod$method
  cat('Number of observations: ', n, '\n')
  cat('Quantile: ', quan, '\n')
  cat('Priors', method, '\n')
  cat('Quantile regression coefficients (Median): \n')
  print(coef(mod)$betatauMedian)
  cat('Quantile regression coefficients (Mean): \n')
  print(coef(mod)$betatauMean)
}

##' ll
##'
##' Log Likelihood given beta, gamma, sigma, alpha, mdzero
##'
##' @param beta
##' @param gamma
##' @param sigma baseline normal(0, sd = sigma)
##' @param alpha PT parameter
##' @param mdzero
##' @param maxm
##' @param y
##' @param X
##' @return log likelihood
##' @author Minzhao Liu, Mike Daniels
##' @useDynLib bqrpt
ll <- function(beta, gamma, sigma, alpha, mdzero, maxm, y, X){
  res <- (y - X%*%beta)/(X%*%gamma)
  n <- length(y)
  ans <- 0
  whicho <- whichn <- rep(0, n)
  ans <- .Fortran('loglik_unippt',
                  nsubject = as.integer(n),
                  mdzero = as.integer(mdzero),
                  maxm = as.integer(maxm),
                  alpha = as.double(alpha),
                  mu = as.double(0),
                  sigma = as.double(sigma^2),
                  b = as.double(res),
                  whicho=as.integer(whicho),
                  whichn=as.integer(whichn),
                  ans = as.double(ans))$ans
  return(ans)
}

##' PostQuantiles
##'
##' calculate the posterior quantiles
##'
##' @param beta
##' @param gamma
##' @param sigma
##' @param alpha
##' @param mdzero
##' @param maxm
##' @param y
##' @param X
##' @return postquantile
##' @author Minzhao Liu, Mike Daniels
##' @useDynLib bqrpt
PostQuantile <- function(beta, gamma, sigma, alpha, mdzero, maxm, y, X, quan){
  n <- length(y)
  res <- (y - X%*%beta)/(X%*%gamma)
  nquan <- length(quan)
  ans <- rep(0, nquan)
  ans <- .Fortran('postquantile',
                  n = as.integer(n),
                  res = as.double(res),
                  sigma2 = as.double(sigma^2),
                  alpha = as.double(alpha),
                  maxm = as.integer(maxm),
                  quan = as.double(quan),
                  ans = as.double(ans),
                  nquan = as.integer(nquan))$ans
  return(ans)
}

##' PostDensity
##'
##' postdensity
##'
##' @param grid
##' @param beta
##' @param gamma
##' @param sigma
##' @param alpha
##' @param mdzero
##' @param maxm
##' @param y
##' @param X
##' @return postdensity
##' @author Minzhao Liu, Mike Daniels
PostDensity <- function(grid, beta, gamma, sigma, alpha, mdzero, maxm, y, X){
  n <- length(y)
  res <- as.vector((y - X%*%beta)/(X%*%gamma))
  ngrid <- length(grid)
  like <- rep(0, ngrid)
  whicho <- whichn <- rep(0, n)

  like <- .Fortran('postdensity',
                   grid = as.double(grid),
                   maxm = as.integer(maxm),
                   mdzero = as.integer(mdzero),
                   n = as.integer(n),
                   alpha = as.double(alpha),
                   mu = as.double(0),
                   sigma2 = as.double(sigma^2),
                   res = as.double(res),
                   whicho = as.integer(whicho),
                   whichn = as.integer(whichn),
                   like = as.double(like),
                   ngrid = as.integer(ngrid))$like
  return(like)

}
