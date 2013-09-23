##' Bayesian Quantile Regression with Polya Tree
##' for univariate outcome.
##'
##' Fit a quantile regression for univariate response
##' with block Gibbs sampling method
##'
##' @title Quantile Regression with Polya Tree
##' for univariate outcome
##' @param y [n] response
##' @param x [n, p] covariates matrix
##' @param mcmc mcmc parameters
##' @param prior priors
##' @param quan quantile requested. Multiple quantiles allowed.
##' @param method \itemize{
##' \item [normal] normal prior
##' \item [ss] spike-slab
##' }
##' @return an object of class \code{HeterPTlm}
##' @author Minzhao Liu, Mike Daniels
##' @export
HeterPTlmMH <- function(y, X, mcmc, prior = NULL, quan = 0.5,
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
  tuneinit <- mcmc$tuneinit
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
  tunebeta <- tuneinit[1:p]
  tunegamma <- tuneinit[(p+1):(p*2)]
  tunesigma <- tuneinit[2*p+1]
  tunealpha <- tuneinit[2*p+2]
  attgamma <- accgamma <- attbeta <- accbeta <- rep(0, p)
  attsigma <- attalpha <- accsigma <- accalpha <- 0
  propv <- sqrt(diag(solve(t(X)%*%X)))

#################################################

  ## first
  loglikeo <- ll(beta, gamma, sigma, alpha, mdzero, maxm, y, X)

  start <- proc.time()[3]

  ## MCMC roll
  for (iscan in 1:nscan) {

    ## beta
    attbeta <- attbeta + 1
    betac <- beta
    betac <- rnorm(p, beta, tunebeta*propv)

    ## gamma
    attgamma <- attgamma + 1
    gammac <- gamma
    gammac <- rnorm(p, gamma, tunegamma*propv)
    gammac[1] <- 1
    while (any(X%*%gammac < 0)) {
      gammac <- rnorm(p, gamma, tunegamma)
      gammac[1] <- 1
    }

    ## sigma
    attsigma <- attsigma + 1
    theta <- log(sigma)
    thetac <- rnorm(1, theta, tunesigma)
    logcgksigmac <- -theta
    logcgksigmao <- -thetac
    sigmac <- exp(thetac)

    ## alpha
    attalpha <- attalpha + 1
    theta <- log(alpha)
    thetac <- rnorm(1, theta, tunealpha)
    logcgkalphac <- -theta
    logcgkalphao <- -thetac
    alphac <- exp(thetac)

    ## deltabeta and deltagamma
    if (method == 'ss') {
      deltabeta <- rbinom(p, 1, pibeta*dnorm(beta, betapm, betapv)/(pibeta*dnorm(beta, betapm, betapv) + (1-pibeta)*dnorm(beta, 0, betapv/1000)))
      for (i in 2:p){
        deltagamma[i] <- rbinom(1, 1, pigamma[i]*dnorm(gamma[i], gammapm[i], gammapv[i])/(pigamma[i]*dnorm(gamma[i], gammapm[i], gammapv[i]) + (1-pigamma[i])*dnorm(gamma[i], 0, gammapv[i]/1000)))
      }

      ## pibeta and pigamma
      pibeta <- rbeta(p, 1 + deltabeta, 1 + (1 - deltabeta))
      pigamma <- rbeta(p, 1 + deltagamma, 1 + 1 - deltagamma)
    }

    ## log likelihood for new candidate
    loglikec <- ll(betac, gammac, sigmac, alphac, mdzero, maxm, y, X)

    ## priors
    logpriorc <- logprioro <- 0

    if (method == 'normal') {
      logpriorc <- logpriorc + sum(dnorm(betac, betapm, betapv, log = T))
      logprioro <- logprioro + sum(dnorm(beta, betapm, betapv, log = T))
      logpriorc <- logpriorc + sum(dnorm(gammac, gammapm, gammapv, log = T))
      logprioro <- logprioro + sum(dnorm(gamma, gammapm, gammapv, log = T))

    } else if (method == 'ss') {
      logpriorc <- logpriorc + sum(ifelse(deltabeta == 0, log((1 - pibeta)*dnorm(betac, 0, betapv/1000)), log(pibeta) + dnorm(betac, betapm, betapv, log = T)))
      logprioro <- logprioro + sum(ifelse(deltabeta == 0, log((1 - pibeta)*dnorm(beta, 0, betapv/1000)), log(pibeta) + dnorm(beta, betapm, betapv, log = T)))
      logpriorc <- logpriorc + sum(ifelse(deltagamma==0, log((1-pigamma)*dnorm(gammac,0,gammapv/1000)), log(pigamma) + dnorm(gammac, gammapm, gammapv, log = T)))
      logprioro <- logprioro + sum(ifelse(deltagamma==0,log((1-pigamma)*dnorm(gamma, 0, gammapv/1000)), log(pigamma) + dnorm(gamma, gammapm, gammapv, log = T)))
    }

    logpriorc <- logpriorc + dgamma(sigmac, a/2, b/2, log = T)
    logprioro <- logprioro + dgamma(sigma, a/2, b/2, log = T)
    logpriorc <- logpriorc + dgamma(alphac, a/2, b/2, log = T)
    logprioro <- logprioro + dgamma(alpha, a/2, b/2, log = T)

    ## additional likelihood for gamma
    loglikeaddc <- -sum(log(X%*%gammac))
    loglikeaddo <- -sum(log(X%*%gamma))

    ## ratio
    ratio <- loglikec + logpriorc - loglikeo - logprioro +
      logcgksigmac - logcgksigmao + logcgkalphac - logcgkalphao + loglikeaddc - loglikeaddo

    if (log(runif(1)) <= ratio) {
      accbeta <- accbeta + 1
      accgamma <- accgamma + 1
      accsigma <- accsigma + 1
      accalpha <- accalpha + 1
      beta <- betac
      gamma <- gammac
      sigma <- sigmac
      alpha <- alphac
      loglikeo <- loglikec
    }

    ## TUNE
    if (attbeta[1] >= 100 & iscan < nburn) {
      tunerate <- exp(min(0.01, 1/sqrt(iscan/100)))
      tunegamma <- tunegamma*ifelse(accgamma/attgamma > arate,
                                    tunerate, 1/tunerate)
      tunebeta <- tunebeta*ifelse(accbeta/attbeta > arate,
                                    tunerate, 1/tunerate)
      tunesigma <- tunesigma*ifelse(accsigma/attsigma > arate, tunerate, 1/tunerate)
      tunealpha <- tunealpha*ifelse(accalpha/attalpha > arate, tunerate, 1/tunerate)
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
