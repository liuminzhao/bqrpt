##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Multivariate Bayesian Linear quantile regression with
##' Mixture of Polya tree priors , here d=2 
##' @param y 
##' @param X 
##' @param nsub 
##' @param mcmc 
##' @param prior 
##' @param quan 
##' @return class HeterPTlm
##' @author Minzhao Liu
HeterPTlmm <- function(y, X, nsub, mcmc, prior, quan){

  dyn.load('heterptlmm.so')
  
  # DATA
  nrec <- length(y)
  q <- nrec/nsub
  p <- dim(X)[2]

  # PT
  if (is.null(prior$maxm)) {
    maxm <- floor(log(nrec)/log(2)) }  else maxm <- prior$maxm
  if (is.null(prior$mdzero)) {mdzero <- 1}  else  mdzero <- prior$mdzero

  # PRIOR
  betapm <- prior$betapm
  betapv <- prior$betapv
  gammapm <- prior$gammapm
  gammapv <- prior$gammapv
  tau <- prior$tau
  a0b0 <- prior$a0b0

  # MCMC
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nburn <- mcmc$nburn
  ndisp <- mcmc$ndisp
  arate <- mcmc$arate

  # QUAN
  quan <- quan
  nquan <- length(quan)
  
  # SAVE
  betasave<- matrix(0, nsave, p)
  gammasave <- matrix(0, nsave, p)
  sigmasave <- matrix(0, nsave, 3) # sigma1, sigma2, rho
  alphasave <- rep(0, nsave)
    #   quansave not used yet
  
  # GRID
  # not used yet
  
  # INITIAL
  beta <- as.vector(solve(t(X)%*%X)%*%t(X)%*%y)
  gamma <- c(1, rep(0, p-1))
  sigmavec <- c(1, 1, 0)
  Sigma <- diag(2)
  alpha <- 1
  v <- as.vector((y-X%*%beta)/(X%*%gamma))
  propv <- solve(t(X)%*%X)
  # WORKING
  whicho <- whichn <- rep(0, nsub)
  
  # DEBUG, TUNING
   # not used yet

  ####################################
  # MCMC USING FORTRAN

  foo <- .Fortran("heterptlmm",
              maxm=as.integer(maxm),
              mdzero=as.integer(mdzero),
              nrec=as.integer(nrec),
              nsub=as.integer(nsub),
              p=as.integer(p),
              q=as.integer(q),
              x=as.double(X),
              y=as.double(y),
              betapm=as.double(betapm),
              betapv=as.double(betapv),
              tau=as.double(tau),
              a0b0=as.double(a0b0),
              gammapm=as.double(gammapm),
              gammapv=as.double(gammapv),
              nburn=as.integer(nburn),
              nskip=as.integer(nskip),
              nsave=as.integer(nsave),
              ndisp=as.integer(ndisp),
              betasave=as.double(betasave),
              gammasave=as.double(gammasave),
              sigmasave=as.double(sigmasave),
              alphasave=as.double(alphasave),
              beta=as.double(beta),
              gamma=as.double(gamma),
              alpha=as.double(alpha),
              Sigma=as.double(Sigma),
              propv=as.double(propv),
              arate=as.double(arate)
              )

  ####################################

  betasave <- matrix(foo$betasave, nsave, p)
  gammasave <- matrix(foo$gammasave, nsave, p)
  alphasave <- foo$alphasave
  sigmasave <- matrix(foo$sigmasave, nsave, 3)

  coef <- list(beta=apply(betasave, 2, median),
               gamma=apply(gammasave,2,median),
               alpha=median(alphasave),
               sigma=apply(sigmasave, 2, median)
               )

  z <- list(coef=coef,
            betasave=betasave,
            gammasave=gammasave,
            alphasave=alphasave,
            sigmasave=sigmasave
            )

  class(z) <- "HeterPTlmm"

  return(z)
  
}

