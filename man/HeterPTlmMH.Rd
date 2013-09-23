\name{HeterPTlmMH}
\alias{HeterPTlmMH}
\title{Quantile Regression with Polya Tree
for univariate outcome}
\usage{
  HeterPTlmMH(y, X, mcmc, prior = NULL, quan = 0.5,
    method = "normal", den = FALSE)
}
\arguments{
  \item{y}{[n] response}

  \item{x}{[n, p] covariates matrix}

  \item{mcmc}{mcmc parameters}

  \item{prior}{priors}

  \item{quan}{quantile requested. Multiple quantiles
  allowed.}

  \item{method}{\itemize{ \item [normal] normal prior \item
  [ss] spike-slab }}
}
\value{
  an object of class \code{HeterPTlm}
}
\description{
  Bayesian Quantile Regression with Polya Tree for
  univariate outcome.
}
\details{
  Fit a quantile regression for univariate response with
  block Gibbs sampling method
}
\author{
  Minzhao Liu, Mike Daniels
}

