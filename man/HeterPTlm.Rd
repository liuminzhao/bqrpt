\name{HeterPTlm}
\alias{coef.HeterPTlm}
\alias{HeterPTlm}
\alias{plot.HeterPTlm}
\alias{summary.HeterPTlm}
\title{Quantile Regression with Polya Tree
for univariate outcome}
\usage{
  HeterPTlm(y, X, mcmc, prior, quan = 0.5, method = "ss",
    den = FALSE)

  \method{coef}{HeterPTlm} (mod, ...)

  \method{plot}{HeterPTlm} (obj, ...)

  \method{summary}{HeterPTlm} (mod, ...)
}
\arguments{
  \item{y}{[n] response}

  \item{x}{[n, p] covariates matrix}

  \item{mcmc}{mcmc parameters}

  \item{prior}{priors}

  \item{quan}{quantile requested. Multiple quantiles
  allowed.}

  \item{method}{\itemize{ \item [ss] spike-slab }}
}
\value{
  an object of class \code{HeterPTlm}
}
\description{
  Bayesian Quantile Regression with Polya Tree for
  univariate outcome.
}
\details{
  Fit a quantile regression for univariate response.
}
\author{
  Minzhao Liu, Mike Daniels
}

