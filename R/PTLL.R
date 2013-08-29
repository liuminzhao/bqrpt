# PT log lik given y, maxm, mu, sigma, mdzero, alpha
PTLL <- function(y, maxm=floor(log(length(y))/log(2)),
                 mu = 0, sigma = 1, mdzero = 0, alpha = 1){

  n <- length(y)
  dyn.load('heterptlm.so')

  whicho <- whichn <- rep(0, n)

  tmp <- .Fortran("loglik_unippt",
                  as.integer(n),
                  as.integer(mdzero),
                  as.integer(maxm),
                  as.double(alpha),
                  as.double(mu),
                  as.double(sigma^2),
                  as.double(y),
                  as.integer(whicho),
                  as.integer(whichn),
                  ans = as.double(0)
                  )

  tmp$ans
}

PTLL2 <- function(y, maxm=floor(log(length(y))/log(2)),
                  mu = 0, sigma = 1, mdzero = 0, alpha = 1){

  n <- length(y)
  dyn.load('heterptlm.so')

  whicho <- whichn <- rep(0, n)

  tmp <- .Fortran("loglik_unippt",
                  as.integer(n),
                  as.integer(mdzero),
                  as.integer(maxm),
                  as.double(alpha),
                  as.double(0),
                  as.double(sigma^2),
                  as.double(y - mu),
                  as.integer(whicho),
                  as.integer(whichn),
                  ans = as.double(0)
                  )

  tmp$ans
}
