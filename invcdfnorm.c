#include <R.h>
#include <Rmath.h>
       
double F77_SUB(invcdfnorm)(double *p, double *mu, double *sigma, int *lower_tail, int *log_p)
{
  return qnorm(*p, *mu, *sigma, *lower_tail, *log_p);
}

double F77_SUB(dnrm)(double *x, double *mu, double *sigma, int *give_log)
{
  return dnorm(*x, *mu, *sigma, *give_log);
}


double F77_SUB(dgamma2)(double *x, double *shape, double *scale, int *give_log)
{
  return dgamma(*x, *shape, *scale, *give_log);
}


void F77_SUB(rndstart)(void) { GetRNGstate(); }

void F77_SUB(rndend)(void) { PutRNGstate(); }     
double F77_SUB(myrnorm)(double *mu, double *sigma) { return rnorm(*mu, *sigma) ;}
double F77_SUB(myrunif)(double *a, double *b) { return runif(*a, *b) ;}

int F77_SUB(ihmssf)(int *i, int *j, int *n){
  int ans;
  int i1;
  int j1;

  i1=*i-1;
  j1=*j-1;

  if (i <= j) 
    {
      ans=*n*i1-*i*i1/2+*j ;
    }
  else 
    {
      ans=*n*j1-*j*j1/2+*i;
    }

  return ans;
}

