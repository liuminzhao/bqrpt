      subroutine myptdensity1(maxm, n, y, 
     &     mu, sigma,
     &     nburn, nsave, nskip, arate,  
     &     grid, ngrid, f, alpha,
     &     musave, sigmasave, alphasave, fsave,
     &     ratesave, tunesave)

c$$$Time-stamp: <liuminzhao 11/05/2012 14:30:59>

      implicit none

C     PT PARAMETER
      integer maxm, mdzero

C     DATA
      integer n
      real*8 y(n), v(n), vc(n)

C     BASELINE
      real*8 mu, sigma, alpha, muc, sigmac, alphac, sigma2, sigma2c
      real*8 sd, sd2, theta, thetac, alpha2, alpha2c

C     MCMC PARAMETER
      integer nburn, nsave, nskip
      real*8 arate

C     GRID PARAMETER
      integer ngrid
      real*8 grid(ngrid)

C     OUTPUT 
      real*8 f(ngrid)

C     SAVE
      real*8 musave(nsave), sigmasave(nsave), alphasave(nsave)
      real*8 fsave(nsave, ngrid)
      real*8 ratesave(nburn, 3), tunesave(nburn, 3)

C     TMP
      integer i, j, k, l, m
      real*8 loglikeo, loglikec, logprioro, logpriorc, logcgkc, logcgko
      integer whicho(n), whichn(n)
      integer iscan, nscan, skipcount, isave
      real*8 tau(2)

C     FUNCTION
      real*8 rnorm, ratio, runif, myrnorm, myrunif, dnrm, dgamma2

C     TUNE 
      integer att(3), acc(3)
      real*8 tune(3)

C     INITIAL
      mdzero = 0

      mu = 0.d0
      sigma = 1.d0
      sigma2 = 1.d0

      alpha = 1.d0

      tau(1) = 0.1
      tau(2) = 0.1

      do i=1,3
         att(i) = 0
         acc(i) = 0
         tune(i) = 0.1
      end do
      
C     MCMC

      isave = 0
      skipcount = 0
      nscan = nburn + (nskip + 1) * nsave

C     FIRST 

      call loglik_unippt(n,mdzero, maxm, alpha, mu, sigma2, y, 
     &     whicho, whichn, loglikeo)

C     ROLLING

      do iscan = 1, nscan

C     MU
         att(1) = att(1) + 1
c         muc = myrnorm(mu, 0.1d0)
         muc = myrnorm(mu, tune(1))
         do i = 1, n
            v(i) = y(i) - mu
            vc(i) = y(i) - muc
         end do

         logprioro = dnrm(mu, 0.d0, 30.d0, 1)
         logpriorc = dnrm(muc, 0.d0, 30.d0, 1)

         loglikec = 0.d0

         call loglik_unippt(n,mdzero, maxm, alpha, 0.d0, sigma2, vc,
     &        whicho, whichn, loglikec)

c$$$         print*, mu
c$$$         print*, loglikeo
c$$$         print*, muc
c$$$         print*, loglikec 

         ratio=loglikec + logpriorc - loglikeo- logprioro

         if ( ratio .gt. 0 .or. log(dble(myrunif(0.d0, 1.d0))).lt. 
     &        ratio) then 
            loglikeo = loglikec
            mu = muc 
            do i = 1, n
               v(i) = vc(i)
            end do
            acc(1) = acc(1) + 1
         end if

C     SIGMA
         att(2) = att(2) + 1
         theta = log(sigma) 
         thetac = myrnorm(theta, tune(2))
         logcgkc = -theta
         logcgko = -thetac
         sigmac = exp(thetac)

         loglikec = 0.d0
         call loglik_unippt(n, mdzero, maxm, alpha, 0.d0, sigmac**2 , v, 
     &        whicho, whichn, loglikec)

         logpriorc=-tau(1)*thetac - tau(2)*exp(-2*thetac)/2
         logprioro=-tau(1)*theta-tau(2)*exp(-2*theta)/2
         ratio=loglikec + logpriorc -loglikeo -logprioro+
     &        logcgkc - logcgko

         if (ratio .gt. 0 .or. log(dble(myrunif(0.d0, 1.d0))).lt. ratio) 
     &        then 
            loglikeo=loglikec
            sigma = sigmac
            sigma2 = sigma**2
            acc(2) = acc(2) + 1
         end if

C     ALPHA

         att(3) = att(3) + 1
         theta = log(alpha) 
         thetac = myrnorm(theta, tune(3))
         logcgkc = -theta
         logcgko = -thetac
         alphac = exp(thetac)

         loglikec = 0.d0
         call loglik_unippt(n, mdzero, maxm, alphac, 0.d0, sigma2 , v, 
     &        whicho, whichn, loglikec)

         logpriorc=dgamma2(alphac, 1.d0, 1.d0,1)
         logprioro=dgamma2(alpha, 1.d0, 1.d0,1)

c         logpriorc=-tau(1)*thetac - tau(2)*exp(-2*thetac)/2
c         logprioro=-tau(1)*theta-tau(2)*exp(-2*theta)/2
         ratio=loglikec + logpriorc -loglikeo -logprioro+
     &        logcgkc - logcgko

         if (ratio .gt. 0 .or. log(dble(myrunif(0.d0, 1.d0))).lt. ratio) 
     &        then 
            loglikeo=loglikec
            alpha = alphac
c            print*, alpha
            acc(3) = acc(3) + 1
         end if

C     TUNING
         if ((att(1).ge.50).and.(iscan.le. nburn)) then 
            do i=1, 3
               if (dble(acc(i))/dble(att(i)) .gt. arate) then 
                  tune(i)=tune(i) + 
     &                 min(0.1d0,dble(10)/sqrt(dble(iscan)))
                 else
                  tune(i)=tune(i)-
     &                  min(0.1d0,dble(10)/sqrt(dble(iscan)))
               end if
               
               if (tune(i).gt. dble(10)) then 
                  tune(i)=10
               end if

               if (tune(i) .lt. 0.01d0) then 
                  tune(i)=0.01
               end if

               ratesave(iscan, i) = dble(acc(i))/dble(att(i))
               tunesave(iscan, i) = tune(i)
               acc(i)=0
               att(i)=0
               
            end do
         end if

C     SAVE

         if (iscan.gt. nburn) then 
            skipcount = skipcount+1
            if (skipcount.gt.nskip) then 

               isave=isave+1

               musave(isave) = mu
               sigmasave(isave) = sigma
               alphasave(isave) = alpha

C     density estimate
               do i=1,ngrid
                  loglikec=0.d0
                  call gridupptprior(grid(i),maxm,mdzero ,n,
     &                 alpha, mu, sigma2, y,
     &                 whicho,whichn,loglikec)
                  f(i)=f(i)+exp(loglikec)  
                  fsave(isave, i) = exp(loglikec)
               end do

               skipcount=0

            end if
         end if
      end do

      do i=1, ngrid
         f(i)=f(i)/dble(nsave)
      end do
      
      return
      end


c=======================================================================
      subroutine loglik_unippt(nsubject,mdzero,maxm,alpha,mu,sigma,b,
     &                        whicho,whichn,logliko)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for 
c     the baseline parameters in a random effect model using
c     in a marginal univariate partially specified PT.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none

c++++ Input
      integer mdzero,maxm,nsubject
      real*8 alpha
      real*8 mu,sigma 
      real*8 b(nsubject)

c++++ Working External
      integer whicho(nsubject),whichn(nsubject)

c++++ Working Internal
      integer countero,countern
      integer i,j,je2,k,k1,k2,l,nint,parti
      integer ok
      real*8 dnrm,invcdfnorm
      real*8 prob
      real*8 quan
      real*8 tmp1,tmp2

c++++ Output
      real*8 logliko

      logliko=0.d0
      do i=1,nsubject
         call rchkusr()   
c+++++++ first observation
         if(i.eq.1)then
              logliko=dnrm(b(1),mu,sqrt(sigma),1)
c+++++++ following observations
           else
              quan=mu
              countero=0
              if(b(i).le.quan) then
                  parti=1
                  do l=1,i-1
                     if(b(l).le.quan)then
                        countero=countero+1
                        whicho(countero)=l
                     end if   
                  end do
               else
                  parti=2
                  do l=1,i-1
                     if(b(l).gt.quan)then
                        countero=countero+1
                        whicho(countero)=l
                     end if   
                  end do
              end if   
              if(mdzero.ne.0)then 
                 logliko=logliko+
     &                log(2.d0*alpha+dble(2*countero))-
     &                log(2.d0*alpha+dble(i-1))
              end if 

              if(countero.eq.0)go to 1
              ok=1
              j=2
              do while(ok.eq.1.and.j.le.maxm)
                 nint=2**j
                 je2=j**2
                 prob=1.d0/dble(nint)
              
                 k1=2*(parti-1)+1
                 k2=2*(parti-1)+2
                 quan=invcdfnorm(dble(k1)*prob,mu,
     &                           sqrt(sigma),1,0)
                 if(b(i).le.quan)then
                   parti=k1
                   k=k1
                  else
                   parti=k2
                    k=k2
                 end if  
                 countern=0
                 if(k.eq.1)then
                    do l=1,countero
                       if(b(whicho(l)).le.quan)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                       end if   
                    end do
                  else if(k.eq.nint)then
                    quan=invcdfnorm(dble(k-1)*prob,mu,
     &                              sqrt(sigma),1,0) 
                    do l=1,countero
                       if(b(whicho(l)).gt.quan)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                       end if   
                    end do
                  else
                    tmp1=invcdfnorm(dble(k-1)*prob,mu,
     &                              sqrt(sigma),1,0)
                    tmp2=invcdfnorm(dble(k  )*prob,mu,
     &                              sqrt(sigma),1,0)

                    if(tmp1.ge.tmp2)then
                       call rexit("Error in the limits")
                    end if  
                 
                    do l=1,countero
                       if(b(whicho(l)).gt.tmp1.and.
     &                    b(whicho(l)).le.tmp2)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                       end if   
                    end do
                 end if
              
                 logliko=logliko+
     &                 log(2.d0*alpha*dble(je2)+dble(2*countern))-
     &                 log(2.d0*alpha*dble(je2)+dble(  countero))

                 if(countern.eq.0)then
                    ok=0
                  else  
                    countero=countern
                    do l=1,countern
                       whicho(l)=whichn(l)
                    end do
                    j=j+1
                 end if   
              end do
1             continue
              logliko=logliko+dnrm(b(i),mu,sqrt(sigma),1)
         end if
      end do
      return
      end

C     =========================================
      subroutine gridupptprior(theta,maxm,mdzero,nsubject,alpha,mu,
     &                         sigma,b,
     &                         whicho,whichn,logprioro)
c=======================================================================
c     This subroutine evaluate the log-contional prior distribution,
c     arising in a marginal univariate partially specified PT, 
c     for a value in a grid 'theta'.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none

c++++ Input
      integer maxm,mdzero,nsubject
      real*8 alpha,mu,sigma,b(nsubject),theta

c++++ Working Externals
      integer whicho(nsubject),whichn(nsubject)

c++++ Working Internals
      integer countero,countern
      integer j,je2,k,k1,k2,l,nint,parti
      integer ok
      real*8 dnrm,invcdfnorm
      real*8 prob
      real*8 quan
      real*8 tmp1,tmp2

c++++ Output
      real*8 logprioro
      
      logprioro=0.d0
      
      quan=mu
      countero=0
      if(theta.le.quan) then
          parti=1
          do l=1,nsubject
             if(b(l).le.quan)then
                countero=countero+1
                whicho(countero)=l
             end if   
          end do
        else
          parti=2
          do l=1,nsubject
             if(b(l).gt.quan)then
                countero=countero+1
                whicho(countero)=l
             end if   
          end do
      end if  

      if(mdzero.ne.0)then 
         logprioro=logprioro+
     &    log(2.d0*alpha+dble(2*countero))-
     &    log(2.d0*alpha+dble(nsubject))
      end if 

      if(countero.eq.0)go to 1

      ok=1
      j=2
      do while(ok.eq.1.and.j.le.maxm)
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
        
         k1=2*(parti-1)+1
         k2=2*(parti-1)+2
         quan=invcdfnorm(dble(k1)*prob,mu,sqrt(sigma),1,0)
      
         if(theta.le.quan)then
           parti=k1
           k=k1
          else
           parti=k2
           k=k2
         end if  
         
         countern=0

         if(k.eq.1)then
            do l=1,countero
               if(b(whicho(l)).le.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else if(k.eq.nint)then
            quan=invcdfnorm(dble(k-1)/dble(nint),mu,
     &                      sqrt(sigma),1,0) 
            do l=1,countero
               if(b(whicho(l)).gt.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else
            tmp1=invcdfnorm(dble(k-1)/dble(nint),mu,
     &                      sqrt(sigma),1,0)
            tmp2=invcdfnorm(dble(k  )/dble(nint),mu,
     &                      sqrt(sigma),1,0)

            if(tmp1.ge.tmp2)then
              call rexit("Error in the limits")
            end if  
         
            do l=1,countero
               if(b(whicho(l)).gt.tmp1.and.
     &            b(whicho(l)).le.tmp2)then
                 countern=countern+1
                 whichn(countern)=whicho(l)
               end if
            end do
         end if
        
         logprioro=logprioro+
     &       log(2.d0)+                
     &       log(     alpha*dble(je2)+dble(  countern))-
     &       log(2.d0*alpha*dble(je2)+dble(  countero))

         if(countern.eq.0)then
             ok=0
           else  
             countero=countern
             do l=1,countern
                whicho(l)=whichn(l)
             end do
             j=j+1
         end if   
      end do
1     continue

      logprioro=logprioro+dnrm(theta,mu,sqrt(sigma),1)
      
      return
      end

C=======================================================
C     POST QUANTILE FOR EACH INTERATION IN MCMC
C     Time-stamp: <liuminzhao 03/13/2012 00:33:09>
C     2012/03/12 
C=======================================================

      subroutine postquantile(n, v, sigma2, alpha, maxm, quan, quansave,
     &     nquan)
      implicit none

C     output and input
      integer n, maxm, nquan, bign(nquan)
      real*8 v(n), sigma2, alpha, quan(nquan), quansave(nquan)

C     working space
      real*8 p(2**maxm), interval(2**maxm-1), mu, tmp
      integer ncount(2**maxm), i,j, count, nint
      integer nmatrix(maxm, 2**maxm)
      
C     FUNCTION 
      real*8 invcdfnorm

C     INITIAL
      mu=0
      nint=2**maxm
      do i=1, nint
         p(i)=0
         ncount(i)=0
      end do 

      do i=1, maxm
         do j=1, nint
            nmatrix(i,j)=0
         end do
      end do

C     begin
      
C     first get interval

      do i=1, nint-1
         interval(i) = invcdfnorm(dble(i)/dble(nint), mu, sqrt(sigma2),
     &        1,0)
      end do

C     get ncount 
      
      do i=1, nint
         count=0
         if ( i.eq.1) then
            do j=1, n
               if (v(j) .lt. interval(i)) then
                  count = count + 1
               end if 
            end do
         else if (i.eq.nint) then 
            do j=1,n
               if (v(j) .gt. interval(i-1)) then
                  count = count + 1
               end if
            end do
         else 
            do j=1,n
               if (v(j).gt.interval(i-1).and.v(j).lt.interval(i)) then
                  count=count+1
               end if
            end do
         end if
         ncount(i) = count
      end do

C     get nmatrix 
      do i=1, nint
         nmatrix(maxm, i)=ncount(i)
      end do

      do i=(maxm-1), 1,-1
         do j=1, 2**i
            nmatrix(i,j) = nmatrix(i+1, 2*j-1)+nmatrix(i+1, 2*j)
         end do
      end do

C     get p

      do j=1, nint
         tmp=1
         do i=maxm, 2, -1
            tmp=tmp*(alpha*dble(i**2)+nmatrix(i, floor(dble(j-1)/
     &           dble(2**(maxm-i)))+1))/(2.d0*alpha*dble(i**2) + 
     &           nmatrix(i-1,floor(dble(j-1)/dble(2**(maxm-i+1)))+1))
         end do
         tmp=tmp*dble(0.5)
         p(j)=tmp
      end do

C     get N and quansave
      do i=1, nquan
         tmp=0
         do j=1, nint
            tmp=tmp+p(j)
            if (tmp.gt. quan(i)) exit 
         end do 
         bign(i) = j
         quansave(i) = invcdfnorm(dble(quan(i)-tmp+dble(j)*p(j))/
     &        dble(2**maxm)/p(j), mu, sqrt(sigma2), 1,0)
      end do

      return
      end 
