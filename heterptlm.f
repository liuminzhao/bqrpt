c==========================================================
      subroutine heterptlm(maxm, mdzero, nrec, p, x, y, betapm, betapv,
     &     tau, betasave, gammasave, gammapm,gammapv, alpha, beta,gamma,  
     &     nsave, sigma2, v, a0b0, mcmc, whicho, whichn, 
     &     sigmasave, alphasave, f, ngrid,grid, quan, quansave, nquan,
     &     ratesave, tunesave, hetersave,
     &     propv)
c===========================================================
c$$$
C          Time-stamp: <liuminzhao 03/26/2012 13:37:32>
c$$$      
c$$$      2012/03/26 set back to 0 for acc* when tuning 
C         add ratesave(100, 2p-1), tunesave(nburn, 2p-1)
c$$$      2012/03/13 add postquantile
c$$$      2012/03/09 mcmc fortran for myheterptlm13.R 
c$$$      
c===========================================================
      implicit none

C     PT parameter
      integer maxm , mdzero

C     obs
      integer nrec, p, ngrid, nquan
      real*8 x(nrec,p) , y(nrec), quan(nquan)

C     prior
      real*8 betapm(p), betapv(p,p)
      real*8 tau(2)

C     mcmc
      integer mcmc(4), nburn, nskip, nsave,ndisp

C     store
      real*8 betasave(nsave, p), gammasave(nsave, p)
      real*8 sigmasave(nsave), alphasave(nsave),quansave(nsave,nquan)
      real*8 gammapv(p,p), gammapm(p), f(ngrid), grid(ngrid)

C     current
      real*8 alpha, beta(p), gamma(p),  sigma2, v(nrec)

C     external working 
      integer whicho(nrec), whichn(nrec)
      real*8 betac(p) , gammac(p)
      real*8 vc(nrec)
      real*8 a0b0(2)

C     internal working
      integer i, iscan, isave, j, nscan, skipcount,k , discount
      real*8 alphac
      real*8 logcgkn, logcgko, loglikec, loglikeo, logpriorc, logprioro
      real*8 loglikaddc, loglikaddo
      real*8 mu,   rnorm, sd,sdc, theta, thetac, ratio

      real*8 tmp1, tmp2, tmpquan(nquan)
      real runif 
      integer ratecount

C     tunning
      real*8 tune1(p), tune2(p), tune3, tune4 
      integer att1(p), att2(p), att3, att4
      integer acc1(p), acc2(p), acc3, acc4
      real*8 propv(p,p)

C     cpu
      real*8 sec00, sec0, sec1, sec

C     FUNCTION
      real*8 myrnorm,myrunif, dnrm, dgamma2

C     DEBUG
      real*8 ratesave(200, 2*p+2)
      real*8 tunesave(mcmc(1), 2*p+2)
      real*8 hetersave(mcmc(1), p)

C     initial 
      nburn=mcmc(1)
      nskip=mcmc(2)
      nsave=mcmc(3)
      ndisp=mcmc(4)

      sd=sqrt(sigma2)

C      mdzero=0

      do i=1, p
         tune1(i)=0.1
         tune2(i)=0.1
         acc1(i)=0
         att1(i)=0
         acc2(i)=0
         att2(i)=0
      end do 
      tune3=0.1
      tune4=0.1
      att3=0
      att4=0
      acc3=0
      acc4=0

      do i=1, nquan
         tmpquan(i)=0
      end do

      mu=0.d0
      ratecount=1

C----------------------------------------------
C     start mcmc 
C---------------------------------------------      

      isave=0
      skipcount=0
      discount=0
      nscan = nburn + (nskip + 1)*(nsave)
      
      call cpu_time(sec0)
      sec00 = 0.d0 

C     FIRST 
      loglikeo = 0.d0

      do i = 1, nrec
         tmp1=0.d0
         tmp2=0.d0
         do j = 1,p
            tmp1= tmp1+x(i,j)*beta(j)
            tmp2= tmp2+x(i,j)*gamma(j)
         end do
         v(i) = (y(i)-tmp1)/tmp2
      end do

      call loglik_unippt(nrec,mdzero, maxm, alpha, mu, sigma2, v, 
     &     whicho, whichn, loglikeo)

C=======================================
C     ROLLING 
C=======================================
      
      do iscan=1, nscan
         
C     ===================================
C     BETA 
         do i=1, p
            att1(i)=att1(i)+1
            do j=1,p
               betac(j) = beta(j)
            end do
            betac(i) = myrnorm(beta(i), tune1(i)*sqrt(propv(i,i)))

            do k = 1, nrec
               tmp1=0.d0
               tmp2=0.d0
               do j = 1,p
                  tmp1= tmp1+x(k,j)*betac(j)
                  tmp2= tmp2+x(k,j)*gamma(j)
               end do
               vc(k) = (y(k)-tmp1)/tmp2
            end do

C     log prior
            logpriorc=dnrm(betac(i), betapm(i), betapv(i,i), 1)
            logprioro=dnrm(beta(i), betapm(i), betapv(i,i), 1)

C     log likelihood

            loglikec=0.d0
            call loglik_unippt(nrec,mdzero, maxm, alpha, mu, sigma2, vc,
     &           whicho, whichn, loglikec)

C     acceptance step
            ratio=loglikec + logpriorc - loglikeo- logprioro

            if (log(dble(myrunif(0.d0, 1.d0))).lt. ratio) then 
               loglikeo=loglikec
               beta(i)=betac(i)
               do k=1, nrec
                  v(k) = vc(k)
               end do
               acc1(i)=acc1(i)+1
            end if

         end do
C     ================================================
c$$$
c$$$C     GAMMA 
         do i=2, p
            att2(i)=att2(i)+1
            do j=1,p
               gammac(j) = gamma(j)
            end do
            gammac(i) = myrnorm(gamma(i), tune2(i)*sqrt(propv(i,i)))
            gammac(1)=1
            do k = 1, nrec
               tmp1=0.d0
               tmp2=0.d0
               do j = 1,p
                  tmp1= tmp1+x(k,j)*beta(j)
                  tmp2= tmp2+x(k,j)*gammac(j)
               end do
               vc(k) = (y(k)-tmp1)/tmp2
            end do
C     =============================================
C     check if gammac makes X'gamma less than 0
            if (iscan .lt. nburn) then 
               tmp1=dble(100)
               do j=1, nrec
                  tmp2=0.d0
                  do k=1,p
                     tmp2=tmp2+x(j,k)*gammac(k)
                  end do
                  if (tmp2 .lt. tmp1) then 
                     tmp1=tmp2
                  end if
               end do
               hetersave(iscan, i)=tmp1
            end if 
C     =============================================

C     log prior
            logpriorc=dnrm(gammac(i), gammapm(i), gammapv(i,i), 1)
            logprioro=dnrm(gamma(i), gammapm(i), gammapv(i,i), 1)

C     log likelihood

            loglikec=0.d0
            call loglik_unippt(nrec,mdzero, maxm, alpha, mu, sigma2, vc,
     &           whicho, whichn, loglikec)

C     additional loglike 
            loglikaddc=0.d0
            loglikaddo=0.d0

            do k=1, nrec
               tmp1=0.d0
               tmp2=0.d0
               do j=1,p
                  tmp1=tmp1+x(k,j)*gammac(j)
                  tmp2=tmp2+x(k,j)*gamma(j)
               end do
               loglikaddc=loglikaddc+log(tmp1)
               loglikaddo=loglikaddo+log(tmp2)
            end do
            
C     acceptance step
            ratio=loglikec + logpriorc - loglikeo- logprioro - 
     &           loglikaddc + loglikaddo

            if (log(dble(myrunif(0.d0, 1.d0))).lt. ratio) then 
               loglikeo=loglikec
               gamma(i)=gammac(i)
               gamma(1)=1
               do k=1, nrec
                  v(k) = vc(k)
               end do
               acc2(i)=acc2(i)+1
            end if
         end do
C     =============================================================
C     SIGMA 

         att3=att3+1
         theta=log(sd)
         thetac=myrnorm(theta, tune3)
         logcgkn=-theta
         logcgko=-thetac
         sdc=exp(thetac)
         
C     likelihood
         loglikec=0.d0
         
         call loglik_unippt(nrec,mdzero, maxm, alpha, mu, sdc**2, v, 
     &        whicho, whichn, loglikec)
         
         
C     log prior
         logpriorc=-tau(1)*thetac - tau(2)*exp(-2*thetac)/2
         logprioro=-tau(1)*theta-tau(2)*exp(-2*theta)/2
C         logprioro=0.d0
C         logpriorc=0.d0

C     acceptance
         ratio=loglikec + logpriorc -loglikeo -logprioro+
     &        logcgkn - logcgko
         
         if (log(dble(myrunif(0.d0, 1.d0))).lt. ratio) then 
            loglikeo=loglikec
            sd=sdc
            sigma2=sd**2
            acc3=acc3+1
         end if
         
C     ===========================================================
C     precision : alpha 
         att4=att4+1
         theta=log(alpha)
         thetac=myrnorm(theta,tune4)
         logcgkn=-theta
         logcgko=-thetac
         alphac=exp(thetac)
         
C     likelihood
         loglikec=0.d0
         
         call loglik_unippt(nrec,mdzero, maxm, alphac, mu, sigma2, v, 
     &        whicho, whichn, loglikec)
         
C     log prior
         logpriorc=dgamma2(alphac, a0b0(1), a0b0(2),1)
         logprioro=dgamma2(alpha, a0b0(1), a0b0(2),1)
C     acceptance
         ratio=loglikec + logpriorc - loglikeo- logprioro +
     &        logcgkn - logcgko
         
         if (log(dble(myrunif(0.d0, 1.d0))).lt. ratio) then 
            loglikeo=loglikec
            alpha=alphac
            acc4=acc4+1
         end if
C     ===================================================
C     TUNING 

         if ((att1(1).ge.100).and.(iscan.le. nburn)) then 
            do i=1, p
C               if (dble(acc1(i))/dble(att1(i)) .gt. 0.25d0) then 
               if (dble(acc1(i))/dble(att1(i)) .gt. 0.4d0) then 
                  tune1(i)=tune1(i) + 
     &                 min(0.1d0,dble(10)/sqrt(dble(iscan)))
                 else
                  tune1(i)=tune1(i)-
     &                  min(0.1d0,dble(10)/sqrt(dble(iscan)))
               end if
               
               if (tune1(i).gt. dble(10)) then 
                  tune1(i)=100
               end if

               if (tune1(i) .lt. 0.01d0) then 
                  tune1(i)=0.01
               end if

C               if (dble(acc2(i))/dble(att2(i)) .gt. 0.25d0) then 
               if (dble(acc2(i))/dble(att2(i)) .gt. 0.4d0) then 
                  tune2(i)=tune2(i) + 
     &                 min(0.1d0,dble(10)/sqrt(dble(iscan)))
                 else
                  tune2(i)=tune2(i)-
     &                  min(0.1d0,dble(10)/sqrt(dble(iscan)))
               end if  


               if (tune2(i).gt. dble(10)) then 
                  tune2(i)=100
               end if

               if (tune2(i) .lt. 0.01d0) then 
                  tune2(i)=0.01
               end if

C     SET UP TO 0


               ratesave(ratecount, i)=dble(acc1(i))/dble(att1(i))
               ratesave(ratecount, 3+i)=dble(acc2(i))/dble(att2(i))
               acc1(i)=0
               att1(i)=0
               acc2(i)=0
               att2(i)=0

            end do 

C            if (dble(acc3)/dble(att3) .gt. 0.25d0) then 
            if (dble(acc3)/dble(att3) .gt. 0.4d0) then 
               tune3=tune3 + 
     &              min(0.01d0,dble(10)/sqrt(dble(iscan)))
            else
               tune3=tune3-
     &              min(0.01d0,dble(10)/sqrt(dble(iscan)))
            end if  
            
            
            if (tune3.gt. dble(10)) then 
               tune3=10
            end if
            
            if (tune3 .lt. 0.01d0) then 
               tune3=0.01
            end if

C            if (dble(acc4)/dble(att4) .gt. 0.25d0) then 
            if (dble(acc4)/dble(att4) .gt. 0.4d0) then 
               tune4=tune4 + 
C     &              min(0.01d0,dble(10)/sqrt(dble(iscan)))
     &              min(0.1d0,dble(10)/sqrt(dble(iscan)))
            else
               tune4=tune4-
C     &              min(0.01d0,dble(10)/sqrt(dble(iscan)))
     &              min(0.1d0,dble(10)/sqrt(dble(iscan)))
            end if  
            
            
            if (tune4.gt. dble(10)) then 
               tune4=10
            end if
            
            if (tune4 .lt. 0.01d0) then 
               tune4=0.01
            end if

            ratesave(ratecount, 7)=dble(acc3)/dble(att3)
            ratesave(ratecount, 8)=dble(acc4)/dble(att4)
            acc3=0
            att3=0
            acc4=0
            att4=0

            ratecount=ratecount+1
         end if


C ========================================================
C     save 

C     TUNESAVE

C         print*, tune1
         if (iscan .lt. nburn) then 
            tunesave(iscan, 1)=tune1(1)
            tunesave(iscan, 2)=tune1(2)
            tunesave(iscan, 3)=tune1(3)
            tunesave(iscan, 4)=tune2(1)
            tunesave(iscan, 5)=tune2(2)
            tunesave(iscan, 6)=tune2(3)
            tunesave(iscan, 7)=tune3
            tunesave(iscan, 8)=tune4
            
         end if 

         if (iscan.gt. nburn) then 
            skipcount = skipcount+1
            if (skipcount.gt.nskip) then 

C               print* , beta , gamma 

               isave=isave+1
               discount=discount+1
C               print*, isave 

               do j=1,p
                  betasave(isave,j)=beta(j)
                  gammasave(isave, j)=gamma(j)
               end do

               sigmasave(isave)=sigma2
               alphasave(isave)=alpha

C     density estimate
               do i=1,ngrid
                  loglikec=0.d0
                  call gridupptprior(grid(i),maxm,mdzero ,nrec,
     &                 alpha,mu,sigma2,v,
     &                 whicho,whichn,loglikec)
                  f(i)=f(i)+exp(loglikec)  
               end do

C     postquantile
               call postquantile(nrec, v, sigma2, alpha, maxm,
     &              quan, tmpquan, nquan)
               
               do i=1, nquan
                  quansave(isave, i)= tmpquan(i)
               end do
               
               skipcount=0
               if (discount .ge. ndisp) then
                  call cpu_time(sec1)
                  sec00=sec00 + (sec1-sec0)
                  sec=sec00
                  sec0 = sec1
                  print*,  isave, " out of " , nsave,
     &                 " for time ", floor(sec)
                  discount=0
               end if
            end if
         end if
      end do
C     +===============================================================

C     final MH rate
      ratesave(100, 1)=dble(acc1(1))/dble(att1(1))
      ratesave(100, 2)=dble(acc1(2))/dble(att1(2))
      ratesave(100, 3)=dble(acc1(3))/dble(att1(3))
      ratesave(100, 4)=0.d0
      ratesave(100, 5)=dble(acc2(2))/dble(att2(2))
      ratesave(100, 6)=dble(acc2(3))/dble(att2(3))
      ratesave(100, 7)=dble(acc3)/dble(att3)
      ratesave(100, 8)=dble(acc4)/dble(att4)
      
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
