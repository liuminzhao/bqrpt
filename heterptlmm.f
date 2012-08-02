c===============================================
      subroutine heterptlmm(maxm, mdzero, 
     &     nrec, nsub, p, q, x, y,
     &     betapm, betapv, tau, a0b0, gammapm, gammapv, 
     &     nburn, nskip, nsave, ndisp,
     &     betasave, gammasave, sigmasave, alphasave, 
c     &     whicho, whichn, betac, gammac, sigmac, vc, alphac, 
     &     beta, gamma, alpha, Sigma, propv, arate,
     &     ngrid, grid, f,
     &     nquan, qtile, quansave)
c$$$      Time-stamp: <2012/07/26>
c##################################################
      
      implicit none

C     PT
      integer maxm, mdzero

C     DATA
      integer nrec, nsub, p, q
      real*8 x(nrec, p), y(nsub, q), propv(p,p)

C     PRIOR
      real*8 betapm(p), betapv(p,p), gammapm(p), gammapv(p,p)
      real*8 tau(2), a0b0(2)

C     MCMC
      integer nburn, nskip, nsave, ndisp

C     SAVE
      real*8 betasave(p, q, nsave), gammasave(p, q, nsave)
      real*8 sigmasave(nsave, 3), alphasave(nsave)

C     WORKING
      integer whicho(nsub), whichn(nsub)
      real*8 betac(p,q), gammac(p,q), vc(nrec)
      real*8 alphac, sigmavecc(3)
      real*8 b(nsub, q), bz(nsub, q), bc(nsub,q)
      real*8 Sigma(2,2), Sigmac(2,2)

      integer i, iscan, isave, j, nscan, skipcount, k, dispcount, l
      real*8 loglikec, loglikeo, logpriorc, logprioro
      real*8 logcgkn, logcgko , loglikaddo, loglikaddc

      real*8 ratio, mu(q)

c     TUNING
      integer att1(q), att2(q), att3(3), att4, att5
      integer acc1(q), acc2(q), acc3(3), acc4, acc5
      real*8 tune1(q), tune2(q), tune3(3), tune4, tune5
      real*8 arate

C     CURRENT
      real*8 beta(p, q), gamma(p, q), alpha, sigmavec(3), v(nrec)

C     TEMP : workmr2c is is candidate for workmr2 (rnorm(q,q))
      integer parti(q)
      real*8 tmp1, tmp2, tmp(q), tmpnsub(nsub)
      real*8 workmr1(q,q), ortho(q,q), workmr2(q,q), workmr2c(q,q)
      real*8 orthoc(q,q)
      real*8 detlogl, detloglc
      real*8 linf(q), lsup(q)
      real*8 theta, thetac 

C     cpu
      real*8 sec00, sec0, sec1, sec

c     grid
      integer ngrid
      real*8 grid(ngrid)
      real*8 f(ngrid, q)

c     quantile
      integer nquan
      real*8 qtile(nquan)
      real*8 quansave(nsave, nquan*q)
      real*8 tmpquan(nquan)

C     FUNCTION
      real*8 myrnorm, dnrm, myrunif, dgamma2

C     INITIAL 
      do i=1, q
         do j=1, q
            workmr2(i,j)=myrnorm(0.d0, 1.d0)
         end do
      end do

      do i=1,q
         mu(i)=0.d0
      end do

      Sigma(1,1)=1.d0
      Sigma(1,2)=0.d0
      Sigma(2,1)=0.d0
      Sigma(2,2)=1.d0

      do i=1,q
         parti(i)=0
         linf(i)=0.d0
         lsup(i)=0.d0
         tmp(i)=0.d0
      end do

      do i=1, nsub
         whicho(i)=0
         whichn(i)=0
      end do

      do i=1, nsub
         do j=1, q
            bz(i,j)=0.d0
         end do
      end do

      do i=1, q
         tune1(i)=0.1
         tune2(i)=0.1
         acc1(i)=0
         att1(i)=0
         acc2(i)=0
         att2(i)=0
      end do 

      do i=1,3
         tune3(i)=0.1
         att3(i)=0
         acc3(i)=0
      end do

      tune4=0.1
      tune5=0.1

      att4=0
      att5=0

      acc4=0
      acc5=0

      sigmavec(1)=sqrt(Sigma(1,1))
      sigmavec(2)=sqrt(Sigma(2,2))
      sigmavec(3)=Sigma(1,2)/sigmavec(1)/sigmavec(2)

      call cpu_time(sec0)
      sec00 = 0.d0

C=========================================
c$$$  START MCMC
c=========================================
      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*nsave

C     FIRST : f(e1, ..., en)
      loglikeo=0.d0
      
C     residual
      do i=1, nsub
         do k=1, q
            tmp1=0.d0
            tmp2=0.d0
            do j = 1,p
               tmp1= tmp1+x(i,j)*beta(j, k)
               tmp2= tmp2+x(i,j)*gamma(j, k)
            end do
            b(i,k) = (y(i,k)-tmp1)/tmp2
         end do
      end do

C     residual matrix
c$$$      do i=1, nsub
c$$$         do j=1, q
c$$$            b(i, j) = v(q*(i-1)+j)
c$$$         end do
c$$$      end do

C     get O
C      call rhaar2(workmr1, workmr2, q, ortho)
      ortho(1,1)=1.d0
      ortho(1,2)=0.d0
      ortho(2,1)=0.d0
      ortho(2,2)=1.d0

      loglikeo=0.d0

C     log|Sigma|
      detlogl=log(Sigma(1,1)*Sigma(2,2)-Sigma(1,2)*Sigma(2,1))

C     first loglikPT 

c$$$      print*, 'before call loglik_mpt'
c$$$      print*, maxm
c$$$      print*, q
c$$$      print*, nsub
c$$$      print*, parti
c$$$      print*, whicho
c$$$      print*, whichn
c      print*, b(1,1)
c      print*, b(1,2)
c$$$      print*, bz
c$$$      print*, alpha
c$$$      print*, detlogl
c$$$      print*, linf
c$$$      print*, lsup
c$$$      print*, mu
c$$$      print*, Sigma
c$$$      print*, tmp
c$$$      print*, mdzero
c$$$      print*, loglikeo

      call loglik_mpt(maxm, q, nsub, parti,
     &     whicho, whichn, b, bz, alpha, detlogl,
     &     linf, lsup, mu, Sigma, 
     &     tmp, ortho, mdzero, loglikeo)

c      print*, loglikeo

CCCCCCCCCCCC  roll cccccccccccccccccccccc

      do iscan=1, nscan
         
c=========================================
c     beta
C=========================================

c         print*, 'begin beta'
         do k=1, q
            att1(k)=att1(k) + 1
            do i=1, p
               do j=1, q
                  betac(i,j)=beta(i,j)
               end do
            end do
            
            do i=1, p
               betac(i,k) = myrnorm(beta(i,k),tune1(k)*sqrt(propv(i,i)))
            end do

c            print*, betac

            do i = 1, nsub
               do j=1,q
                  tmp1=0.d0
                  tmp2=0.d0
                  do l = 1,p
                     tmp1= tmp1+x(i,l)*betac(l, j)
                     tmp2= tmp2+x(i,l)*gamma(l, j)
                  end do
                  bc(i, j) = (y(i,j)-tmp1)/tmp2
               end do
            end do

c     log prior
            logprioro=0.d0
            logpriorc=0.d0
            
            do i=1, p
               logpriorc= logpriorc + dnrm(betac(i,k), betapm(i), 
     &              betapv(i,i),1)
               logprioro= logprioro + dnrm(beta(i,k), betapm(i), 
     &              betapv(i,i),1)
            end do

c     log likelihood
            loglikec=0.d0

c$$$            do k=1, nsub
c$$$               do j=1, q
c$$$                  bc(k, j) = vc(q*(k-1)+j)
c$$$               end do
c$$$            end do

            call loglik_mpt(maxm, q, nsub, parti,
     &           whicho, whichn, bc, bz, alpha, detlogl,
     &           linf, lsup, mu, Sigma, 
     &           tmp, ortho, mdzero, loglikec)

c     acceptance step
            ratio=loglikec + logpriorc - logprioro - loglikeo

c$$$            print*, loglikec
c$$$            print*, loglikeo
c$$$
c$$$            print*, logpriorc
c$$$            print*, logprioro
c$$$
c$$$            print*, ratio

            if (log(dble(myrunif(0.d0, 1.d0))).lt. ratio) then 
               loglikeo=loglikec
               do i=1,p
                  beta(i, k)=betac(i,k)
               end do
c$$$               do k=1, nrec
c$$$                  v(k) = vc(k)
c$$$               end do
               do i=1, nsub
                  do j=1, q
                     b(i, j) = bc(i,j)
                  end do
               end do
               acc1(k)=acc1(k)+1
            end if
            
         end do

c========================================
c     GAMMA 
C========================================

c         print*, 'begin gamma'

         do k=1, q
            att2(k)=att2(k)+1
            do i=1,p
               do j=1,q
                  gammac(i,j) = gamma(i,j)
               end do
            end do

            do i=2,p
               gammac(i,k)=myrnorm(gamma(i,k),tune2(k)*sqrt(propv(i,i)))
            end do
c            gammac(1)=1.d0
            
C     ==================================
C     take gammac only when min(x'gammac) >0 

            tmp1=dble(100)
            do i=1, nsub
               do j=1, q
                  tmp2=0.d0
                  do l=1,p
                     tmp2=tmp2+x(i,l)*gammac(l,j)
                  end do
                  if (tmp2 .lt. tmp1) then 
                     tmp1=tmp2
                  end if
               end do
            end do

            if (tmp1 .gt. 0.d0) then 
               do i = 1, nsub
                  do j=1,q
                     tmp1=0.d0
                     tmp2=0.d0
                     do l = 1,p
                        tmp1= tmp1+x(i,l)*beta(l,j)
                        tmp2= tmp2+x(i,l)*gammac(l,j)
                     end do
                     bc(i,j) = (y(i,j)-tmp1)/tmp2
                  end do
               end do

c$$$               do k=1, nsub
c$$$                  do j=1, q
c$$$                     bc(k, j) = vc(q*(k-1)+j)
c$$$                  end do
c$$$               end do
               
C     log prior
               logprioro=0.d0
               logpriorc=0.d0
               
               do i=1, p
                  logpriorc=logpriorc+ dnrm(gammac(i,k), gammapm(i), 
     &                 gammapv(i,i), 1)
                  logprioro=logprioro+ dnrm(gamma(i,k), gammapm(i), 
     &                 gammapv(i,i), 1)
               end do

C     log likelihood
               
               loglikec=0.d0

               call loglik_mpt(maxm, q, nsub, parti,
     &              whicho, whichn, bc, bz, alpha, detlogl,
     &              linf, lsup, mu, Sigma, 
     &              tmp, ortho, mdzero, loglikec)
               
C     additional loglike 
               loglikaddc=0.d0
               loglikaddo=0.d0
               
               do i=1, nsub
                  do j=1, q
                     tmp1=0.d0
                     tmp2=0.d0
                     do l=1,p
                        tmp1=tmp1+x(i,l)*gammac(l,j)
                        tmp2=tmp2+x(i,l)*gamma(l,j)
                     end do
                     loglikaddc=loglikaddc+log(tmp1)
                     loglikaddo=loglikaddo+log(tmp2)
                  end do
               end do
               
C     acceptance step
               ratio=loglikec + logpriorc - loglikeo- logprioro - 
     &              loglikaddc + loglikaddo

               if (log(dble(myrunif(0.d0, 1.d0))).lt. ratio) then 
                  loglikeo=loglikec
                  do i=1,p
                     gamma(i,k)=gammac(i,k)
                  end do
c                  gamma(1)=1
c$$$                  do k=1, nrec
c$$$                     v(k) = vc(k)
c$$$                  end do
                  do i=1, nsub
                     do j=1, q
                        b(i, j) = bc(i,j)
                     end do
                  end do
                  acc2(k)=acc2(k)+1
               end if
            end if               
         end do

C========================================
C     Sigma 
c========================================

c         print*, 'begin sigma'

c     sigmavec(1)
         att3(1)=att3(1)+1
         theta=log(sigmavec(1))
         thetac=myrnorm(theta, tune3(1))
         logcgkn=-theta
         logcgko=-thetac
         Sigmac(1,1)=exp(2*thetac)
         Sigmac(2,2)=sigmavec(2)**2
         Sigmac(1,2)=exp(thetac)*sigmavec(2)*sigmavec(3)
         Sigmac(2,1)=exp(thetac)*sigmavec(2)*sigmavec(3)

         detloglc=log(Sigmac(1,1)*Sigmac(2,2)-Sigmac(1,2)**2)
         
c     likelihood
         loglikec=0.d0

         call loglik_mpt(maxm, q, nsub, parti,
     &        whicho, whichn, b, bz, alpha, detloglc,
     &        linf, lsup, mu, Sigmac, 
     &        tmp, ortho, mdzero, loglikec)

c     log prior
         logpriorc=-tau(1)*thetac-tau(2)*exp(-2*thetac)/2
         logprioro=-tau(1)*theta-tau(2)*exp(-2*theta)/2

c     acceptance
         ratio=loglikec + logpriorc -loglikeo -logprioro+
     &        logcgkn - logcgko
         
         if (log(dble(myrunif(0.d0, 1.d0))).lt. ratio) then 
            loglikeo=loglikec
            detlogl=detloglc
            sigmavec(1)=exp(thetac)
            Sigma(1,1)=Sigmac(1,1)
            Sigma(1,2)=Sigmac(1,2)
            Sigma(2,1)=Sigmac(2,1)
            acc3(1)=acc3(1)+1
         end if

c     sigmavec(2)
         att3(2)=att3(2)+1
         theta=log(sigmavec(2))
         thetac=myrnorm(theta, tune3(2))
         logcgkn=-theta
         logcgko=-thetac
         Sigmac(1,1)=sigmavec(1)**2
         Sigmac(2,2)=exp(2*thetac)
         Sigmac(1,2)=exp(thetac)*sigmavec(1)*sigmavec(3)
         Sigmac(2,1)=exp(thetac)*sigmavec(1)*sigmavec(3)

         detloglc=log(Sigmac(1,1)*Sigmac(2,2)-Sigmac(1,2)**2)
         
c     likelihood
         loglikec=0.d0

         call loglik_mpt(maxm, q, nsub, parti,
     &        whicho, whichn, b, bz, alpha, detloglc,
     &        linf, lsup, mu, Sigmac, 
     &        tmp, ortho, mdzero, loglikec)

c     log prior
         logpriorc=-tau(1)*thetac-tau(2)*exp(-2*thetac)/2
         logprioro=-tau(1)*theta-tau(2)*exp(-2*theta)/2

c     acceptance
         ratio=loglikec + logpriorc -loglikeo -logprioro+
     &        logcgkn - logcgko
         
         if (log(dble(myrunif(0.d0, 1.d0))).lt. ratio) then 
            loglikeo=loglikec
            detlogl=detloglc
            sigmavec(2)=exp(thetac)
            Sigma(2,2)=Sigmac(2,2)
            Sigma(1,2)=Sigmac(1,2)
            Sigma(2,1)=Sigmac(2,1)
            acc3(2)=acc3(2)+1
         end if

c     sigmavec(3): rho
         att3(3)=att3(3)+1
         theta=sigmavec(3)
         thetac=myrnorm(theta, tune3(3))
         thetac=max(thetac, -0.99)
         thetac=min(thetac, 0.99)
         Sigmac(1,1)=sigmavec(1)**2
         Sigmac(2,2)=sigmavec(2)**2
         Sigmac(1,2)=thetac*sigmavec(1)*sigmavec(2)
         Sigmac(2,1)=thetac*sigmavec(1)*sigmavec(2)

         detloglc=log(Sigmac(1,1)*Sigmac(2,2)-Sigmac(1,2)**2)
         
c     likelihood
         loglikec=0.d0

         call loglik_mpt(maxm, q, nsub, parti,
     &        whicho, whichn, b, bz, alpha, detloglc,
     &        linf, lsup, mu, Sigmac, 
     &        tmp, ortho, mdzero, loglikec)

c     log prior
c         logpriorc=-tau(1)*thetac-tau(2)*exp(-2*thetac)/2
c         logprioro=-tau(1)*theta-tau(2)*exp(-2*theta)/2

c     acceptance
c         ratio=loglikec + logpriorc -loglikeo -logprioro+
c     &        logcgkn - logcgko
         ratio=loglikec-loglikeo
         
         if (log(dble(myrunif(0.d0, 1.d0))).lt. ratio) then 
            loglikeo=loglikec
            detlogl=detloglc
            sigmavec(3)=thetac
            Sigma(1,2)=Sigmac(1,2)
            Sigma(2,1)=Sigmac(2,1)
            acc3(3)=acc3(3)+1
         end if

c========================================
c     Alpha 
c========================================

c         print*, 'begin alpha'

         att4=att4+1
         theta=log(alpha)
         thetac=myrnorm(theta,tune4)
         logcgkn=-theta
         logcgko=-thetac
         alphac=exp(thetac)
         
C     likelihood
         loglikec=0.d0

         call loglik_mpt(maxm, q, nsub, parti,
     &        whicho, whichn, b, bz, alphac, detlogl,
     &        linf, lsup, mu, Sigma, 
     &        tmp, ortho, mdzero, loglikec)
         
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

c=========================================
c     Ortho 
c========================================

c         print*, 'begin ortho'

         att5=att5+1
         do i=1,q
            do j=1,q   
               workmr2c(i,j)=myrnorm(workmr2(i,j),tune5)
            end do
         end do
         call rhaar2(workmr1,workmr2,q,orthoc)
         
         call loglik_mpt(maxm, q, nsub, parti,
     &        whicho, whichn, b, bz, alpha, detlogl,
     &        linf, lsup, mu, Sigma, 
     &        tmp, orthoc, mdzero, loglikec)

c         print*, 'prob here' 

c     acceptance
         ratio=loglikec-loglikeo

         if(log(dble(myrunif(0.d0, 1.d0))).lt.ratio)then
            acc5=acc5+1
            do i=1,q
               do j=1,q
                  ortho(i,j)=orthoc(i,j)
                  workmr2(i,j)=workmr2c(i,j)
               end do
            end do
            loglikeo=loglikec
         end if

c=========================================
c     Tuning
c ======================================
         
         if ((att1(1).ge.100).and.(iscan.le.nburn)) then 
            
c            print*, 'tune1'

            do i=1, q
               if (dble(acc1(i))/dble(att1(i)) .gt. arate) then 
                  tune1(i)=tune1(i) + 
     &                 min(0.1d0,dble(10)/sqrt(dble(iscan)))
                 else
                  tune1(i)=tune1(i)-
     &                  min(0.1d0,dble(10)/sqrt(dble(iscan)))
               end if
               
               if (tune1(i).gt. dble(100)) then 
                  tune1(i)=100.d0
               end if

               if (tune1(i) .lt. 0.01d0) then 
                  tune1(i)=0.01
               end if

               if (dble(acc2(i))/dble(att2(i)) .gt. arate) then 
                  tune2(i)=tune2(i) + 
     &                 min(0.1d0,dble(10)/sqrt(dble(iscan)))
                 else
                  tune2(i)=tune2(i)-
     &                  min(0.1d0,dble(10)/sqrt(dble(iscan)))
               end if  

               if (tune2(i).gt. dble(100)) then 
                  tune2(i)=100.d0
               end if

               if (tune2(i) .lt. 0.01d0) then 
                  tune2(i)=0.01
               end if

C     SET UP TO 0

c               ratesave(ratecount, i)=dble(acc1(i))/dble(att1(i))
c               ratesave(ratecount, 3+i)=dble(acc2(i))/dble(att2(i))
               acc1(i)=0
               att1(i)=0
               acc2(i)=0
               att2(i)=0

            end do 

c            print*, 'tune3'

            do i=1, 3
               if (dble(acc3(i))/dble(att3(i)) .gt. arate) then 
                  tune3(i)=tune3(i) + 
     &                 min(0.1d0,dble(10)/sqrt(dble(iscan)))
                 else
                  tune3(i)=tune3(i)-
     &                  min(0.1d0,dble(10)/sqrt(dble(iscan)))
               end if
               
               if (tune3(i).gt. dble(100)) then 
                  tune3(i)=100.d0
               end if

               if (tune3(i) .lt. 0.01d0) then 
                  tune3(i)=0.01
               end if

C     SET UP TO 0

               acc3(i)=0
               att3(i)=0
               
            end do 

c            print*, 'tune4'
         
            if (dble(acc4)/dble(att4) .gt. arate) then 
               tune4=tune4 + 
     &              min(0.1d0,dble(10)/sqrt(dble(iscan)))
            else
               tune4=tune4-
     &              min(0.1d0,dble(10)/sqrt(dble(iscan)))
            end if  
            
            
            if (tune4.gt. dble(10)) then 
               tune4=10.d0
            end if
            
            if (tune4 .lt. 0.01d0) then 
               tune4=0.01
            end if

c            ratesave(ratecount, 7)=dble(acc3)/dble(att3)
c            ratesave(ratecount, 8)=dble(acc4)/dble(att4)

            acc4=0
            att4=0

c            ratecount=ratecount+1

c            print*, 'tune5'

            if (dble(acc5)/dble(att5) .gt. arate) then 
               tune5=tune5 + 
     &              min(0.1d0,dble(10)/sqrt(dble(iscan)))
            else
               tune5=tune5-
     &              min(0.1d0,dble(10)/sqrt(dble(iscan)))
            end if  
            
            
            if (tune5.gt. dble(10)) then 
               tune5=10.d0
            end if
            
            if (tune5 .lt. 0.01d0) then 
               tune5=0.01
            end if

c            ratesave(ratecount, 7)=dble(acc3)/dble(att3)
c            ratesave(ratecount, 8)=dble(acc5)/dble(att5)

            acc5=0
            att5=0

         end if


c=================================
c     Save
c====================================

         if (iscan.gt. nburn) then 
c            print*, 'begin save'
            skipcount = skipcount+1
            if (skipcount.gt.nskip) then 

               isave=isave+1
               dispcount=dispcount+1


               do j=1,p
                  do k=1,q
                     betasave(j,k, isave)=beta(j,k)
                     gammasave( j,k, isave)=gamma(j,k)
                  end do
               end do

               do j=1,3
                  sigmasave(isave,j)=sigmavec(j)
               end do

               alphasave(isave)=alpha

c     marginal density
               do j=1, q
                  do k=1, nsub
                     tmpnsub(k)=b(k, j)
                  end do

                  do i= 1, ngrid
                     loglikec=0.d0
                     call gridupptprior(grid(i), maxm, mdzero, nsub,
     &                    alpha, mu(j), Sigma(j,j), tmpnsub,
     &                    whicho, whichn, loglikec)
                     f(i, j)=f(i,j)+exp(loglikec)
                  end do
               end do

c     postquantile 
               do j=1, q
                  do k=1, nsub
                     tmpnsub(k)= b(k, j)
                  end do
                  
                  call postquantile(nsub, tmpnsub, Sigma(j,j), alpha,
     &                 maxm, qtile, tmpquan, nquan)
                  
                  do i=1, nquan
                     quansave(isave, (j-1)*nquan+i)=tmpquan(i)
                  end do
               end do

               skipcount=0
               if (dispcount .ge. ndisp) then
                  call cpu_time(sec1)
                  sec00=sec00 + (sec1-sec0)
                  sec=sec00
                  sec0 = sec1
                  print*,  isave, " out of " , nsave,
     &                 " for time ", floor(sec)
c                 print*, beta
c                  print*, gamma


                  dispcount=0
               end if
            end if
         end if
         
      end do

      do j=1, q
         do i=1, ngrid
            f(i, j) = f(i,j)/dble(nsave)
         end do
      end do
      
c      print*, 'before return'
c      print*, betasave(1000,3,2)

      return
      end

c======================================================================
      subroutine loglik_mpt(maxm, q, nsub, parti, whicho, whichn, 
     &     b, bz, alpha, detlogl, linf, lsup, 
     &     mu, Sigma, tmp, ortho, mdzero, loglik)

C======================================================================
c     given v, mu, Sigma, Ortho, 
C     Calculate f(v1, ..., vn; mu, Sigma, O)
c======================================================================
      implicit none

c     Input
      integer maxm, q, nsub, mdzero
      integer parti(q)
      integer whicho(nsub), whichn(nsub)
      real*8 b(nsub, q), bz(nsub, q), alpha, detlogl
      real*8 linf(q), lsup(q)
      real*8 mu(q), Sigma(q,q), ortho(q,q)
      real*8 tmp(q)

c     output
      real*8 loglik

c     function
      integer ihmssf

c     working
      integer countero, countern, final
      integer i,j,je2,k,k1,k2,l,nint,ok
      real*8 dnrm, invcdfnorm, prob, quan, tmp1
      real*8 workmhr(q*(q+1)/2), workchol(q,q), sigmainv(q,q), lo(q,q)

C     routine
c     inital
    
      loglik=0.d0

      do i=1, q
         do j=1,q 
            workchol(i,j)=0.d0
            sigmainv(i,j)=0.d0
            lo(i,j)=0.d0
         end do
      end do

      do i=1,q
         tmp(i)=0.d0
      end do

c      print*, ortho
c      print*, 'before call cholesky'
c      print*, Sigma

c     get (LO)^-1, where Simga=LL'
      call cholesky(q, Sigma, workmhr)

c      print*, 'after call cholesky'

      do i=1,q
         do j=1,i
            workchol(i,j)=workmhr(ihmssf(i,j,q))
         end do
      end do

      do i=1,q
         do j=1,q
            tmp1=0.d0
            do k=1, q
               tmp1=tmp1+workchol(i,k)*ortho(k,j)
            end do
            lo(i,j)=tmp1
         end do
      end do

c     print*, "get LO"
c      print*, lo

c     so far , inverse sigmainv (2,2) by hand from LO
      tmp1=lo(1,1)*lo(2,2)-lo(1,2)*lo(2,1)
      sigmainv(1,1)=1.d0/tmp1*lo(2,2)
      sigmainv(1,2)=-1.d0/tmp1*lo(2,1)
      sigmainv(2,1)=-1.d0/tmp1*lo(1,2)
      sigmainv(2,2)=1.d0/tmp1*lo(1,1)
      
c      print*, sigmainv
c      print*, "before loglipt_mucan"

      call loglikpt_mucan(maxm, q, nsub, parti,
     &     whicho, whichn, b, bz, alpha, detlogl,
     &     linf, lsup, mu, sigmainv, 
     &     tmp, mdzero, loglik)

      return
      end

c======================================================================            
      subroutine rhaar2(xwork,x,n,q)
c======================================================================      
c     subroutine to generate a orthogonal random matrix with
c     haar distribution given a matrix of normally distributed
c     elements. This can be used for random walks in x.
c
c     A.J.V., 2008
c======================================================================      
      implicit none
      integer n
      real*8 x(n,n),q(n,n)
      real*8 xwork(n,n)
  
      integer i,j,k
      real*8 ck,dk,scale,sigma,sums,tau   

      do i=1,n
         do j=1,n
            xwork(i,j)=x(i,j)
            q(i,j)=0.d0
         end do
         q(i,i)=1.d0
      end do

      do k=1,n-1
         scale=0.d0
         do i=1,k
            scale=max(scale,abs(xwork(i,k)))
         end do
         if(scale.ne.0.d0)then
            sums=0.d0  
            do i=k,n
               sums=sums+xwork(i,k)**2
            end do  
            sigma=sign(sqrt(sums),xwork(k,k)) 
            xwork(k,k)=xwork(k,k)+sigma
            ck=sigma*xwork(k,k) 

            dk=-sign(1.d0,sigma)

            do j=k+1,n
               sums=0.d0
               do i=k,n
                  sums=sums+xwork(i,k)*xwork(i,j) 
               end do
               tau=sums/ck
               do i=k,n
                  xwork(i,j)=dk*(xwork(i,j)-tau*xwork(i,k))
               end do
            end do

            do i=1,n
               sums=0.d0
               do j=k,n
                  sums=sums+q(i,j)*xwork(j,k)
               end do
               tau=sums/ck
               do j=k,n
                  q(i,j)=dk*(q(i,j)-tau*xwork(j,k)) 
               end do
            end do 

         end if  
      end do

      if(xwork(n,n).lt.0.d0)then
         do i=1,n
            q(i,n)=-q(i,n)
         end do
      end if  

      return
      end 

c=======================================================================      
      subroutine cholesky(n,a,l)
c=======================================================================      
c     Subroutine to do a Double precision Half stored CHOLesky
c     decomposition.  Calculate Cholesky decomposition of symmetric,
c     positive definite matrix A which is LOWER HALF-STORED.  The matrix
c     l is the output.
c     A.J.V., 2005
      implicit none
      integer n,i,ii,j,jj,k,kk
      real*8 a(n,n),l(n*(n+1)/2)
      real*8 aii,scal,rtemp

      jj=0
      do i=1,n
         do j=i,n
            jj=jj+1
            l(jj)=a(i,j)
         end do
      end do   
      
      ii = 1
      do i = 1,n-1
         aii = l(ii)
         if (aii.le.0.d0) then
            call rexit("matrix is not pd in chol subroutine")
            return
         end if
         aii = sqrt(aii)
         l(ii) = aii
         scal = 1.d0/aii
         do j = ii+1,ii+n-i
            l(j) = scal*l(j)
         end do

         jj = 1
         do j = ii+1,ii+n-i
            if (l(j).ne.0.d0) then
               rtemp = -1.d0 * l(j)
               kk = ii + jj + n - i
               do k = j,ii+n-i
                  l(kk) = l(kk) + l(k)*rtemp
                  kk = kk + 1
               end do
            end if
            jj = jj + n - i - j + ii + 1
         end do
         ii = ii + n - i + 1
      end do
      aii = l(ii)
      if (aii.le.0.d0) then
            call rexit("matrix is not pd in chol subroutine")
          return 
      end if
      aii = sqrt(aii)
      l(ii) = aii
      return
      end


c=======================================================================                  
      subroutine loglikpt_mucan(m,nrand,nsubject,parti,
     &                          whicho,whichn,b,bzc,cpar,detlogl,
     &                          linf,lsup,muc,sigmainv,
     &                          vec,fixed,loglikc)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the candidate 
c     value of the baseline mean in a marginal Multivariate PT.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 22-06-2008.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      real*8 b(nsubject,nrand),bzc(nsubject,nrand),cpar,detlogl
      real*8 linf(nrand),lsup(nrand)
      real*8 muc(nrand),sigmainv(nrand,nrand)
      real*8 vec(nrand)

c-----Output
      real*8 loglikc

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      real*8 dnrm,invcdfnorm,prob,quan,tmp1

c-----Routine

      loglikc=0.d0

      do i=1,nsubject
         do j=1,nrand
            tmp1=0.d0
            do k=1,nrand
               tmp1=tmp1+sigmainv(j,k)*(b(i,k)-muc(k))
            end do
            vec(j)=tmp1
         end do
         
         do j=1,nrand
            bzc(i,j)=vec(j)
         end do

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ first subject
         if(i.eq.1)then
            loglikc=-0.5d0*detlogl
            do j=1,nrand
               loglikc=loglikc+dnrm(bzc(i,j),0.d0, 1.d0, 1)
            end do   

c+++++++ following subjects
          else

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,0.d0,1.d0,1,0)

            countero=0
            
            do j=1,nrand
               if(bzc(i,j).le.quan)then
                  linf(j)=-999999.d0
                  lsup(j)=quan
                  parti(j)=1
                else
                  linf(j)=quan
                  lsup(j)= 999999.d0
                  parti(j)=2
               end if
            end do
           
            do l=1,i-1
               final=1
               do j=1,nrand
                  if(bzc(l,j).gt.lsup(j).or.bzc(l,j).lt.linf(j))then
                    final=0
                  end if
               end do
               
               if(final.eq.1)then
                  countero=countero+1
                  whicho(countero)=l
               end if   
            end do

            if(fixed.ne.1)then 
              loglikc=loglikc+
     &         log((2.d0**nrand)*cpar+dble(2.d0**nrand)*dble(countero))-
     &         log((2.d0**nrand)*cpar+dble(i-1))
            end if 

            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1.and.j.le.m)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)

               do k=1,nrand
                  k1=2*(parti(k)-1)+1
                  k2=2*(parti(k)-1)+2
                  quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                  if(bzc(i,k).le.quan)then
                    parti(k)=k1 
                    lsup(k)=quan
                   else 
                    parti(k)=k2
                    linf(k)=quan
                  end if
               end do                 
               
               countern=0
               do l=1,countero
                  final=1
                  do k=1,nrand
                     if(bzc(whicho(l),k).gt.lsup(k).or.
     &                  bzc(whicho(l),k).lt.linf(k)    )then
                        final=0 
                     end if   
                  end do
                  
                  if(final.eq.1)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                  end if
               end do

               loglikc=loglikc+
     &           log((2.d0**nrand)*cpar*dble(je2)+
     &               dble(2.d0**nrand)*dble(countern))-
     &           log((2.d0**nrand)*cpar*dble(je2)+dble(countero))

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

1           continue

            loglikc=loglikc-0.5d0*detlogl
            do j=1,nrand
               loglikc=loglikc+dnrm(bzc(i,j),0.d0, 1.d0, 1)
            end do   
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
