      subroutine myptdensitybi(maxm, nsub, q, y, 
     &     mu, Sigma,
     &     nburn, nsave, nskip, arate,  
     &     grid, ngrid, f, alpha,
     &     musave, sigmasave, alphasave, fsave,
     &     ratesave, tunesave)

c$$$Time-stamp: <liuminzhao 11/08/2012 15:07:57>

      implicit none

C     PT PARAMETER
      integer maxm, mdzero

C     DATA
      integer nsub, q
      real*8 y(nsub, q), v(nsub, q), vc(nsub, q), vz(nsub, q)

C     BASELINE
      real*8 mu(2), Sigma(q, q), alpha, muc(2), Sigmac(2, 2)
      real*8 alphac, sigma2, sigma2c
      real*8 sd, sd2, theta, thetac, alpha2, alpha2c
      real*8 sigmavec(3), sigmavecc(3)

C     MCMC PARAMETER
      integer nburn, nsave, nskip
      real*8 arate

C     GRID PARAMETER
      integer ngrid
      real*8 grid(ngrid)

C     OUTPUT 
      real*8 f(ngrid, q)

C     SAVE
      real*8 musave(nsave, 2), sigmasave(nsave, 3), alphasave(nsave)
      real*8 fsave(nsave, ngrid, q)
      real*8 ratesave(nburn, 7), tunesave(nburn, 7)

C     TMP
      integer i, j, k, l, m
      real*8 loglikeo, loglikec, logprioro, logpriorc, logcgkc, logcgko
      integer whicho(nsub), whichn(nsub)
      integer iscan, nscan, skipcount, isave
      real*8 tau(2)

      integer parti(q)
      real*8 ortho(q, q), orthoc(q, q)
      real*8 workmr1(q, q), workmr2(q, q), workmr2c(q, q)
      real*8 detlogl, detloglc
      real*8 linf(q), lsup(q)
      real*8 tmp1, tmp2, tmp(q), tmpnsub(nsub)

C     FUNCTION
      real*8 rnorm, ratio, runif, myrnorm, myrunif, dnrm, dgamma2

C     TUNE 
      integer att(7), acc(7)
      real*8 tune(7)

C     INITIAL
      mdzero = 1

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
      end do

      do i=1, nsub
         whicho(i)=0
         whichn(i)=0
      end do

      do i=1, nsub
         do j=1, q
            vz(i,j)=0.d0
         end do
      end do

      alpha = 1.d0

      tau(1) = 0.1
      tau(2) = 0.1

      do i=1,7
         att(i) = 0
         acc(i) = 0
         tune(i) = 0.1
      end do
      
      sigmavec(1)=sqrt(Sigma(1,1))
      sigmavec(2)=sqrt(Sigma(2,2))
      sigmavec(3)=Sigma(1,2)/sigmavec(1)/sigmavec(2)

C     MCMC

      isave = 0
      skipcount = 0
      nscan = nburn + (nskip + 1) * nsave

C     FIRST 
c      print*, "first"

      ortho(1,1)=1.d0
      ortho(1,2)=0.d0
      ortho(2,1)=0.d0
      ortho(2,2)=1.d0

      loglikeo=0.d0

C     log|Sigma|
      detlogl=log(Sigma(1,1)*Sigma(2,2)-Sigma(1,2)*Sigma(2,1))

      call loglik_mpt(maxm, q, nsub, parti,
     &     whicho, whichn, y, vz, alpha, detlogl,
     &     linf, lsup, mu, Sigma, 
     &     tmp, ortho, mdzero, loglikeo)

C     ROLLING

      do iscan = 1, nscan

C     MU
c         print*, mu

         do l = 1, 2
            att(l) = att(l) + 1
            muc(1) = mu(1)
            muc(2) = mu(2)
            muc(l) = myrnorm(mu(l), tune(l))

            logprioro = dnrm(mu(l), 0.d0, 30.d0, 1)
            logpriorc = dnrm(muc(l), 0.d0, 30.d0, 1)

            loglikec = 0.d0

            call loglik_mpt(maxm, q, nsub, parti,
     &           whicho, whichn, y, vz, alpha, detlogl,
     &           linf, lsup, muc, Sigma, 
     &           tmp, ortho, mdzero, loglikec)

            ratio=loglikec + logpriorc - loglikeo- logprioro

            if ( ratio .gt. 0 .or. log(dble(myrunif(0.d0, 1.d0))).lt. 
     &           ratio) then 
               loglikeo = loglikec
               mu(l) = muc(l)
               acc(l) = acc(l) + 1
            end if
         end do

C     SIGMA
c         print*, "sigma"
c         print*, Sigma

c     Sigma(1,1)
         att(3) = att(3) + 1
         theta = log(sigmavec(1)) 
         thetac = myrnorm(theta, tune(3))
         logcgkc = -theta
         logcgko = -thetac
         Sigmac(1,1)=exp(2*thetac)
         Sigmac(2,2)=sigmavec(2)**2
         Sigmac(1,2)=exp(thetac)*sigmavec(2)*sigmavec(3)
         Sigmac(2,1)=exp(thetac)*sigmavec(2)*sigmavec(3)

         detloglc=log(Sigmac(1,1)*Sigmac(2,2)-Sigmac(1,2)**2)

         loglikec = 0.d0

c         print*, Sigmac

         call loglik_mpt(maxm, q, nsub, parti,
     &        whicho, whichn, y, vz, alpha, detloglc,
     &        linf, lsup, mu, Sigmac, 
     &        tmp, ortho, mdzero, loglikec)

         logpriorc=dgamma2(exp(thetac), 1.d0, 1.d0,1)
         logprioro=dgamma2(exp(theta), 1.d0, 1.d0,1)

c         logpriorc=-tau(1)*thetac - tau(2)*exp(-2*thetac)/2
c         logprioro=-tau(1)*theta-tau(2)*exp(-2*theta)/2
         ratio=loglikec + logpriorc -loglikeo -logprioro+
     &        logcgkc - logcgko

         if (ratio .gt. 0 .or. log(dble(myrunif(0.d0, 1.d0))).lt. ratio) 
     &        then 
            loglikeo=loglikec
            detlogl=detloglc
            sigmavec(1)=exp(thetac)
            Sigma(1,1)=Sigmac(1,1)
            Sigma(1,2)=Sigmac(1,2)
            Sigma(2,1)=Sigmac(2,1)
            acc(3) = acc(3) + 1
         end if

c     Sigma(2,2)


         att(4) = att(4) + 1
         theta = log(sigmavec(2)) 
         thetac = myrnorm(theta, tune(4))
         logcgkc = -theta
         logcgko = -thetac
         Sigmac(1,1)=sigmavec(1)**2
         Sigmac(2,2)=exp(2*thetac)
         Sigmac(1,2)=exp(thetac)*sigmavec(1)*sigmavec(3)
         Sigmac(2,1)=exp(thetac)*sigmavec(1)*sigmavec(3)

         detloglc=log(Sigmac(1,1)*Sigmac(2,2)-Sigmac(1,2)**2)

         loglikec = 0.d0

c         print*, Sigmac

         call loglik_mpt(maxm, q, nsub, parti,
     &        whicho, whichn, y, vz, alpha, detloglc,
     &        linf, lsup, mu, Sigmac, 
     &        tmp, ortho, mdzero, loglikec)


         logpriorc=dgamma2(exp(thetac), 1.d0, 1.d0,1)
         logprioro=dgamma2(exp(theta), 1.d0, 1.d0,1)

c         logpriorc=-tau(1)*thetac - tau(2)*exp(-2*thetac)/2
c         logprioro=-tau(1)*theta-tau(2)*exp(-2*theta)/2

         ratio=loglikec + logpriorc -loglikeo -logprioro+
     &        logcgkc - logcgko

         if (ratio .gt. 0 .or. log(dble(myrunif(0.d0, 1.d0))).lt. ratio) 
     &        then 
            loglikeo=loglikec
            detlogl=detloglc
            sigmavec(2)=exp(thetac)
            Sigma(2,2)=Sigmac(2,2)
            Sigma(1,2)=Sigmac(1,2)
            Sigma(2,1)=Sigmac(2,1)
            acc(4) = acc(4) + 1
         end if

c     Sigma(1,2)

c$$$         att(5) = att(5) + 1
c$$$         theta=sigmavec(3)
c$$$         thetac=myrnorm(theta, tune(3))
c$$$         thetac=max(thetac, -0.99)
c$$$         thetac=min(thetac, 0.99)
c$$$         Sigmac(1,1)=sigmavec(1)**2
c$$$         Sigmac(2,2)=sigmavec(2)**2
c$$$         Sigmac(1,2)=thetac*sigmavec(1)*sigmavec(2)
c$$$         Sigmac(2,1)=thetac*sigmavec(1)*sigmavec(2)
c$$$
c$$$         detloglc=log(Sigmac(1,1)*Sigmac(2,2)-Sigmac(1,2)**2)
c$$$
c$$$         loglikec = 0.d0
c$$$
c$$$         call loglik_mpt(maxm, q, nsub, parti,
c$$$     &        whicho, whichn, y, vz, alpha, detloglc,
c$$$     &        linf, lsup, mu, Sigmac, 
c$$$     &        tmp, ortho, mdzero, loglikec)
c$$$
c$$$         ratio=loglikec  -loglikeo
c$$$
c$$$         if (ratio .gt. 0 .or. log(dble(myrunif(0.d0, 1.d0))).lt. ratio) 
c$$$     &        then 
c$$$            loglikeo=loglikec
c$$$            detlogl=detloglc
c$$$            sigmavec(3)= thetac
c$$$            Sigma(1,2)=Sigmac(1,2)
c$$$            Sigma(2,1)=Sigmac(2,1)
c$$$            acc(5) = acc(5) + 1
c$$$         end if


C     ALPHA
c         print*, "alpha"

         att(6) = att(6) + 1
         theta = log(alpha) 
         thetac = myrnorm(theta, tune(6))
         logcgkc = -theta
         logcgko = -thetac
         alphac = exp(thetac)

         loglikec = 0.d0

         call loglik_mpt(maxm, q, nsub, parti,
     &        whicho, whichn, y, vz, alphac, detlogl,
     &        linf, lsup, mu, Sigma, 
     &        tmp, ortho, mdzero, loglikec)


         logpriorc=dgamma2(alphac, 1.d0, 1.d0, 1)
         logprioro=dgamma2(alpha, 1.d0, 1.d0, 1)

         ratio=loglikec + logpriorc -loglikeo -logprioro+
     &        logcgkc - logcgko

         if (ratio .gt. 0 .or. log(dble(myrunif(0.d0, 1.d0))).lt. ratio) 
     &        then 
            loglikeo=loglikec
            alpha = alphac
            acc(6) = acc(6) + 1
         end if

c=========================================
c     Ortho 
c========================================

c         print*, "ortho"

         att(7) = att(7) + 1
c$$$         do i=1,q
c$$$            do j=1,q   
c$$$               workmr2c(i,j)=myrnorm(workmr2(i,j),tune(7))
c$$$            end do
c$$$         end do
c$$$         call rhaar2(workmr1,workmr2,q,orthoc)
c$$$
c$$$         call loglik_mpt(maxm, q, nsub, parti,
c$$$     &        whicho, whichn, y, vz, alpha, detlogl,
c$$$     &        linf, lsup, mu, Sigma, 
c$$$     &        tmp, orthoc, mdzero, loglikec)
c$$$
c$$$c     acceptance
c$$$         ratio=loglikec-loglikeo
c$$$
c$$$         if(ratio .gt. 0 .or. log(dble(myrunif(0.d0, 1.d0)))
c$$$     &        .lt.ratio) then
c$$$            acc(7) = acc(7) + 1
c$$$            do i=1,q
c$$$               do j=1,q
c$$$                  ortho(i,j)=orthoc(i,j)
c$$$                  workmr2(i,j)=workmr2c(i,j)
c$$$               end do
c$$$            end do
c$$$            loglikeo=loglikec
c$$$         end if


C     TUNING
         if ((att(1).ge.50).and.(iscan.le. nburn)) then 
            do i=1, 7
               if (dble(acc(i))/dble(att(i)) .gt. arate) then 
                  tune(i)=tune(i) + 
     &                 min(0.01d0,dble(10)/sqrt(dble(iscan)))
                 else
                  tune(i)=tune(i)-
     &                  min(0.01d0,dble(10)/sqrt(dble(iscan)))
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
               
               do i = 1, 2
                  musave(isave, i) = mu(i)
               end do

               do i = 1, 3
                  sigmasave(isave, i) = sigmavec(i)
               end do

               alphasave(isave) = alpha

C     density estimate

               do j = 1, 2
                  do i = 1, nsub
                     tmpnsub(i) = y(i, j)
                  end do
                  
                  do i=1,ngrid
                     loglikec=0.d0
                     call gridupptprior(grid(i),maxm,mdzero ,nsub,
     &                    alpha, mu(j), Sigma(j, j), tmpnsub,
     &                    whicho,whichn,loglikec)
                     f(i, j) = f(i, j) + exp(loglikec)  
                     fsave(isave, i, j) = exp(loglikec)
                  end do
               end do

               skipcount=0
               
            end if
         end if
      end do
      
      do j = 1, 2
         do i=1, ngrid
            f(i, j)=f(i, j)/dble(nsave)
         end do
      end do
      
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

c$$$      if(mdzero.ne.0)then 
c$$$         logprioro=logprioro+
c$$$     &    log(2.d0*alpha+dble(2*countero))-
c$$$     &    log(2.d0*alpha+dble(nsubject))
c$$$      end if 

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
