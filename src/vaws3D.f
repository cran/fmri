CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   extract element of 3D array
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function getlwght(lwght,dw1,dw2,dw3,j1,j2,j3)
      integer dw1,dw2,dw3,j1,j2,j3
      real*8 lwght(dw1,dw2,dw3)
      getlwght=lwght(j1,j2,j3)
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute Location Kernel (Compact support only, based on x^2
C                                   ignores scaling)
C
C          Kern=1     Plateau
C          Kern=2     Epanechnicov
C          Kern=3     Gaussian
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lkern(kern,xsq)
      implicit logical (a-z)
      integer kern
      real*8 xsq
      IF (xsq.ge.1) THEN
         lkern=0.d0
      ELSE IF (kern.eq.1) THEN
         IF(xsq.le.0.5d0) THEN
            lkern=1.d0
         ELSE
            lkern=2.d0*(1.d0-xsq)
         END IF
      ELSE IF (kern.eq.2) THEN
         lkern=1.d0-xsq
      ELSE IF (kern.eq.3) THEN
         lkern=exp(-xsq*8.d0)
      ELSE IF (kern.eq.4) THEN
         lkern=exp(-xsq*18.d0)
      ELSE
C        use Epanechnikov
         lkern=1.d0-xsq
      ENDIF
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute aws-weights  w_{ij}
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awswghts(n1,n2,n3,j1,j2,j3,dv0,thi,theta,
     1                    vwghts,skern,spf,spmin,spmax,bii,wj)
      implicit logical (a-z)
      integer n1,n2,n3,j1,j2,j3,dv0,skern
      real*8 thi(dv0),theta(n1,n2,n3,dv0),vwghts(dv0),spf,spmin,spmax,
     1       bii,wj,wjin
      integer k
      real*8 sij,z
      wjin=wj
      sij=0.d0
C  compute distance in sij
      DO k=1,dv0
         z=thi(k)-theta(j1,j2,j3,k)
         sij=sij+z*z*vwghts(k)
      END DO
      sij=bii*sij
      IF (sij.gt.spmax) THEN
         wj=0.d0
      ELSE IF (skern.eq.1) THEN
C  skern == "Plateau"
         wj=wj*min(1.d0,1.d0-spf*(sij-spmin))
      ELSE IF (skern.eq.2) THEN
C  skern == "Triangle"
         wj=wj*(1.d0-sij)
      ELSE
C  skern == "Exp"
         IF (sij.gt.spmin) wj=wj*exp(-spf*(sij-spmin))
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chaws2(y,si2,mask,wlse,n1,n2,n3,dv,dv0,hakt,lambda,
     1                  theta,bi,thn,kern,skern,spmin,spmax,lwght,
     2                  wght,vwghts,swjy,thi)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      integer n1,n2,n3,kern,skern,dv,dv0
      logical aws,wlse,mask(n1,n2,n3)
      real*8 y(n1,n2,n3,dv),theta(n1,n2,n3,dv0),bi(n1,n2,n3),
     1       thn(n1,n2,n3,dv),lambda,spmax,wght(2),si2(n1,n2,n3),
     1       hakt,lwght(*),spmin,vwghts(dv0),thi(dv0),getlwght
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n
      real*8 bii,swj,swjy(dv),wj,hakt2,spf,si2j,si2i,swjv
      external getlwght
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               if(mask(i1,i2,i3)) THEN
                  DO k=1,dv
                     thn(i1,i2,i3,k)=0.d0
                  END DO
                  CYCLE
               END IF
               si2i=si2(i1,i2,i3)
               bii=bi(i1,i2,i3)/lambda
C   scaling of sij outside the loop
               swj=0.d0
               swjv=0.d0
               DO k=1,dv
                  swjy(k)=0.d0
               END DO
               DO k=1,dv0
                  thi(k)=theta(i1,i2,i3,k)
               END DO
               DO jw3=1,dlw3
                  j3=jw3-clw3+i3
                  if(j3.lt.1.or.j3.gt.n3) CYCLE
                  DO jw2=1,dlw2
                     j2=jw2-clw2+i2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     DO jw1=1,dlw1
C  first stochastic term
                        j1=jw1-clw1+i1
                        if(j1.lt.1.or.j1.gt.n1) CYCLE
                        IF(mask(j1,j2,j3)) CYCLE
                        wj=getlwght(lwght,dlw1,dlw2,dlw3,jw1,jw2,jw3)
                        if(wj.le.0.d0) CYCLE
                        si2j=si2(j1,j2,j3)
                        IF (aws) THEN
                           call awswghts(n1,n2,n3,j1,j2,j3,dv0,thi,
     1                     theta,vwghts,skern,spf,spmin,spmax,bii,wj)
                        END IF
                        if(wlse) THEN 
                           wj=wj*si2j
                        ELSE
                           swjv=swjv+wj/si2j
                        END IF
                        swj=swj+wj
                        DO k=1,dv
                           swjy(k)=swjy(k)+wj*y(j1,j2,j3,k)
                        END DO
                     END DO
                  END DO
               END DO
               DO k=1,dv
                  thn(i1,i2,i3,k)=swjy(k)/swj
               END DO
               IF(wlse) THEN
                  bi(i1,i2,i3)=swj
               ELSE
                  bi(i1,i2,i3)=swj*swj/swjv
               END IF
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chawsv(y,res,si2,mask,wlse,n1,n2,n3,n4,dv,dv0,hakt,
     1                  lambda,theta,bi,resnew,thn,kern,skern,spmin,
     2                  spmax,lwght,wght,vwghts,swjy,thi,resi)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      integer n1,n2,n3,n4,kern,skern,dv,dv0
      logical aws,wlse,mask(n1,n2,n3)
      real*8 res(n4,n1,n2,n3),y(n1,n2,n3,dv),theta(n1,n2,n3,dv0),
     1       bi(n1,n2,n3),thn(n1,n2,n3,dv),lambda,spmax,wght(2),
     1       si2(n1,n2,n3),hakt,lwght(*),spmin,vwghts(dv0),thi(dv0),
     1       resi(n4),getlwght,resnew(n4,n1,n2,n3)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n
      real*8 bii,swj,swjy(dv),wj,hakt2,spf,si2j,si2i,swjv,sresisq
      external getlwght
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               if(mask(i1,i2,i3)) THEN
                  DO k=1,dv
                     thn(i1,i2,i3,k)=0.d0
                  END DO
                  CYCLE
               END IF
               si2i=si2(i1,i2,i3)
               bii=bi(i1,i2,i3)/lambda
C   scaling of sij outside the loop
               swj=0.d0
               swjv=0.d0
               DO k=1,dv
                  swjy(k)=0.d0
               END DO
               DO k=1,n4
                  resi(k)=0.d0
               END DO
               DO k=1,dv0
                  thi(k)=theta(i1,i2,i3,k)
               END DO
               DO jw3=1,dlw3
                  j3=jw3-clw3+i3
                  if(j3.lt.1.or.j3.gt.n3) CYCLE
                  DO jw2=1,dlw2
                     j2=jw2-clw2+i2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     DO jw1=1,dlw1
C  first stochastic term
                        j1=jw1-clw1+i1
                        if(j1.lt.1.or.j1.gt.n1) CYCLE
                        IF(mask(j1,j2,j3)) CYCLE
                        wj=getlwght(lwght,dlw1,dlw2,dlw3,jw1,jw2,jw3)
                        if(wj.le.0.d0) CYCLE
                        si2j=si2(j1,j2,j3)
                        IF (aws) THEN
                           call awswghts(n1,n2,n3,j1,j2,j3,dv0,thi,
     1                     theta,vwghts,skern,spf,spmin,spmax,bii,wj)
                        END IF
                        if(wlse) THEN 
                           wj=wj*si2j
                        ELSE
                           swjv=swjv+wj/si2j
                        END IF
                        swj=swj+wj
                        DO k=1,dv
                           swjy(k)=swjy(k)+wj*y(j1,j2,j3,k)
                        END DO
                        call daxpy(n4,wj,res(1,j1,j2,j3),1,resi,1)
                     END DO
                  END DO
               END DO
               DO k=1,dv
                  thn(i1,i2,i3,k)=swjy(k)/swj
               END DO
               sresisq=0.d0
               DO k=1,n4
                  resnew(k,i1,i2,i3)=resi(k)/swj
                  sresisq=sresisq+resi(k)*resi(k)
               END DO
               bi(i1,i2,i3)=swj*swj/sresisq*(n4-1.d0)
C  thats the inverse variance
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ihaws2(y,si2,mask,wlse,n1,n2,n3,dv,dv0,hakt,lambda,
     1                  theta,bi,thn,kern,skern,spmin,spmax,lwght,
     2                  wght,vwghts,swjy,thi)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi       (input)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      integer dv,dv0,n1,n2,n3,kern,skern
      logical aws,wlse,mask(n1,n2,n3)
      real*8 theta(n1,n2,n3,dv0),bi(n1,n2,n3),y(dv,n1,n2,n3),
     1       lambda,spmax,wght(2),si2(n1,n2,n3),thn(dv,n1,n2,n3),
     1       hakt,lwght(*),spmin,vwghts(dv0),thi(dv0),getlwght
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n
      real*8 bii,swj,swjy(dv),wj,hakt2,spf,si2j,si2i
      external getlwght
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               if(mask(i1,i2,i3)) THEN
                  DO k=1,dv
                     thn(k,i1,i2,i3)=0.d0
                  END DO
                  CYCLE
               END IF
               si2i=si2(i1,i2,i3)
               bii=bi(i1,i2,i3)/lambda
C   scaling of sij outside the loop
               swj=0.d0
               DO k=1,dv
                  swjy(k)=0.d0
               END DO
               DO k=1,dv0
                  thi(k)=theta(i1,i2,i3,k)
               END DO
               DO jw3=1,dlw3
                  j3=jw3-clw3+i3
                  if(j3.lt.1.or.j3.gt.n3) CYCLE
                  DO jw2=1,dlw2
                     j2=jw2-clw2+i2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     DO jw1=1,dlw1
C  first stochastic term
                        j1=jw1-clw1+i1
                        if(j1.lt.1.or.j1.gt.n1) CYCLE
                        IF(mask(j1,j2,j3)) CYCLE
                        wj=getlwght(lwght,dlw1,dlw2,dlw3,jw1,jw2,jw3)
                        if(wj.le.0.d0) CYCLE
                        si2j=si2(j1,j2,j3)
                        IF (aws) THEN
                           call awswghts(n1,n2,n3,j1,j2,j3,dv0,thi,
     1                     theta,vwghts,skern,spf,spmin,spmax,bii,wj)
                        END IF
                        if(wlse) THEN 
                           wj=wj*si2j
                        END IF
                        swj=swj+wj
                        call daxpy(dv,wj,y(1,j1,j2,j3),1,swjy,1)
                     END DO
                  END DO
               END DO
               DO k=1,dv
                  thn(k,i1,i2,i3)=swjy(k)/swj
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smooth3d(y,si2,mask,wlse,n1,n2,n3,dv,hakt,
     1                    thn,kern,lwght,wght,swjy)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      integer n1,n2,n3,kern,dv
      logical wlse,mask(n1,n2,n3)
      real*8 y(n1,n2,n3,dv),thn(n1,n2,n3,dv),wght(2),
     1       si2(n1,n2,n3),hakt,lwght(*),getlwght
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n
      real*8 swj,swjy(dv),wj,hakt2
      external getlwght
      hakt2=hakt*hakt
      ih1=hakt
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               if(mask(i1,i2,i3)) THEN
                  DO k=1,dv
                     thn(i1,i2,i3,k)=0.d0
                  END DO
                  CYCLE
               END IF
C   scaling of sij outside the loop
               swj=0.d0
               DO k=1,dv
                  swjy(k)=0.d0
               END DO
               DO jw3=1,dlw3
                  j3=jw3-clw3+i3
                  if(j3.lt.1.or.j3.gt.n3) CYCLE
                  DO jw2=1,dlw2
                     j2=jw2-clw2+i2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     DO jw1=1,dlw1
C  first stochastic term
                        j1=jw1-clw1+i1
                        if(j1.lt.1.or.j1.gt.n1) CYCLE
                        IF(mask(j1,j2,j3)) CYCLE
                        wj=getlwght(lwght,dlw1,dlw2,dlw3,jw1,jw2,jw3)
                        if(wj.le.0.d0) CYCLE
                        if(wlse) THEN 
                           wj=wj*si2(j1,j2,j3)
                        END IF
                        swj=swj+wj
                        DO k=1,dv
                           swjy(k)=swjy(k)+wj*y(j1,j2,j3,k)
                        END DO
                     END DO
                  END DO
               END DO
               DO k=1,dv
                  thn(i1,i2,i3,k)=swjy(k)/swj
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
C
C   calculate location weights in lwght
C
      subroutine locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      implicit logical (a-z)
      integer dlw1,dlw2,dlw3,kern
      real*8 wght(2),hakt2,lwght(dlw1,dlw2,dlw3),lkern
      external lkern
      real*8 z1,z2,z3
      integer j1,j2,j3,clw1,clw2,clw3,ih1,ih2
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
      DO j3=1,dlw3
         Do j2=1,dlw2
            DO j1=1,dlw1
               lwght(j1,j2,j3)=0.d0
            END DO
         END DO
         z3=(clw3-j3)*wght(2)
         z3=z3*z3
         ih2=sqrt(hakt2-z3)/wght(1)
         DO j2=clw2-ih2,clw2+ih2
            z2=(clw2-j2)*wght(1)
            z2=z3+z2*z2
            ih1=sqrt(hakt2-z2)
            DO j1=clw1-ih1,clw1+ih1
               z1=clw1-j1
               lwght(j1,j2,j3)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      RETURN
      END
C
C
C
