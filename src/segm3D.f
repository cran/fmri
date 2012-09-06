CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine segm3d(y,res,si2,mask,wlse,n1,n2,n3,nt,df,hakt,
     1                  lambda,theta,bi,thn,lwght,wght,swres,pval,
     3                  segm,delta,thresh,fov,vq,vest0i,varest,
     4                  restrict)
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
      integer n1,n2,n3,nt,kern,segm(n1,n2,n3)
      logical aws,wlse,mask(n1,n2,n3),restrict
      real*8 y(n1,n2,n3),theta(n1,n2,n3),bi(n1,n2,n3),delta,thresh,
     1      thn(n1,n2,n3),lambda,wght(2),si2(n1,n2,n3),pval(n1,n2,n3),
     1      hakt,lwght(*),thi,getlwght,swres(nt),fov,vq(n1,n2,n3),
     1      varest(n1,n2,n3),res(nt,n1,n2,n3),vest0i(n1,n2,n3),df,
     1      kv(n1,n2,n3)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,segmi
      real*8 bii,swj,swjy,wj,hakt2,spf,si2j,si2i,vqi,
     1       varesti,fpchisq,ti,thij,sij,z,si,swr,z1,
     2       a,b,dn,pvali
      external getlwght,fpchisq
      kern=1
      hakt2=hakt*hakt
      spf=4.d0/3.d0
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
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
      IF(hakt.gt.1.25) THEN
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               if(mask(i1,i2,i3)) THEN
                  thn(i1,i2,i3)=0.d0
                  CYCLE
               END IF
               vqi=vq(i1,i2,i3)
               thi=theta(i1,i2,i3)
               si2i=vest0i(i1,i2,i3)
               varesti=varest(i1,i2,i3)
               dn=varesti/si2i*fov
               call getdfnab(df,dn,a,b)
C   this should be more conservative using actual variance reduction instead of theoretical
               ti=max(0.d0,abs(thi)-delta)
               IF(a*ti/sqrt(varesti/vqi)-b.gt.thresh) THEN
                   pval(i1,i2,i3)=0.d0
               ELSE
                   pval(i1,i2,i3)=1.d0
               END IF
            END DO
         END DO
      END DO
      END IF 
C   scaling of sij outside the loop
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               if(mask(i1,i2,i3)) THEN
                  thn(i1,i2,i3)=0.d0
                  CYCLE
               END IF
               segmi=segm(i1,i2,i3)
               pvali=pval(i1,i2,i3)
               vqi=vq(i1,i2,i3)
               si2i=vest0i(i1,i2,i3)
               bii=bi(i1,i2,i3)/lambda
               swj=0.d0
               swjy=0.d0
               DO k=1,nt
                  swres(k)=0.d0
               END DO
               thi=theta(i1,i2,i3)
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
                           thij=thi-theta(j1,j2,j3)
                           sij=thij*thij*bii
                           if(restrict) THEN
C restrict smoothing within segmented areas
                           if(abs(segmi).eq.1) THEN
                           if(segmi*segm(j1,j2,j3).gt.0) THEN
C
C   allow for nonadaptive smoothing if values are in the same segment
C   since pvali << 1 adaptation is significantly reduced
C
                                 sij=pvali*sij
                              ELSE
C
C   no smoothing if voxel i is classified as 1 or -1 and 
C                   voxel i and j are in different segments
C
                                 CYCLE
                              END IF
                           END IF
                           END IF
C endif for restrict smoothing within segmented areas
                           IF(sij.gt.1.d0) CYCLE
                        IF(sij.gt.0.25d0) wj=wj*(1.d0-spf*(sij-0.25d0))
                        END IF
                        if(wlse)  wj=wj*si2j
                        swj=swj+wj
                        swjy=swjy+wj*y(j1,j2,j3)
C  weighted sum of residuals
                        call daxpy(nt,wj,res(1,j1,j2,j3),1,swres,1)
                     END DO
                  END DO
               END DO
               z=0.d0
               z1=0.d0
C
C   now calculate variance of estimates from smoothed residuals
C
               DO k=1,nt
                  swr=swres(k)/swj
                  z1=z1+swr
                  z=z+swr*swr
               END DO
               thi=swjy/swj
               z1=z1/nt
               si = (z/nt - z1*z1)
               if(restrict) THEN
C  smoothing restricted within segmented ares
               if(segmi.eq.1) THEN
                  if(thi.lt.theta(i1,i2,i3)) THEN
                     thi = theta(i1,i2,i3)
                  ELSE
                     varest(i1,i2,i3)=si
                     bi(i1,i2,i3)=si2i/si*si2(i1,i2,i3)
                  END IF
               END IF
               if(segmi.eq.-1) THEN
                  if(thi.gt.theta(i1,i2,i3)) THEN
                     thi = theta(i1,i2,i3)
                  ELSE
                     varest(i1,i2,i3)=si
                     bi(i1,i2,i3)=si2i/si*si2(i1,i2,i3)
                  END IF
               END IF 
               END IF
C  end if for smoothing restricted within segmented ares
               thn(i1,i2,i3)=thi
C               si = si/nt
               if(restrict.and.segmi.ne.0) CYCLE
               varest(i1,i2,i3)=si
               bi(i1,i2,i3)=si2i/si*si2(i1,i2,i3)
C   keep the detected segment
               dn=si/si2i*fov
               call getdfnab(df,dn,a,b)
C 
C   note that a and b refer to  1/a_n and b_n/a_n
C
C   this should be more conservative using actual variance reduction instead of theoretical
               si=sqrt(si/vqi)
               kv(i1,i2,i3)=a*(thi-delta)/si-b
C   thats the SD of thi
               if(a*(thi+delta)/si+b.lt.-thresh) THEN
                  segm(i1,i2,i3)=-1
               ELSE IF (a*(thi-delta)/si-b.gt.thresh) THEN
                  segm(i1,i2,i3)=1
               END IF
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine getdfnab(df,n,a,b)
C
C   this function computes approximations for constants a=1/a_n and b=b_n/a_n
C   such that for the maximum T_n of n r.v. from student t_df  
C   a_n T_n +b_n  ~ \Phi_\df    (asymp. extreme value distribution for t_df)
C
C   approximation formulaes obtained from samples of 100000 Extremes 
C   df \in 10:264   n \in  100 : 20000   
C   max. approx error < 0.002 for a   and < 0.006   for b
C
      implicit logical (a-z)
      real*8 a,b,df,n
      real*8 dfinv,ldf,dfq,ninvh,x1,x2,x3,x4,x5,x6,lna,lnb,ln
      dfinv=1.d0/(df-1.d0)
      dfq=sqrt(sqrt(df))
      ldf=log(df+1.d1)
      ninvh=exp(-.2d0*log(n+8.d0))
      ln=log(n)
      lna=exp(0.01*log(ln))
      lnb=exp(1.85*log(ln))
      x1=1.d0/(df+lnb)
      x2=df/(df+lna)
      x3=df/lna
      x4=1.d0/(df+lna)
      x5=df/(df+lnb)
      x6=df/lnb
      a=-7.754959d1-1.210474d1*dfinv+2.996616d-2*dfq-3.428621d-2*ldf-
     -   3.702433d-3*ln+2.354546*lna+1.055396d-4*lnb-4.086823*x1+
     +   7.527587d1*x2-1.180929d-5*x3+9.221818d+1*x4-2.438421d-2*x5+
     +   6.206907d-5*x6
      b= 1.28747d3-2.657702d1*dfinv-1.229952d-1*dfq+1.483752d-1*ldf+
     +   1.909285d-2*ln-2.286336d+1*lna-3.68511d-4*lnb+2.490921*x1-
     -   1.265445d3*x2+5.276941d-5*x3-1.255524d3*x4-1.510486d-1*x5-
     -   7.066693d-4*x6-2.079471d-1*ninvh-5.203836d-3*ninvh*dfq
      RETURN
      END
