CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate segmentation
C   for fmri data
C   called in segm3D file (segm.r)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine segm3d(y,res,si2,pos,wlse,n1,n2,n3,nt,df,hakt,
     1                  lambda,theta,bi,thn,lwght,wght,swres,pval,
     3                  segm,delta,thresh,fov,vq,vest0i,varest,
     4                  restrict)
C
C   y        observed values of regression function
C   res      residuals (scaled with spm$resscale)
C   si2      variance estimate^{-1} (1/arg to function segm3D)
C   mask     brain mask
C   wlse     (arg weighted to function segm3D)
C   n1,n2,n3    design dimensions
C   nt       length of fmri series
C   df       degrees of freedom for estimates in spm (skalar)
C   hakt     actual bandwidth
C   lambda   lambda (for t(,df))
C   theta    estimates from last step   (input)
C   bi       \sum  Wi (wi depends on si2 in case of wlse)  (output)
C   thn      new estimates (output)
C   lwght    array for non-adaptive weights (auxiliary)
C   wght     rations of voxel dimensions
C   swres    vector for sum_j(w_{ij} res_j) (auxiliary)
C   pval     array of p-values  (auxiliary for penalized smoothing within segments)
C   segm     array of segment assignments (output)
C   delta    half width for hypothesis (input)
C   thresh   threshold for local test in segment assignments
C   fov      sum(mask) (input)
C   vq       array of rescaled variance varesti/si2 (input)
C   vest0i   array of residual variances (input)
C   varest   array of smoothed residual variances (input)
C   restrict penalize smoothing within segments (using pval)
C
      implicit none
      integer n1,n2,n3,nt,kern,segm(*)
      logical aws
      integer wlse,pos(n1,n2,n3),restrict
      double precision y(*),theta(*),bi(*),delta,thn(*),lambda,
     1      wght(2),si2(*),pval(*),hakt,lwght(*),thi,swres(nt),fov,
     2      vq(*),varest(*),res(nt,*),vest0i(*),df,thresh
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,iindp,jindp,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,dlw12,k,segmi
      double precision bii,swj,swjy,wj,hakt2,spf,si2j,si2i,vqi,
     1       varesti,fpchisq,ti,thij,sij,z,si,swr,z1,
     2       a,b,dn,pvali
      external fpchisq
      kern=1
      hakt2=hakt*hakt
      spf=4.d0/3.d0
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/wght(2))
      ih2=FLOOR(hakt/wght(1))
      ih1=FLOOR(hakt)
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      dlw12=dlw1*dlw2
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
               iindp = pos(i1,i2,i3)
               if(iindp.eq.0) CYCLE
               vqi=vq(iindp)
               thi=theta(iindp)
               si2i=vest0i(iindp)
               varesti=varest(iindp)
               dn=varesti/si2i*fov
               call getdfnab(df,dn,a,b)
C   this should be more conservative using actual variance reduction instead of theoretical
               ti=max(0.d0,abs(thi)-delta)
               IF(a*ti/sqrt(varesti/vqi)-b.gt.thresh) THEN
                   pval(iindp)=0.d0
               ELSE
                   pval(iindp)=1.d0
               END IF
            END DO
         END DO
      END DO
      END IF
C   scaling of sij outside the loop
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               iindp=pos(i1,i2,i3)
               if(iindp.eq.0) CYCLE
               segmi=segm(iindp)
               pvali=pval(iindp)
               vqi=vq(iindp)
               si2i=vest0i(iindp)
               bii=bi(iindp)/lambda
               swj=0.d0
               swjy=0.d0
               DO k=1,nt
                  swres(k)=0.d0
               END DO
               thi=theta(iindp)
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
                        jindp=pos(j1,j2,j3)
                        IF(jindp.eq.0) CYCLE
                        wj=lwght(jw1+(jw2-1)*dlw1+(jw3-1)*dlw12)
                        if(wj.le.0.d0) CYCLE
                        si2j=si2(jindp)
                        IF (aws) THEN
                           thij=thi-theta(jindp)
                           sij=thij*thij*bii
                           if(restrict.ne.0) THEN
C restrict smoothing within segmented areas
                           if(abs(segmi).eq.1) THEN
                           if(segmi*segm(jindp).gt.0) THEN
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
                        if(wlse.ne.0)  wj=wj*si2j
                        swj=swj+wj
                        swjy=swjy+wj*y(jindp)
C  weighted sum of residuals
                        call daxpy(nt,wj,res(1,jindp),1,swres,1)
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
               if(restrict.ne.0) THEN
C  smoothing restricted within segmented ares
               if(segmi.eq.1) THEN
                  if(thi.lt.theta(iindp)) THEN
                     thi = theta(iindp)
                  ELSE
                     varest(iindp)=si
                     bi(iindp)=si2i/si*si2(iindp)
                  END IF
               END IF
               if(segmi.eq.-1) THEN
                  if(thi.gt.theta(iindp)) THEN
                     thi = theta(iindp)
                  ELSE
                     varest(iindp)=si
                     bi(iindp)=si2i/si*si2(iindp)
                  END IF
               END IF
               END IF
C  end if for smoothing restricted within segmented ares
               thn(iindp)=thi
C               si = si/nt
               if(restrict*segmi.ne.0) CYCLE
               varest(iindp)=si
               bi(iindp)=si2i/si*si2(iindp)
C   keep the detected segment
               dn=si/si2i*fov
               call getdfnab(df,dn,a,b)
C
C   note that a and b refer to  1/a_n and b_n/a_n
C
C   this should be more conservative using actual variance reduction instead of theoretical
               si=sqrt(si/vqi)
C   thats the SD of thi
               if(a*(thi+delta)/si+b.lt.-thresh) THEN
                  segm(iindp)=-1
               ELSE IF (a*(thi-delta)/si-b.gt.thresh) THEN
                  segm(iindp)=1
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
      implicit none
      double precision a,b,df,n
      double precision dfinv,ldf,dfq,ninvh,x1,x2,x3,x4,x5,x6,lna,lnb,ln
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
C
C   calculate location weights in lwght
C
      subroutine locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      implicit none
      integer dlw1,dlw2,dlw3,kern
      double precision wght(2),hakt2,lwght(dlw1,dlw2,dlw3),lkern
      external lkern
      double precision z1,z2,z3
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
         ih2=FLOOR(sqrt(hakt2-z3)/wght(1))
         DO j2=clw2-ih2,clw2+ih2
            IF(j2.lt.1.or.j2.gt.dlw2) CYCLE
            z2=(clw2-j2)*wght(1)
            z2=z3+z2*z2
            ih1=FLOOR(sqrt(hakt2-z2))
            DO j1=clw1-ih1,clw1+ih1
               IF(j1.lt.1.or.j1.gt.dlw1) CYCLE
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
      double precision function lkern(kern,xsq)
      implicit none
      integer kern
      double precision xsq
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
