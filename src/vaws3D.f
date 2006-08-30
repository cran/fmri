CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute Location Kernel (Compact support only, based on x^2
C                                   ignores scaling)
C
C          Kern=1     Uniform
C          Kern=2     Epanechnicov
C          Kern=3     Biweight
C          Kern=4     Triweight
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lkern(kern,xsq)
      implicit logical (a-z)
      integer kern
      real*8 xsq,z
      IF (xsq.ge.1) THEN
         lkern=0.d0
      ELSE IF (kern.eq.1) THEN
         lkern=1.d0
      ELSE IF (kern.eq.2) THEN
         lkern=1.d0-xsq
      ELSE IF (kern.eq.3) THEN
         z=1.d0-xsq
         lkern=z*z
      ELSE IF (kern.eq.4) THEN
         z=1.d0-xsq
         lkern=z*z*z
      ELSE IF (kern.eq.5) THEN
         lkern=dexp(-xsq*8.d0)
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
C  skern == "Triangle"
         wj=wj*(1.d0-sij)
      ELSE
C  skern == "Exp"
         IF (sij.gt.spmin) wj=wj*dexp(-spf*(sij-spmin))
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chaws2(y,si2,n1,n2,n3,dv,dv0,hakt,lambda,theta,bi,
     1    ai,kern,skern,spmin,spmax,lwght,wght,vwghts,swjy,thi,narm)
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
      logical aws,narm
      real*8 y(n1,n2,n3,dv),theta(n1,n2,n3,dv0),bi(n1,n2,n3),
     1       ai(n1,n2,n3,dv),lambda,spmax,wght(2),si2(n1,n2,n3),
     1       hakt,lwght(1),spmin,vwghts(dv0),thi(dv0)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n,ij(6)
      real*8 bii,swj,swjy(dv),wj,hakt2,spf,si2j,si2i
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
      dlw1=min0(2*n1-1,2*ih1+1)
      dlw2=min0(2*n2-1,2*ih2+1)
      dlw3=min0(2*n3-1,2*ih3+1)
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
	        si2i=si2(i1,i2,i3)
                if(narm.and.si2i.lt.1d-18) CYCLE
C    si2j.lt.1d-18   indicates that we have an NA in (j1,j2,j3)
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
		  jwind3=(jw3-1)*dlw1*dlw2
                  DO jw2=1,dlw2
	             j2=jw2-clw2+i2
	             if(j2.lt.1.or.j2.gt.n2) CYCLE
		     jwind2=jwind3+(jw2-1)*dlw1
                     DO jw1=1,dlw1
C  first stochastic term
                        wj=lwght(jw1+jwind2)
			if(wj.le.0.d0) CYCLE
	                j1=jw1-clw1+i1
	                if(j1.lt.1.or.j1.gt.n1) CYCLE
			si2j=si2(j1,j2,j3)
                        if(narm.and.si2j.lt.1d-18) CYCLE
C    si2j.lt.1d-18   indicates that we have an NA in (j1,j2,j3)
                        IF (aws) THEN
                           call awswghts(n1,n2,n3,j1,j2,j3,dv0,thi,
     1                     theta,vwghts,skern,spf,spmin,spmax,bii,wj)
                        END IF
                        swj=swj+wj*si2j
			DO k=1,dv
                           swjy(k)=swjy(k)+wj*y(j1,j2,j3,k)*si2j
			END DO
                     END DO
                  END DO
               END DO
	       DO k=1,dv
                  ai(i1,i2,i3,k)=swjy(k)
	       END DO
               bi(i1,i2,i3)=swj
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
      subroutine chawsvr(si2,n1,n2,n3,dv,dv0,hakt,lambda,theta,bi,
     1   var,vred,kern,skern,spmin,spmax,lwght,gwght,swght,dgw,wght,
     2   vwghts,thi,narm)
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
      integer n1,n2,n3,kern,skern,dv,dv0,dgw(3)
      logical aws,narm
      real*8 theta(n1,n2,n3,dv0),bi(n1,n2,n3),lambda,var(n1,n2,n3),
     1     swght(n1,n2,n3),hakt,lwght(1),spmin,vwghts(dv0),thi(dv0),
     2     si2(n1,n2,n3),spmax,wght(2),gwght(1),vred(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n,
     2        dgw1,dgw2,dgw3,j1a,j1e,j2a,j2e,j3a,j3e,
     3        clw10,clw20,clw30
      real*8 bii,swj,swj2,wj,hakt2,swj2vr,spf,si2i,si2j
      integer ngw
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      dgw1=dgw(1)
      dgw2=dgw(2)
      dgw3=dgw(3)
      ngw=dgw1*dgw2*dgw3
      clw30=hakt/wght(2)
      clw20=hakt/wght(1)
      clw10=hakt
      n=n1*n2*n3
      dlw1=min0(2*n1-1,2*clw10+1)
      dlw2=min0(2*n2-1,2*clw20+1)
      dlw3=min0(2*n3-1,2*clw30+1)
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      clw10=clw1-1
      clw20=clw2-1
      clw30=clw3-1
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	        si2i=si2(i1,i2,i3)
                if(narm.and.si2i.lt.1d-18) CYCLE
C    si2j.lt.1d-18   indicates that we have an NA in (j1,j2,j3)
               bii=bi(i1,i2,i3)/lambda
C   scaling of sij outside the loop
               swj=0.d0
	       swj2=0.d0
	       DO k=1,dv0
	          thi(k)=theta(i1,i2,i3,k)
	       END DO
C   fill swght with zeros where needed
               j3a=max0(i3-clw30,1)
               j3e=min0(clw30+i3,n3)
               j2a=max0(i2-clw20,1)
               j2e=min0(clw20+i2,n2)
               j1a=max0(i1-clw10,1)
               j1e=min0(clw10+i1,n1)
               DO j1=j1a,j1e
                  DO j2=j2a,j2e
                     DO j3=j3a,j3e
		        swght(j1,j2,j3)=0.d0
                     END DO
                  END DO
	       END DO
               DO jw3=1,dlw3
	          j3=jw3-clw3+i3
	          if(j3.lt.1.or.j3.gt.n3) CYCLE
		  jwind3=(jw3-1)*dlw1*dlw2
                  DO jw2=1,dlw2
	             j2=jw2-clw2+i2
	             if(j2.lt.1.or.j2.gt.n2) CYCLE
		     jwind2=jwind3+(jw2-1)*dlw1
                     DO jw1=1,dlw1
                        wj=lwght(jw1+jwind2)
			if(wj.le.0.d0) CYCLE
	                j1=jw1-clw1+i1
	                if(j1.lt.1.or.j1.gt.n1) CYCLE
	                si2j=si2(j1,j2,j3)
                        if(narm.and.si2j.lt.1d-18) CYCLE
C    si2j.lt.1d-18   indicates that we have an NA in (j1,j2,j3)
                        IF (aws) THEN
                           call awswghts(n1,n2,n3,j1,j2,j3,dv0,thi,
     1                     theta,vwghts,skern,spf,spmin,spmax,bii,wj)
                        END IF
                        swght(j1,j2,j3)=wj*si2j
                     END DO
                  END DO
               END DO
C
C     now the convolution
C
               call conv3D1(n1,n2,n3,dgw1,dgw2,dgw3,
     1            j1a,j1e,j2a,j2e,j3a,j3e,gwght,swght,si2,swj2,swj2vr)
	       var(i1,i2,i3)=swj2
	       vred(i1,i2,i3)=swj2vr
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Calculate    \tilde{Q}_{h,aws} = \sum_j \tilde{w}_ij^2 *si2_j   in bi2
C                Q_{h,aws} = \sum_j \tilde{w}_ij^2                  in vi2
C                Q_{h,loc} = \sum_j K_h(i,j)^2                      in vi20
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chawsvr1(si2,n1,n2,n3,dv0,hakt,lambda,theta,bi2,bi0,
     1    vi2,vi20,kern,skern,spmin,spmax,lwght,wght,vwghts,thi,narm)
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
      logical aws,narm
      real*8 theta(n1,n2,n3,dv0),bi2(n1,n2,n3),vi2(n1,n2,n3),
     1       vi20(n1,n2,n3),lambda,spmax,wght(2),si2(n1,n2,n3),
     2       hakt,lwght(1),spmin,vwghts(dv0),thi(dv0),bi0(n1,n2,n3)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n
      real*8 bii,swj,swj0,swj00,slwj0,wj,hakt2,spf,si2i,si2j
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      aws=lambda.lt.1d40
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      n=n1*n2*n3
      dlw1=min0(2*n1-1,2*ih1+1)
      dlw2=min0(2*n2-1,2*ih2+1)
      dlw3=min0(2*n3-1,2*ih3+1)
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
	        si2i=si2(i1,i2,i3)
                if(narm.and.si2i.lt.1d-18) CYCLE
               bii=bi2(i1,i2,i3)/lambda
C   scaling of sij outside the loop
               swj=0.d0
               swj0=0.d0
               swj00=0.d0
	       slwj0=0.d0
	       DO k=1,dv0
	          thi(k)=theta(i1,i2,i3,k)
	       END DO
               DO jw3=1,dlw3
	          j3=jw3-clw3+i3
	          if(j3.lt.1.or.j3.gt.n3) CYCLE
		  jwind3=(jw3-1)*dlw1*dlw2
                  DO jw2=1,dlw2
	             j2=jw2-clw2+i2
	             if(j2.lt.1.or.j2.gt.n2) CYCLE
		     jwind2=jwind3+(jw2-1)*dlw1
                     DO jw1=1,dlw1
C  first stochastic term
                        wj=lwght(jw1+jwind2)
			if(wj.le.0.d0) CYCLE
	                j1=jw1-clw1+i1
	                if(j1.lt.1.or.j1.gt.n1) CYCLE
	                si2j=si2(j1,j2,j3)
                        if(narm.and.si2j.lt.1d-18) CYCLE
C    si2j.lt.1d-18   indicates that we have an NA in (j1,j2,j3)
 			swj00=swj00+wj*wj
			wj=wj*si2(j1,j2,j3)
			slwj0=slwj0+wj
                        IF (aws) THEN
                           call awswghts(n1,n2,n3,j1,j2,j3,dv0,thi,
     1                     theta,vwghts,skern,spf,spmin,spmax,bii,wj)
                        END IF
			wj=wj*wj
                        swj=swj+wj/si2j  
			swj0=swj0+wj
                     END DO
                  END DO
               END DO
               bi2(i1,i2,i3)=swj
               bi0(i1,i2,i3)=slwj0
               vi2(i1,i2,i3)=swj0
               vi20(i1,i2,i3)=swj00
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Calculate  Q_{hxg}  = \sum_l [\sum_j K_h(i,j) K_g(j,l)]^2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chawsvr2(n1,n2,n3,hakt,vred,kern,lwght,gwght,
     1                    swght,dgw,wght)
C   
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      integer n1,n2,n3,kern,dgw(3)
      real*8 wght(2),gwght(1),swght(n1,n2,n3),
     1       hakt,lwght(1),vred(n1,n2,n3)
      integer i1,i2,i3,dlw1,dlw2,dlw3,n,
     1        dgw1,dgw2,dgw3,j1a,j1e,j2a,j2e,j3a,j3e,
     2        clw10,clw20,clw30,cn1,cn2,cn3,i1c,i2c,i3c,cr1,cr2,cr3
      real*8 swj2,hakt2,swj2vr,varc,vredc
      integer ngw
      hakt2=hakt*hakt
      dgw1=dgw(1)
      dgw2=dgw(2)
      dgw3=dgw(3)
      ngw=dgw1*dgw2*dgw3
      clw30=hakt/wght(2)
      clw20=hakt/wght(1)
      clw10=hakt
      n=n1*n2*n3
      dlw1=min0(2*n1-1,2*clw10+1)
      dlw2=min0(2*n2-1,2*clw20+1)
      dlw3=min0(2*n3-1,2*clw30+1)
      clw10=(dlw1-1)/2
      clw20=(dlw2-1)/2
      clw30=(dlw3-1)/2
      cr1=clw10+(dgw1-1)/2
      cr2=clw20+(dgw2-1)/2
      cr3=clw30+(dgw3-1)/2
      cn1=(n1+1)/2
      cn2=(n2+1)/2
      cn3=(n3+1)/2
C   
C    location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
C
C
C     get correct values for the center
C
C
      i3=cn3
      i2=cn3
      i1=cn1
      call fillwgh0(n1,n2,n3,i1,i2,i3,dlw1,dlw2,dlw3,lwght,
     1              swght,wght,hakt2)      
C
C     now the convolution
C
      j3a=max0(i3-clw30,1)
      j3e=min0(clw30+i3,n3)
      j2a=max0(i2-clw20,1)
      j2e=min0(clw20+i2,n2)
      j1a=max0(i1-clw10,1)
      j1e=min0(clw10+i1,n1)
      call conv3D0(n1,n2,n3,dgw1,dgw2,dgw3,
     1       j1a,j1e,j2a,j2e,j3a,j3e,gwght,swght,swj2)
      vredc=swj2
      call rchkusr()
C
C
C     calculate everything up to symmetry
C
C      
      DO i3=1,cn3
         DO i2=1,cn2
            DO i1=1,cn1
               if(i1.gt.cr1.and.i2.gt.cr2.and.i3.gt.cr3) THEN
C   we are in the center and can use vred
	          vred(i1,i2,i3)=vredc
	       ELSE
                  call fillwgh0(n1,n2,n3,i1,i2,i3,dlw1,dlw2,dlw3,
     1                          lwght,swght,wght,hakt2)      
C
C     now the convolution
C
                  j3a=max0(i3-clw30,1)
                  j3e=min0(clw30+i3,n3)
                  j2a=max0(i2-clw20,1)
                  j2e=min0(clw20+i2,n2)
                  j1a=max0(i1-clw10,1)
                  j1e=min0(clw10+i1,n1)
                  call conv3D0(n1,n2,n3,dgw1,dgw2,dgw3,j1a,j1e,
     1                   j2a,j2e,j3a,j3e,gwght,swght,swj2)
	          vred(i1,i2,i3)=swj2
                  call rchkusr()
	       END IF
            END DO
         END DO
      END DO
C
C
C      Now use symmetry
C
C
      DO i1=1,cn1
         i1c=n1+1-i1
         DO i2=1,cn2
            i2c=n2+1-i2
	    DO i3=1,cn3
               i3c=n3+1-i3	       
	       swj2=vred(i1,i2,i3)
	       vred(i1c,i2,i3)=swj2
	       vred(i1c,i2c,i3)=swj2
	       vred(i1c,i2,i3c)=swj2
	       vred(i1c,i2c,i3c)=swj2
	       vred(i1,i2c,i3)=swj2
	       vred(i1,i2c,i3c)=swj2
	       vred(i1,i2,i3c)=swj2
	    END DO
	 END DO
      END DO	 
      RETURN
      END
C
C            Arrange cube of weights
C
      subroutine fillwght(n1,n2,n3,i1,i2,i3,dlw1,dlw2,dlw3,si2,lwght,
     1                    swght,wght,hakt2)
      implicit logical (a-z)
      integer n1,n2,n3,i1,i2,i3,dlw1,dlw2,dlw3
      real*8 si2(n1,n2,n3),swght(n1,n2,n3),lwght(dlw1,dlw2,dlw3),
     1       wght(2),hakt2
      integer clw1,clw2,clw3,j1a,j1e,j2a,j2e,j3a,j3e,clw10,clw20,
     1        clw30,j1,j2,j3,jw1,jw2,jw3
      real*8 wj
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
      clw10=clw1-1
      clw20=clw2-1
      clw30=clw3-1
      j3a=max0(i3-clw30,1)
      j3e=min0(clw30+i3,n3)
      j2a=max0(i2-clw20,1)
      j2e=min0(clw20+i2,n2)
      j1a=max0(i1-clw10,1)
      j1e=min0(clw10+i1,n1)
      DO j1=j1a,j1e
         DO j2=j2a,j2e
            DO j3=j3a,j3e
	       swght(j1,j2,j3)=0.d0
            END DO
         END DO
      END DO
      DO jw3=1,dlw3
	 j3=jw3-clw3+i3
	 if(j3.lt.1.or.j3.gt.n3) CYCLE
         DO jw2=1,dlw2
	    j2=jw2-clw2+i2
	    if(j2.lt.1.or.j2.gt.n2) CYCLE
            DO jw1=1,dlw1
	       wj=lwght(jw1,jw2,jw3)
	       if(wj.le.0.d0) CYCLE
C  first stochastic term
	       j1=jw1-clw1+i1
	       if(j1.lt.1.or.j1.gt.n1) CYCLE
               swght(j1,j2,j3)=wj*si2(j1,j2,j3)
            END DO
         END DO
      END DO
C
C     swght contains location weights for point i1,i2,i3
C
      RETURN
      END
C
C            Arrange cube of weights (homogeneous case)
C
      subroutine fillwgh0(n1,n2,n3,i1,i2,i3,dlw1,dlw2,dlw3,lwght,
     1                    swght,wght,hakt2)
      implicit logical (a-z)
      integer n1,n2,n3,i1,i2,i3,dlw1,dlw2,dlw3
      real*8 swght(n1,n2,n3),lwght(dlw1,dlw2,dlw3),wght(2),hakt2
      integer clw1,clw2,clw3,j1a,j1e,j2a,j2e,j3a,j3e,clw10,clw20,
     1        clw30,j1,j2,j3,jw1,jw2,jw3
      real*8 wj
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
      clw10=clw1-1
      clw20=clw2-1
      clw30=clw3-1
      j3a=max0(i3-clw30,1)
      j3e=min0(clw30+i3,n3)
      j2a=max0(i2-clw20,1)
      j2e=min0(clw20+i2,n2)
      j1a=max0(i1-clw10,1)
      j1e=min0(clw10+i1,n1)
      DO j1=j1a,j1e
         DO j2=j2a,j2e
            DO j3=j3a,j3e
	       swght(j1,j2,j3)=0.d0
            END DO
         END DO
      END DO
      DO jw3=1,dlw3
	 j3=jw3-clw3+i3
	 if(j3.lt.1.or.j3.gt.n3) CYCLE
         DO jw2=1,dlw2
	    j2=jw2-clw2+i2
	    if(j2.lt.1.or.j2.gt.n2) CYCLE
            DO jw1=1,dlw1
	       wj=lwght(jw1,jw2,jw3)
	       if(wj.le.0.d0) CYCLE
C  first stochastic term
	       j1=jw1-clw1+i1
	       if(j1.lt.1.or.j1.gt.n1) CYCLE
               swght(j1,j2,j3)=wj
            END DO
         END DO
      END DO
C
C     swght contains location weights for point i1,i2,i3
C
      RETURN
      END
      subroutine conv3D1(n1,n2,n3,dgw1,dgw2,dgw3,
     1             j1a,j1e,j2a,j2e,j3a,j3e,gwght,swght,si2,sw2,sw2v)
      implicit logical (a-z)
      integer n1,n2,n3,dgw1,dgw2,dgw3,j1a,j1e,j2a,j2e,j3a,j3e
      real*8 gwght(dgw1,dgw2,dgw3),swght(n1,n2,n3),sw2,sw2v,
     1       si2(n1,n2,n3)
      integer cgw10,cgw20,cgw30,l1,l2,l3,m1,m2,m3,k1,k2,k3,
     1       cgw1,cgw2,cgw3,m10,m20,m30
      real*8 sw
      cgw10=(dgw1-1)/2
      cgw20=(dgw2-1)/2
      cgw30=(dgw3-1)/2
      cgw1=cgw10+1
      cgw2=cgw20+1
      cgw3=cgw30+1
      sw2=0.d0
      sw2v=0.d0
      DO l3=j3a-cgw30,j3e+cgw30
	 if(l3.lt.1.or.l3.gt.n3) CYCLE
	 DO l2=j2a-cgw20,j2e+cgw20
	    if(l2.lt.1.or.l2.gt.n2) CYCLE
	    DO l1=j1a-cgw10,j1e+cgw10
	       if(l1.lt.1.or.l1.gt.n1) CYCLE
	       sw=0.d0
	       DO m10=-cgw10,cgw10
		  k1=m10+l1
	          if(k1.lt.j1a.or.k1.gt.j1e) CYCLE
                  m1=m10+cgw1
	          DO m20=-cgw20,cgw20
		     k2=m20+l2
	             if(k2.lt.j2a.or.k2.gt.j2e) CYCLE
                     m2=m20+cgw2
		     DO m30=-cgw30,cgw30
			k3=m30+l3
			if(k3.lt.j3a.or.k3.gt.j3e) CYCLE
                        m3=m30+cgw3
			sw=sw+swght(k1,k2,k3)*gwght(m1,m2,m3)
		     END DO
		  END DO
	       END DO
	       sw2=sw2+sw*sw/si2(l1,l2,l3)
	       sw2v=sw2v+sw*sw
	    END DO
         END DO
      END DO
C
C    sw2 contains  sum_l (sum_j w_ij v_jl)^2 sigma2_l
C    sw2v contains  sum_l (sum_j w_ij v_jl)^2 
C    with  w_ij in swght and vjl in gwght 
C
      RETURN
      END
      subroutine conv3D0(n1,n2,n3,dgw1,dgw2,dgw3,
     1             j1a,j1e,j2a,j2e,j3a,j3e,gwght,swght,sw2)
      implicit logical (a-z)
      integer n1,n2,n3,dgw1,dgw2,dgw3,j1a,j1e,j2a,j2e,j3a,j3e
      real*8 gwght(dgw1,dgw2,dgw3),swght(n1,n2,n3),sw2
      integer cgw10,cgw20,cgw30,l1,l2,l3,m1,m2,m3,k1,k2,k3,
     1       cgw1,cgw2,cgw3,m10,m20,m30
      real*8 sw
      cgw10=(dgw1-1)/2
      cgw20=(dgw2-1)/2
      cgw30=(dgw3-1)/2
      cgw1=cgw10+1
      cgw2=cgw20+1
      cgw3=cgw30+1
      sw2=0.d0
      DO l3=j3a-cgw30,j3e+cgw30
	 if(l3.lt.1.or.l3.gt.n3) CYCLE
	 DO l2=j2a-cgw20,j2e+cgw20
	    if(l2.lt.1.or.l2.gt.n2) CYCLE
	    DO l1=j1a-cgw10,j1e+cgw10
	       if(l1.lt.1.or.l1.gt.n1) CYCLE
	       sw=0.d0
	       DO m10=-cgw10,cgw10
		  k1=m10+l1
	          if(k1.lt.j1a.or.k1.gt.j1e) CYCLE
                  m1=m10+cgw1
	          DO m20=-cgw20,cgw20
		     k2=m20+l2
	             if(k2.lt.j2a.or.k2.gt.j2e) CYCLE
                     m2=m20+cgw2
		     DO m30=-cgw30,cgw30
			k3=m30+l3
			if(k3.lt.j3a.or.k3.gt.j3e) CYCLE
                        m3=m30+cgw3
			sw=sw+swght(k1,k2,k3)*gwght(m1,m2,m3)
		     END DO
		  END DO
	       END DO
	       sw2=sw2+sw*sw
	    END DO
         END DO
      END DO
C
C    sw2 contains  sum_l (sum_j w_ij v_jl)^2 
C    with  w_ij in swght and vjl in gwght 
C
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
         ih2=dsqrt(hakt2-z3)/wght(1)
         DO j2=clw2-ih2,clw2+ih2
            z2=(clw2-j2)*wght(1)
            z2=z3+z2*z2
            ih1=dsqrt(hakt2-z2)
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
      subroutine nqg(gwght,gwght2,dg1,dg2,dg3,n1,n2,n3,qg,ng)
      implicit logical (a-z)
      integer dg1,dg2,dg3,n1,n2,n3
      real*8 gwght(dg1,dg2,dg3),gwght2(dg1,dg2,dg3),qg(n1,n2,n3),
     1       ng(n1,n2,n3)
      integer cg1,cg2,cg3,cn1,cn2,cn3,i1,i2,i3,kk,ll,mm
      real*8 zn,zq,zrown,zrowq,z
      cg1=dg1/2
      cg2=dg2/2
      cg3=dg3/2
      zn=0.d0
      zq=0.d0
      DO i1=1,dg1
         DO i2=1,dg2
	    DO i3=1,dg3
	       zn=zn+gwght(i1,i2,i3)
	       zq=zq+gwght2(i1,i2,i3)
	    END DO
	 END DO
      END DO
      cn1=(n1+1)/2
      cn2=(n2+1)/2
      cn3=(n3+1)/2
      DO i1=1,cn1
         DO i2=1,cn2
	    DO i3=1,cn3
	    ng(i1,i2,i3)=zn
	    qg(i1,i2,i3)=zq
	    END DO
	 END DO
      END DO
C   now handle boundary surfaces
C      call intpr("boundary",8,cg1,1)
      call rchkusr()
      ll=cg1
      DO WHILE(ll.gt.0) 
         zrown=0.d0
	 zrowq=0.d0
	 DO i2=1,dg2
	    DO i3=1,dg3
	       zrown=zrown+gwght(cg1+1-ll,i2,i3)
	       zrowq=zrowq+gwght2(cg1+1-ll,i2,i3)
	    END DO
	 END DO
         DO i1=1,ll
            DO i2=1,cn2
               DO i3=1,cn3
	          ng(i1,i2,i3)=ng(i1,i2,i3)-zrown 
	          qg(i1,i2,i3)=qg(i1,i2,i3)-zrowq 
	       END DO
	    END DO
	 END DO
	 ll=ll-1
      END DO
C      call intpr("boundary",8,cg2,1)
      call rchkusr()
      ll=cg2
      DO WHILE(ll.gt.0) 
         zrown=0.d0
	 zrowq=0.d0
	 DO i1=1,dg1
	    DO i3=1,dg3
	       zrown=zrown+gwght(i1,cg2+1-ll,i3)
	       zrowq=zrowq+gwght2(i1,cg2+1-ll,i3)
	    END DO
	 END DO
         DO i2=1,ll
            DO i1=1,cn1
               DO i3=1,cn3
	          ng(i1,i2,i3)=ng(i1,i2,i3)-zrown 
	          qg(i1,i2,i3)=qg(i1,i2,i3)-zrowq 
	       END DO
	    END DO
	 END DO
	 ll=ll-1
      END DO
C      call intpr("boundary",8,cg3,1)
      call rchkusr()
      ll=cg3
      DO WHILE(ll.gt.0) 
         zrown=0.d0
	 zrowq=0.d0
	 DO i1=1,dg1
	    DO i2=1,dg2
	       zrown=zrown+gwght(i1,i2,cg3+1-ll)
	       zrowq=zrowq+gwght2(i1,i2,cg3+1-ll)
	    END DO
	 END DO
         DO i3=1,ll
            DO i1=1,cn1
               DO i2=1,cn2
	          ng(i1,i2,i3)=ng(i1,i2,i3)-zrown 
	          qg(i1,i2,i3)=qg(i1,i2,i3)-zrowq 
	       END DO
	    END DO
	 END DO
	 ll=ll-1
      END DO
C  now edges  
C      call intpr("edges",5,cg3,1)
      ll=cg1
      mm=cg2
      call rchkusr()
      DO WHILE(ll.gt.0)
         DO WHILE(mm.gt.0)
            zrown=0.d0
            zrowq=0.d0
	    DO i3=1,dg3
	       zrown=zrown+gwght(cg1+1-ll,cg2+1-mm,i3)
	       zrowq=zrowq+gwght2(cg1+1-ll,cg2+1-mm,i3)
	    END DO
            DO i1=1,ll
               DO i2=1,mm
                  DO i3=1,cn3
	             ng(i1,i2,i3)=ng(i1,i2,i3)+zrown 
	             qg(i1,i2,i3)=qg(i1,i2,i3)+zrowq 
	          END DO
	       END DO
	    END DO
	    mm=mm-1
	 END DO
	 ll=ll-1
      END DO
C      call intpr("edges",5,cg2,1)
      ll=cg1
      mm=cg3
      call rchkusr()
      DO WHILE(ll.gt.0)
         DO WHILE(mm.gt.0)
            zrown=0.d0
            zrowq=0.d0
	    DO i2=1,dg2
	       zrown=zrown+gwght(cg1+1-ll,i2,cg3+1-mm)
	       zrowq=zrowq+gwght2(cg1+1-ll,i2,cg3+1-mm)
	    END DO
            DO i1=1,ll
               DO i3=1,mm
                  DO i2=1,cn2
	             ng(i1,i2,i3)=ng(i1,i2,i3)+zrown 
	             qg(i1,i2,i3)=qg(i1,i2,i3)+zrowq 
	          END DO
	       END DO
	    END DO
	    mm=mm-1
	 END DO
	 ll=ll-1
      END DO
C      call intpr("edges",5,cg3,1)
      call rchkusr()
      ll=cg2
      mm=cg3
      DO WHILE(ll.gt.0)
         DO WHILE(mm.gt.0)
            zrown=0.d0
            zrowq=0.d0
	    DO i1=1,dg1
	       zrown=zrown+gwght(i1,cg2+1-ll,cg3+1-mm)
	       zrowq=zrowq+gwght2(i1,cg2+1-ll,cg3+1-mm)
	    END DO
            DO i2=1,ll
               DO i3=1,mm
                  DO i1=1,cn1
	             ng(i1,i2,i3)=ng(i1,i2,i3)+zrown 
	             qg(i1,i2,i3)=qg(i1,i2,i3)+zrowq 
	          END DO
	       END DO
	    END DO
	    mm=mm-1
	 END DO
	 ll=ll-1
      END DO
      call rchkusr()
C  now the corner
C      call intpr("corner",6,cg3,1)
      kk=cg1
      ll=cg2
      mm=cg3
      DO WHILE(kk.gt.0)
         DO WHILE(ll.gt.0)
            DO WHILE(mm.gt.0)
               DO i2=1,ll
                  DO i3=1,mm
                     DO i1=1,kk
	                ng(i1,i2,i3)=ng(i1,i2,i3)-
     1                               gwght(cg1+1-kk,cg2+1-ll,cg3+1-mm)
	                qg(i1,i2,i3)=qg(i1,i2,i3)-
     1                              gwght2(cg1+1-kk,cg2+1-ll,cg3+1-mm)
	             END DO
	          END DO
	       END DO
	       mm=mm-1
	    END DO
	    ll=ll-1
	 END DO
	 kk=kk-1
      END DO
      call rchkusr()      
C  now symmetries
      DO i1=1,cn1
         DO i2=1,cn2
	    DO i3=1,cn3
	       z=ng(i1,i2,i3)
	       ng(n1+1-i1,i2,i3)=z
	       ng(n1+1-i1,n2+1-i2,i3)=z
	       ng(n1+1-i1,n2+1-i2,n3+1-i3)=z
	       ng(n1+1-i1,i2,n3+1-i3)=z
	       ng(i1,n2+1-i2,i3)=z
	       ng(i1,n2+1-i2,n3+1-i3)=z
	       ng(i1,i2,n3+1-i3)=z
	       z=qg(i1,i2,i3)
	       qg(n1+1-i1,i2,i3)=z
	       qg(n1+1-i1,n2+1-i2,i3)=z
	       qg(n1+1-i1,n2+1-i2,n3+1-i3)=z
	       qg(n1+1-i1,i2,n3+1-i3)=z
	       qg(i1,n2+1-i2,i3)=z
	       qg(i1,n2+1-i2,n3+1-i3)=z
	       qg(i1,i2,n3+1-i3)=z
	    END DO
	 END DO
      END DO
      RETURN
      END
