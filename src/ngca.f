      subroutine fd1(z,s2,fv,dv)
      real*8 z,s2,fv,dv
      real*8 z2,zexp,z2s2
      z2=z*z
      z2s2=z2/s2
      zexp=exp(-z2s2/2.d0)
      fv=z*z2*zexp
      dv=z2*(3.d0-z2s2)*zexp
      RETURN 
      END
      subroutine fd2(z,s,fv,dv)
      real*8 z,s,fv,dv
      real*8 ztanh
      ztanh=tanh(s*z)
      fv=ztanh
      dv=s*(1.d0-ztanh*ztanh)
      RETURN 
      END
      subroutine fd3(z,s,fv,dv)
      real*8 z,s,fv,dv
      real*8 sz
      sz=s*z
      fv=sin(sz)
      dv=s*cos(sz)
      RETURN 
      END
      subroutine fd4(z,s,fv,dv)
      real*8 z,s,fv,dv
      real*8 sz
      sz=s*z
      fv=cos(sz)
      dv=-s*sin(sz)
      RETURN 
      END
      subroutine fastica(y,omega,d,n,l,ifun,t,beta,v,normv,s,omegak,
     1                   delta)
      implicit logical (a-z)
      integer d,n,l,t,ifun(l)
      real*8 y(d,n),omega(d,l),v(d,l),beta(d),s(l),normv(l),omegak(d),
     1       delta
      real*8 dnrm2,ddot
      integer i,j,k,o,ifunk,icount
      real*8 nk,z,fw,fwd1,nbeta,zv,zdelta,sk,ninv,eps
      external dnrm2,ddot,dscal
      call dblepr("delta",5,delta,1)
      eps=1.d-6
      nk=0.d0
      ninv=1.d0/n
C this is just to avoid some warnings 
      DO k=1,l
         sk=s(k)
         ifunk=ifun(k)
         call dcopy(d,omega(1,k),1,omegak,1) 
C         call dblepr("omega0",6,omegak,d)
         DO j=1,t
C Start Loop t
            nk=0.d0
            DO o=1,d
               beta(o)=0.d0
            END DO
            DO i=1,n
               z=ddot(d,omegak,1,y(1,i),1)
C    <omega_k,y_i>
C            call dblepr("z",1,z,1)
               select case (ifunk)
               case(1) 
                  call fd1(z,sk,fw,fwd1)
               case(2) 
                  call fd2(z,sk,fw,fwd1)
               case(3) 
                  call fd3(z,sk,fw,fwd1)
               case(4) 
                  call fd4(z,sk,fw,fwd1)
               case default
                  fw=0.d0
                  fwd1=1.d0
               end select
               DO o=1,d
                  z=y(o,i)*fw-omegak(o)*fwd1
                  beta(o)=beta(o)+z
                  nk=nk+z*z
               END DO
            END DO
            call dscal(d,ninv,beta,1)
C   beta <- beta/n
            nbeta=dnrm2(d,beta,1)
C   nbeta <- ||beta||_2
C            call dblepr("nbeta",5,nbeta,1)
            IF(nbeta.gt.eps) THEN
               zdelta = 0.d0
               DO o=1,d
                  z=beta(o)/nbeta
                  zdelta =zdelta+omegak(o)*z
                  omegak(o)=z
               END DO
C            call dblepr("omegak",6,omegak,d)
            ELSE
               EXIT
C   keep omega
            END IF
           zdelta=acos(abs(zdelta))
C                call dblepr("zdelta",6,zdelta,1)
           IF(zdelta.lt.delta) THEN
C               call intpr("stopped at time",15,j,1)
C               call dblepr("zdelta",6,zdelta,1)
               EXIT
            END IF
C   omega does not change
C End Loop t
         END DO
         nk=nk/n-nbeta*nbeta
         nk=sqrt(n/nk)
         zv=0.d0
         DO o=1,d
            z=beta(o)*nk
            v(o,k)=z
            zv=zv+z*z
         END DO
         IF(zv.gt.1d-8) THEN
            normv(k)=sqrt(zv)
         ELSE
            normv(k)=0.d0
         END IF
         call rchkusr()
         icount=k/100
         if(icount*100.eq.k) call intpr("testfunction nr.",15,k,1)
      END DO
      RETURN
      END
      subroutine smtime(x,n1,n2,n3,nt,mask,h,xnew,w,lw)
      implicit logical (a-z)
      integer n1,n2,n3,nt,lw
      logical mask(n1,n2,n3)
      real*8 x(n1,n2,n3,nt),xnew(n1,n2,n3,nt),w(lw),h
      integer i1,i2,i3,j,jn,ja,je,clw,clw1
      real*8 z,sw
      clw=lw/2
      clw1=clw+1
C  Use Epanechnikov kernel
      w(clw1)=1.d0
      DO jn=1,clw
         z=jn/h
         z=1.d0-z*z
         w(clw1-jn)=z
         w(clw1+jn)=z
      END DO
      DO j=1,nt
         ja = max(1,j-clw)
         je = min(nt,j+clw)
         sw=0.d0
         DO jn=ja,je
            sw=sw+w(jn-j+clw1)
         END DO
         DO i1=1,n1
            DO i2=1,n2
               DO i3=1,n3
                  IF(.not.mask(i1,i2,i3)) CYCLE
                  z=0.d0
                  DO jn=ja,je
                     z=z+w(jn-j+clw1)*x(i1,i2,i3,jn)
                  END DO
                  xnew(i1,i2,i3,j)=z/sw
               END DO
            END DO
         END DO
         call rchkusr()
      END DO
      RETURN
      END
      subroutine smspace(x,n1,n2,n3,nt,mask,h,xnew,vext,w,lw1,lw2,lw3)
      implicit logical (a-z)
      integer n1,n2,n3,nt,lw1,lw2,lw3
      logical mask(n1,n2,n3)
      real*8 x(n1,n2,n3,nt),xnew(n1,n2,n3,nt),vext(3),w(lw1,lw2,lw3),h
      integer j1,j2,j3,jt,jn1,jn2,jn3,ja1,je1,ja2,je2,ja3,je3,
     1        clw1,clw2,clw3,clw11,clw21,clw31,jj1,jj2
      real*8 z,z1,z2,z3,sw
      clw1=lw1/2
      clw2=lw2/2
      clw3=lw3/2
      clw11=clw1+1
      clw21=clw2+1
      clw31=clw3+1
C  Use Epanechnikov kernel
      w(clw11,clw21,clw31)=1.d0
      DO j1=1,clw1
         z1=j1/h*vext(1)
         z1=z1*z1
         DO j2=1,clw2
            z2=j2/h*vext(2)
            z2=z2*z2+z1
            DO j3=1,clw3
               z3=j3/h*vext(3)
               z=max(0.d0,1.d0-z2-z2*z3)
               w(clw11-j1,clw21-j2,clw31-j3)=z
               w(clw11-j1,clw21-j2,clw31+j3)=z
               w(clw11-j1,clw21+j2,clw31-j3)=z
               w(clw11-j1,clw21+j2,clw31+j3)=z
               w(clw11+j1,clw21-j2,clw31-j3)=z
               w(clw11+j1,clw21-j2,clw31+j3)=z
               w(clw11+j1,clw21+j2,clw31-j3)=z
               w(clw11+j1,clw21+j2,clw31+j3)=z
            END DO
         END DO
      END DO
      DO j1=1,n1
         DO j2=1,n2
            DO j3=1,n3
               if(.not.mask(j1,j2,j3)) CYCLE
               sw=0.d0
               ja1=max(1,j1-clw1)
               je1=min(n1,j1+clw1)
               ja2=max(1,j2-clw2)
               je2=min(n2,j2+clw2)
               ja3=max(1,j3-clw3)
               je3=min(n3,j3+clw3)
               DO jn1=ja1,je1
                  jj1=jn1-j1+clw11
                  DO jn2=ja2,je2
                     jj2=jn2-j2+clw21
                     DO jn3=ja3,je3
                        if(.not.mask(jn1,jn2,jn3)) CYCLE
                        sw=sw+w(jj1,jj2,jn3-j3+clw31)
                     END DO
                  END DO
               END DO
               DO jt=1,nt
                  z=0.d0
                  DO jn1=ja1,je1
                     jj1=jn1-j1+clw11
                     DO jn2=ja2,je2
                        jj2=jn2-j2+clw21
                        DO jn3=ja3,je3
                        if(.not.mask(jn1,jn2,jn3)) CYCLE
                        z=z+x(jn1,jn2,jn3,jt)*w(jj1,jj2,jn3-j3+clw31)
                        END DO
                     END DO
                  END DO
                  xnew(j1,j2,j3,jt)=z/sw
               END DO
            END DO
            call rchkusr()
         END DO
      END DO
      RETURN
      END
