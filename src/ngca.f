      subroutine fd1(z,s2,fv,dv)
      real*8 z,s2,fv,dv
      real*8 z2,zexp
      z2=z*z
      zexp=exp(-z2/2.d0/s2)
      fv=z*z2*zexp
      dv=z2*(3.d0-z2/s2)*zexp
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
      subroutine fastica(y,omega,d,n,l,ifun,t,beta,v,normv,s)
      implicit logical (a-z)
      integer d,n,l,t,ifun(l)
      real*8 y(d,n),omega(d,l),v(d,l),beta(d),s(l),normv(l)
      integer i,j,k,o
      real*8 nk,z,fw,fwd1,nbeta,srnbeta,zv,zdelta
      nk=0.d0
      nbeta=0.d0
C this is just to avoid some warnings 
      DO k=1,l
         DO j=1,t
            nk=0.d0
            DO o=1,d
               beta(o)=0.d0
            END DO
            DO i=1,n
               z=0.d0
               DO o=1,d
                  z=z+omega(o,k)*y(o,i)
               END DO
               select case (ifun(k))
               case(1) 
                  call fd1(z,s(k),fw,fwd1)
               case(2) 
                  call fd2(z,s(k),fw,fwd1)
               case(3) 
                  call fd3(z,s(k),fw,fwd1)
               case(4) 
                  call fd4(z,s(k),fw,fwd1)
               case default
                  fw=0.d0
                  fwd1=1.d0
               end select
               DO o=1,d
                  z=y(o,i)*fw-omega(o,k)*fwd1
                  beta(o)=beta(o)+z
                  nk=nk+z*z
               END DO
            END DO
            DO o=1,d
               beta(o)=beta(o)/n
            END DO
            nbeta=0.d0
            DO o=1,d
               z=beta(o)
               nbeta=nbeta+z*z
            END DO
            IF(nbeta.gt.1d-6) THEN
               srnbeta=sqrt(nbeta)
               zdelta = 0.d0
               DO o=1,d
                  z=beta(o)/srnbeta
                  zdelta =zdelta+abs(omega(o,k)-z)
                  omega(o,k)=z
               END DO
            ELSE
               EXIT
C   keep omega
            END IF
            IF(zdelta.lt.1.d-5)  EXIT
         END DO
         nk=nk/n-nbeta
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
