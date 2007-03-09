      real*8 function f1(z,s)
      real*8 z,s,z2
      z2=z*z
      f1=1.66d0*z*z2*dexp(-z2/2.d0/s/s)
      RETURN
      END
      real*8 function f1d1(z,s)
      real*8 z,s,z2,s2
      z2=z*z
      s2=s*s
      f1d1=1.66d0*z2*(3.d0-z2/s2)*dexp(-z2/2.d0/s2)
      RETURN
      END
      real*8 function f2(z,s)
      real*8 z,s
      f2=2.d0*dtanh(s*z)
      RETURN
      END
      real*8 function f2d1(z,s)
      real*8 z,s,chsz
      chsz=dcosh(z*s)
      f2d1=2.d0*s/chsz/chsz
      RETURN
      END
      real*8 function f3(z,s)
      real*8 z,s
      f3=dsin(s*z)
      RETURN
      END
      real*8 function f3d1(z,s)
      real*8 z,s
      f3d1=s*dcos(s*z)
      RETURN
      END
      real*8 function f4(z,s)
      real*8 z,s
      f4=dcos(s*z)
      RETURN
      END
      real*8 function f4d1(z,s)
      real*8 z,s
      f4d1=-s*dsin(s*z)
      RETURN
      END
      Subroutine fastica(y,omega,d,n,l,t,beta,v,normv,s)
      integer d,n,l,t
      real*8 y(d,n),omega(d,l,4),v(d,l,4),beta(d),s(l,4),normv(l,4)
      real*8 f1,f1d1,f2,f2d1,f3,f3d1,f4,f4d1
      integer i,j,k,o,ifunct
      real*8 nk,z,fw,fwd1,nbeta,srnbeta,zv
      nk=0.d0
      nbeta=0.d0
C this is just to avoid some warnings 
      DO ifunct=1,4
      DO k=1,l
         DO j=1,t
            nk=0.d0
            DO o=1,d
               beta(o)=0.d0
            END DO
            DO i=1,n
               z=0.d0
               DO o=1,d
                  z=z+omega(o,k,ifunct)*y(o,i)
               END DO
               select case (ifunct)
               case(1) 
                  fw=f1(z,s(k,ifunct))
                  fwd1=f1d1(z,s(k,ifunct))
               case(2) 
                  fw=f2(z,s(k,ifunct))
                  fwd1=f2d1(z,s(k,ifunct))
               case(3) 
                  fw=f3(z,s(k,ifunct))
                  fwd1=f3d1(z,s(k,ifunct))
               case(4) 
                  fw=f4(z,s(k,ifunct))
                  fwd1=f4d1(z,s(k,ifunct))
               case default
                  fw=0.d0
                  fwd1=1.d0
               end select
               DO o=1,d
                  z=y(o,i)*fw-omega(o,k,ifunct)*fwd1
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
            srnbeta=dsqrt(nbeta)
            DO o=1,d
               omega(o,k,ifunct)=beta(o)/srnbeta
            END DO
         END DO
         nk=nk/n-nbeta
         nk=dsqrt(n/nk)
         zv=0.d0
         DO o=1,d
            z=beta(o)*nk
            v(o,k,ifunct)=z
            zv=zv+z*z
         END DO
         IF(nbeta.gt.1d-4) THEN
            normv(k,ifunct)=dsqrt(zv)
         ELSE
            normv(k,ifunct)=0.d0
         END IF
         call rchkusr()
      END DO
      END DO      
      RETURN
      END
 