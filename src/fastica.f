      subroutine fastICAw(X,T,n,tol,ifun,parifun,maxit,w,work,iwork)
      implicit logical (a-z)
      integer n,T,ifun,maxit
      real*8 X(T,n),w(T,T),tol,work(1),parifun
      logical convgd
C   length(work): 14*T^2+T
C   length(iwork): 8*T
      integer i,j,iter,indS,indU,indV,indwork,indend
      real*8 z,z1,z2
C     Step 1 Orthonormalization of w
      indS = 1
      indU = indS+T*T
      indV = IndU+T*T
      indwork = IndV+T*T
      indend = Indwork+11*T*T
      iter=0
      convgd=.FALSE.
      DO while(iter.lt.maxit.and..not.convgd)
         call orthonor(w,T,work(indS),work(indU),work(indV),
     1              work(indwork),iwork)
C     step 2 generate new guess for w
         convgd=.TRUE.
         DO i=1,T
            call updtw(w,T,X,n,work,ifun,parifun)
            z1=0.d0
            z2=0.d0
            DO j=1,T
               z=work(j)
               z1=z1+w(i,j)*z
               z2=z2+z*z
               w(i,j)=z
            END DO
            if(abs(abs(z1/sqrt(z2))-1).gt.tol) convgd=.FALSE.
         END DO
      END DO
      call orthonor(w,T,work(indS),work(indU),work(indV),
     1              work(indwork),iwork)
      return
      end
      subroutine orthonor(A,m,S,U,VT,work,iwork)
      implicit logical (a-z)
      integer m,iwork(1)
      real*8 A(m,m),S(m),U(m,m),VT(m,m),work(1)
      integer lwork,info
      lwork = 11*m*m
      call DGESDD("A",m,m,A,m,S,U,m,VT,m,work,lwork,IWORK,INFO)
      if(info.ne.0) call intpr("orthonor:info",13,info,1)
      if(info.eq.0.and.work(1).gt.lwork) THEN
         call intpr("value of m",10,m,1)
         call dblepr("optimal lwork",13,work,1)
         call intpr("actual lwork",12,lwork,1)
      END IF
      call DGEMM("N","N",m,m,m,1.d0,U,m,U,m,0.d0,A,m)
      return
      end
      subroutine updtw(w,i,m,x,n,sw,ifun,parifun)
      implicit logical (a-z)
      integer m,n,i,ifun
      real*8 w(m,m),x(m,n),sw(m),parifun
      integer j,k
      real*8 z,gval,hval
      DO k=1,m
         sw(k)=0.d0
      END DO 
      DO j=1,n
         z=0.d0
         DO k=1,m
            z=z+w(i,k)*x(k,j)
         END DO
         call ghval(z,ifun,parifun,gval,hval)
         DO k=1,m
            sw(k)=x(j,j)*gval+w(i,j)*hval
         END DO
      END DO
      RETURN
      END
      subroutine ghval(z,ifun,parifun,gval,hval)
      implicit logical (a-z)
      real*8 z,gval,hval,parifun
      integer ifun
      real*8 y,ya,yb,ys,ye
      SELECT CASE (ifun)
         CASE DEFAULT
            call intpr("Illegal score function number",29,ifun,1)
            RETURN
         CASE (1) 
C  this is G=1/parifun log cosh(parifun*z) 
            gval=tanh(parifun*z)
            y=cosh(parifun*z)
            hval=parifun/y/y
         CASE (2)
C  this is G=-exp(-z^2/2) 
            y = exp(-0.5d0*z*z)
            gval=z*y
            hval=(1.d0-z*z)*y
         CASE (3)
C  this is G=z^4
            y = z*z
            gval=z*y
            hval=3.d0*y
         CASE (4)
C  this is  log(exp(parifun*z) + exp(-z)) a skewed variant of case 1     
            ya=exp(parifun*z)
            yb=exp(-z)
            y=parifun*ya-yb
            ys=ya+yb
            gval=y/ys
            hval=-gval*gval+(parifun*parifun*ya+yb)/ys
         CASE (5)
C  this is G=-exp(-z^2/2*(1-parifun*tanh(z)) a skewed variant of case 2
C  parifun in (-1,1) determines strength of skewness
            y=tanh(z)
            ys=1-y*y
            yb=(1-parifun*y)
            ye=exp(-z*z/2.d0*yb)
            ya=parifun/2.d0*z*z*ys-z*yb
            gval=ye*(parifun/2.d0*z*z*ys - z*yb)
            hval=ye*(-1.d0+parifun*(2.d0*z*ys-z*z*ys*y+y)+ya*ya)
      END SELECT
      RETURN
      END
