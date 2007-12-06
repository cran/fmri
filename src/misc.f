      subroutine mcorrlag(res,mask,n1,n2,n3,nv,scorr,lag)

      implicit logical(a-z)
      integer n1,n2,n3,nv,lag(3)
      real*8 scorr,res(n1,n2,n3,nv)
      logical mask(n1,n2,n3)
      real*8 z2,y2,resi,resip1,vrm,vrmp1,zk,zcorr,z
      integer i1,i2,i3,i4,l1,l2,l3,k
      zk=nv
      l1=lag(1)
      l2=lag(2)
      l3=lag(3)
      z=0.d0
      k=0.d0
C  correlation in x
      do i1=1,n1-l1
         do i2=1,n2-l2
            do i3=1,n3-l3
         if (.not.(mask(i1,i2,i3).and.mask(i1+l1,i2+l2,i3+l3))) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1+l1,i2+l2,i3+l3,i4)
                  z2=z2+resi*resi
                  y2=y2+resip1*resip1
                  zcorr=zcorr+resi*resip1
               enddo
               vrm=z2/zk
               vrmp1=y2/zk
               vrm=vrm*vrmp1
               if(vrm.gt.1e-10) THEN
                  z=z+zcorr/zk/sqrt(vrm)
                  k=k+1
               end if
            enddo
         enddo
      enddo
      scorr=z/k
      return
      end
      subroutine sweepm(res,mask,n1,n2,n3,nv)

      implicit logical(a-z)
      integer n1,n2,n3,nv
      real*8 res(n1,n2,n3,nv)
      logical mask(n1,n2,n3)
      integer i1,i2,i3,k
      real*8 z
      Do i1=1,n1
         Do i2=1,n2
            Do i3=1,n3
               if (.not.mask(i1,i2,i3)) CYCLE
               z=0.d0
               DO k=1,nv
                  z=z+res(i1,i2,i3,k)
               END DO
               z=z/nv
               DO k=1,nv
                  res(i1,i2,i3,k)=res(i1,i2,i3,k)-z
               END DO
            END DO
         END DO
      END DO
      return
      end
      subroutine mcorr(res,mask,n1,n2,n3,nv,scorr,l1,l2,l3)

      implicit logical(a-z)
      integer n1,n2,n3,nv,l1,l2,l3,lag(3)
      real*8 scorr(l1,l2,l3),res(n1,n2,n3,nv)
      logical mask(n1,n2,n3)
      integer i1,i2,i3
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call mcorrlag(res,mask,n1,n2,n3,nv,scorr(i1,i2,i3),lag)
               call rchkusr()  
            END DO
         END DO
      END DO
      return
      end
      subroutine imcorrl(res,mask,n1,n2,n3,nv,scorr,lag)

      implicit logical(a-z)
      integer n1,n2,n3,nv,lag(3),res(n1,n2,n3,nv)
      real*8 scorr
      logical mask(n1,n2,n3)
      real*8 z2,y2,resi,resip1,vrm,vrmp1,zk,zcorr,z
      integer i1,i2,i3,i4,l1,l2,l3,k
      zk=nv
      l1=lag(1)
      l2=lag(2)
      l3=lag(3)
      z=0.d0
      k=0.d0
C  correlation in x
      do i1=1,n1-l1
         do i2=1,n2-l2
            do i3=1,n3-l3
         if (.not.(mask(i1,i2,i3).and.mask(i1+l1,i2+l2,i3+l3))) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1+l1,i2+l2,i3+l3,i4)
                  z2=z2+resi*resi
                  y2=y2+resip1*resip1
                  zcorr=zcorr+resi*resip1
               enddo
               vrm=z2/zk
               vrmp1=y2/zk
               vrm=vrm*vrmp1
               if(vrm.gt.1e-10) THEN
                  z=z+zcorr/zk/sqrt(vrm)
                  k=k+1
               end if
            enddo
         enddo
      enddo
      scorr=z/k
      return
      end
      subroutine imcorr(res,mask,n1,n2,n3,nv,scorr,l1,l2,l3)

      implicit logical(a-z)
      integer n1,n2,n3,nv,l1,l2,l3,lag(3),res(n1,n2,n3,nv)
      real*8 scorr(l1,l2,l3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call imcorrl(res,mask,n1,n2,n3,nv,scorr(i1,i2,i3),lag)
               call rchkusr()  
            END DO
         END DO
      END DO
      return
      end
      subroutine thcorr(w,n1,n2,n3,scorr,l1,l2,l3)

      implicit logical(a-z)
      integer n1,n2,n3,l1,l2,l3,lag(3)
      real*8 scorr(l1,l2,l3),w(n1,n2,n3)
      integer i1,i2,i3
      real*8 z,zcorr
      z=0.d0
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               z=z+w(i1,i2,i3)*w(i1,i2,i3)
            END DO
         END DO
      END DO
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call thcorlag(w,n1,n2,n3,zcorr,lag)
               scorr(i1,i2,i3)=zcorr/z
               call rchkusr()  
            END DO
         END DO
      END DO
      return
      end
      subroutine thcorlag(w,n1,n2,n3,scorr,lag)

      implicit logical(a-z)
      integer n1,n2,n3,lag(3)
      real*8 scorr,w(n1,n2,n3)
      integer i1,i2,i3,c1,c2,c3,j1,j2,j3,l1,l2,l3
      real*8 z
      c1=(n1-1)/2
      c2=(n2-1)/2
      c3=(n3-1)/2
      z=0.d0
      Do i1=-c1,c1
         j1=i1+c1+1 
         l1=lag(1)-i1+c1+1
         if(l1.lt.1.or.l1.gt.n1) CYCLE
         DO i2=-c2,c2
            j2=i2+c2+1 
            l2=lag(2)-i2+c2+1
            if(l2.lt.1.or.l2.gt.n2) CYCLE
            DO i3=-c3,c3
               j3=i3+c3+1
               l3=lag(3)-i3+c3+1
               if(l3.lt.1.or.l3.gt.n3) CYCLE
               z=z+w(j1,j2,j3)*w(l1,l2,l3)
            END DO
         END DO
      END DO
      scorr=z
      return
      end

      subroutine ivar(res,resscale,mask,n1,n2,n3,nv,var)

      implicit logical(a-z)
      integer n1,n2,n3,nv,res(n1,n2,n3,nv)
      real*8 resscale,var(n1,n2,n3)
      logical mask(n1,n2,n3)
      real*8 z2,zk,resi,ressc2
      integer i1,i2,i3,i4
      zk=nv
      ressc2=resscale*resscale
      do i1=1,n1
         do i2=1,n2
            do i3=1,n3
               var(i1,i2,i3)=1.d20
               if (.not.mask(i1,i2,i3)) CYCLE
               z2=0.d0
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  z2=z2+resi*resi
               enddo
               var(i1,i2,i3)=z2/(zk-1.d0)*ressc2
            enddo
         enddo
      enddo
      return
      end
