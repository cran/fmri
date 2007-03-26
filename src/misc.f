      subroutine corr(res,mask,n1,n2,n3,nv,scorr)

      implicit logical(a-z)
      integer n1,n2,n3,nv
      real*8 scorr(3),res(n1,n2,n3,nv)
      logical mask(n1,n2,n3)
      real*8 z,z2,resi,vrm
      integer i1,i2,i3,i4,k
      
      scorr(1)=0.d0
      scorr(2)=0.d0
      scorr(3)=0.d0
      z=0.d0
      z2=0.d0
      k=0

      do i1=1,n1
         do i2=1,n2
            do i3=1,n3
               if (.not.mask(i1,i2,i3)) CYCLE
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  if (resi.eq.0.d0) CYCLE
                  z=z+resi
                  z2=z2+resi*resi
                  k=k+1
               enddo
            enddo
         enddo
      enddo
      if (k.gt.0) then 
         z=z/k
         vrm=z2/k-(z*z)
         k=0
         z=0.d0
         do i1=1,n1-1
            do i2=1,n2
               do i3=1,n3
                  if (.not.(mask(i1,i2,i3).and.mask(i1+1,i2,i3))) CYCLE
                  do i4=1,nv
                     resi=res(i1,i2,i3,i4)*res(i1+1,i2,i3,i4)
                     if (resi.eq.0.d0) CYCLE
                     z=z+resi
                     k=k+1
                  enddo
               enddo
            enddo
         enddo
         if (k.gt.0) scorr(1)=z/k/vrm
         
         k=0
         z=0.d0
         do i1=1,n1
            do i2=1,n2-1
               do i3=1,n3
                  if (.not.(mask(i1,i2,i3).and.mask(i1,i2+1,i3))) CYCLE
                  do i4=1,nv
                     resi=res(i1,i2,i3,i4)*res(i1,i2+1,i3,i4)
                     if (resi.eq.0.d0) CYCLE
                     z=z+resi
                     k=k+1
                  enddo
               enddo
            enddo
         enddo
         if (k.gt.0) scorr(2)=z/k/vrm
         
         k=0
         z=0.d0
         do i1=1,n1
            do i2=1,n2
               do i3=1,n3-1
                  if (.not.(mask(i1,i2,i3).and.mask(i1,i2,i3+1))) CYCLE
                  do i4=1,nv
                     resi=res(i1,i2,i3,i4)*res(i1,i2,i3+1,i4)
                     if (resi.eq.0.d0) CYCLE
                     z=z+resi
                     k=k+1
                  enddo
               enddo
            enddo
         enddo
         if (k.gt.0) scorr(3)=z/k/vrm
      endif

      return
      end
