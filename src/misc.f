
      subroutine thcorr(w,n1,n2,n3,scorr,l1,l2,l3)

      implicit none
      integer n1,n2,n3,l1,l2,l3,lag(3)
      double precision scorr(l1,l2,l3),w(n1,n2,n3)
      integer i1,i2,i3
      double precision z,zcorr
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

      implicit none
      integer n1,n2,n3,lag(3)
      double precision scorr,w(n1,n2,n3)
      integer i1,i2,i3,c1,c2,c3,j1,j2,j3,l1,l2,l3
      double precision z
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

      subroutine lconnect(segm,n1,n2,n3,i1,i2,i3,ind1,ind2,ind3,
     1                   checked,mask)
C
C   assumes that we search for a connected region in segm==.TRUE.
C   that contains seed voxel (i1,i2,i3)
C   result: mask == .TRUE. if voxel is connected to seed
      implicit none
      integer n1,n2,n3,i1,i2,i3,ind1(*),ind2(*),ind3(*)
      logical final
      integer checked(*),mask(n1,n2,n3),segm(n1,n2,n3)
      integer j1,j2,j3,k,l1,l2,l3,lind,lind0,n
C     first find pixel close to (i1,i2) with segm(j1,j2)=0
      n=n1*n2*n3
      DO j1=1,n1
         DO j2=1,n2
            DO j3=1,n3
               mask(j1,j2,j3)=0
            END DO
         END DO
      END DO
      if(segm(i1,i2,i3).eq.0) THEN
         final=.FALSE.
         DO k=1,n1
            DO l1=-k,k
               DO l2=-k,k
                  DO l3=-k,k
                     if(max(abs(l1),abs(l2),abs(l3)).ne.k) CYCLE
                     j1=i1+l1
                     if(j1.lt.1.or.j1.gt.n1) CYCLE
                     j2=i2+l2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     j3=i3+l3
                     if(j3.lt.1.or.j3.gt.n3) CYCLE
                     if(segm(j1,j2,j3).ne.0) THEN
                        final=.TRUE.
                        i1=j1
                        i2=j2
                        i3=j3
                     END IF
                     if(final) EXIT
                  END DO
                  if(final) EXIT
               END DO
               if(final) EXIT
            END DO
            if(final) EXIT
         END DO
      END IF
      mask(i1,i2,i3)=1
      ind1(1)=i1
      ind2(1)=i2
      ind3(1)=i3
      lind=1
      lind0=1
      DO k=1,n1*n2*n3
         checked(k)=0
      END DO
      final=.FALSE.
      DO while(.not.final)
         DO k=1,lind0
            if(checked(k).ne.0) CYCLE
            DO l1=-1,1
               DO l2=-1,1
                  DO l3=-1,1
                     if(l1.eq.0.and.l2.eq.0.and.i3.eq.0) CYCLE
                     j1=ind1(k)+l1
                     if(j1.lt.1.or.j1.gt.n1) CYCLE
                     j2=ind2(k)+l2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     j3=ind3(k)+l3
                     if(j3.lt.1.or.j3.gt.n3) CYCLE
                if(segm(j1,j2,j3).ne.0.and.mask(j1,j2,j3).eq.0) THEN
                        mask(j1,j2,j3)=1
                        lind=lind+1
                        if(lind.gt.n) THEN
               call intpr("lconnect: lind exeeds maximum of",32,n,1)
                            return
                        END IF
                        ind1(lind)=j1
                        ind2(lind)=j2
                        ind3(lind)=j3
                     END IF
                  END DO
               END DO
            END DO
         END DO
         if(lind.eq.lind0) THEN
            final=.TRUE.
         ELSE
            lind0=lind
         END IF
      END DO
      RETURN
      END
