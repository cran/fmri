      subroutine getvofh(bw,kern,wght,vol)
      implicit logical(a-z)
      integer kern
      double precision bw,wght(2),vol,sofw3D
      external sofw3D
      vol=sofw3D(bw,kern,wght)
      RETURN
      END
      double precision function sofw3D(bw,kern,wght)
      implicit logical(a-z)
      integer kern
      double precision bw,wght(2)
      integer j1,j2,j3,dlw1,dlw2,dlw3,clw1,clw2,clw3,ih1,ih2,ih3
      double precision sw,sw2,h2,lkern,z1,z2,z3,z
      external lkern
      h2=bw*bw
C
C   first calculate location weights
C
      ih3=FLOOR(bw/wght(2))
      ih2=FLOOR(bw/wght(1))
      ih1=FLOOR(bw)
      dlw1=2*ih1+1
      dlw2=2*ih2+1
      dlw3=2*ih3+1
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
      sw=0.d0
      sw2=0.d0
      DO j3=1,dlw3
         z3=(clw3-j3)*wght(2)
         z3=z3*z3
         ih2=FLOOR(sqrt(h2-z3)/wght(1))
         DO j2=clw2-ih2,clw2+ih2
            z2=(clw2-j2)*wght(1)
            z2=z3+z2*z2
            ih1=FLOOR(sqrt(h2-z2))
            DO j1=clw1-ih1,clw1+ih1
               z1=clw1-j1
               z=lkern(kern,(z1*z1+z2)/h2)
               sw=sw+z
               sw2=sw2+z*z
            END DO
         END DO
      END DO
      sofw3D=sw*sw/sw2
      RETURN
      END
      subroutine sofw3Df(bw,kern,wght,fw)
      implicit logical(a-z)
      integer kern
      double precision bw,wght(2),fw,sofw3D
      external sofw3D
      fw=sofw3D(bw,kern,wght)
      RETURN
      END
      
      double precision function sofw3D0(bw,kern,wght)
      implicit logical(a-z)
      integer kern
      double precision bw,wght(2)
      integer j1,j2,j3,dlw1,dlw2,dlw3,clw1,clw2,clw3,ih1,ih2,ih3
      double precision sw,h2,lkern,z1,z2,z3,z
      external lkern
      h2=bw*bw
C
C   first calculate location weights
C
      ih3=FLOOR(bw/wght(2))
      ih2=FLOOR(bw/wght(1))
      ih1=FLOOR(bw)
      dlw1=2*ih1+1
      dlw2=2*ih2+1
      dlw3=2*ih3+1
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
      sw=0.d0
      DO j3=1,dlw3
         z3=(clw3-j3)*wght(2)
         z3=z3*z3
         ih2=FLOOR(sqrt(h2-z3)/wght(1))
         DO j2=clw2-ih2,clw2+ih2
            z2=(clw2-j2)*wght(1)
            z2=z3+z2*z2
            ih1=FLOOR(sqrt(h2-z2))
            DO j1=clw1-ih1,clw1+ih1
               z1=clw1-j1
               z=lkern(kern,(z1*z1+z2)/h2)
               sw=sw+z
            END DO
         END DO
      END DO
      sofw3D0=sw
      RETURN
      END
      subroutine ni2var(bw,kern,wght,quot)
      implicit logical(a-z)
      integer kern
      double precision bw,wght(2),quot
      double precision sofw3D0,sofw3D
      external sofw3D0,sofw3D
      quot=sofw3D0(bw,kern,wght)/sofw3D(bw,kern,wght)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   determine sum of location weights for a given geometry a(3) and given 
C   bandwidth
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Algorithmus zur Nullstellenbestimmung einer monotonen Funktion auf(0,\infty)
      subroutine gethani(x,y,kern,value,wght,eps,bw)
      implicit logical(a-z)
      integer kern
      double precision x,y,value,wght(2),eps,bw
      double precision fw1,fw2,fw3,z
      double precision sofw3D
      external sofw3D
      if(x.ge.y) RETURN
      fw1=sofw3D(x,kern,wght)
      fw2=sofw3D(y,kern,wght)
      DO WHILE(fw1.gt.value)
         x=x*x/y
         fw1=sofw3D(x,kern,wght)
      END DO
      DO WHILE(fw2.le.value)
         y=y*y/x
         fw2=sofw3D(y,kern,wght)
      END DO
      DO WHILE(min(fw2/value,value/fw1).gt.1.d0+eps.and.y-x.gt.1e-6)
C         z=x+(value-fw1)/(fw2-fw1)*(y-x)
         z=(x+y)/2.d0
         fw3=sofw3D(z,kern,wght)
         if(fw3.le.value) THEN
            x=z
            fw1=fw3
         ENDIF
         if(fw3.ge.value) THEN
            y=z
            fw2=fw3
         ENDIF
               call rchkusr()
      END DO
      if(fw2/value.gt.value/fw1) THEN
          bw=x+(value-fw1)/(fw2-fw1)*(y-x)
      ELSE
          bw=y-(fw2-value)/(fw2-fw1)*(y-x)
      ENDIF
      RETURN
      END  

      subroutine mcorrlag(res,mask,n1,n2,n3,nv,scorr,lag)

      implicit logical(a-z)
      integer n1,n2,n3,nv,lag(3)
      double precision scorr,res(nv,n1,n2,n3)
      logical mask(n1,n2,n3)
      double precision z2,y2,resi,resip1,vrm,vrmp1,zk,zcorr,z
      integer i1,i2,i3,i4,l1,l2,l3,k
      zk=nv
      l1=lag(1)
      l2=lag(2)
      l3=lag(3)
      z=0.d0
      k=0
C  correlation in x
      do i1=1,n1-l1
         do i2=1,n2-l2
            do i3=1,n3-l3
         if (.not.(mask(i1,i2,i3).and.mask(i1+l1,i2+l2,i3+l3))) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i4,i1,i2,i3)
                  resip1=res(i4,i1+l1,i2+l2,i3+l3)
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
      subroutine mcorlag1(res,mask,indm,nvox,n1,n2,n3,nv,scorr,lag)
      implicit logical(a-z)
      integer n1,n2,n3,nv,lag(3),nvox,indm(nvox)
      double precision scorr,res(nv,nvox)
      logical mask(n1,n2,n3)
      double precision z2,y2,resi,resip1,vrm,vrmp1,zk,zcorr,z
      integer i1,i2,i3,i4,l1,l2,l3,k,i,j,m,m0
      zk=nv
      l1=lag(1)
      l2=lag(2)
      l3=lag(3)
      z=0.d0
      k=0
C  correlation in x
      do j=1,nvox
         i=indm(j)
         i3=i/n1/n2+1
         i=i-(i3-1)*n1*n2
         if(i.eq.0) THEN
            i=n1*n2
            i3=i3-1
         ENDIF
         i2=i/n1+1
         i1=i-(i2-1)*n1
         if(i1.eq.0) THEN
            i1=n1
            i2=i2-1
         ENDIF
         if(i3.gt.n3-l3) CYCLE
         if(i2.gt.n2-l2) CYCLE
         if(i1.gt.n1-l1) CYCLE
         if (.not.mask(i1+l1,i2+l2,i3+l3)) CYCLE
         z2=0.d0
         y2=0.d0
         zcorr=0.d0
C   get position in indm corresponding to voxel i1+l1,i2+l2,i3+l3
         i=i1+l1+(i2+l2-1)*n1+(i3+l3-1)*n1*n2
         DO m=1,nvox
            m0=m
            if(i.eq.indm(m)) EXIT
         END DO
         do i4=1,nv
            resi=res(i4,j)
            resip1=res(i4,m0)
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
      scorr=z/k
      return
      end
      subroutine sweepm(res,mask,n1,n2,n3,nv)

      implicit logical(a-z)
      integer n1,n2,n3,nv
      double precision res(nv,n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3,k
      double precision z
      Do i1=1,n1
         Do i2=1,n2
            Do i3=1,n3
               if (.not.mask(i1,i2,i3)) CYCLE
               z=0.d0
               DO k=1,nv
                  z=z+res(k,i1,i2,i3)
               END DO
               z=z/nv
               DO k=1,nv
                  res(k,i1,i2,i3)=res(k,i1,i2,i3)-z
               END DO
            END DO
         END DO
      END DO
      return
      end
      subroutine sweepm0(res,n,nv)

      implicit logical(a-z)
      integer n,nv
      double precision res(nv,n)
      integer i,k
      double precision z
      Do i=1,n
         z=0.d0
         DO k=1,nv
            z=z+res(k,i)
         END DO
         z=z/nv
         DO k=1,nv
            res(k,i)=res(k,i)-z
         END DO
      END DO
      return
      end
      subroutine mean3D(res,n1,n2,n3,nv,mres)

      implicit logical(a-z)
      integer n1,n2,n3,nv
      double precision res(nv,n1,n2,n3),mres(n1,n2,n3)
      integer i1,i2,i3,k
      double precision z
      Do i1=1,n1
         Do i2=1,n2
            Do i3=1,n3
               z=0.d0
               DO k=1,nv
                  z=z+res(k,i1,i2,i3)
               END DO
               mres(i1,i2,i3)=z/nv
            END DO
         END DO
      END DO
      return
      end
      subroutine mcorr(res,mask,n1,n2,n3,nv,scorr,l1,l2,l3)

      implicit logical(a-z)
      integer n1,n2,n3,nv,l1,l2,l3,lag(3)
      double precision scorr(l1,l2,l3),res(nv,n1,n2,n3)
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
      subroutine mcorr1(res,mask,indm,nvox,n1,n2,n3,nv,scorr,l1,l2,l3)
      implicit logical(a-z)
      integer n1,n2,n3,nv,l1,l2,l3,lag(3),nvox,indm(nvox)
      double precision scorr(l1,l2,l3),res(nv,n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call mcorlag1(res,mask,indm,nvox,n1,n2,n3,nv,
     1                       scorr(i1,i2,i3),lag)
               call rchkusr()  
            END DO
         END DO
      END DO
      return
      end
      subroutine imcorrl(res,mask,n1,n2,n3,nv,scorr,lag)

      implicit logical(a-z)
      integer n1,n2,n3,nv,lag(3)
      double precision scorr,res(nv,n1,n2,n3)
      logical mask(n1,n2,n3)
      double precision z2,y2,resi,resip1,vrm,vrmp1,zk,zcorr,z
      integer i1,i2,i3,i4,l1,l2,l3,k
      zk=nv
      l1=lag(1)
      l2=lag(2)
      l3=lag(3)
      z=0.d0
      k=0
C  correlation in x
      do i1=1,n1-l1
         do i2=1,n2-l2
            do i3=1,n3-l3
         if (.not.(mask(i1,i2,i3).and.mask(i1+l1,i2+l2,i3+l3))) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i4,i1,i2,i3)
                  resip1=res(i4,i1+l1,i2+l2,i3+l3)
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
      integer n1,n2,n3,nv,l1,l2,l3,lag(3)
      double precision scorr(l1,l2,l3),res(nv,n1,n2,n3)
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

      implicit logical(a-z)
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

      subroutine ivar(res,resscale,mask,n1,n2,n3,nv,var)
C
C   compute variance estimates !!! (not the inverse)
C
      implicit logical(a-z)
      integer n1,n2,n3,nv
      double precision resscale,var(n1,n2,n3),res(nv,n1,n2,n3)
      logical mask(n1,n2,n3)
      double precision z2,zk,resi,ressc2,z1
      integer i1,i2,i3,i4
      zk=nv
      ressc2=resscale*resscale
      do i1=1,n1
         do i2=1,n2
            do i3=1,n3
               var(i1,i2,i3)=1.d20
               if (.not.mask(i1,i2,i3)) CYCLE
               z2=0.d0
               z1=0.d0
               do i4=1,nv
                  resi=res(i4,i1,i2,i3)
                  z1=z1+resi
                  z2=z2+resi*resi
               enddo
               z1 = z1/zk
               z2 = z2/zk
               var(i1,i2,i3)=(z2-z1*z1)*ressc2
            enddo
         enddo
      enddo
      return
      end
      subroutine lconnect(segm,n1,n2,n3,i1,i2,i3,ind1,ind2,ind3,
     1                   checked,mask)
C
C   assumes that we search for a connected region in segm==.TRUE.
C   that contains seed voxel (i1,i2,i3)
C   result: mask == .TRUE. if voxel is connected to seed
      implicit logical (a-z)
      integer n1,n2,n3,i1,i2,i3,ind1(*),ind2(*),ind3(*)
      logical final,checked(*),mask(n1,n2,n3),segm(n1,n2,n3)
      integer j1,j2,j3,k,l1,l2,l3,lind,lind0,n
C     first find pixel close to (i1,i2) with segm(j1,j2)=0
      n=n1*n2*n3
      DO j1=1,n1
         DO j2=1,n2
            DO j3=1,n3
               mask(j1,j2,j3)=.FALSE.
            END DO
         END DO
      END DO
      if(.not.segm(i1,i2,i3)) THEN
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
                     if(segm(j1,j2,j3)) THEN
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
      mask(i1,i2,i3)=.TRUE.
      ind1(1)=i1
      ind2(1)=i2
      ind3(1)=i3
      lind=1
      lind0=1
      DO k=1,n1*n2*n3
         checked(k)=.FALSE.
      END DO
      final=.FALSE.
      DO while(.not.final)
         DO k=1,lind0
            if(checked(k)) CYCLE
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
                     if(segm(j1,j2,j3).and..not.mask(j1,j2,j3)) THEN
                        mask(j1,j2,j3)=.TRUE.
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
      subroutine countsgm(s,n1,n2,n3,anz)
      implicit logical (a-z)
      integer n1,n2,n3,anz(7)
      logical s(n1,n2,n3)
      integer i1,i2,i3,m
      DO i1=2,n1-1
         DO i2=2,n2-1
            DO i3=2,n3-1
               if(.not.s(i1,i2,i3)) CYCLE
               m=1
C  test neighbors
               if(s(i1-1,i2,i3)) m=m+1
               if(s(i1+1,i2,i3)) m=m+1
               if(s(i1,i2-1,i3)) m=m+1
               if(s(i1,i2+1,i3)) m=m+1
               if(s(i1,i2,i3-1)) m=m+1
               if(s(i1,i2,i3+1)) m=m+1
               anz(m)=anz(m)+1
            END DO
         END DO
      END DO
      RETURN
      END
      
      
