      subroutine slight(stat,mask,n1,n2,n3,slght,nsl,slstat)
C
C  compute mean over searchlights
C
      implicit none
      integer n1,n2,n3,nsl,slght(3,nsl)
      integer mask(n1,n2,n3)
      double precision stat(n1,n2,n3),slstat(n1,n2,n3)
      integer i1,i2,i3,anz,k,j1,j2,j3
      double precision z
      do i1=1,n1
        do i2=1,n2
          do i3=1,n3
             if(mask(i1,i2,i3).eq.0) CYCLE
             z=0.d0
             anz=0
             do k=1,nsl
               j1=i1+slght(1,k)
               if(j1.lt.1.or.j1.gt.n1) CYCLE
               j2=i2+slght(2,k)
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               j3=i3+slght(3,k)
               if(j3.lt.1.or.j3.gt.n3) CYCLE
               if(mask(j1,j2,j3).eq.0) CYCLE
               anz=anz+1
               z=z+stat(j1,j2,j3)
             end do
             slstat(i1,i2,i3)=z/anz
          END DO
        END DO
      END DO
      return
      END
      subroutine getslpv(stat,n,p,kv,nsim,pval)
C
C  compute pvalues using empirical distribution
C
      implicit none
      integer n,nsim
      double precision stat(n),p(nsim),kv(nsim),pval(n)
      integer i,j
      double precision si,kvmax
      kvmax=kv(nsim)
      Do i=1,n
         si=stat(i)
         pval(i)=1
         if(si.ge.kvmax) THEN
            pval(i)=1-p(nsim)
            CYCLE
         END if
         DO j=nsim-1,1,-1
            if(si.ge.kv(j)) THEN
               pval(i)=1-p(j)
               EXIT
            END IF
         END DO
      end do
      return
      END

      subroutine extrpatt(beta,voxel,n1,n2,n3,nb,sl,nsl,pattern,nvox)
      implicit none
      integer n1,n2,n3,nb,nsl,nvox,sl(3,nsl)
      integer voxel(n1,n2,n3)
      double precision beta(n1,n2,n3,nb), pattern(nb,nsl,nvox)
      integer i1,i2,i3,ib,ivox,k,j1,j2,j3
      ivox = 0
      do i3=1,n3
        do i2=1,n2
          do i1=1,n1
             if(voxel(i1,i2,i3).eq.0) CYCLE
             ivox=ivox+1
             do k=1,nsl
               j1=i1+sl(1,k)
               if(j1.lt.1.or.j1.gt.n1) CYCLE
               j2=i2+sl(2,k)
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               j3=i3+sl(3,k)
               if(j3.lt.1.or.j3.gt.n3) CYCLE
               do ib=1,nb
                  pattern(ib,k,ivox)=beta(i1,i2,i3,ib)
               end do
             end do
          END DO
        END DO
      END DO
      return
      END
