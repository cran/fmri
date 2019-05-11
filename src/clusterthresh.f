
      subroutine ccluster(x,n1,n2,n3,z)
        implicit none
        integer :: n1,n2,n3,x(n1,n2,n3),z(n1,n2,n3)
C   x - cluster size, initially 1 if above threshold, zero else
C   z - cluster number
        integer :: i1,i2,i3,k,l,nk,nl
C initialize
        k=0
        DO i1 = 1,n1
          DO i2 = 1,n2
            DO i3 = 1,n3
              if(x(i1,i2,i3).eq.0) CYCLE
              k=k+1
              z(i1,i2,i3)=k
            END DO
          END DO
        END DO
C determine clusters
        DO i1 = 1,n1
          DO i2 = 1,n2
            DO i3 = 1,n3
              if(x(i1,i2,i3).eq.0) CYCLE
              if(i1.lt.n1.and.x(i1+1,i2,i3).gt.0) THEN
                 if(z(i1,i2,i3).ne.z(i1+1,i2,i3)) THEN
                   k=z(i1,i2,i3)
                   l=z(i1+1,i2,i3)
                   nl=x(i1+1,i2,i3)
                   nk=x(i1,i2,i3)
                   call jcluster(x,z,n1*n2*n3,k,l,nk,nl)
                END IF
              END IF
              if(i2.lt.n2.and.x(i1,i2+1,i3).gt.0) THEN
                 if(z(i1,i2,i3).ne.z(i1,i2+1,i3)) THEN
                   k=z(i1,i2,i3)
                   l=z(i1,i2+1,i3)
                   nl=x(i1,i2+1,i3)
                   nk=x(i1,i2,i3)
                   call jcluster(x,z,n1*n2*n3,k,l,nk,nl)
                END IF
              END IF
              if(i3.lt.n3.and.x(i1,i2,i3+1).gt.0) THEN
                 if(z(i1,i2,i3).ne.z(i1,i2,i3+1)) THEN
                   k=z(i1,i2,i3)
                   l=z(i1,i2,i3+1)
                   nl=x(i1,i2,i3+1)
                   nk=x(i1,i2,i3)
                   call jcluster(x,z,n1*n2*n3,k,l,nk,nl)
                END IF
              END IF
            END DO
          END DO
        END DO
      return
      END


      subroutine jcluster(y,z,n,k,l,nk,nl)
C  join clusters with numbers k and l
        implicit none
        integer :: n,k,l,nk,nl,y(n),z(n)
        integer :: i
        DO i=1,n
          if(z(i).eq.k) THEN
            y(i)=y(i)+nl
          END IF
          if(z(i).eq.l) THEN
            y(i)=y(i)+nk
            z(i)=k
          END IF
        END DO
        return
      END
