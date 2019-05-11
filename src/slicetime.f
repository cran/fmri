      function sinc(x)
        implicit none
        double precision :: sinc
        double precision :: x, xpi
        double precision, parameter :: pi = acos(-1d0)
        xpi = x*pi
        sinc = 1.d0
        if(x /= 0) sinc = sin(xpi)/xpi
      end function sinc

      subroutine sincfilter(t,nt,x,nx,ft,wr)
        implicit none
        integer :: nt, nx, wr
        double precision :: t(nt), x(nx), ft(nt), meanx, minx, maxx
        double precision :: sinc
        integer :: i, j, wrm1
        double precision :: y
        external sinc
        wrm1 = wr-1
        meanx = 0.d0
        minx = x(1)
        maxx = minx
        DO i = 1, nx
           meanx = meanx + x(i)
           minx = min(x(i),minx)
           maxx = max(x(i),maxx)
        END DO
        meanx = meanx/nx
        DO i = 1, nt
           y=0.d0
           DO j = 0, wrm1
             y = y + (x(1)-meanx)*sinc(t(i)+j)
           END DO
           DO j = 1, nx
             y = y + (x(j)-meanx)*sinc(t(i)-j)
           END DO
           DO j = nx, nx+wr
             y = y + (x(nx)-meanx)*sinc(t(i)-j)
           END DO
           y = y + meanx
           y = max(y,minx)
           y = min(y,maxx)
           ft(i) = y
        END DO
        return
      END

      subroutine slicetim(x,nt,n1,n2,n3,y,t,sliceord)
        implicit none
        integer :: n1,n2,n3,nt,sliceord(n3)
        double precision :: x(nt,n1,n2,n3), y(nt,n1,n2,n3), t(nt)
        integer :: i1,i2,i3,it,wr
        double precision :: rn3,dt
        wr=8
        rn3 = n3
        DO i3=1,n3
           dt=sliceord(i3)-1
           DO it=1,nt
              t(it)=it-dt/rn3
           END DO
           DO i2=1,n2
              DO i1=1,n1
              call sincfilter(t,nt,x(1,i1,i2,i3),nt,y(1,i1,i2,i3),wr)
              END DO
           END DO
        END DO
        return
      END
