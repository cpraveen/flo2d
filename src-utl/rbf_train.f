      subroutine rbf_train(nwp, rw, drw, wt)
      implicit none
      integer nwp
      double precision    rw  (2,*),
     1        drw (2,*),
     2        wt  (2,*)

      integer i, j, n, ipvt(nwp+2+1)
      double precision    dx, dy, ds, z(nwp+2+1), rcond
      double precision    a(nwp+2+1, nwp+2+1), b(nwp+2+1)

      print*,'Constructing RBF approximation'

      n      = nwp+2+1
      a(:,:) = 0.0d0

c     construct rbf coefficient matrix
      do i=1,nwp
         a(i,i) = 0.0d0
         do j=i+1,nwp
            dx = rw(1,i) - rw(1,j)
            dy = rw(2,i) - rw(2,j)
            ds = dsqrt(dx*dx + dy*dy)
            a(i,j) = ds*ds*dlog(ds)
         enddo
         a(i,nwp+1) = 1.0d0
         a(i,nwp+2) = rw(1,i)
         a(i,nwp+3) = rw(2,i)
      enddo

c     Copy lower diagonal part due to symmetry
      do i=1,nwp
         do j=1,i-1
            a(i,j) = a(j,i)
         enddo
      enddo

c     Extra equations due to linear polynomial
      do i=1,nwp
         a(nwp+1,i) = 1.0d0
         a(nwp+2,i) = rw(1,i)
         a(nwp+3,i) = rw(2,i)
      enddo

c     Perform LU Decomposition
      call dgeco(a, n, n, ipvt, rcond, z)

      write(*,'("    Condition number = ",e10.4)') 1.0d0/rcond

c     Solve for x-coordinate
      do i=1,nwp
         b(i) = drw(1,i)
      enddo
      b(nwp+1) = 0.0d0
      b(nwp+2) = 0.0d0
      b(nwp+3) = 0.0d0
      call dgesl(a, n, n, ipvt, b, 0)
      do i=1,nwp+3
         wt(1,i) = b(i)
      enddo

c     Solve for y-coordinate
      do i=1,nwp
         b(i) = drw(2,i)
      enddo
      b(nwp+1) = 0.0d0
      b(nwp+2) = 0.0d0
      b(nwp+3) = 0.0d0
      call dgesl(a, n, n, ipvt, b, 0)
      do i=1,nwp+3
         wt(2,i) = b(i)
      enddo

      return
      end
