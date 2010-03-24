      program main
      implicit none
      include 'param.h'
      include 'deform.h'

      real(dp) :: l1, h1, alpha1, x1, y1
      real(dp) :: l3, h3, alpha3, x3, y3
      integer  :: i, c
      integer  :: n1, n2, n3, nout, nwp, nstep
      real(dp),allocatable :: coord(:,:)
      integer ,allocatable :: elem(:,:)
      integer :: betype(nbemax), bdedge(2,nbemax)
      real(dp) :: xx, yy, x, y, theta, dtheta, drmax=0.0d0

      real(dp), allocatable :: rw(:,:), drw(:,:), wt(:,:)

      call read_input

      allocate( coord(2,np) )
      allocate( elem(3,nt) )

c     Number of points on first,second,third ellipese, and on outer
      n1   = 200
      n2   = 400
      n3   = 200
      nout = 86

      PI = 4.0d0*datan(1.0d0)

c     Semi-major/minor axes of first and third ellipse
      l1 = 2.0d0
      h1 = 0.5d0
      l3 = 3.5d0
      h3 = 0.5d0

      el1_nbeg = 1
      el1_nend = el1_nbeg + n1 - 1
      el2_nbeg = el1_nend + 1
      el2_nend = el2_nbeg + n2 - 1
      el3_nbeg = el2_nend + 1
      el3_nend = el3_nbeg + n3 - 1
      out_nbeg = el3_nend + 1
      out_nend = out_nbeg + nout - 1

c     Read design variables
      open(10, file="control.dat", status="old")
      read(10,*) x1, y1, alpha1, x3, y3, alpha3
      close(10)

      write(*,*) x1, y1, alpha1
      write(*,*) x3, y3, alpha3

      alpha1 = alpha1*PI/180.0d0
      alpha3 = alpha3*PI/180.0d0

      call read_grid(coord, elem, betype, bdedge)

c     Total number of moving/fixed points
      nwp = n1 + n2 + n3 + nout
      allocate( rw(2,nwp) )
      allocate(drw(2,nwp) )
      allocate( wt(2,nwp+3) )

c     Counter
      c = 0

c     First ellipse
      theta = 0.0d0
      dtheta= 2.0d0*PI/n1
      do i=el1_nbeg,el1_nend
         theta = theta - dtheta
         c     = c + 1
         xx    = 0.5d0*l1*dcos(theta)
         yy    = 0.5d0*h1*dsin(theta)
         x     = dcos(alpha1)*xx + dsin(alpha1)*yy + x1
         y     =-dsin(alpha1)*xx + dcos(alpha1)*yy + y1
         drw(1,c) = x - coord(1,i)
         drw(2,c) = y - coord(2,i)
      enddo

c     Second ellipse
      do i=el2_nbeg,el2_nend
         c     = c + 1
         drw(1,c) = 0.0d0
         drw(2,c) = 0.0d0
      enddo

c     Third ellipse
      theta = 0.0d0
      dtheta= 2.0d0*PI/n3
      do i=el3_nbeg,el3_nend
         theta = theta - dtheta
         c     = c + 1
         xx    = 0.5d0*l3*dcos(theta)
         yy    = 0.5d0*h3*dsin(theta)
         x     = dcos(alpha3)*xx + dsin(alpha3)*yy + x3
         y     =-dsin(alpha3)*xx + dcos(alpha3)*yy + y3
         drw(1,c) = x - coord(1,i)
         drw(2,c) = y - coord(2,i)
      enddo

c     Outer boundary
      do i=out_nbeg,out_nend
         c        = c + 1
         drw(1,c) = 0.0d0
         drw(2,c) = 0.0d0
      enddo

      if(c.ne.nwp) stop "c is not equal to nwp"

      do i=1,nwp
         drmax = max(drmax, dsqrt(drw(1,i)**2 + drw(2,i)**2))
      enddo
      print*,'Maximum movement =',drmax

c     Copy boundary points from coord into rw
      call set_bdpts(coord, nwp, rw)

c     Train rbf
      call rbf_train(nwp, rw, drw, wt)

      call divergence(np, coord, nwp, rw, wt, nstep)
      print*,'Number of steps =',nstep

c     Break deformation into smaller steps
      drw = drw/nstep

      do i=1,nstep

         call set_bdpts(coord, nwp, rw)

c        Train rbf
         call rbf_train(nwp, rw, drw, wt)

c        RBF deformation
         call rbf_eval(np, coord, nwp, rw, wt)

      enddo

c     Write out new grid
      print*,'Overwriting grid file grid.fm'
      open(20, file="grid.fm")
      write(20,*) np, nt, nbe
      do i=1,np
         write(20,'(i8,2e24.15)') i, coord(:,i)
      enddo
      do i=1,nt
         write(20,*) i, elem(:,i)
      enddo
      do i=1,nbe
         write(20,*) i, bdedge(:,i), betype(i)
      enddo
      close(20)

      stop
      end
