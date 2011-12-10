!-----------------------------------------------------------------------------
!.....Read grid data from a file
!.....Currently supports only triangular elements
!-----------------------------------------------------------------------------
      subroutine read_grid(coord, elem, betype, bdedge)
      implicit none
      include 'param.h'
      real(dp) :: coord(2,*)
      integer  :: elem(3,*), betype(*), bdedge(2,*)

      integer  :: ngrid, ip, i, j

      print*,'Reading grid from file ',trim(gridfile)

      ngrid = 10
      open(ngrid, file=gridfile, status="old")
      rewind(ngrid)
      read(ngrid,*) np, nt, nbe
      write(*, '( " Number of points         :", i8)') np
      write(*, '( " Number of triangles      :", i8)') nt
      write(*, '( " Number of boundary edges :", i8)') nbe

!     Read point coordinates
      do ip=1,np
         read(ngrid,*) i, coord(1,ip), coord(2,ip)
      enddo

!     Read triangles
      do ip=1,nt
         read(ngrid,*) i, elem(1, ip), elem(2, ip), elem(3, ip)
      enddo

!     Read boundary edges
      do i=1,nbe
         read(ngrid,*) j, bdedge(1,i), bdedge(2,i), betype(i)
      enddo

      close(ngrid)

!     Use this for finding sensitivity wrt AOA: rotate grid instead of
!     rotating the free-stream velocity vector.
!     call rotate_grid(coord)

!     Find bounding box
      xmin = 1000000.0d0
      ymin = 1000000.0d0
      xmax =-1000000.0d0
      ymax =-1000000.0d0
      do ip=1,np
         xmin = dmin1(xmin, coord(1,ip))
         ymin = dmin1(ymin, coord(2,ip))

         xmax = dmax1(xmax, coord(1,ip))
         ymax = dmax1(ymax, coord(2,ip))
      enddo

      write(*,'(" Bounding box:")')
      write(*,'(10x, "xmin =", f8.3)') xmin
      write(*,'(10x, "xmax =", f8.3)') xmax
      write(*,'(10x, "ymin =", f8.3)') ymin
      write(*,'(10x, "ymax =", f8.3)') ymax

      return
      end
