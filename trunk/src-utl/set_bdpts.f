      subroutine set_bdpts(coord, nwp, rw)
      implicit none
      include 'param.h'
      include 'deform.h'

      integer  :: nwp
      real(dp) :: coord(2,*), rw(2,*)

      integer  :: i, c

c     Counter
      c = 0

c     First ellipse
      do i=el1_nbeg,el1_nend
         c     = c + 1
         rw (1,c) = coord(1,i)
         rw (2,c) = coord(2,i)
      enddo

c     Second ellipse
      do i=el2_nbeg,el2_nend
         c     = c + 1
         rw (1,c) = coord(1,i)
         rw (2,c) = coord(2,i)
      enddo

c     Third ellipse
      do i=el3_nbeg,el3_nend
         c     = c + 1
         rw (1,c) = coord(1,i)
         rw (2,c) = coord(2,i)
      enddo

c     Outer boundary
      do i=out_nbeg,out_nend
         c        = c + 1
         rw (1,c) = coord(1,i)
         rw (2,c) = coord(2,i)
      enddo

      if(c.ne.nwp) stop "set_bdpts: c is not equal to nwp"

      end
