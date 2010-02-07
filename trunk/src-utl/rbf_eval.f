      subroutine rbf_eval(np, coord, nwp, rw, wt)
      implicit none
      integer np, nwp
      double precision    coord(2,np),
     1        rw (2,*),
     2        wt (2,*)

      integer i, j, k
      double precision    dx, dy, dx1, dy1, dr1, rbf

      do i=1,np
         dx = 0.0d0
         dy = 0.0d0
         do k=1,nwp
            dx1 = coord(1,i) - rw(1,k)
            dy1 = coord(2,i) - rw(2,k)
            dr1 = dsqrt(dx1*dx1 + dy1*dy1)
            if(dr1.eq.0.0d0)then
               rbf = 0.0d0
            else
               rbf = dr1*dr1*dlog(dr1)
            endif
            dx  = dx + wt(1,k)*rbf
            dy  = dy + wt(2,k)*rbf
         enddo
         dx       = dx + wt(1,nwp+1) + wt(1,nwp+2)*coord(1,i) +
     1              wt(1,nwp+3)*coord(2,i)
         dy       = dy + wt(2,nwp+1) + wt(2,nwp+2)*coord(1,i) +
     1              wt(2,nwp+3)*coord(2,i)
         coord(1,i) = coord(1,i) + dx
         coord(2,i) = coord(2,i) + dy
      enddo

      end
