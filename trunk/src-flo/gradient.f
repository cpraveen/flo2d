C------------------------------------------------------------------------------
C Compute gradients at vertices used for viscous flux
C qx(1,.) = u_x, qx(2,.) = v_x, qx(3,.) = (pressure/density/(gamma-1))_x
C------------------------------------------------------------------------------
      subroutine gradient(elem, edge, spts, coord, varea, qc, qv, 
     1                    qx, qy)
      implicit none
      include 'param.h'
      integer  :: elem(3,*), edge(2,*), spts(*)
      real(dp) :: coord(2,*), varea(*), qc(nvar,*), qv(nvar,*), qx(3,*),
     1            qy(3,*)

      integer  :: i, j, v1, v2, v3, e1, e2
      real(dp) :: const, fact, qtemp(nvar)

      do i=1,np
         qx(1,i) = 0.0d0
         qx(2,i) = 0.0d0
         qx(3,i) = 0.0d0
         qy(1,i) = 0.0d0
         qy(2,i) = 0.0d0
         qy(3,i) = 0.0d0
      enddo

      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)
         call gradint(coord(1,v1), coord(1,v2), coord(1,v3),
     +                qx(1,v1), qx(1,v2), qx(1,v3),
     +                qy(1,v1), qy(1,v2), qy(1,v3), qc(1,i))
      enddo

      const = 2.0d0/3.0d0
      do i=1,np
         qx(1,i) = const*qx(1,i)
         qx(2,i) = const*qx(2,i)
         qx(3,i) = const*qx(3,i)
         qy(1,i) = const*qy(1,i)
         qy(2,i) = const*qy(2,i)
         qy(3,i) = const*qy(3,i)
      enddo

C Loop over all boundary edges
      do i=nin+1,ne
         v1 = edge(1,i)
         v2 = edge(2,i)
         call gradbnd(coord(1,v1), coord(1,v2), qv(1,v1), qv(1,v2),
     +               qx(1,v1), qx(1,v2), qy(1,v1), qy(1,v2))
      enddo

C Divide by the area
      do i=1,np
         fact    = 1.0d0/varea(i)
         qx(1,i) = qx(1,i)*fact
         qx(2,i) = qx(2,i)*fact
         qx(3,i) = qx(3,i)*fact
         qy(1,i) = qy(1,i)*fact
         qy(2,i) = qy(2,i)*fact
         qy(3,i) = qy(3,i)*fact
      enddo

      end
