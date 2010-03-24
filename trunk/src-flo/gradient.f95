!------------------------------------------------------------------------------
! Compute gradients at vertices used for viscous flux
! qx(1,.) = u_x, qx(2,.) = v_x, qx(3,.) = (pressure/density/(gamma-1))_x
!------------------------------------------------------------------------------
      subroutine gradient(elem, edge, spts, coord, varea, qc, qv, &
                          qx, qy)
      implicit none
      include 'param.h'
      integer  :: elem(3,*), edge(2,*), spts(*)
      real(dp) :: coord(2,*), varea(*), qc(nvar,*), qv(nvar,*), qx(3,*), &
                  qy(3,*)

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
         call gradint(coord(1,v1), coord(1,v2), coord(1,v3), &
                      qx(1,v1), qx(1,v2), qx(1,v3), &
                      qy(1,v1), qy(1,v2), qy(1,v3), qc(1,i))
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

! Loop over all boundary edges
      do i=nin+1,ne
         v1 = edge(1,i)
         v2 = edge(2,i)
         call gradbnd(coord(1,v1), coord(1,v2), qv(1,v1), qv(1,v2), &
                     qx(1,v1), qx(1,v2), qy(1,v1), qy(1,v2))
      enddo

! Divide by the area
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
!------------------------------------------------------------------------------
! Contribution to gradient from interior edges
!------------------------------------------------------------------------------
      subroutine gradint(x1, x2, x3, qx1, qx2, qx3, qy1, qy2, qy3, qc)
      implicit none
      include 'common.h'
      real(dp) :: x1(2), x2(2), x3(2), qx1(3), qx2(3), &
                  qx3(3), qy1(3), qy2(3), qy3(3), qc(nvar)

      real(dp) :: u, v, T, dx1, dx2, dx3, dy1, dy2, dy3

      u      = qc(2)
      v      = qc(3)
      T      = qc(4)/qc(1)/GAMMA1

      dx1    = x3(1) - x2(1)
      dy1    = x3(2) - x2(2)

      dx2    = x1(1) - x3(1)
      dy2    = x1(2) - x3(2)

      dx3    = x2(1) - x1(1)
      dy3    = x2(2) - x1(2)

      qx1(1) = qx1(1) + u*dy1
      qx1(2) = qx1(2) + v*dy1
      qx1(3) = qx1(3) + T*dy1

      qy1(1) = qy1(1) - u*dx1
      qy1(2) = qy1(2) - v*dx1
      qy1(3) = qy1(3) - T*dx1

      qx2(1) = qx2(1) + u*dy2
      qx2(2) = qx2(2) + v*dy2
      qx2(3) = qx2(3) + T*dy2

      qy2(1) = qy2(1) - u*dx2
      qy2(2) = qy2(2) - v*dx2
      qy2(3) = qy2(3) - T*dx2

      qx3(1) = qx3(1) + u*dy3
      qx3(2) = qx3(2) + v*dy3
      qx3(3) = qx3(3) + T*dy3

      qy3(1) = qy3(1) - u*dx3
      qy3(2) = qy3(2) - v*dx3
      qy3(3) = qy3(3) - T*dx3

      end
!------------------------------------------------------------------------------
! Contribution to gradient from boundary edges
!------------------------------------------------------------------------------
      subroutine gradbnd(x1, x2, qv1, qv2, qx1, qx2, qy1, qy2)
      implicit none
      include 'common.h'
      real(dp) :: x1(2), x2(2), qv1(nvar), qv2(nvar), qx1(3), &
                  qx2(3), qy1(3), qy2(3)

      real(dp) :: dx, dy, u1, v1, u2, v2, T1, T2, const

      const  = 1.0d0/9.0d0

      dx     = const*( x2(1) - x1(1) )
      dy     = const*( x2(2) - x1(2) )

      u1     = qv1(2)
      v1     = qv1(3)
      T1     = qv1(4)/qv1(1)/GAMMA1

      u2     = qv2(2)
      v2     = qv2(3)
      T2     = qv2(4)/qv2(1)/GAMMA1

      qx1(1) = qx1(1) + (4.0d0*u1 + 2.0d0*u2)*dy
      qx1(2) = qx1(2) + (4.0d0*v1 + 2.0d0*v2)*dy
      qx1(3) = qx1(3) + (4.0d0*T1 + 2.0d0*T2)*dy

      qy1(1) = qy1(1) - (4.0d0*u1 + 2.0d0*u2)*dx
      qy1(2) = qy1(2) - (4.0d0*v1 + 2.0d0*v2)*dx
      qy1(3) = qy1(3) - (4.0d0*T1 + 2.0d0*T2)*dx

      qx2(1) = qx2(1) + (2.0d0*u1 + 4.0d0*u2)*dy
      qx2(2) = qx2(2) + (2.0d0*v1 + 4.0d0*v2)*dy
      qx2(3) = qx2(3) + (2.0d0*T1 + 4.0d0*T2)*dy

      qy2(1) = qy2(1) - (2.0d0*u1 + 4.0d0*u2)*dx
      qy2(2) = qy2(2) - (2.0d0*v1 + 4.0d0*v2)*dx
      qy2(3) = qy2(3) - (2.0d0*T1 + 4.0d0*T2)*dx

      end
