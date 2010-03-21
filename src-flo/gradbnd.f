C------------------------------------------------------------------------------
C Contribution to gradient from boundary edges
C------------------------------------------------------------------------------
      subroutine gradbnd(x1, x2, qv1, qv2, qx1, qx2, qy1, qy2)
      implicit none
      include 'common.h'
      real(dp) :: x1(2), x2(2), qv1(nvar), qv2(nvar), qx1(3),
     +            qx2(3), qy1(3), qy2(3)

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
