C------------------------------------------------------------------------------
C Contribution to gradient from interior edges
C------------------------------------------------------------------------------
      subroutine gradint(x1, x2, x3, qx1, qx2, qx3, qy1, qy2, qy3, 
     +                   a1, a2, a3, qc)
      implicit none
      include 'common.h'
      real(dp) :: x1(2), x2(2), x3(2), qx1(3), qx2(3),
     +            qx3(3), qy1(3), qy2(3), qy3(3),
     +            a1, a2, a3, qc(nvar)

      real(dp) :: u, v, T, dx1, dx2, dx3, dy1, dy2, dy3, a

      u      = qc(2)
      v      = qc(3)
      T      = qc(4)/qc(1)/GAMMA1

      dx1    = x3(1) - x2(1)
      dy1    = x3(2) - x2(2)

      dx2    = x1(1) - x3(1)
      dy2    = x1(2) - x3(2)

      dx3    = x2(1) - x1(1)
      dy3    = x2(2) - x1(2)

      a      = 0.5d0*(-dx3*dy2 + dx2*dy3)
      a1     = a1 + a
      a2     = a2 + a
      a3     = a3 + a

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
