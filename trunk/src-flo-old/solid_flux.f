      subroutine solid_flux(x1, x2, qc, res)
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), qc(nvar), res(nvar)

      double precision fx, fy

      fx=  qc(4)*(x2(2) - x1(2))
      fy= -qc(4)*(x2(1) - x1(1))

      res(2) = res(2) + fx
      res(3) = res(3) + fy

      return
      end
