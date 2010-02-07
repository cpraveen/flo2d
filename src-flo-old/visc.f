C------------------------------------------------------------------------------
C Add viscous fluxes to the FV residual
C------------------------------------------------------------------------------
      subroutine viscflux(edge, tedge, coord, qv, qx, qy, res)
      implicit none
      include 'param.h'
      integer          edge(2,nemax), tedge(2,nemax)
      double precision coord(2,npmax), qv(nvar,npmax), qx(3,npmax),
     +                 qy(3,npmax), res(nvar,ntmax)

      integer          i, e1, e2, c1, c2
      double precision resdummy(nvar)

C Viscous flux for interior edges
      do i=1,nin
         e1 = edge(1,i)
         e2 = edge(2,i)
         c1 = tedge(1,i)
         c2 = tedge(2,i)
         call viscres(coord(1,e1), coord(1,e2), qv(1,e1), qv(1,e2),
     +                qx(1,e1), qx(1,e2), qy(1,e1), qy(1,e2),
     +                res(1,c1), res(1,c2))
      enddo

C Viscous flux for boundary edges
      do i=1,nvar
         resdummy(i) = 0.0d0
      enddo

      do i=nin+1,ne
         e1 = edge(1,i)
         e2 = edge(2,i)
         c1 = tedge(1,i)
         call viscres(coord(1,e1), coord(1,e2), qv(1,e1), qv(1,e2),
     +                qx(1,e1), qx(1,e2), qy(1,e1), qy(1,e2),
     +                res(1,c1), resdummy)
      enddo

      return
      end

C------------------------------------------------------------------------------
C For an edge, Add contribution of viscous fluxes to the residual vector
C------------------------------------------------------------------------------
      subroutine viscres(x1, x2, qv1, qv2, qx1, qx2, qy1, qy2, 
     +                   res1, res2)
      implicit none
      include 'common.h'
      include 'visc.h'
      double precision x1(2), x2(2), qv1(nvar), qv2(nvar), qx1(3),
     +                 qx2(3), qy1(3), qy2(3), res1(nvar), res2(nvar)

      double precision u, v, txx1, txx2, tyy1, tyy2, 
     +                 txy1, txy2, txx, txy, tyy, ex1, ey1, ex2, ey2,
     +                 ex, ey, xmom, ymom, ener, mul1, mul2, mut1, mut2,
     +                 mu1, mu2, dx, dy, gp1, gp2, n2b3, sutherland
      external         sutherland


      n2b3 = 2.0d0/3.0d0

      mul1 = sutherland(qv1)
      mul2 = sutherland(qv2)

      mu1  = mul1
      mu2  = mul2

      txx1 = n2b3*mu1*(2.0d0*qx1(1) - qy1(2))
      tyy1 = n2b3*mu1*(2.0d0*qy1(2) - qx1(1))
      txy1 =      mu1*(      qy1(1) + qx1(2))

      txx2 = n2b3*mu2*(2.0d0*qx2(1) - qy2(2))
      tyy2 = n2b3*mu2*(2.0d0*qy2(2) - qx2(1))
      txy2 =      mu2*(      qy2(1) + qx2(2))

      txx  = 0.5d0*(txx1 + txx2)
      txy  = 0.5d0*(txy1 + txy2)
      tyy  = 0.5d0*(tyy1 + tyy2)

      gp1  = gamma*(mul1/prandtl)
      ex1  = -gp1*qx1(3)
      ey1  = -gp1*qy1(3)

      gp2  = gamma*(mul2/prandtl)
      ex2  = -gp2*qx2(3)
      ey2  = -gp2*qy2(3)

      ex   = 0.5d0*(ex1 + ex2)
      ey   = 0.5d0*(ey1 + ey2)

      dx   = x2(1) - x1(1)
      dy   = x2(2) - x1(2)

      u    = 0.5d0*(qv1(2) + qv2(2))
      v    = 0.5d0*(qv1(3) + qv2(3))

      xmom = txx*dy - txy*dx
      ymom = txy*dy - tyy*dx
      ener = (txx*u + txy*v - ex)*dy - (txy*u + tyy*v - ey)*dx

      res1(2) = res1(2) - xmom
      res1(3) = res1(3) - ymom
      res1(4) = res1(4) - ener

      res2(2) = res2(2) + xmom
      res2(3) = res2(3) + ymom
      res2(4) = res2(4) + ener

      return
      end
