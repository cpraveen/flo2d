!-----------------------------------------------------------------------------
! convert conserved variable to primitive variable
!-----------------------------------------------------------------------------
      subroutine con2prim(con, prim)
      implicit none
      include 'common.h'
      real(dp) :: prim(nvar), con(nvar)

      prim(1) = con(1)
      prim(2) = con(2)/con(1)
      prim(3) = con(3)/con(1)
      prim(4) = GAMMA1*(con(4) - 0.5d0*(con(2)**2 + con(3)**2)/con(1))

      end
