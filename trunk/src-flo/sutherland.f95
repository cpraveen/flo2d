!------------------------------------------------------------------------------
! Sutherland viscosity formula
!------------------------------------------------------------------------------
      function sutherland(q)
      implicit none
      include 'param.h'
      real(dp) :: q(nvar), sutherland

      real(dp) :: temp, num, den

      temp = q(4)/q(1)/GAS_CONST
      num  = T_inf + SCONST
      den  = temp  + SCONST
      sutherland = (temp/T_inf)**1.5d0*(num/den)/Rey

      end
