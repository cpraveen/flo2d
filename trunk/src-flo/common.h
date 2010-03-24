      integer, parameter :: dp = kind(1.0d0)

      integer, parameter :: no=0
      integer, parameter :: yes=1

! nvar = number of variables, fixed at 5, which includes turbulent
! viscosity
      integer, parameter :: nvar=4

!     GAMMA = ratio of specific heats
!     GAMMA1= GAMMA - 1
!     GAS_CONST = gas constant, this can be set to 1.0
!     M_PI = value of pi
      real(dp) :: GAMMA, GAMMA1, GAS_CONST, M_PI
      common/const/GAMMA, GAMMA1, GAS_CONST, M_PI

