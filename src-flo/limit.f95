#if defined LIMITED

! Limiter function
      function limit(dl, dr)
      implicit none
      include 'common.h'
      real(dp) :: dl, dr, limit

      real(dp) :: da, TOL, R1, R2, R3, R4

      TOL   = 1.0d-10

      da    = 0.5d0*( dl + dr )
      R1    = (dl - dr)**2 + TOL
      R2    = dl*dl + dr*dr + 2.0d0*abs(dl*dr) + TOL
      R4    = R1/R2
      limit = (1.0d0 - R4)*da

      return
      end

#elif defined UNLIMITED

! No limiting; just average
      function limit(dl, dr)
      implicit none
      include 'common.h'
      real(dp) :: dl, dr, limit

      limit = 0.5d0*(dl + dr)

      return
      end

#endif
