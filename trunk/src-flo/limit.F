#if defined LIMITED

C Limiter function
      double precision function limit(dl, dr)
      implicit none
      double precision dl, dr
      intrinsic        dmax1

      double precision da, TOL, R1, R2, R3

      TOL   = 1.0d-10

      da    = 0.5d0*( dl + dr )
      R1    = dabs(dl - dr)
      R2    = dmax1( dabs(dl) + dabs(dr), TOL )
      R3    = (R1/R2)**3
      limit = (1.0d0 - R3)*da

      return
      end

#elif defined UNLIMITED

C No limiting; just average
      double precision function limit(dl, dr)
      implicit none
      double precision dl, dr

      limit = 0.5d0*(dl + dr)

      return
      end

#endif
