#if defined LIMITED

C Limiter function
      double precision function limit(dl, dr)
      implicit none
      double precision dl, dr
      intrinsic        dmax1

      double precision da, TOL, R1, R2, R3, R4

      TOL   = 1.0d-10

      da    = 0.5d0*( dl + dr )
      R1    = (dl - dr)**2 + TOL
      R2    = dl*dl + dr*dr + 2.0d0*dabs(dl*dr) + TOL
      R4    = R1/R2
      limit = (1.0d0 - R4)*da

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
