      subroutine pres_match(edge, coord, qv)
      implicit none
      include 'param.h'
      real(dp) :: coord(2,*), qv(nvar,*)
      integer  :: edge(2,*)

      real(dp) :: ptarg(np)
      real(dp) :: dx, dy, ds, f1, f2, obj
      integer  :: i, j, p1, p2

      open(30, file='ptarg.dat', status='old')
      read(30,*)(ptarg(j),j=1,np)
      close(30)

      obj = 0.0d0
      do i=nsw1,nsw2
         p1 = edge(1,i)
         p2 = edge(2,i)

         dx = coord(1,p2) - coord(1,p1)
         dy = coord(2,p2) - coord(2,p1)
         ds = sqrt(dx**2 + dy**2)

         f1 = 0.5d0*(qv(4,p1) - ptarg(p1))**2 
         f2 = 0.5d0*(qv(4,p2) - ptarg(p2))**2 
         obj= obj + 0.5d0*(f1 + f2)*ds
      enddo

      print*,'Pressure matching function =',obj

      end
