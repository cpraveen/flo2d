      subroutine divergence(np, coord, nwp, rw, wt, nstep)
      implicit none
      integer np, nwp, nstep
      double precision    coord(2,*),
     1        rw (2,*),
     2        wt (2,*)

      integer i, j, k, nstep1, nstep2
      double precision dx1, dy1, dr1, lg
      double precision div(np), div_min, div_max
      double precision omg(np), omg_min, omg_max
      double precision shr(np), shr_min, shr_max

      div_min = +1.0d+20
      div_max = -1.0d-20
      omg_min = +1.0d+20
      omg_max = -1.0d-20
      shr_min = +1.0d+20
      shr_max = -1.0d-20
      do i=1,np
         div(i) = 0.0d0
         omg(i) = 0.0d0
         shr(i) = 0.0d0
         do k=1,nwp
            dx1 = coord(1,i) - rw(1,k)
            dy1 = coord(2,i) - rw(2,k)
            dr1 = dsqrt(dx1*dx1 + dy1*dy1)
            if(dr1.eq.0.0d0)then
               lg = 0.0d0
            else
               lg = dlog(dr1)
            endif
            div(i) = div(i) + wt(1,k)*dx1*(1.0d0 + 2.0d0*lg) +
     1                        wt(2,k)*dy1*(1.0d0 + 2.0d0*lg)
            omg(i) = omg(i) - wt(1,k)*dy1*(1.0d0 + 2.0d0*lg) +
     1                        wt(2,k)*dx1*(1.0d0 + 2.0d0*lg)
            shr(i) = shr(i) + wt(1,k)*dy1*(1.0d0 + 2.0d0*lg) +
     1                        wt(2,k)*dx1*(1.0d0 + 2.0d0*lg)
         enddo
         div(i)  = div(i) + wt(1,nwp+2) + wt(2,nwp+3)
         div_min = min(div_min, div(i))
         div_max = max(div_max, div(i))
         omg(i)  = omg(i) - wt(1,nwp+3) + wt(2,nwp+2)
         omg_min = min(omg_min, omg(i))
         omg_max = max(omg_max, omg(i))
         shr(i)  = shr(i) + wt(1,nwp+3) + wt(2,nwp+2)
         shr_min = min(shr_min, shr(i))
         shr_max = max(shr_max, shr(i))
      enddo

      print*,'Divergence min =',div_min
      print*,'Divergence max =',div_max
      print*,'Vorticity  min =',omg_min
      print*,'Vorticity  max =',omg_max
      print*,'Shear      min =',shr_min
      print*,'Shear      max =',shr_max

      nstep1 = max( abs(div_min), abs(div_max) ) / 0.1d0
      nstep1 = min(nstep1, 15)
      nstep1 = max(nstep1, 1)

      nstep2 = max( abs(omg_min), abs(omg_max) ) / 0.05d0
      nstep2 = min(nstep2, 15)
      nstep2 = max(nstep2, 1)

      nstep  = max(nstep1, nstep2)

      end
