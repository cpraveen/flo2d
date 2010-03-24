!----------------------------------------------------------------------------
! Add elements of matrix for computing averaging weights
!----------------------------------------------------------------------------
      subroutine afact1(x1, x2, x3, sax1, sax2, sax3, say1, say2, &
                        say3, sax21, sax22, sax23, say21, say22,  &
                        say23, saxy1, saxy2, saxy3)
      implicit none
      include 'common.h'
      real(dp) :: x1(2), x2(2), x3(2), sax1, sax2, sax3, say1, &
                  say2, say3, sax21, sax22, sax23, say21, say22,  &
                  say23, saxy1, saxy2, saxy3

      real(dp) :: xt, yt, dx1, dy1, dx2, dy2, dx3, dy3

      xt    = (x1(1) + x2(1) + x3(1))/3.0d0
      yt    = (x1(2) + x2(2) + x3(2))/3.0d0

      dx1   = xt    - x1(1)
      dy1   = yt    - x1(2)
      sax1  = sax1  + dx1
      say1  = say1  + dy1
      sax21 = sax21 + dx1**2
      say21 = say21 + dy1**2
      saxy1 = saxy1 + dx1*dy1

      dx2   = xt    - x2(1)
      dy2   = yt    - x2(2)
      sax2  = sax2  + dx2
      say2  = say2  + dy2
      sax22 = sax22 + dx2**2
      say22 = say22 + dy2**2
      saxy2 = saxy2 + dx2*dy2

      dx3   = xt    - x3(1)
      dy3   = yt    - x3(2)
      sax3  = sax3  + dx3
      say3  = say3  + dy3
      sax23 = sax23 + dx3**2
      say23 = say23 + dy3**2
      saxy3 = saxy3 + dx3*dy3

      end
!----------------------------------------------------------------------------
! Add elements of matrix for computing averaging weights with boundary
! correction
!----------------------------------------------------------------------------
      subroutine afact2(x1, x2, x3, xp, nx, ny, sax, say, sax2, &
                        say2, saxy)
      implicit none
      include 'common.h'
      real(dp) :: x1(2), x2(2), x3(2), xp(2), nx, ny, sax, say, &
                  sax2, say2, saxy

      real(dp) :: xt, yt, xn, xg, yg, dx1, dy1

      xt   = (x1(1) + x2(1) + x3(1))/3.0d0
      yt   = (x1(2) + x2(2) + x3(2))/3.0d0

      xn   = (xt - xp(1))*nx + (yt - xp(2))*ny
      xg   = xt - 2.0d0*xn*nx
      yg   = yt - 2.0d0*xn*ny

      dx1  = xg - xp(1)
      dy1  = yg - xp(2)

      sax  = sax  + dx1
      say  = say  + dy1
      sax2 = sax2 + dx1**2
      say2 = say2 + dy1**2
      saxy = saxy + dx1*dy1

      end
!----------------------------------------------------------------------------
! Computes the averaging weights
! afact(1) and afact(2) must be set to zero before calling this
! function. If det==0 then afact(1)=afact(2)=0. This could happen if
! there is a "corner triangle".
!----------------------------------------------------------------------------
      subroutine afact3(sax, say, sax2, say2, saxy, afact)
      implicit none
      include 'common.h'
      real(dp) :: sax, say, sax2, say2, saxy, afact(3)
      real(dp) :: det

      det      = sax2*say2 - saxy**2
      if(det .ne. 0.0d0)then
         afact(1) = afact(1) + (say2*sax - saxy*say)/det
         afact(2) = afact(2) + (sax2*say - saxy*sax)/det
      endif

      end
!----------------------------------------------------------------------------
! Computes denominator in vertex averaging formula
!----------------------------------------------------------------------------
      subroutine afact4(x1, x2, x3, af1, af2, af3)
      implicit none
      include 'common.h'
      real(dp) :: x1(2), x2(2), x3(2), af1(3), af2(3), af3(3)

      real(dp) :: xt, yt, dx1, dx2, dx3, dy1, dy2, dy3, w1, w2, w3

      xt = ( x1(1) + x2(1) + x3(1) )/3.0d0
      yt = ( x1(2) + x2(2) + x3(2) )/3.0d0

      dx1= xt - x1(1)
      dy1= yt - x1(2)
      w1 = 1.0d0 - (af1(1)*dx1 + af1(2)*dy1)

      dx2= xt - x2(1)
      dy2= yt - x2(2)
      w2 = 1.0d0 - (af2(1)*dx2 + af2(2)*dy2)

      dx3= xt - x3(1)
      dy3= yt - x3(2)
      w3 = 1.0d0 - (af3(1)*dx3 + af3(2)*dy3)

      af1(3) = af1(3) + w1
      af2(3) = af2(3) + w2
      af3(3) = af3(3) + w3

      end
!----------------------------------------------------------------------------
! Check min and max range of weights
!----------------------------------------------------------------------------
      subroutine checkweights(elem, coord, afact)
      implicit none
      include 'common.h'
      include 'size.h'
      integer  :: elem(3,*)
      real(dp) :: coord(2,*), afact(3,*)

      integer  :: i, v1, v2, v3
      real(dp) :: dx1, dx2, dx3, dy1, dy2, dy3, &
                  xt, yt, w1, w2, w3, wmin, wmax

      wmin = 1.0d20
      wmax =-1.0d20
      open(20, file='wt.dat')
      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)

         xt = (coord(1,v1) + coord(1,v2) + coord(1,v3))/3.0d0
         yt = (coord(2,v1) + coord(2,v2) + coord(2,v3))/3.0d0

         dx1      = xt    - coord(1,v1)
         dy1      = yt    - coord(2,v1)
         w1       = 1.0d0 - (afact(1,v1)*dx1 + afact(2,v1)*dy1)

         dx2      = xt    - coord(1,v2)
         dy2      = yt    - coord(2,v2)
         w2       = 1.0d0 - (afact(1,v2)*dx2 + afact(2,v2)*dy2)

         dx3      = xt    - coord(1,v3)
         dy3      = yt    - coord(2,v3)
         w3       = 1.0d0 - (afact(1,v3)*dx3 + afact(2,v3)*dy3)

         w1       = w1/afact(3,v1)
         w2       = w2/afact(3,v2)
         w3       = w3/afact(3,v3)

         wmin     = dmin1(wmin, w1)
         wmin     = dmin1(wmin, w2)
         wmin     = dmin1(wmin, w3)

         wmax     = dmax1(wmax, w1)
         wmax     = dmax1(wmax, w2)
         wmax     = dmax1(wmax, w3)

         if(w1 .le. 0.0d0) then
            write(20,*) coord(1,v1), coord(2,v1)
            write(20,*) xt, yt
            write(20,*)
         endif

         if(w2 .le. 0.0d0) then
            write(20,*) coord(1,v2), coord(2,v2)
            write(20,*) xt, yt
            write(20,*)
         endif

         if(w3 .le. 0.0d0) then
            write(20,*) coord(1,v3), coord(2,v3)
            write(20,*) xt, yt
            write(20,*)
         endif

      enddo
      close(20)

      write(*,'(2x, "Minimum weight           =", f8.4)') wmin
      write(*,'(2x, "Maximum weight           =", f8.4)') wmax

      end
