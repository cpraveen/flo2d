!-----------------------------------------------------------------------------
      subroutine geometric(elem, edge, tedge, esue, vedge, spts, &
                           ptype, coord, drmin, tarea, varea, af)
      implicit none
      include 'param.h'

      integer  :: ptype(*), elem(3,*), edge(2,*), tedge(2,*), spts(*), &
                  vedge(2,*), esue(3,*)
      real(dp) :: coord(2,*), tarea(*), varea(*), drmin(*), af(3,*)

!     local variables
      integer  :: esup1(mesup*np), esup2(np+1), psup1(mpsup*np), &
                  psup2(np+1), betype(nbe), bdedge(2,nbe)

!     Read grid from file
      call read_grid(coord, elem, betype, bdedge)

!     Set point types
      call set_ptype(ptype, spts, betype, bdedge)

!     Make sure elements are oriented ccw
      call tri_orient(elem, coord)

!     Find elements surrounding a point
      call el_surr_point(elem, esup1, esup2)

!     Find points surrounding a point
      call pt_surr_pt(esup1, esup2, elem, psup1, psup2)

!     Create edges
      call create_edge(psup1, psup2, edge)

!     Find element adjoining each edge
      call el_surr_edge(esup1, esup2, elem, edge, tedge, vedge)

!     Put edge list in a particular order
      call order_edges(edge, tedge, vedge, betype, bdedge)

!     Element surrounding element
      call el_surr_el(elem, edge, tedge, esue)

!     Smooth the grid using Laplacian smoothing
!     call smooth(ptype, psup1, psup2, elem, esup1, esup2, coord)

!     Renumber triangles
      if(timemode .ne. 1)then
         call renumber(elem, esue, tedge)
      endif

!     Calculate triangle areas
      call tri_area(coord, elem, tarea, varea)

!     Length scale for time-step calculation
      call dtlength(coord, tarea, elem, drmin)

!     Write grid in gnuplot format for visualization
      call write_grid(coord, edge, tedge)

!     Generate .gnu files for visualization
      call prep_gnuplot

!     Area averaging factors
      call avgfact(ptype, elem, edge, coord, tarea, af)

      end

!-----------------------------------------------------------------------------
! Rotate grid by aoa
!-----------------------------------------------------------------------------
      subroutine rotate_grid(x)
      implicit none
      include 'param.h'
      real(dp) :: x(2,*)

      integer  :: i
      real(dp) :: x1, y1

      print*,'========Translating and rotating the grid================'

! Translate grid to ref point of vortex model
      do i=1,np
         x(1,i) = x(1,i) - xref
         x(2,i) = x(2,i) - yref
      enddo

! Rotate grid
      aoa = aoa_deg*M_PI/180.0d0
      do i=1,np
         x1 = x(1,i)
         y1 = x(2,i)
         x(1,i) = cos(aoa)*x1 + sin(aoa)*y1
         x(2,i) = cos(aoa)*y1 - sin(aoa)*x1
      enddo

      aoa     = 0.0d0
      aoa_deg = 0.0d0
      xref    = 0.0d0
      yref    = 0.0d0

      end
!-----------------------------------------------------------------------------
! Generate gnuplot script file for plotting grid
!-----------------------------------------------------------------------------
      subroutine prep_gnuplot
      implicit none
      include 'param.h'
      integer gnu

      gnu = 10
      open(unit=gnu, file='grid.gnu')
      write(gnu,*)"set xrange[",xmin,":",xmax,"]"
      write(gnu,*)"set yrange[",ymin,":",ymax,"]"
      write(gnu,*)"set size ratio -1"
      write(gnu,*)"set nokey"
      write(gnu,*)"p 'BD.DAT' w l"
      write(gnu,*)"pause 5"
      write(gnu,*)"p 'GRID.DAT' w l"
      write(gnu,*)"pause 5"
      close(gnu)

      end
!-----------------------------------------------------------------------------
!.....Check whether ordering of triangle is counter-clockwise
!.....Otherwise correct it
!-----------------------------------------------------------------------------
      subroutine tri_orient(elem, coord)
      implicit none
      include 'param.h'
      integer  :: elem(3,*)
      real(dp) :: coord(2,*)

      integer  :: cw, ccw, tmp, i, p1, p2, p3
      real(dp) :: dx1, dy1, dx2, dy2, cross

      cw = 0
      ccw= 0

      do i=1,nt
         p1    = elem(1,i)
         p2    = elem(2,i)
         p3    = elem(3,i)

         dx1   = coord(1,p2) - coord(1,p1)
         dy1   = coord(2,p2) - coord(2,p1)

         dx2   = coord(1,p3) - coord(1,p2)
         dy2   = coord(2,p3) - coord(2,p2)

         cross = dx1*dy2 - dx2*dy1

         if(cross .eq. 0.0d0)then
            print*,'Fatal: triangle',i,' is degenerate'
            stop
         endif

         if(cross .lt. 0.0d0)then
            cw        = cw + 1
            tmp       = elem(2,i)
            elem(2,i) = elem(3,i)
            elem(3,i) = tmp
         else
            ccw       = ccw + 1
         endif
      enddo

      write(*,'(" No of cw  triangles       =", i8)') cw
      write(*,'(" No of cww triangles       =", i8)') ccw

      end
!-----------------------------------------------------------------------------
!.....Calculate element and control volume areas for median cell
!-----------------------------------------------------------------------------
      subroutine tri_area(coord, elem, tarea, varea)
      implicit none
      include 'param.h'
      integer  :: elem(3,*)
      real(dp) :: coord(2,*), tarea(*), varea(*)

      real(dp) :: dx1, dy1, dx2, dy2
      integer  :: i, n1, n2, n3

      write(*,'(" Finding triangle areas")')

      do i=1,np
         varea(i) = 0.0d0
      enddo

      maxelarea = 0.0d0
      minelarea = 1.0d8
      do i=1,nt
         n1 = elem(1,i)
         n2 = elem(2,i)
         n3 = elem(3,i)

! Triangle area
         dx1= coord(1,n2) - coord(1,n1)
         dy1= coord(2,n2) - coord(2,n1)

         dx2= coord(1,n3) - coord(1,n1)
         dy2= coord(2,n3) - coord(2,n1)

         tarea(i) = 0.5d0*( dx1*dy2 - dx2*dy1 )
         maxelarea = max(maxelarea, tarea(i))
         minelarea = min(minelarea, tarea(i))

         varea(n1) = varea(n1) + tarea(i)
         varea(n2) = varea(n2) + tarea(i)
         varea(n3) = varea(n3) + tarea(i)

      enddo

      do i=1,np
         varea(i) = varea(i)*4.0d0/9.0d0
      enddo

      write(*,'(2x,"Minimum triangle area    =",e12.4)')minelarea
      write(*,'(2x,"Maximum triangle area    =",e12.4)')maxelarea

      end
!-----------------------------------------------------------------------------
!.....Finds elements surrounding a point
!.....Taken from Lohner
!.....esup1 stores the elements
!.....ordering is such that the elements surrounding point ipoin are stored in
!.....locations esup2(ipoin)+1 to esup2(ipoin+1)
!-----------------------------------------------------------------------------
      subroutine el_surr_point(elem, esup1, esup2)
      implicit none
      include 'param.h'
      integer :: esup1(*), esup2(*), elem(3,*)

      integer :: i, ie, inode, ipoi1, ipoin, istor

      do i=1,np+1
            esup2(i) = 0
      enddo

      do ie=1,nt
         do inode=1,3
            ipoi1        = elem(inode, ie) + 1
            esup2(ipoi1) = esup2(ipoi1) + 1
         enddo
      enddo

      do ipoin=2, np+1
         esup2(ipoin) = esup2(ipoin) + esup2(ipoin-1)
      enddo

      do ie=1, nt
         do inode=1,3
            ipoin        = elem(inode, ie)
            istor        = esup2(ipoin) + 1
            esup2(ipoin) = istor
            esup1(istor) = ie
         enddo
      enddo

      do ipoin=np+1, 2, -1
         esup2(ipoin) = esup2(ipoin-1)
      enddo

      esup2(1) = 0

      end
!-----------------------------------------------------------------------------
!.....Finds points surrounding a point
!.....Taken from Lohner
!.....psup1 contains the points
!.....Neighbours of ipoin are between psup2(ipoin)+1 to psup2(ipoin+1)
!-----------------------------------------------------------------------------
      subroutine pt_surr_pt(esup1, esup2, elem, psup1, psup2)
      implicit none
      include 'param.h'
      integer :: esup1(*), esup2(*)
      integer :: psup1(*), psup2(*), elem(3,*)

      integer :: ipoin, jpoin, inode, istor, ie, iesup, lpoin(np)

      do ipoin=1,np
         lpoin(ipoin) = 0
      enddo

      psup2(1) = 0
      istor     = 0

      do ipoin=1,np
         do iesup=esup2(ipoin)+1, esup2(ipoin+1)
            ie = esup1(iesup)
            do inode=1,3
               jpoin = elem(inode, ie)
               if(jpoin.ne.ipoin .and. lpoin(jpoin).ne.ipoin) then
                  istor = istor + 1
                  psup1(istor) = jpoin
                  lpoin(jpoin) = ipoin
               endif
            enddo
         enddo
         psup2(ipoin+1) = istor
      enddo

      end
!-----------------------------------------------------------------------------
! Create edges of triangles
!-----------------------------------------------------------------------------
      subroutine create_edge(psup1, psup2, edge)
      implicit none
      include 'param.h'
      integer :: psup1(*), psup2(*), edge(2,*)

      integer :: ipoin, ip, neigh

      ne = 0

      do ipoin = 1,np
         do ip=psup2(ipoin)+1, psup2(ipoin+1)
            neigh = psup1(ip)
            if(neigh.gt.ipoin) then
               ne          = ne + 1
               if(ne.gt.nemax) stop "create_edge: increase nemax"
               edge(1,ne) = ipoin
               edge(2,ne) = neigh
            endif
         enddo
      enddo

      write(*,'(" Number of edges           =", i8)') ne

      end
!-----------------------------------------------------------------------------
! Sort edges based on their type
!-----------------------------------------------------------------------------
      subroutine order_edges(edge, tedge, vedge, betype, bdedge)
      implicit none
      include 'param.h'
      integer :: edge(2,*), tedge(2,*), vedge(2,*), betype(*), &
                 bdedge(2,*)

      integer :: swedge(6,nbe), inedge(6,ne), ffedge(6,nbe)
      integer :: n1, n2, i, j, ibe

      write(*,'(" Sorting edges based on their type ...")')

      nsw = 0
      nin = 0
      nff = 0
      nif = 0
      nof = 0

! Put edges into different arrays based on type
      do i=1,ne
         n1 = edge(1,i)
         n2 = edge(2,i)

         if(tedge(2,i).gt.0)then ! It is an interior edge
            nin = nin + 1
            inedge(1,nin) = edge(1,i)
            inedge(2,nin) = edge(2,i)
            inedge(3,nin) = tedge(1,i)
            inedge(4,nin) = tedge(2,i)
            inedge(5,nin) = vedge(1,i)
            inedge(6,nin) = vedge(2,i)
         else ! It is a boundary edge
            ibe = 0
            j   = 0
            do while(ibe.eq.0 .and. j.le.nbe)
               j = j + 1
               if(n1.eq.bdedge(1,j) .and. n2.eq.bdedge(2,j))then
                  ibe = j
               endif
            enddo
            if(ibe.eq.0)then
               print*,'Fatal error: could not locate boundary edge'
               print*,'Edge number =',i
               print*,'n1, n2      =',n1,n2
               stop
            endif
            if(betype(j).eq.solid)then
               nsw = nsw + 1
               swedge(1,nsw) = edge(1,i)
               swedge(2,nsw) = edge(2,i)
               swedge(3,nsw) = tedge(1,i)
               swedge(4,nsw) = -solid
               swedge(5,nsw) = vedge(1,i)
               swedge(6,nsw) = vedge(2,i)
            elseif(betype(j).eq.farfield)then
               nff = nff + 1
               ffedge(1,nff) = edge(1,i)
               ffedge(2,nff) = edge(2,i)
               ffedge(3,nff) = tedge(1,i)
               ffedge(4,nff) = -farfield
               ffedge(5,nff) = vedge(1,i)
               ffedge(6,nff) = vedge(2,i)
            else
               print*,'Fatal error: unknown edge type =',betype(j)
               print*,'Boundary edge number =',j
            endif
         endif

      enddo

! Now put them back into the common array "edge"
! Order is important
      ne = 0

      nin1 = ne + 1
      do i=1,nin
         ne = ne + 1
         edge(1,ne) = inedge(1,i)
         edge(2,ne) = inedge(2,i)
         tedge(1,ne)= inedge(3,i)
         tedge(2,ne)= inedge(4,i)
         vedge(1,ne)= inedge(5,i)
         vedge(2,ne)= inedge(6,i)
      enddo
      nin2 = ne

      nsw1 = ne + 1
      do i=1,nsw
         ne = ne + 1
         edge(1,ne) = swedge(1,i)
         edge(2,ne) = swedge(2,i)
         tedge(1,ne)= swedge(3,i)
         tedge(2,ne)= swedge(4,i)
         vedge(1,ne)= swedge(5,i)
         vedge(2,ne)= swedge(6,i)
      enddo
      nsw2 = ne

      nff1 = ne + 1
      do i=1,nff
         ne = ne + 1
         edge(1,ne) = ffedge(1,i)
         edge(2,ne) = ffedge(2,i)
         tedge(1,ne)= ffedge(3,i)
         tedge(2,ne)= ffedge(4,i)
         vedge(1,ne)= ffedge(5,i)
         vedge(2,ne)= ffedge(6,i)
      enddo
      nff2 = ne

      if(nin+nsw+nff .ne. ne)then
         print*,'order_edges: Error -> number of edges dont add up'
         stop
      endif

      write(*,'(2x, "Number of interior edges =", i8)')nin
      write(*,'(2x, "Number of solid    edges =", i8)')nsw
      write(*,'(2x, "Number of farfield edges =", i8)')nff
      write(*,'(2x, "Total no of edges        =", i8)')ne

      end
!-----------------------------------------------------------------------------
!.....For each edge find the elements to its right and left
!-----------------------------------------------------------------------------
      subroutine el_surr_edge(esup1, esup2, elem, edge, tedge, vedge)
      implicit none
      include 'param.h'
      integer :: esup1(*), esup2(*), tedge(2,*), &
                 vedge(2,*), elem(3,*), edge(2,*)

      integer :: i, jt, n1, n2, el, tmp

      do i=1,ne
            tedge(1,i) = 0
            tedge(2,i) = 0
            n1 = edge(1,i)
            n2 = edge(2,i)
            do jt=esup2(n1)+1, esup2(n1+1)
                  el = esup1(jt)
                  if( (n1.eq.elem(1,el) .and. n2.eq.elem(2,el)) .or. &
                      (n1.eq.elem(2,el) .and. n2.eq.elem(3,el)) .or. &
                      (n1.eq.elem(3,el) .and. n2.eq.elem(1,el)) ) &
                        tedge(1,i) = el
                  if( (n2.eq.elem(1,el) .and. n1.eq.elem(2,el)) .or. &
                      (n2.eq.elem(2,el) .and. n1.eq.elem(3,el)) .or. &
                      (n2.eq.elem(3,el) .and. n1.eq.elem(1,el)) ) &
                        tedge(2,i) = el
            enddo

            if(tedge(1,i) .eq. 0)then
                  tedge(1,i) = tedge(2,i)
                  tedge(2,i) = 0
                  tmp          = edge(1,i)
                  edge(1,i)    = edge(2,i)
                  edge(2,i)    = tmp
                  if(tedge(1,i) .eq. 0)then
                     print*,'Fatal: No edge neighbour'
                     stop
                  endif
            endif
      enddo

! For each edge, find vertex of adjacent triangle which is opposite to
! the edge. This is used for limited reconstruction of inviscid fluxes
      do i=1,ne
         vedge(1,i) = 0
         vedge(2,i) = 0
         n1 = edge(1,i)
         n2 = edge(2,i)
         if(tedge(1,i) .ne. 0)then
            el = tedge(1,i)
            if(elem(1,el) .ne. n1 .and. &
               elem(1,el) .ne. n2) vedge(1,i) = elem(1,el)
            if(elem(2,el) .ne. n1 .and. &
               elem(2,el) .ne. n2) vedge(1,i) = elem(2,el)
            if(elem(3,el) .ne. n1 .and. &
               elem(3,el) .ne. n2) vedge(1,i) = elem(3,el)
         endif
         
         if(tedge(2,i) .ne. 0)then
            el = tedge(2,i)
            if(elem(1,el) .ne. n1 .and. &
               elem(1,el) .ne. n2) vedge(2,i) = elem(1,el)
            if(elem(2,el) .ne. n1 .and. &
               elem(2,el) .ne. n2) vedge(2,i) = elem(2,el)
            if(elem(3,el) .ne. n1 .and. &
               elem(3,el) .ne. n2) vedge(2,i) = elem(3,el)
         endif
      enddo

      end
!-----------------------------------------------------------------------------
! Find triangles surrounding a triangle, ie, face neighbours
!-----------------------------------------------------------------------------
      subroutine el_surr_el(elem, edge, tedge, esue)
      implicit none
      include 'param.h'
      integer :: elem(3,*), edge(2,*), tedge(2,*), esue(3,*)

      integer :: it, ie, e1, e2, c1, c2, p1, p2, p3, q1, q2, q3

      do it=1,nt
         esue(1,it) = 0
         esue(2,it) = 0
         esue(3,it) = 0
      enddo

      do ie=1,ne
         e1 = edge(1,ie)
         e2 = edge(2,ie)

         c1 = tedge(1,ie)
         c2 = tedge(2,ie)

         p1 = elem(1,c1)
         p2 = elem(2,c1)
         p3 = elem(3,c1)

         ! c2 is a neighbour of c1
         ! For boundary edge, c2 = -betype
         if    (p1 .eq. e1)then
            esue(1,c1) = c2
         elseif(p2 .eq. e1)then
            esue(2,c1) = c2
         elseif(p3 .eq. e1)then
            esue(3,c1) = c2
         else
            print*,'esue: Fatal error'
            stop
         endif

         if(c2 .gt. 0)then

            q1 = elem(1,c2)
            q2 = elem(2,c2)
            q3 = elem(3,c2)
            if    (q1 .eq. e2)then
               esue(1,c2) = c1
            elseif(q2 .eq. e2)then
               esue(2,c2) = c1
            elseif(q3 .eq. e2)then
               esue(3,c2) = c1
            else
               print*,'esue: Fatal error'
               stop
            endif

         endif

      enddo

      end

!-----------------------------------------------------------------------------
! Renumber cells for lusgs
!-----------------------------------------------------------------------------
      subroutine renumber(elem, esue, tedge)
      implicit none
      include 'size.h'
      integer :: elem(3,*), esue(3,*), tedge(2,*)

      integer :: oldnum(nt), newnum(nt), telem(3,nt), tesue(3,nt)
      integer :: i, j, it, t1, t2, t3, tcount

      print*,'Renumbering cells using Cuthill-McKee...'

      do i=1,nt
         oldnum(i) = 0
         newnum(i) = 0
      enddo

      tcount    = 1
      oldnum(1) = 1
      newnum(1) = 1

      do i=1,nt
         it = oldnum(i)
         if(it .eq. 0)then
            print*,'renumber: Fatal error. it is zero for i=',i
            stop
         endif

         t1 = esue(1,it)
         t2 = esue(2,it)
         t3 = esue(3,it)

         if(t1 .gt. 0 .and. newnum(t1) .eq. 0)then
            tcount         = tcount + 1
            oldnum(tcount) = t1
            newnum(t1)     = tcount
         endif

         if(t2 .gt. 0 .and. newnum(t2).eq.0)then
            tcount         = tcount + 1
            oldnum(tcount) = t2
            newnum(t2)     = tcount
         endif

         if(t3 .gt. 0 .and. newnum(t3) .eq. 0)then
            tcount         = tcount + 1
            oldnum(tcount) = t3
            newnum(t3)     = tcount
         endif
      enddo

      if(nt .ne. tcount)then
         print*,'renumber: tcount does not match nt.'
         print*,'          Possible bug'
         stop
      endif

! Now rearrange some other element data structure
      do i=1,nt
         telem(1,i) = elem(1,i)
         telem(2,i) = elem(2,i)
         telem(3,i) = elem(3,i)
         tesue(1,i) = esue(1,i)
         tesue(2,i) = esue(2,i)
         tesue(3,i) = esue(3,i)
      enddo

      do i=1,nt
         it        = oldnum(i)
         elem(1,i) = telem(1,it)
         elem(2,i) = telem(2,it)
         elem(3,i) = telem(3,it)
         do j=1,3
            t1 = tesue(j,it)
            if(t1 .gt. 0)then
               esue(j,i) = newnum(t1)
            else
               esue(j,i) = t1
            endif
         enddo
      enddo

      do i=1,ne
         t1 = tedge(1,i)
         t2 = tedge(2,i)
         if(t1 .gt. 0)then
            tedge(1,i) = newnum(t1)
         else
            tedge(1,i) = t1
         endif
         if(t2 .gt. 0)then
            tedge(2,i) = newnum(t2)
         else
            tedge(2,i) = t2
         endif
      enddo

      end
!-----------------------------------------------------------------------------
!.....Calculate length used for time step
!.....For each triangle find the minimum altitude, or
!.....area/perimeter
!-----------------------------------------------------------------------------
      subroutine dtlength(coord, tarea, elem, drmin)
      implicit none
      include 'param.h'
      integer  :: elem(3,*)
      real(dp) :: coord(2,*), tarea(*), drmin(*)

      integer  :: it, n1, n2, n3
      real(dp) :: dx1, dx2, dx3, dy1, dy2, dy3, dr1, dr2, dr3, &
                  h1, h2, h3, perim

      do it=1,nt
         n1  = elem(1,it)
         n2  = elem(2,it)
         n3  = elem(3,it)

         dx1 = coord(1,n2) - coord(1,n3)
         dy1 = coord(2,n2) - coord(2,n3)
         dr1 = sqrt(dx1**2 + dy1**2)
         h1  = 2.0d0*tarea(it)/dr1

         dx2 = coord(1,n3) - coord(1,n1)
         dy2 = coord(2,n3) - coord(2,n1)
         dr2 = sqrt(dx2**2 + dy2**2)
         h2  = 2.0d0*tarea(it)/dr2

         dx3 = coord(1,n1) - coord(1,n2)
         dy3 = coord(2,n1) - coord(2,n2)
         dr3 = sqrt(dx3**2 + dy3**2)
         h3  = 2.0d0*tarea(it)/dr3

         perim     = dr1 + dr2 + dr3
         drmin(it) = tarea(it)/perim
      enddo

      end
!-----------------------------------------------------------------------------
! Write grid into a file for visualization with gnuplot
!-----------------------------------------------------------------------------
      subroutine write_grid(coord, edge, tedge)
      implicit none
      include 'param.h'
      real(dp) :: coord(2,*)
      integer  :: edge(2,*), tedge(2,*)

      integer  :: gfile, i, n1, n2

!     Write boundary edges to BD.DAT
      open(unit=10, file='BD.DAT')
      do i=1,ne
      n1 = edge(1,i)
      n2 = edge(2,i)
      if( tedge(1,i)*tedge(2,i) .le. 0)then
         write(10,*)coord(1,n1), coord(2,n1)
         write(10,*)coord(1,n2), coord(2,n2)
         write(10,*)
      endif
      enddo
      close(10)

!     Write grid into file GRID.DAT
      gfile=15
      open(unit=gfile, file='GRID.DAT')

      do i=1,ne
         n1 = edge(1,i)
         n2 = edge(2,i)
         write(gfile,*) coord(1,n1), coord(2,n1)
         write(gfile,*) coord(1,n2), coord(2,n2)
         write(gfile,*)
      enddo

      close(gfile)

!     call system('gnuplot -noraise grid.gnu &')

      end
!----------------------------------------------------------------------
!     Set point type, whether point is interior or boundary
!----------------------------------------------------------------------
      subroutine set_ptype(ptype, spts, betype, bdedge)
      implicit none
      include "param.h"
      integer :: ptype(*), spts(*), betype(*), bdedge(2,*)

      integer :: i, p1, p2

      do i=1,np
         ptype(i) = interior
      enddo

!     Set all boundary points to type "boundary"
      do i=1,nbe
         p1 = bdedge(1,i)
         p2 = bdedge(2,i)

         ptype(p1) = boundary
         ptype(p2) = boundary
      enddo

!     Now set all solid points to type "solid"
      do i=1,nbe
         p1 = bdedge(1,i)
         p2 = bdedge(2,i)

         if(betype(i).eq.solid)then
            ptype(p1) = solid
            ptype(p2) = solid
         endif
      enddo

!     Put all solid points into array spts
      nsp = 0
      do i=1,np
         if(ptype(i).eq.solid)then
            nsp = nsp + 1
            spts(nsp) = i
         endif
      enddo

!     Count number of interior points: just for checking
      nip = 0
      do i=1,np
         if(ptype(i).eq.interior) nip = nip + 1
      enddo
      nbp = np - nip

      write(*,'(" No of interior points     =", i8)') nip
      write(*,'(" No of boundary points     =", i8)') nbp
      write(*,'(" No of solid    points     =", i8)') nsp

      end
