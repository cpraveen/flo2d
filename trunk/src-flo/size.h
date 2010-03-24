! ntmax = number of triangles (elements in general)
! nemax = number of edges
      integer ntmax, nemax
      common/maxdims/ntmax,nemax

      integer nspmax, nfpmax, nopmax, nbpmax, nbemax
      parameter(nspmax= 2000, &
                nfpmax= 2000, &
                nopmax= 2000, &
                nbpmax= 2000, &
                nbemax= 2000)

! Actual number of points, elements, edges and boundary edges
!     np = number of points
!     nt = number of triangles/elements
!     ne = number of edges
!     nbe= number of boundary edges
!     nsp= number of solid boundary points
!     nfp= number of farfield boundary points
!     nop= number of outflow boundary points
!     nbp= number of boundary points, must equal nsp+nfp+nop
!     nip= number of interior points
      integer np,nt,ne,nip,nsp,nfp,nop,nbp,nbe
      common/dims/np,nt,ne,nip,nsp,nfp,nop,nbp,nbe

      integer nsw, nin, nff, nif, nof, &
              nsw1, nin1, nff1, nif1, nof1, &
              nsw2, nin2, nff2, nif2, nof2
      common/nedge/nsw, nin, nff, nif, nof, &
                   nsw1, nin1, nff1, nif1, nof1, &
                   nsw2, nin2, nff2, nif2, nof2

! Maximum elements surrounding a point
      integer mesup
      parameter(mesup=10)

! Maximum elements surrounding a point
      integer mesubp
      parameter(mesubp=8)

! Maximum points surrounding a point
      integer mpsup
      parameter(mpsup=10)
