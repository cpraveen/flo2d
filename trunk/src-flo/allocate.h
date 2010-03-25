      allocate( elem(3,nt) )
      allocate( edge(2,nemax) )
      allocate( tedge(2,nemax) )
      allocate( vedge(2,nemax) )
      allocate( esue(3,nt) )
      allocate( ptype(np) )
      allocate( coord(2,np) )
      allocate( qc(nvar,nt) )
      allocate( qcd(nvar,nt) )
      allocate( qcold(nvar,nt) )
      allocate( dt(nt) )
      allocate( af(3,np) )
      allocate( qv(nvar,np) )
      allocate( qvd(nvar,np) )
      allocate( tarea(nt) )
      allocate( varea(np) )
      allocate( drmin(nt) )
      allocate( res(nvar,nt) )
      allocate( resd(nvar,nt) )
      allocate( qx(3,np) )
      allocate( qy(3,np) )
      allocate( qxd(3,np) )
      allocate( qyd(3,np) )


!     Variables for gmres, see gmres.h
      if(timemode .eq. 3)then
         ! Sizes of these arrays not clear, following seems to work
         allocate( a_jlu (2*4*nvar*nvar*nt) )
         allocate( a_alu (2*4*nvar*nvar*nt) )
         allocate( a_levs(2*4*nvar*nvar*nt) )

         !Sizes of these arrays of perfectly correct
         allocate( a_jmat(4*nvar*nvar*nt) )
         allocate( a_ia  (nvar*nt+1) )
         allocate( a_ja  (4*nvar*nvar*nt) )
         allocate( a_ju  (nvar*nt+1) )
         allocate( a_jw  (3*nvar*nt) )
         allocate( a_w   (nvar*nt) )
         
         ! Set pointers
         jlu  => a_jlu
         alu  => a_alu
         levs => a_levs
         jmat => a_jmat
         ia   => a_ia
         ja   => a_ja
         ju   => a_ju
         jw   => a_jw
         w    => a_w
      endif
