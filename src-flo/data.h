
      integer          ja(4*nvar*nvar*ntmax), ia(nvar*ntmax+1)
      integer          jlu(4*nvar*nvar*ntmax), ju(nvar*ntmax+1)
      integer          levs(nvar*ntmax)
      integer          jw(3*nvar*ntmax)
      integer          lfil
      parameter(lfil=2)
      double precision jmat(4*nvar*nvar*ntmax), alu(4*nvar*nvar*ntmax)
      double precision w(nvar*ntmax)

      common/variables/ja,ia,jlu,ju,levs,jw,jmat,alu,w
