      double precision dmat(nvar,nvar,npmax)

      integer          ja(4*nvar*nvar*ntmax), ia(nvar*ntmax+1)
      integer          jlu(4*nvar*nvar*ntmax), ju(nvar*ntmax+1)
      double precision jmat(4*nvar*nvar*ntmax), alu(4*nvar*nvar*ntmax)

      common/variables/dmat,ja,ia,jlu,ju,jmat,alu
