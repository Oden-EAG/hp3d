c----------------------------------------------------------------------
c
c   routine name       - solelm
c
c----------------------------------------------------------------------
c
c   latest revision    - Apr 14
c
c   purpose            - routine calculates unconstrained solution dof
c                        for a 3D element
c
c   arguments :
c     in:
c           Mdle       - middle node of an element
c     out:
c           ZdofH      - the element unconstrained H1 dof
c           ZdofE      - the element unconstrained H(curl) dof
c           ZdofV      - the element unconstrained H(div) dof
c           ZdofQ      - the element unconstrained L2 dof
c
c----------------------------------------------------------------------
c
      subroutine solelm(Mdle, ZdofH,ZdofE,ZdofV,ZdofQ)
c
      use data_structure3D
#include "syscom.blk"
c
      dimension ZdofH(MAXEQNH,MAXbrickH)
      dimension ZdofE(MAXEQNE,MAXbrickE)
      dimension ZdofV(MAXEQNV,MAXbrickV)
      dimension ZdofQ(MAXEQNQ,MAXbrickQ)
c
c  ...element order of approximation
      dimension norder(19)
c
c  ...modified element nodes and corresponding number of dof
      dimension nodm  (MAXNODM),ndofmH(MAXNODM),
     .          ndofmE(MAXNODM),ndofmV(MAXNODM)
c
      dimension nrconH(MAXbrickH),nrconE(MAXbrickE),nrconV(MAXbrickV),
     .          nacH(NACDIM,MAXbrickH),constrH(NACDIM,MAXbrickH),
     .          nacE(NACDIM,MAXbrickE),constrE(NACDIM,MAXbrickE),
     .          nacV(NACDIM,MAXbrickV),constrV(NACDIM,MAXbrickV)
c
c  ...modified element dof
      dimension zvalH(MAXEQNH,2*MAXbrickH)
      dimension zvalE(MAXEQNE,2*MAXbrickE)
      dimension zvalV(MAXEQNV,2*MAXbrickV)
c
c---------------------------------------------------------------------
c
      iprint=0
c
c  ...determine element order of approximation
      call find_order(Mdle, norder)
c
      if (iprint.ge.1) then
        write(*,7000) Mdle
 7000   format('solelm: Mdle   = ',i6)
        call print_order(NODES(Mdle)%type,norder)
      endif
c
c  ...determine number of local dof
      call celndof(NODES(Mdle)%type,norder,
     .             nrdoflH,nrdoflE,nrdoflV,nrdoflQ)
c
c  ...determine constraints' coefficients
      call logic(Mdle,2,
     .           nodm,ndofmH,ndofmE,ndofmV,nrnodm,
     .           nrconH,nacH,constrH,
     .           nrconE,nacE,constrE,
     .           nrconV,nacV,constrV)
c
c---------------------------------------------------------------------
c
c  ...initiate dof's (use ZdofQ for simplicity)
      zvalH=ZERO ; zvalE=ZERO ; zvalV=ZERO ; ZdofQ=ZERO
c
c  ...initiate counter for some mystic reason
      kH=0; kE=0; kV=0; kQ=0
c
c  ...copy dof's into the local arrays
      do j=1,nrnodm
        call dof_out(nodm(j), kH,kE,kV,kQ,zvalH,zvalE,zvalV,zdofQ)
      enddo
c
c  ...initiate the dof
      ZdofH(1:MAXEQNH,1:nrdoflH) = ZERO
      ZdofE(1:MAXEQNE,1:nrdoflE) = ZERO
      ZdofV(1:MAXEQNV,1:nrdoflV) = ZERO
c
c
c-----------------------------------------------------------------------
c     H1 DOF'S
c-----------------------------------------------------------------------
c
c  ...loop through the local dof
      do k=1,nrdoflH
c
        if (iprint.eq.1) then
          write(*,*)'k = ',k
          write(*,7001) (nacH(l,k),l=1,nrconH(k))
 7001     format(10i10)
          write(*,7002) (constrH(l,k),l=1,nrconH(k))
 7002     format(10(2x,f8.3))
          do ivar=1,1
            write(*,7003) ivar,(zvalH(ivar,nacH(l,k)),l=1,nrconH(k))
#if C_MODE
 7003       format('ivar = ',i3,2x,10(2e12.5,2x))
#else
 7003       format('ivar = ',i3,2x,20e12.5)
#endif
          enddo
        endif
c
c  .....accumulate for the values
        do kp=1,nrconH(k)
          l = nacH(kp,k)
          do ivar=1,NRHVAR*NRCOMS
            ZdofH(ivar,k) = ZdofH(ivar,k) + constrH(kp,k)*zvalH(ivar,l)
          enddo
        enddo
        if (iprint.eq.1) then
ccc           do icomp=1,NRCOMS
           do icomp=1,1
             ibeg = (icomp-1)*NRHVAR+1
             iend = ibeg+NRHVAR-1
             write(*,7004) icomp,(ZdofH(ivar,k),ivar=ibeg,iend)
 7004        format('icomp = ',i2,' ZdofH = ',4(2e12.5,2x))
           enddo
ccc           call pause
        endif
c
c  ...loop through local dof
      enddo
c
c
c-----------------------------------------------------------------------
c     H(curl) DOF'S
c-----------------------------------------------------------------------
c
c  ...loop through the local dof
      do kE=1,nrdoflE
        if (iprint.eq.2) then
          write(*,*) 'kE =',kE
          write(*,7001) (nacE(l,kE),l=1,nrconE(kE))
          write(*,7002) (constrE(l,kE),l=1,nrconE(kE))
          do ivar=1,NREVAR*NRCOMS
            write(*,7002) (zvalE(ivar,nacE(l,kE)),l=1,nrconE(kE))
          enddo
        endif
c
c  .....accumulate for the values
        do kp=1,nrconE(kE)
          l = nacE(kp,kE)
          do ivar=1,NREVAR*NRCOMS
            ZdofE(ivar,kE) = ZdofE(ivar,kE)
     .                      + constrE(kp,kE)*zvalE(ivar,l)
          enddo
        enddo
        if (iprint.eq.2) then
           write(*,7005) kE,(ZdofE(ivar,kE),ivar=1,nvarE)
 7005      format(' solelm: kE,ZdofE = ',i2,8e12.5)
ccc           call pause
        endif
c
c  ...loop through local H(curl) dof
      enddo
c
c
c-----------------------------------------------------------------------
c     H(div) DOF'S
c-----------------------------------------------------------------------
c
c  ...loop through the local dof
      do kV=1,nrdoflV
        if (iprint.eq.1) then
          write(*,*) 'kV=',kV
          write(*,7001) (nacV(l,kV),l=1,nrconV(kV))
          write(*,7002) (constrV(l,kV),l=1,nrconV(kV))
          do ivar=1,NRVVAR*NRCOMS
            write(*,7002) (zvalV(ivar,nacV(l,kV)),l=1,nrconV(kV))
          enddo
        endif
c
c  .....accumulate for the values
        do kp=1,nrconV(kV)
          l = nacV(kp,kV)
          do ivar=1,NRVVAR*NRCOMS
            ZdofV(ivar,kV) = ZdofV(ivar,kV)
     .                      + constrV(kp,kV)*zvalV(ivar,l)
          enddo
        enddo
        if (iprint.eq.3) then
           write(*,7006) kV,ZdofV(1:NRVVAR*NRCOMS,kV)
 7006      format('solelm: kV,ZdofV = ',i2,8e12.5)
ccc           call pause
        endif
c
c  ...loop through local H(div) dof
      enddo
c
c-----------------------------------------------------------------------
c     L2 DOF'S , nothing to do...
c-----------------------------------------------------------------------
c
      if (iprint.eq.4) then
        do kQ=1,nrdoflQ
          nvarQ = NRQVAR*NRCOMS
          write(*,7007) kQ,ZdofQ(1:nvarQ,kQ)
 7007     format('solelm: ZdofQ(1:nvarQ,',i2,') = ',10e12.5)
        enddo
        call pause
      endif
c
      if (iprint.ne.0) call pause
c
c
      end subroutine solelm
