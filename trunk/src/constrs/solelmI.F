c----------------------------------------------------------------------
c
c   routine name       - solelmI
c
c----------------------------------------------------------------------
c
c   latest revision    - Aug 18
c
c   purpose            - routine calculates unconstrained INTERFACE
c                        solution dof for a 3D element
c
c   arguments :
c     in:
c           Mdle       - middle node of an element
c     out:
c           ZdofH      - the element unconstrained interface H1 dof
c           ZdofE      - the element unconstrained interface H(curl) dof
c           ZdofV      - the element unconstrained interface H(div) dof
c
c----------------------------------------------------------------------
c
      subroutine solelmI(Mdle, ZdofH,ZdofE,ZdofV)
c
      use data_structure3D
#include "syscom.blk"
c
      dimension ZdofH(MAXEQNH,MAXbrickH-MAXmdlbH)
      dimension ZdofE(MAXEQNE,MAXbrickE-MAXmdlbE)
      dimension ZdofV(MAXEQNV,MAXbrickV-MAXmdlbV)
c
c  ...element type
      integer :: ntype
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
c  ...determine element order of approximation
      call find_order(Mdle, norder)
      ntype = NODES(Mdle)%ntype
      select case(ntype)
        case(MDLB); norder(19) = 111
        case(MDLP); norder(15) = 11
        case(MDLN); norder(11) = 1
        case(MDLD); norder(14) = 1
      end select
c
c  ...determine number of local dof
      call celndof(ntype,norder,
     .             nrdoflH,nrdoflE,nrdoflV,nrdoflQ)
c
c  ...determine constraints' coefficients
      call logic(Mdle,2,
     .           nodm,ndofmH,ndofmE,ndofmV,nrnodm,
     .           nrconH,nacH,constrH,
     .           nrconE,nacE,constrE,
     .           nrconV,nacV,constrV)
c
c  ...eliminate the middle node from the list of modified element nodes
      nrnodm = nrnodm-1
c
c---------------------------------------------------------------------
c
c  ...initiate dof's
      zvalH=ZERO ; zvalE=ZERO ; zvalV=ZERO
c
c  ...initiate counters
      kH=0; kE=0; kV=0; kQ=0
c
c  ...copy dof's into the local arrays
      do j=1,nrnodm
        call dof_out(nodm(j), kH,kE,kV,kQ,zvalH,zvalE,zvalV,zvoid)
      enddo
      if (kQ.ne.0) then
        write(*,*) 'solelmI: INCONSISTENCY, kQ = ',kQ
        stop
      endif
c
c      write(*,*) 'solelmI:'
c      write(*,5010) kH,kE,kV
c 5010 format('kH = ',I2,', kE = ',I2,', kV = ',I2)
c
c      write(*,5020) nrdoflH,nrdoflE,nrdoflV
c 5020 format('nrdoflH = ',I2,', nrdoflE = ',I2,', nrdoflV = ',I2)
c
c      kH = MAXbrickH-MAXmdlbH
c      kE = MAXbrickE-MAXmdlbE
c      kV = MAXbrickV-MAXmdlbV
c      write(*,5030) kH, kE, kV
c 5030 format('MAXbrickH-MAXmdlbH = ',I3,',
c     .        MAXbrickE-MAXmdlbE = ',I3,', MAXbrickV-MAXmdlbV = ',I3)
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
      do kH=1,nrdoflH
c
c  .....accumulate for the values
        do kp=1,nrconH(kH)
          l = nacH(kp,kH)
          do ivar=1,NRHVAR*NRCOMS
            ZdofH(ivar,kH) = ZdofH(ivar,kH)
     .                     + constrH(kp,kH)*zvalH(ivar,l)
          enddo
        enddo
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
c
c  .....accumulate for the values
        do kp=1,nrconE(kE)
          l = nacE(kp,kE)
          do ivar=1,NREVAR*NRCOMS
            ZdofE(ivar,kE) = ZdofE(ivar,kE)
     .                     + constrE(kp,kE)*zvalE(ivar,l)
          enddo
        enddo
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
c
c  .....accumulate for the values
        do kp=1,nrconV(kV)
          l = nacV(kp,kV)
          do ivar=1,NRVVAR*NRCOMS
            ZdofV(ivar,kV) = ZdofV(ivar,kV)
     .                     + constrV(kp,kV)*zvalV(ivar,l)
          enddo
        enddo
c
c  ...loop through local H(div) dof
      enddo
c
c
      end subroutine solelmI
