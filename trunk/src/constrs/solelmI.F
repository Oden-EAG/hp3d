c----------------------------------------------------------------------
c
c   routine name       - solelmI
c
c----------------------------------------------------------------------
c
c   latest revision    - Sep 2023
c
c   purpose            - routine calculates unconstrained INTERFACE
c                        solution dof for a 3D element
c   remark: this routine must be OMP thread-safe
c
c   arguments :
c     in:
c           Mdle       - middle node of an element
c           Ncoms      - solution component set: 1,...,NRCOMS
c     out:
c           ZdofH      - the element unconstrained interface H1 dof
c           ZdofE      - the element unconstrained interface H(curl) dof
c           ZdofV      - the element unconstrained interface H(div) dof
c
c----------------------------------------------------------------------
#include "typedefs.h"
c
      subroutine solelmI_coms(Mdle,Ncoms, ZdofH,ZdofE,ZdofV)
c
      use data_structure3D
      implicit none
c
      integer, intent(in)  :: Mdle
      integer, intent(in)  :: Ncoms
      VTYPE  , intent(out) :: ZdofH(MAXEQNH,MAXbrickH-MAXmdlbH)
      VTYPE  , intent(out) :: ZdofE(MAXEQNE,MAXbrickE-MAXmdlbE)
      VTYPE  , intent(out) :: ZdofV(MAXEQNV,MAXbrickV-MAXmdlbV)
c
c  ...element type
      integer :: ntype
c
c  ...element order of approximation
      integer :: norder(19)
c
c  ...modified element nodes and corresponding number of dof
      integer :: nodm  (MAXNODM),ndofmH(MAXNODM),
     .           ndofmE(MAXNODM),ndofmV(MAXNODM)
c
      integer :: nrconH(MAXbrickH),nacH(NACDIM,MAXbrickH),
     .           nrconE(MAXbrickE),nacE(NACDIM,MAXbrickE),
     .           nrconV(MAXbrickV),nacV(NACDIM,MAXbrickV)
c
      real(8) :: constrH(NACDIM,MAXbrickH),
     .           constrE(NACDIM,MAXbrickE),
     .           constrV(NACDIM,MAXbrickV)
c
c  ...modified element dof
      VTYPE :: zvalH(MAXEQNH,2*MAXbrickH)
      VTYPE :: zvalE(MAXEQNE,2*MAXbrickE)
      VTYPE :: zvalV(MAXEQNV,2*MAXbrickV)
      VTYPE :: zvoid
c
      integer :: kH,kE,kV,kQ,ivar,j,kp,l
      integer :: nrnodm,nrdoflH,nrdoflE,nrdoflV,nrdoflQ
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
      ZdofH=ZERO ; ZdofE=ZERO ; ZdofV=ZERO
c
c  ...initiate counters (needed for dof_out)
      kH=0; kE=0; kV=0; kQ=0
c
c  ...copy dof's into the local arrays
      do j=1,nrnodm
        call dof_out(nodm(j),Ncoms, kH,kE,kV,kQ,zvalH,zvalE,zvalV,zvoid)
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
          do ivar=1,NRHVAR*NRRHS
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
          do ivar=1,NREVAR*NRRHS
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
          do ivar=1,NRVVAR*NRRHS
            ZdofV(ivar,kV) = ZdofV(ivar,kV)
     .                     + constrV(kp,kV)*zvalV(ivar,l)
          enddo
        enddo
c
c  ...loop through local H(div) dof
      enddo
c
c
      end subroutine solelmI_coms
c
c
c-----------------------------------------------------------------------
c  ...see routine solelmI_coms for description
c     last rev: Sep 2023
      subroutine solelmI(Mdle, ZdofH,ZdofE,ZdofV)
c
      use data_structure3D
      implicit none
c
      integer, intent(in)  :: Mdle
      VTYPE  , intent(out) :: ZdofH(MAXEQNH,MAXbrickH-MAXmdlbH)
      VTYPE  , intent(out) :: ZdofE(MAXEQNE,MAXbrickE-MAXmdlbE)
      VTYPE  , intent(out) :: ZdofV(MAXEQNV,MAXbrickV-MAXmdlbV)
c
c  ...by default, request solution component set N_COMS
      call solelmI_coms(Mdle,N_COMS, ZdofH,ZdofE,ZdofV)
c
      end subroutine solelmI
