!----------------------------------------------------------------------
!
!   routine name       - solelmI
!
!----------------------------------------------------------------------
!
!   latest revision    - Sep 2023
!
!   purpose            - routine calculates unconstrained INTERFACE
!                        solution dof for a 3D element
!   remark: this routine must be OMP thread-safe
!
!   arguments :
!     in:
!           Mdle       - middle node of an element
!           Ncoms      - solution component set: 1,...,NRCOMS
!     out:
!           ZdofH      - the element unconstrained interface H1 dof
!           ZdofE      - the element unconstrained interface H(curl) dof
!           ZdofV      - the element unconstrained interface H(div) dof
!
!----------------------------------------------------------------------
#include "typedefs.h"
!
   subroutine solelmI_coms(Mdle,Ncoms, ZdofH,ZdofE,ZdofV)
!
      use data_structure3D
      implicit none
!
      integer, intent(in)  :: Mdle
      integer, intent(in)  :: Ncoms
      VTYPE  , intent(out) :: ZdofH(MAXEQNH,MAXbrickH-MAXmdlbH)
      VTYPE  , intent(out) :: ZdofE(MAXEQNE,MAXbrickE-MAXmdlbE)
      VTYPE  , intent(out) :: ZdofV(MAXEQNV,MAXbrickV-MAXmdlbV)
!
!  ...element type
      integer :: ntype
!
!  ...element order of approximation
      integer :: norder(19)
!
!  ...modified element nodes and corresponding number of dof
      integer :: nodm  (MAXNODM),ndofmH(MAXNODM),  &
                 ndofmE(MAXNODM),ndofmV(MAXNODM)
!
      integer :: nrconH(MAXbrickH),nacH(NACDIM,MAXbrickH),  &
                 nrconE(MAXbrickE),nacE(NACDIM,MAXbrickE),  &
                 nrconV(MAXbrickV),nacV(NACDIM,MAXbrickV)
!
      real(8) :: constrH(NACDIM,MAXbrickH),  &
                 constrE(NACDIM,MAXbrickE),  &
                 constrV(NACDIM,MAXbrickV)
!
!  ...modified element dof
      VTYPE :: zvalH(MAXEQNH,2*MAXbrickH)
      VTYPE :: zvalE(MAXEQNE,2*MAXbrickE)
      VTYPE :: zvalV(MAXEQNV,2*MAXbrickV)
      VTYPE :: zvoid
!
      integer :: kH,kE,kV,kQ,ivar,j,kp,l
      integer :: nrnodm,nrdoflH,nrdoflE,nrdoflV,nrdoflQ
      integer :: nvarH,nvarE,nvarV
!
!---------------------------------------------------------------------
!
!  ...determine element order of approximation
      call find_order(Mdle, norder)
      ntype = NODES(Mdle)%ntype
      select case(ntype)
         case(MDLB); norder(19) = 111
         case(MDLP); norder(15) = 11
         case(MDLN); norder(11) = 1
         case(MDLD); norder(14) = 1
      end select
!
!  ...determine number of local dof
      call celndof(ntype,norder, &
                   nrdoflH,nrdoflE,nrdoflV,nrdoflQ)
!
!  ...determine constraints' coefficients
      call logic(Mdle,2,                           &
                 nodm,ndofmH,ndofmE,ndofmV,nrnodm, &
                 nrconH,nacH,constrH,              &
                 nrconE,nacE,constrE,              &
                 nrconV,nacV,constrV)
!
!  ...eliminate the middle node from the list of modified element nodes
      nrnodm = nrnodm-1
!
!---------------------------------------------------------------------
!
!  ...initiate dof's
      ZdofH=ZERO ; ZdofE=ZERO ; ZdofV=ZERO
!
!  ...initiate counters (needed for dof_out)
      kH=0; kE=0; kV=0; kQ=0
!
!  ...copy dof's into the local arrays
      do j=1,nrnodm
         call dof_out(nodm(j),Ncoms, kH,kE,kV,kQ,zvalH,zvalE,zvalV,zvoid)
      enddo
      if (kQ.ne.0) then
         write(*,*) 'solelmI: INCONSISTENCY, kQ = ',kQ
         stop
      endif
!
!      write(*,*) 'solelmI:'
!      write(*,5010) kH,kE,kV
! 5010 format('kH = ',I2,', kE = ',I2,', kV = ',I2)
!
!      write(*,5020) nrdoflH,nrdoflE,nrdoflV
! 5020 format('nrdoflH = ',I2,', nrdoflE = ',I2,', nrdoflV = ',I2)
!
!      kH = MAXbrickH-MAXmdlbH
!      kE = MAXbrickE-MAXmdlbE
!      kV = MAXbrickV-MAXmdlbV
!      write(*,5030) kH, kE, kV
! 5030 format('MAXbrickH-MAXmdlbH = ',I3,',  &
!              MAXbrickE-MAXmdlbE = ',I3,', MAXbrickV-MAXmdlbV = ',I3)
!
      nvarH = NRHVAR*NRRHS
      nvarE = NREVAR*NRRHS
      nvarV = NRVVAR*NRRHS
!
!-----------------------------------------------------------------------
!     H1 DOF'S
!-----------------------------------------------------------------------
      if (nvarH .eq. 0) goto 20
!
!  ...loop through the local dof
      do kH=1,nrdoflH
!
!  ......accumulate for the values
         do kp=1,nrconH(kH)
            l = nacH(kp,kH)
            do ivar=1,NRHVAR*NRRHS
               ZdofH(ivar,kH) = ZdofH(ivar,kH) + constrH(kp,kH)*zvalH(ivar,l)
            enddo
         enddo
!
!  ...loop through local dof
      enddo
!
   20 continue
!
!-----------------------------------------------------------------------
!     H(curl) DOF'S
!-----------------------------------------------------------------------
      if (nvarE .eq. 0) goto 30
!
!  ...loop through the local dof
      do kE=1,nrdoflE
!
!  ......accumulate for the values
         do kp=1,nrconE(kE)
            l = nacE(kp,kE)
            do ivar=1,NREVAR*NRRHS
               ZdofE(ivar,kE) = ZdofE(ivar,kE) + constrE(kp,kE)*zvalE(ivar,l)
            enddo
         enddo
!
!  ...loop through local H(curl) dof
      enddo
!
   30 continue
!
!-----------------------------------------------------------------------
!     H(div) DOF'S
!-----------------------------------------------------------------------
      if (nvarV .eq. 0) goto 40
!
!  ...loop through the local dof
      do kV=1,nrdoflV
!
!  ......accumulate for the values
         do kp=1,nrconV(kV)
            l = nacV(kp,kV)
            do ivar=1,NRVVAR*NRRHS
               ZdofV(ivar,kV) = ZdofV(ivar,kV) + constrV(kp,kV)*zvalV(ivar,l)
            enddo
         enddo
!
!  ...loop through local H(div) dof
      enddo
!
   40 continue
!
   end subroutine solelmI_coms
!
!
!-----------------------------------------------------------------------
!..see routine solelmI_coms for description
!  last rev: Sep 2023
   subroutine solelmI(Mdle, ZdofH,ZdofE,ZdofV)
!
      use data_structure3D
      implicit none
!
      integer, intent(in)  :: Mdle
      VTYPE  , intent(out) :: ZdofH(MAXEQNH,MAXbrickH-MAXmdlbH)
      VTYPE  , intent(out) :: ZdofE(MAXEQNE,MAXbrickE-MAXmdlbE)
      VTYPE  , intent(out) :: ZdofV(MAXEQNV,MAXbrickV-MAXmdlbV)
!
!  ...by default, request solution component set N_COMS
      call solelmI_coms(Mdle,N_COMS, ZdofH,ZdofE,ZdofV)
!
   end subroutine solelmI
