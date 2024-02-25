!----------------------------------------------------------------------
!
!   routine name       - solelm
!
!----------------------------------------------------------------------
!
!   latest revision    - Sep 2023
!
!   purpose            - routine calculates unconstrained solution dof
!                        for a 3D element
!   remark: this routine must be OMP thread-safe
!
!   arguments :
!     in:
!           Mdle       - middle node of an element
!           Ncoms      - solution component set: 1,...,NRCOMS
!     out:
!           ZdofH      - the element unconstrained H1 dof
!           ZdofE      - the element unconstrained H(curl) dof
!           ZdofV      - the element unconstrained H(div) dof
!           ZdofQ      - the element unconstrained L2 dof
!
!----------------------------------------------------------------------
#include "typedefs.h"
!
   subroutine solelm_coms(Mdle,Ncoms, ZdofH,ZdofE,ZdofV,ZdofQ)
!
      use data_structure3D
!
      implicit none
!
      integer, intent(in)  :: Mdle
      integer, intent(in)  :: Ncoms
      VTYPE  , intent(out) :: ZdofH(MAXEQNH,MAXbrickH)
      VTYPE  , intent(out) :: ZdofE(MAXEQNE,MAXbrickE)
      VTYPE  , intent(out) :: ZdofV(MAXEQNV,MAXbrickV)
      VTYPE  , intent(out) :: ZdofQ(MAXEQNQ,MAXbrickQ)
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
!
!  ...miscellanea
      integer :: nrdoflH,nrdoflE,nrdoflV,nrdoflQ
      integer :: kH,kE,kV,kQ,j,l,kp,ivar,nrnodm
      integer :: nvarH,nvarE,nvarV
!
#if DEBUG_MODE
      integer :: nvarQ
      integer :: iprint
      iprint=0
#endif
!
!---------------------------------------------------------------------
!
!  ...determine element order of approximation
      call find_order(Mdle, norder)
!
#if DEBUG_MODE
      if (iprint.ge.1) then
         write(*,7000) Mdle
         call print_order(NODES(Mdle)%ntype,norder)
      endif
#endif
!
!  ...determine number of local dof
      call celndof(NODES(Mdle)%ntype,norder, &
                   nrdoflH,nrdoflE,nrdoflV,nrdoflQ)
!
!  ...determine constraints' coefficients
      call logic(Mdle,2,                           &
                 nodm,ndofmH,ndofmE,ndofmV,nrnodm, &
                 nrconH,nacH,constrH,              &
                 nrconE,nacE,constrE,              &
                 nrconV,nacV,constrV)
!
!---------------------------------------------------------------------
!
!  ...initiate dof's
      ZdofH=ZERO ; ZdofE=ZERO ; ZdofV=ZERO ; ZdofQ=ZERO
!
!  ...initiate counters (needed for dof_out)
      kH=0; kE=0; kV=0; kQ=0
!
!  ...copy dof's into the local arrays (use ZdofQ for simplicity)
      do j=1,nrnodm
         call dof_out(nodm(j),Ncoms, kH,kE,kV,kQ,zvalH,zvalE,zvalV,ZdofQ)
      enddo
!
#if DEBUG_MODE
   7000 format(' solelm: Mdle = ',i10)
   7001 format(' solelm: nac    = ',10i10)
   7002 format(' solelm: constr = ',10(2x,f8.3))
#if C_MODE
   7003 format(' solelm: ivar,zval = ',i2,', ',10(2e12.5,2x))
   7004 format(' solelm: ZdofH = ',10(2e12.5,2x))
   7005 format(' solelm: ZdofE = ',10(2e12.5,2x))
   7006 format(' solelm: ZdofV = ',10(2e12.5,2x))
   7007 format(' solelm: ZdofQ = ',10(2e12.5,2x))
#else
   7003 format(' solelm: ivar,zval = ',i2,', ',20e12.5)
   7004 format(' solelm: ZdofH = ',20e12.5)
   7005 format(' solelm: ZdofE = ',20e12.5)
   7006 format(' solelm: ZdofV = ',20e12.5)
   7007 format(' solelm: ZdofQ = ',20e12.5)
#endif
#endif
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
#if DEBUG_MODE
         if (iprint.eq.1) then
            write(*,*) 'solelm: kH = ',kH
            write(*,7001) (   nacH(l,kH),l=1,nrconH(kH))
            write(*,7002) (constrH(l,kH),l=1,nrconH(kH))
            do ivar=1,nvarH
               write(*,7003) ivar,(zvalH(ivar,nacH(l,kH)),l=1,nrconH(kH))
            enddo
         endif
#endif
!
!  ......accumulate for the values
         do kp=1,nrconH(kH)
            l = nacH(kp,kH)
            do ivar=1,nvarH
               ZdofH(ivar,kH) = ZdofH(ivar,kH) + constrH(kp,kH)*zvalH(ivar,l)
            enddo
         enddo
!
#if DEBUG_MODE
         if (iprint.eq.1) then
            write(*,7004) ZdofH(1:nvarH,kH)
         endif
#endif
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
#if DEBUG_MODE
         if (iprint.eq.2) then
            write(*,*) 'solelm: kE =',kE
            write(*,7001) (   nacE(l,kE),l=1,nrconE(kE))
            write(*,7002) (constrE(l,kE),l=1,nrconE(kE))
            do ivar=1,nvarE
               write(*,7003) ivar,(zvalE(ivar,nacE(l,kE)),l=1,nrconE(kE))
            enddo
         endif
#endif
!
!  ......accumulate for the values
         do kp=1,nrconE(kE)
            l = nacE(kp,kE)
            do ivar=1,nvarE
               ZdofE(ivar,kE) = ZdofE(ivar,kE) + constrE(kp,kE)*zvalE(ivar,l)
            enddo
         enddo
!
#if DEBUG_MODE
         if (iprint.eq.2) then
            write(*,7005) ZdofE(1:nvarE,kE)
         endif
#endif
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
#if DEBUG_MODE
         if (iprint.eq.3) then
            write(*,*) 'solelm: kV = ',kV
            write(*,7001) (   nacV(l,kV),l=1,nrconV(kV))
            write(*,7002) (constrV(l,kV),l=1,nrconV(kV))
            do ivar=1,nvarV
               write(*,7003) ivar,(zvalV(ivar,nacV(l,kV)),l=1,nrconV(kV))
            enddo
         endif
#endif
!
!  ......accumulate for the values
         do kp=1,nrconV(kV)
            l = nacV(kp,kV)
            do ivar=1,nvarV
               ZdofV(ivar,kV) = ZdofV(ivar,kV) + constrV(kp,kV)*zvalV(ivar,l)
            enddo
         enddo
!
#if DEBUG_MODE
         if (iprint.eq.3) then
            write(*,7006) ZdofV(1:nvarV,kV)
         endif
#endif
!
!  ...loop through local H(div) dof
      enddo
!
   40 continue
!
!-----------------------------------------------------------------------
!     L2 DOF'S , nothing to do...
!-----------------------------------------------------------------------
!
#if DEBUG_MODE
      if (iprint.eq.4) then
         nvarQ = NRQVAR*NRRHS
         if (nvarQ > 0) then
            do kQ=1,nrdoflQ
               write(*,*) 'solelm: kQ = ',kQ
               write(*,7007) ZdofQ(1:nvarQ,kQ)
            enddo
         endif
      endif
#endif
!
   end subroutine solelm_coms
!
!
!-----------------------------------------------------------------------
!..see routine solelm_coms for description
!  last rev: Sep 2023
   subroutine solelm(Mdle, ZdofH,ZdofE,ZdofV,ZdofQ)
!
      use data_structure3D
      implicit none
!
      integer, intent(in)  :: Mdle
      VTYPE  , intent(out) :: ZdofH(MAXEQNH,MAXbrickH)
      VTYPE  , intent(out) :: ZdofE(MAXEQNE,MAXbrickE)
      VTYPE  , intent(out) :: ZdofV(MAXEQNV,MAXbrickV)
      VTYPE  , intent(out) :: ZdofQ(MAXEQNQ,MAXbrickQ)
!
!  ...by default, request solution component set N_COMS
      call solelm_coms(Mdle,N_COMS, ZdofH,ZdofE,ZdofV,ZdofQ)
!
   end subroutine solelm
