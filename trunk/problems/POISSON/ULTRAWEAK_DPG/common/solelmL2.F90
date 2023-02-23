! ----------------------------------------------------------------------

!    routine name       - solelmL2

! ----------------------------------------------------------------------

!    latest revision    - Feb 23

!    purpose            - routine calculates unconstrained L2 solution dof
!                         for a 3D element

!    arguments :
!      in:
!            Mdle       - middle node of an element
!      out:
!            ZdofH      - the element unconstrained H1 dof
! !           ZdofQ      - the element unconstrained L2 dof

! ----------------------------------------------------------------------

      subroutine solelm_L2(Mdle,ZdofQ)

      use data_structure3D
! #include "syscom.blk"
    
      implicit none
    
      integer, intent(in) :: Mdle
      real, intent(out) :: ZdofQ(MAXEQNQ,MAXbrickQ)


      integer :: norder(19)

!  ...modified element nodes and corresponding number of dof
      integer :: nodm  (MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM),ndofmV(MAXNODM)
!
      integer :: nrconH(MAXbrickH),nrconE(MAXbrickE),nrconV(MAXbrickV), &
               nacH(NACDIM,MAXbrickH), nacE(NACDIM,MAXbrickE), nacV(NACDIM,MAXbrickV)

      real(8) :: constrH(NACDIM,MAXbrickH), constrE(NACDIM,MAXbrickE), constrV(NACDIM,MAXbrickV)
      
      integer :: kH,kE,kV,kQ,j
      integer :: nrdoflH,nrdoflE,nrdoflV,nrdoflQ,nrnodm

      real :: zvalH(MAXEQNH,2*MAXbrickH)
      real :: zvalE(MAXEQNE,2*MAXbrickE)
      real :: zvalV(MAXEQNV,2*MAXbrickV)

! ---------------------------------------------------------------------

    integer :: iprint=0
!   ...determine element order of approximation
      call find_order(Mdle, norder)

      if (iprint.ge.1) then
        write(*,7000) Mdle
 7000   format('solelm: Mdle   = ',i6)
        call print_order(NODES(Mdle)%type,norder)
      endif
! c  ...determine number of local dof
      call celndof(NODES(Mdle)%type,norder, &
                  nrdoflH,nrdoflE,nrdoflV,nrdoflQ)

!   ...determine constraints' coefficients
      call logic(Mdle,2, &
                nodm,ndofmH,ndofmE,ndofmV,nrnodm, &
                nrconH,nacH,constrH, &
                nrconE,nacE,constrE, &
                nrconV,nacV,constrV)

!---------------------------------------------------------------------

!  ...initiate dof's (use ZdofQ for simplicity)
      zvalH=ZERO ; zvalE=ZERO ; zvalV=ZERO ; ZdofQ=ZERO
!
!  ...initiate counter for some mystic reason
      kH=0; kE=0; kV=0; kQ=0

!  ...copy dof's into the local arrays
      do j=1,nrnodm
        call dof_out(nodm(j), kH,kE,kV,kQ,zvalH,zvalE,zvalV,zdofQ)
      enddo

!  ...initiate the dof
    !   ZdofH(1:MAXEQNH,1:nrdoflH) = ZERO
    !   ZdofE(1:MAXEQNE,1:nrdoflE) = ZERO
    !   ZdofV(1:MAXEQNV,1:nrdoflV) = ZERO
!----------------------------------------------------------
!     L2 DOF'S , nothing to do...
!-----------------------------------------------------------------------

!       if (iprint.eq.4) then
!         do kQ=1,nrdoflQ
!           nvarQ = NRQVAR*NRCOMS
!           write(*,7007) kQ,ZdofQ(1:nvarQ,kQ)
!  7007     format('solelm: ZdofQ(1:nvarQ,',i2,') = ',10e12.5)
!         enddo
!         call pause
!       endif

!       if (iprint.ne.0) call pause

      end subroutine solelm_L2
