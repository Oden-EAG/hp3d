!-----------------------------------------------------------------------
!> @brief Routine extracts L2 dofs for an element
!> @param[in]  Mdle     - middle node number of the coarse element
!> @param[in]  ZdofQ    - array containing the L2 dofs
!> @date June 2024
!-----------------------------------------------------------------------
subroutine solelm_L2(Mdle,ZdofQ)
!
   use data_structure3D
!
   implicit none
!   
   integer, intent(in)  :: Mdle
   real,    intent(out) :: ZdofQ(MAXEQNQ,MAXbrickQ)
!
!
   integer :: norder(19)
!
!..modified element nodes and corresponding number of dof
   integer :: nodm  (MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM),ndofmV(MAXNODM)
!
   integer :: nrconH(MAXbrickH),nrconE(MAXbrickE),nrconV(MAXbrickV), &
              nacH(NACDIM,MAXbrickH), nacE(NACDIM,MAXbrickE), nacV(NACDIM,MAXbrickV)

   real(8) :: constrH(NACDIM,MAXbrickH), constrE(NACDIM,MAXbrickE), constrV(NACDIM,MAXbrickV)
!      
   integer :: kH,kE,kV,kQ,j
   integer :: nrdoflH,nrdoflE,nrdoflV,nrdoflQ,nrnodm
!
   real :: zvalH(MAXEQNH,2*MAXbrickH)
   real :: zvalE(MAXEQNE,2*MAXbrickE)
   real :: zvalV(MAXEQNV,2*MAXbrickV)
!
! ---------------------------------------------------------------------
!
   integer :: iprint=0
!..determine element order of approximation
   call find_order(Mdle, norder)
!
   if (iprint.ge.1) then
      write(*,7000) Mdle
7000   format('solelm: Mdle   = ',i6)
      call print_order(NODES(Mdle)%ntype,norder)
   endif
!..determine number of local dof
   call celndof(NODES(Mdle)%ntype,norder, &
               nrdoflH,nrdoflE,nrdoflV,nrdoflQ)
!
!..determine constraints' coefficients
   call logic(Mdle,2, &
               nodm,ndofmH,ndofmE,ndofmV,nrnodm, &
               nrconH,nacH,constrH, &
               nrconE,nacE,constrE, &
               nrconV,nacV,constrV)
!---------------------------------------------------------------------
!
!..initiate dof's (use ZdofQ for simplicity)
   zvalH=ZERO ; zvalE=ZERO ; zvalV=ZERO ; ZdofQ=ZERO
!
!..initiate counter for some mystic reason
   kH=0; kE=0; kV=0; kQ=0
!
!..copy dof's into the local arrays
   do j=1,nrnodm
      call dof_out(nodm(j),N_COMS, kH,kE,kV,kQ,zvalH,zvalE,zvalV,zdofQ)
   enddo
!
end subroutine solelm_L2
