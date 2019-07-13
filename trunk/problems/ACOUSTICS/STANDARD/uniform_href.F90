!--------------------------------------------------------------------
!
!     routine name      - uniform_href
!
!--------------------------------------------------------------------
!
!     latest revision:  - July 17
!
!     purpose:          -  performs single uniform h-refinement
!
!     arguments:
!
!     in:
!             Irefine   = 0 if no refinement -> only display
!                       = 1 if uniform refinements
!             Nreflag   = 1 for h-refinements
!                       = 2 for p-refinements
!             Factor    - if element error \ge Factor*max_error
!                         then the element is refined
!     out:
!             Nstop     = 1 if refinement did not increase dof
!
!---------------------------------------------------------------------
!
   subroutine uniform_href(Irefine,Nreflag,Factor, Nstop)
!
   use control
   use data_structure3D
   use environment
   use common_prob_data
   use m_assembly, ONLY: NRDOF_CON, NRDOF_TOT
!
   implicit none
!   
   integer, intent(in)  :: Irefine
   integer, intent(in)  :: Nreflag
   real*8,  intent(in)  :: Factor
   integer, intent(out) :: Nstop
!
   integer, parameter :: max_step = 300
   integer, dimension(max_step), save :: nrdof_tot_mesh
   integer, dimension(max_step), save :: nrdof_con_mesh
   real*8,  dimension(max_step), save :: rate_mesh
   real*8,  dimension(max_step), save :: error_mesh
   real*8,  dimension(max_step), save :: rel_error_mesh
   real*8,  dimension(max_step), save :: rate_error_mesh

!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
   integer :: nrdofField, nrdofDPG
   integer, dimension(NR_PHYSA) :: nflag
   integer, dimension(NRELES)   :: N_elem
   real*8 :: err, rnorm
   real*8 :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
   integer :: i,iprint,ic,mdle,iel,kref
   integer :: idec_href, iref

!
!-----------------------------------------------------------------------
!
!..initialize
   iprint = 0
!
!..increase step if necessary
!..irefineold=0 means no refinement was made in the previous step
!..if first call or if a refinement was made, increase step
   if ((istep.eq.0).or.(irefineold.ne.0)) istep=istep+1
   irefineold=Irefine
!
!..get dof count from solver
! call ndof_DPG !calculates nrdofDPG and nrdofField
   nrdof_tot_mesh(istep) = NRDOF_TOT
   nrdof_con_mesh(istep) = NRDOF_CON
!
!..field variables flag - used inside ndof_DPG and for error calc.
   nflag(1:NR_PHYSA)=1
!
!..create list of mdle nods numbers
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      N_elem(iel) = mdle
   enddo
!   
!..initialize global error
   err = 0.d0; rnorm = 0.d0
!
   ! L2PROJ = .TRUE.
!$OMP PARALLEL DEFAULT(PRIVATE)                              &
!$OMP SHARED(N_elem,NRELES,NEXACT,nflag)  &
!$OMP REDUCTION(+:err,rnorm)           
!$OMP DO SCHEDULE(DYNAMIC)
   do iel=1,NRELES
      if (NEXACT.ge.1) then
         call element_error(N_elem(iel),nflag,errorH,errorE,errorV,errorQ,  &
                                        rnormH,rnormE,rnormV,rnormQ)
         err   = err   + errorH
         rnorm = rnorm + rnormH
      endif
   enddo
!$OMP END DO        
!$OMP END PARALLEL
!
!..update and display convergence history
   if (NEXACT.ge.1) then
      error_mesh(istep) = sqrt(err)
      rel_error_mesh(istep) = sqrt(err/rnorm)
   endif
!
!..compute decrease rate for the residual and error
   select case(istep)
   case(1)
      rate_mesh(istep) = 0.d0
      if (NEXACT.ge.1) rate_error_mesh(istep) = 0.d0
   case default
      if (NEXACT.ge.1) then
         rate_error_mesh(istep) = &
         log(rel_error_mesh(istep-1)/rel_error_mesh(istep))/  &
         log(float(nrdof_tot_mesh(istep-1))/float(nrdof_tot_mesh(istep)))
      endif
   end select
!
!..print out the history of refinements
   write(*,*)
   write(*,*)
   write(*,*) 'HISTORY OF REFINEMENTS'
   if (NEXACT.gt.0) write(*,7006)
 7006  format(' mesh |',' nrdof_tot |','nrdof_con|','|',   &
              ' field error  |','rel field error|','   rate ')
   write(*,*)
    
   do i=1,istep
      if (NEXACT.gt.0) then
      write(*,7004) i,nrdof_tot_mesh(i),nrdof_con_mesh(i) &
                    ,error_mesh(i), &
                    rel_error_mesh(i),rate_error_mesh(i)
 7004 format(2x,i2,'  | ',2(i8,' | '), 2(' | ',e12.5),'  |',f7.2)
      endif
   enddo
!
!




!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
!..use appropriate strategy depening on refinement type
!  NOTE: if Irefine=0 nothing is done, so only relevant info (the
!        code above this line) is displayed.
   select case(Irefine)
!..uniform refinements
   case(IUNIFORM)
      call global_href
      call update_gdof
      if (IBC_PROB .ne. BC_IMPEDANCE) then
         write(*,*) 'Calling update_ddof'
         call update_ddof
      endif   
      ! call verify_orient
      ! call verify_neig
!
   end select
!
   Nstop = 0 ! as expected, elements were refined
! 
! 
   end subroutine uniform_href
