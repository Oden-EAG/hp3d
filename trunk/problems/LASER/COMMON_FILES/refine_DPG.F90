!--------------------------------------------------------------------
!
!     routine name      - refine_DPG
!
!--------------------------------------------------------------------
!
!     last revision:    - Oct 2019
!
!     purpose:          - refines elements assuming problem has
!                         already been solved. If uniform refinements
!                         it refines everything. Otherwise it follows
!                         greedy strategy based on residual to
!                         determine which elements to refine.
!                         It also displays dof, residuals and
!                         relevant errors at each step.
!
!     arguments:
!
!     in:
!             Irefine   = 0 if no refinement -> only display
!                       = 1 if uniform refinements
!                       = 2 if anisotropic refinements (provide kref)
!                       = 3 if adaptive refinements
!             Nreflag   = 1 for h-refinements
!                       = 2 for p-refinements
!                       = 3 for h- or p-refinement, depending
!                           upon the element size
!             Factor    - if element error \ge Factor*max_error
!                         then the element is refined
!             Nflag     - integer vector of length NR_PHYSA to indicate
!                         which PHYS Attribute one should compute the error for
!             PhysNick  - nickname for PHYS Attribute: physNick = HEVQ
!                         e.g., physNick = 1001 means 1H1 and 1L2 variable, etc.
!             Ires      - logical, .true. if compute residual is true
!     out:
!             Nstop     = 1 if refinement did not increase dof
!
!---------------------------------------------------------------------
!
subroutine refine_DPG(Irefine,Nreflag,Factor,Nflag,PhysNick,Ires, Nstop)
!
   use commonParam
   use control
   use data_structure3D
   use environment  , only: QUIET_MODE
   use assembly_sc  , only: NRDOF_CON, NRDOF_TOT
   use parametersDPG, only: NORD_ADD
   use par_mesh     , only: DISTRIBUTED,HOST_MESH
   use mpi_param    , only: ROOT,RANK
   use MPI          , only: MPI_COMM_WORLD,MPI_SUM,MPI_COMM_WORLD,   &
                            MPI_REAL8,MPI_IN_PLACE,MPI_MAX
!
   implicit none
!
   integer, intent(in)  :: Irefine
   integer, intent(in)  :: Nreflag
   real(8), intent(in)  :: Factor
   integer, intent(in)  :: Nflag(NR_PHYSA)
   integer, intent(in)  :: PhysNick
   logical, intent(in)  :: Ires
   integer, intent(out) :: Nstop
!
   integer, parameter :: max_step = 20
   integer, dimension(max_step), save :: nrdof_tot_mesh
   integer, dimension(max_step), save :: nrdof_con_mesh
   real(8), dimension(max_step), save :: residual_mesh
   real(8), dimension(max_step), save :: rate_mesh
   real(8), dimension(max_step), save :: error_mesh
   real(8), dimension(max_step), save :: rel_error_mesh
   real(8), dimension(max_step), save :: rate_error_mesh
!
   real(8) :: elem_resid(NRELES)
   integer :: elem_ref_flag(NRELES)
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
!
   real(8) :: resid_subd,resid_tot,elem_resid_max
   real(8) :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
   real(8) :: error_tot,rnorm_tot,error_subd,rnorm_subd
!
   integer :: i,ic,mdle,iel,kref,iii,subd,count,ierr
   integer :: vref(3), jref(3)
   integer :: nord, nordx, nordy, nordz, nordyz
   integer :: nordx_new, nordy_new, nordz_new, nord_new
   integer :: idec_pref,idec_href, iref, enforce_flag, icpref
!
!..geometry dof (work space for nodcor)
   real(8) :: xnod(3,MAXbrickH)
   real(8) :: maxz,minz,midz
!
!..element type
   character(len=4) :: etype
!
   real(8) :: MPI_Wtime,start_time,end_time
!
!..printing flag
   integer :: iprint = 0
!
!-----------------------------------------------------------------------
!
!..initialize
   elem_resid = 0.d0; elem_ref_flag = 0.d0
!
!..increase step if necessary
!..irefineold=0 means no refinement was made in the previous step
!..if first call or if a refinement was made, increase step
   if ((istep.eq.0).or.(irefineold.ne.0)) istep=istep+1
   irefineold=Irefine
!
!..get dof count from solver
   nrdof_tot_mesh(istep) = NRDOF_TOT
   nrdof_con_mesh(istep) = NRDOF_CON
   if (NRDOF_TOT .eq. 0) then
      istep=istep-1
      goto 90
   endif
!!
!!..fetch active elements
!   if (.not. DISTRIBUTED) then
!      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
!      NRELES_SUBD = NRELES
!   endif
!
!..initialize residual and error
   resid_subd = 0.d0
   error_subd = 0.d0
   rnorm_subd = 0.d0
   elem_resid_max = 0.d0
!
!..start timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..residual/error computation
!
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(errorH,errorE,errorV,errorQ,            &
!$OMP         rnormH,rnormE,rnormV,rnormQ,            &
!$OMP         mdle,subd,etype,xnod,maxz,minz,midz)    &
!$OMP SCHEDULE(DYNAMIC)                               &
!$OMP REDUCTION(+:resid_subd,error_subd,rnorm_subd)   &
!$OMP REDUCTION(MAX:elem_resid_max)
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
      if (DISTRIBUTED .and. (RANK .ne. subd)) cycle
      if (USE_PML) then
         xnod = 0.d0
         call nodcor(Mdle, xnod)
         etype = NODES(Mdle)%type
         select case(etype)
            case('mdlb')
               maxz = maxval(xnod(3,1:8))
               minz = minval(xnod(3,1:8))
            case('mdlp')
               maxz = maxval(xnod(3,1:6))
               minz = minval(xnod(3,1:6))
            case default
               write(*,*) 'refine_DPG: invalid etype param. stop.'
               stop
         end select
         !midz = minz + (maxz-minz)/2.d0
         if (maxz .gt. PML_REGION) cycle
      endif
      if(ires) then
         call elem_residual(mdle, elem_resid(iel),elem_ref_flag(iel))
         resid_subd = resid_subd + elem_resid(iel)
         if (elem_resid(iel) > elem_resid_max) elem_resid_max = elem_resid(iel)
      endif
      if (NEXACT.ge.1) then
         call element_error(mdle,Nflag, errorH,errorE,errorV,errorQ,  &
                                        rnormH,rnormE,rnormV,rnormQ)
         select case(PhysNick)
            case(1)
               error_subd = error_subd + errorQ
               rnorm_subd = rnorm_subd + rnormQ
            case(10)
               error_subd = error_subd + errorV
               rnorm_subd = rnorm_subd + rnormV
            case(100)
               error_subd = error_subd + errorE
               rnorm_subd = rnorm_subd + rnormE
            case(1000)
               error_subd = error_subd + errorH
               rnorm_subd = rnorm_subd + rnormH
            case(1001)
               error_subd = error_subd + errorH + errorQ
               rnorm_subd = rnorm_subd + rnormH + rnormQ
            case default
               error_subd = error_subd + errorH + errorE + errorV + errorQ
               rnorm_subd = rnorm_subd + rnormH + rnormE + rnormV + rnormQ
         end select
      endif
   enddo
!$OMP END PARALLEL DO
!
   resid_tot = 0.d0; error_tot = 0.d0; rnorm_tot = 0.d0
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      count = 1
      call MPI_REDUCE(resid_subd,resid_tot,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(error_subd,error_tot,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(rnorm_subd,rnorm_tot,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid_max,count,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      count = NRELES
      call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid,count,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
   else
      resid_tot = resid_subd
      error_tot = error_subd
      rnorm_tot = rnorm_subd
   endif
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (.not. QUIET_MODE) then
      if (RANK .eq. ROOT) write(*,1010) end_time-start_time
 1010 format(' elem_residual  : ',f12.5,'  seconds')
   endif
!
!..mark elements for adaptive refinement
!   do iel=1,NRELES
!      mdle = ELEM_ORDER(iel)
!      call find_domain(mdle, i)
!      if(elem_resid(iel)>0.75d0*elem_resid_max) then
!         !write(*,1020) 'ndom=',i,', mdle=',mdle,', r =',elem_resid(iel)
!      endif
! 1020 format(A,I2,A,I5,A,es12.5)
!   enddo
!
   if (RANK .ne. ROOT) goto 90
!
!..update and display convergence history
   residual_mesh(istep) = sqrt(resid_tot)
   if (NEXACT.ge.1) then
      error_mesh(istep)     = sqrt(error_tot)
      rel_error_mesh(istep) = sqrt(error_tot/rnorm_tot)
   endif
!
!..compute decrease rate for the residual and error
   select case(istep)
   case(1)
      rate_mesh(istep) = 0.d0
      if (NEXACT.ge.1) rate_error_mesh(istep) = 0.d0
   case default
      rate_mesh(istep) =  &
      log(residual_mesh(istep-1)/residual_mesh(istep))/  &
      log(float(nrdof_tot_mesh(istep-1))/float(nrdof_tot_mesh(istep)))
      if (NEXACT.ge.1) then
         rate_error_mesh(istep) = &
         log(error_mesh(istep-1)/error_mesh(istep))/  &
         log(float(nrdof_tot_mesh(istep-1))/float(nrdof_tot_mesh(istep)))
      endif
   end select
!
!..print out the history of refinements
   write(*,*)
   write(*,*) 'HISTORY OF REFINEMENTS'
   if (NEXACT.eq.0) write(*,7005)
   if (NEXACT.gt.0) write(*,7006)
 7005  format(' mesh |',' nrdof_tot |',' nrdof_con |','    residual   |','   residual rate  ',/)
 7006  format(' mesh |',' nrdof_tot |',' nrdof_con |','    residual   |','   residual rate  |', &
              ' field error  |','rel field error|','   error rate ',/)
!
   do i=1,istep
      if (NEXACT.eq.0) then
         write(*,7003) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i)
 7003    format(2x,i2,'  | ',2(i9,' | '),es12.5,'  |',f7.2)
      else
         write(*,7004) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i), &
                       error_mesh(i),rel_error_mesh(i),rate_error_mesh(i)
 7004    format(2x,i2,'  | ',2(i9,' | '),es12.5,'  |',f7.2,'          ', &
                2(' | ',es12.5),'  |',f7.2)
      endif
      if (i .eq. istep) write(*,*)
   enddo
!
   90 continue
!
!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
!..use appropriate strategy depending on refinement type
   Nstop = 0
   select case(Irefine)
!  ...uniform refinements
      case(IUNIFORM)
!     ...isotropic href
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call global_href
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
         call update_gdof
         call update_Ddof
!   ...anisotropic refinements
      case(IANISOTROPIC)
!     ...anisotropic href (z)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call global_href_aniso(0,1) ! refine in z
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
         call update_gdof
         call update_Ddof
!
      case default; Nstop = 1
!
   end select
!
  2020 format(' global ref : ',f12.5,'  seconds')
  2025 format(' close mesh : ',f12.5,'  seconds')
  2030 format(A,I8,', ',I9)
!
!-----------------------------------------------------------------------
!                              FINALIZE
!-----------------------------------------------------------------------
!
   if (IBCFLAG .eq. 3) then
      call propagate_flag(2,9)
      call propagate_flag(3,9)
   endif
!
end subroutine refine_DPG
!
!
!-----------------------------------------------------------------------
! subroutine: href_solve (uniform refinements and solve)
!-----------------------------------------------------------------------
subroutine href_solve(Nflag,PhysNick, Nstop)
!
   use commonParam
   use control
   use data_structure3D
   use MPI           , only: MPI_COMM_WORLD,MPI_INTEGER
   use mpi_param     , only: RANK,ROOT
   use par_mesh      , only: DISTRIBUTED,HOST_MESH,distr_mesh
   use zoltan_wrapper, only: zoltan_w_set_lb,zoltan_w_partition
!
   implicit none
!
   integer, intent(in)  :: Nflag(NR_PHYSA)
   integer, intent(in)  :: PhysNick
   integer, intent(out) :: Nstop
!
   integer :: nsteps,count,src,ierr,i
   logical :: ires
!
   Nstop = 0
   ires = .true.
!
   if (RANK .eq. ROOT) then
      nsteps=0
      do while (nsteps.le.0)
         write(*,*) 'Provide: number of uniform h-refinements'
         read(*,*) nsteps
      enddo
   else
      write(6,405) '[', RANK, '] : ','Waiting for broadcast from master...'
 405  format(A,I4,A,A)
   endif
   count = 1; src = ROOT
   call MPI_BCAST (nsteps,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
!
!..initial solve
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      call par_mumps_sc('H')
   else
      !call mumps_sc('G')
      call pardiso_sc('H')
   endif
!
!..do refinements and solve
   do i=0,nsteps
!  ...display error and refine if necessary
      if (i.ne.nsteps) then
!     ...Uniform refinement and solve
         call refine_DPG(IUNIFORM,1,0.25d0,Nflag,PhysNick,ires, Nstop)
         if (DISTRIBUTED .and. (.not. HOST_MESH)) then
            !call zoltan_w_set_lb(1)
            !call distr_mesh
            !call print_partition
            call par_mumps_sc('H')
         else
            !call mumps_sc('H')
            call pardiso_sc('H')
         endif
      else
!     ...Last step only display (no refinement)
         call refine_DPG(INOREFINEMENT,1,0.25d0,Nflag,PhysNick,ires, Nstop)
      endif
   enddo
!
end subroutine href_solve


