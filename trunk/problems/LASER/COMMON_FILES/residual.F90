!----------------------------------------------------------------------
!
!    routine            - residual
!
!----------------------------------------------------------------------
!
!     latest revision:  - Sep 2019
!
!     purpose:          - compute and print residual
!
!----------------------------------------------------------------------
subroutine residual()
!
   use data_structure3D
   use commonParam
   use control
   use environment
   use assembly_sc, only: NRDOF_TOT,NRDOF_CON
   use par_mesh   , only: DISTRIBUTED,HOST_MESH
   use mpi_param  , only: ROOT,RANK
   use MPI        , only: MPI_SUM,MPI_COMM_WORLD,MPI_REAL8
!
   implicit none
!
!..workspace for element_error routine
   real(8) :: resid_subd,resid_tot
   integer :: iel,mdle,count,ierr
!
   real(8) :: elem_resid
   integer :: elem_ref_flag
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
!..fetch active elements
   if ((DISTRIBUTED) .and. (.not. HOST_MESH)) then
      if (RANK .eq. ROOT) then
         write(*,*) 'residual: mesh is distributed. computing error in parallel...'
      endif
   else
      if (RANK .ne. ROOT) goto 90
      write(*,*) 'residual: mesh is not distributed (or on host). computing error on host...'
      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      NRELES_SUBD = NRELES
   endif
!
!..initialize residual
   resid_subd = 0.d0
!
!..start timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..residual/error computation
!
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(mdle,etype,xnod,maxz,minz,midz,         &
!$OMP         elem_resid,elem_ref_flag)               &
!$OMP SCHEDULE(DYNAMIC)                               &
!$OMP REDUCTION(+:resid_subd)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
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
      call elem_residual(mdle, elem_resid,elem_ref_flag)
      resid_subd = resid_subd + elem_resid
   enddo
!$OMP END PARALLEL DO
!
   resid_tot = 0.d0
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      count = 1
      call MPI_REDUCE(resid_subd,resid_tot,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
   else
      resid_tot = resid_subd
   endif
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (.not. QUIET_MODE) then
      if (RANK .eq. ROOT) write(*,1010) end_time-start_time
 1010 format(' elem_residual  : ',f12.5,'  seconds')
   endif
!
   if (RANK .eq. ROOT) then
      write(*,7020) NRDOF_TOT,sqrt(resid_tot)
 7020 format('residual: NRDOF_TOT, RESIDUAL = ',i9,',  ',es12.5)
   endif
!
   90 continue
!
end subroutine residual

