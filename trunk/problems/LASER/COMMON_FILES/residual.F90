!----------------------------------------------------------------------
!
!    routine            - residual
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 2019
!
!     purpose:          - compute and print residual
!
!----------------------------------------------------------------------
subroutine residual(res)
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
!..computed total residual
   real(8), intent(out) :: res
!
!..workspace for element_error routine
   real(8) :: resid_subd,resid_tot
   integer :: iel,mdle,count,ierr
!
   real(8) :: elem_resid
   integer :: elem_ref_flag
!
!..timer
   real(8) :: MPI_Wtime,start_time,end_time
!
!..printing flag
   integer :: iprint = 0
!
!-----------------------------------------------------------------------
!
!..start timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..fetch active elements
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
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
!..residual/error computation
!
!$OMP PARALLEL DO                               &
!$OMP PRIVATE(mdle,elem_resid,elem_ref_flag)    &
!$OMP SCHEDULE(DYNAMIC)                         &
!$OMP REDUCTION(+:resid_subd)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      call elem_residual(mdle, elem_resid,elem_ref_flag)
      if (USE_PML .and. Is_pml(mdle)) then
      !  treat PML differently if needed
         resid_subd = resid_subd + elem_resid
      else
         resid_subd = resid_subd + elem_resid
      endif
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
 90 continue
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (.not. QUIET_MODE) then
      if (RANK .eq. ROOT) write(*,1010) end_time-start_time
 1010 format(' elem_residual  : ',f12.5,'  seconds')
   endif
!
   res = 0.d0
   if (RANK .eq. ROOT) then
      res = sqrt(resid_tot)
      write(*,7020) NRDOF_TOT,NRDOF_CON,res
 7020 format(' residual: NRDOF_TOT, NRDOF_CON, RESIDUAL = ',i9,',  ',i9,',  ',es12.5)
   endif
!
end subroutine residual
!
!
