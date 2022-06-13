!----------------------------------------------------------------------
!
!    routine            - exact_error
!
!----------------------------------------------------------------------
!
!     latest revision:  - July 2019
!
!     purpose:          - compute and print exact error
!
!----------------------------------------------------------------------
subroutine exact_error
!
   use data_structure3D
   use common_prob_data
   use environment
   use assembly_sc, only: NRDOF_TOT
   use par_mesh   , only: DISTRIBUTED,HOST_MESH
   use mpi_param  , only: ROOT,RANK
   use MPI        , only: MPI_SUM,MPI_COMM_WORLD,MPI_REAL8
!
   implicit none
!
   integer :: iflag(NR_PHYSA)
!
!..workspace for element_error routine
   real(8) :: errorH,rnormH,errorE,rnormE
   real(8) :: errorV,rnormV,errorQ,rnormQ
   real(8) :: err,rnorm,error_subd,rnorm_subd
   integer :: iel,count,ierr
!
!----------------------------------------------------------------------
!
!..compute the error only for H1 field
   iflag(1:NR_PHYSA) = (/1,0/)
!
!..fetch active elements
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      if (RANK .eq. ROOT) then
         write(*,*) 'exact_error: mesh is distributed. ', &
                                 'computing error in parallel...'
      endif
   else
      if (RANK .ne. ROOT) goto 90
      write(*,*) 'exact_error: mesh is not distributed (or on host). ', &
                              'computing error on host...'
      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      NRELES_SUBD = NRELES
   endif
!
   error_subd = 0.d0; rnorm_subd = 0.d0
!
!$OMP PARALLEL DO                            &
!$OMP PRIVATE(errorH,errorE,errorV,errorQ,   &
!$OMP         rnormH,rnormE,rnormV,rnormQ)   &
!$OMP REDUCTION(+:error_subd,rnorm_subd)
   do iel=1,NRELES_SUBD
      call element_error(ELEM_SUBD(iel),iflag,           &
                         errorH,errorE,errorV,errorQ,    &
                         rnormH,rnormE,rnormV,rnormQ)
      error_subd = error_subd + errorH
      rnorm_subd = rnorm_subd + rnormH
   enddo
!$OMP END PARALLEL DO
!
   err = 0.d0; rnorm = 0.d0
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      count = 1
      call MPI_REDUCE(error_subd,err  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(rnorm_subd,rnorm,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
   else
      err   = error_subd
      rnorm = rnorm_subd
   endif
!
   if (RANK .eq. ROOT) then
      write(*,7020) NRDOF_TOT,sqrt(err),sqrt(rnorm)
 7020 format('exact_error: NRDOF_TOT, H1 ERROR AND NORM = ',i8,3x,2es12.5)
      write(*,7030) NRDOF_TOT,sqrt(err/rnorm)
 7030 format('exact_error: NRDOF_TOT, RELATIVE H1 ERROR = ',i8,3x,es12.5)
   endif
!
   90 continue
!
end subroutine exact_error
