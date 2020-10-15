!--------------------------------------------------------------
!> Purpose : computes norm of the error and the exact solution
!!
!! @param[out]  Err   - norm of the error
!! @param[out]  Rnorm - norm of the exact solution
!--------------------------------------------------------------
subroutine exact_error
  use data_structure3D
  use common_prob_data
  use environment, only : L2PROJ
  ! use assembly_sc, only: NRDOF_TOT,NRDOF_CON
  use par_mesh   , only: DISTRIBUTED,HOST_MESH
  use mpi_param
  use MPI        , only: MPI_BCAST,MPI_COMM_WORLD,MPI_INTEGER
  implicit none
!
  integer, dimension(NR_PHYSA) :: iflag
  integer :: src,count,ierr
!------------------------------------------------------------

  select case(IERROR_PROB)
  case(IERROR_L2)
    L2PROJ = .TRUE.
  case(IERROR_NATURAL)
    L2PROJ = .FALSE.
  case default
    write(*,*) 'exact_error : Error calculation type not supported'
  end select
!
  if (RANK .eq. ROOT) then
 333  write(*,7000)
 7000 format('Declare the attribute to calculate the error of: ',  &
             '   1)Displacement, 2)Stress, 3)Combined')
  read(*,*) IERROR_ATTR
  endif
  if (NUM_PROCS .gt. 1) then
      count = 1; src = ROOT
      call MPI_BCAST (IERROR_ATTR,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
   endif
!
  iflag = 0
  select case(IERROR_ATTR)
  case(DISPLACEMENT)
    iflag(1:2) = (/1,0/)
  case(STRESS)
    iflag(1:2) = (/0,1/)
  case(COMBINED)
    iflag(1:2) = (/1,1/)
  end select
!
!..fetch active elements
  if (DISTRIBUTED .and. (.not. HOST_MESH)) then
    if (RANK .eq. ROOT) then
       write(*,*) 'exact_error: mesh is distributed. computing error in parallel...'
    endif
  else
    if (RANK .ne. ROOT) goto 90
    write(*,*) 'exact_error: mesh is not distributed (or on host). computing error on host...'
    ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
    NRELES_SUBD = NRELES
  endif
!
  call compute_error(iflag,IERROR_ATTR)
!
  90 continue
!
end subroutine exact_error
