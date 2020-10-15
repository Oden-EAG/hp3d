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
!
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
             '   1)Displacement, 2)Stress, 3)Combined, 4) Symm. Lagrange mult.')
  read(*,*) IERROR_ATTR
!
  iflag = 0
  select case(IERROR_ATTR)
  case(DISPLACEMENT)
    iflag = (/0,0,1,0,0/)
  case(STRESS)
    iflag = (/0,0,0,1,0/)
  case(COMBINED)
    iflag = (/0,0,1,1,0/)
  case(LAGRANGE)
    iflag = (/0,0,0,0,1/)
  end select
  endif
  ! 
  if (NUM_PROCS .gt. 1) then
      count = 1; src = ROOT
      call MPI_BCAST (iflag,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
   endif
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
  call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
end subroutine exact_error
