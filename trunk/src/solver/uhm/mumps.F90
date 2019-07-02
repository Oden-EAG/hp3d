!--------------------------------------------------------  
!> Purpose : analyze the sparse system
subroutine uhm_analyze_mumps_1
  use uhm
  implicit none

  real*8 :: &
       t_base, t_after, t_analyze

  ! perform multi-level bisection
  call uhm_timer        (t_base)
  call uhm_build_tree   (UHM_MESH)
  call uhm_timer        (t_after)
  t_analyze = t_after - t_base

  call uhm_lock         (UHM_MESH)

  call uhm_create_matrix_without_buffer(UHM_MESH, UHM_DATATYPE, UHM_N_RHS)
  call uhm_create_leaf_matrix_buffer(UHM_MESH)

  write(*,7000) t_analyze
7000 format('Time analysis UHM      = ', f12.5, ' SEC')

end subroutine uhm_analyze_mumps_1

!--------------------------------------------------------  
!> Purpose : analyze the sparse system
subroutine uhm_analyze_mumps_2
  use uhm
  use control, only: SYMMETRIC
  implicit none
  
  integer :: n_thread, i
  real*8 :: &
       t_base, t_after, t_export, t_analyze

  ! ** set environement
  call uhm_mumps_set_show_n_rhs( UHM_MUMPS, 30);

  
  ! host involve calculation
  call uhm_mumps_set_par(UHM_MUMPS, 1)
  ! 0-unsym, 1-spd, 2-general sym
  call uhm_mumps_set_sym(UHM_MUMPS, UHM_SYM_MUMPS(UHM_ISYM_FLAG))
  call uhm_mumps_set_comm(UHM_MUMPS, -987654)
  
  ! ** initialization
  call uhm_mumps_init(UHM_MUMPS)

  call uhm_timer        (t_base)
  call uhm_mesh_export_matrix_mumps( &
       UHM_MESH, UHM_MUMPS, UHM_N_RHS, UHM_ISYM_FLAG.eq.SYMMETRIC )
  call uhm_timer        (t_after)
  t_export = t_after - t_base

  write(*,7001) t_export
7001 format('Time export to MUMPS = ', f12.5, ' SEC')

  ! ** setting control parameters
  !call uhm_mumps_set_icntl(UHM_MUMPS,1, 6); ! stream for error message : default 6
  !call uhm_mumps_set_icntl(UHM_MUMPS,2, 0); ! stream for warning message : default 0
  !call uhm_mumps_set_icntl(UHM_MUMPS,3, 6); ! stream for output message : default 0
  !call uhm_mumps_set_icntl(UHM_MUMPS,4, 3); ! error level - only message printed
  call uhm_mumps_set_icntl(UHM_MUMPS,5, 0); ! 0-assembled, 1-elemental
  !call uhm_mumps_set_icntl(UHM_MUMPS,6, 7); ! permutation : default  7 automatic
  call uhm_mumps_set_icntl(UHM_MUMPS,7, 5);
  ! ordering : 0-AMD, 1-user, 2-AMF, 3-SCOTCH, 4-PORD, 5-METIS, 6-QAMD, 7-auto

  !call uhm_mumps_set_icntl(UHM_MUMPS,8, 77); ! default scaling           
  !call uhm_mumps_set_icntl(UHM_MUMPS,9, 1);  ! solve Ax = b, otherwise A^T x = b

  !call uhm_mumps_set_icntl(UHM_MUMPS,13,0); ! ScaLAPACK used in root
  !call uhm_mumps_set_icntl(UHM_MUMPS,14,20); ! increment of workspace
  !call uhm_mumps_set_icntl(UHM_MUMPS,18,0); ! input matrix centralized on host

  !call uhm_mumps_set_cntl(UHM_MUMPS,1, 0.01); ! threshold pivot default 0.01


  ! ** analyze
  call uhm_timer        (t_base)
  call uhm_mumps_analyze(UHM_MUMPS)
  call uhm_timer        (t_after)
  t_analyze = t_after - t_base

  write(*,7002) t_analyze
7002 format('Time analysis MUMPS  = ', f12.5, ' SEC')

end subroutine uhm_analyze_mumps_2
