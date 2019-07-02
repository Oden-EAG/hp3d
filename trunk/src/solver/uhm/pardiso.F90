!--------------------------------------------------------  
!> Purpose : analyze the sparse system
subroutine uhm_analyze_pardiso_1
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

end subroutine uhm_analyze_pardiso_1

!--------------------------------------------------------  
!> Purpose : analyze the sparse system
subroutine uhm_analyze_pardiso_2
  use uhm
  use control, only: SYMMETRIC
  implicit none
  
  integer :: n_thread, i
  real*8 :: &
       t_base, t_after, t_export, t_analyze

  call uhm_timer        (t_base)

  call uhm_mesh_export_matrix_pardiso( &
       UHM_MESH, UHM_PARDISO, UHM_N_RHS, UHM_ISYM_FLAG.eq.SYMMETRIC )
  call uhm_timer        (t_after)
  t_export = t_after - t_base

  write(*,7001) t_export
7001 format('Time export to PARDISO = ', f12.5, ' SEC')

  do i=1,64
     call uhm_pardiso_set_iparm(UHM_PARDISO, i, 0)
  end do

  call uhm_get_num_threads     (n_thread)
  
  call uhm_pardiso_set_iparm(UHM_PARDISO, 1, 1)
  call uhm_pardiso_set_iparm(UHM_PARDISO, 2, 2)
  call uhm_pardiso_set_iparm(UHM_PARDISO, 3, n_thread)
  call uhm_pardiso_set_iparm(UHM_PARDISO, 7,0)
  call uhm_pardiso_set_iparm(UHM_PARDISO, 10,13)
  call uhm_pardiso_set_iparm(UHM_PARDISO, 11, 1)
  call uhm_pardiso_set_iparm(UHM_PARDISO, 13, 1)
  call uhm_pardiso_set_iparm(UHM_PARDISO, 21, 1)
  call uhm_pardiso_set_iparm(UHM_PARDISO, 24, 1)
  call uhm_pardiso_set_iparm(UHM_PARDISO, 25, 1)
  call uhm_pardiso_set_iparm(UHM_PARDISO, 32, 0) ! 0-direct, 1-iterative 
  call uhm_pardiso_set_iparm(UHM_PARDISO, 52, 1)

  call uhm_timer          (t_base)
  call uhm_pardiso_init   (UHM_PARDISO)
  call uhm_pardiso_analyze(UHM_PARDISO)
  call uhm_timer          (t_after)
  t_analyze = t_after - t_base

  write(*,7002) t_analyze
7002 format('Time analysis PARDISO  = ', f12.5, ' SEC')

end subroutine uhm_analyze_pardiso_2
