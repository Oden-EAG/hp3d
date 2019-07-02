!--------------------------------------------------------  
!> Purpose : analyze the sparse system
subroutine uhm_analyze_wsmp_1
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

end subroutine uhm_analyze_wsmp_1

!--------------------------------------------------------  
!> Purpose : analyze the sparse system
subroutine uhm_analyze_wsmp_2
  use uhm
  implicit none
  
  integer :: n_thread, i
  real*8 :: &
       t_base, t_after, t_export, t_analyze

  call uhm_wsmp_set_show_n_rhs( UHM_WSMP, 30);
  call uhm_timer        (t_base)
  call uhm_mesh_export_matrix_wsmp( &
       UHM_MESH, UHM_WSMP, UHM_N_RHS)
  call uhm_timer        (t_after)
  t_export = t_after - t_base

  write(*,7001) t_export
7001 format('Time export to WSMP = ', f12.5, ' SEC')

  do i=1,64
     call uhm_wsmp_set_iparm(UHM_WSMP, i, 0)
  end do

  call uhm_get_num_threads     (n_thread)

  call uhm_wsmp_set_iparm(UHM_WSMP,1, 1)
  call uhm_wsmp_set_iparm(UHM_WSMP,2, 2)
  call uhm_wsmp_set_iparm(UHM_WSMP,3, n_thread)
  call uhm_wsmp_set_iparm(UHM_WSMP,7,0)
  call uhm_wsmp_set_iparm(UHM_WSMP,10,13)
  call uhm_wsmp_set_iparm(UHM_WSMP,11, 1)
  call uhm_wsmp_set_iparm(UHM_WSMP,13, 1)
  call uhm_wsmp_set_iparm(UHM_WSMP,21, 1)
  call uhm_wsmp_set_iparm(UHM_WSMP,24, 1)
  call uhm_wsmp_set_iparm(UHM_WSMP,25, 1)
  call uhm_wsmp_set_iparm(UHM_WSMP,52, 1)

  call uhm_timer        (t_base)
  call uhm_wsmp_init(UHM_WSMP)
  call uhm_wsmp_analyze(UHM_WSMP)
  call uhm_timer        (t_after)
  t_analyze = t_after - t_base

  write(*,7002) t_analyze
7002 format('Time analysis WSMP  = ', f12.5, ' SEC')

end subroutine uhm_analyze_wsmp_2
