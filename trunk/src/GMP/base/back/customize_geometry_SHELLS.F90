!-------------------------------------------------------------------------------
subroutine customize_geometry_shells
!-------------------------------------------------------------------------------
  use kinds
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! PARAMETERS FOR add_layer2plane
  integer,                       parameter :: n_plane   = 4
  integer,                       parameter :: nr_bound  = 4
  integer, dimension(nr_bound),  parameter :: ns_bound  = (/1, 2, 3, 5/)
  integer,                       parameter :: n_layers  = 2
  real(DP), dimension(n_layers), parameter :: thickness = (/0.2d0, 0.2d0/)
  integer, dimension(n_layers),  parameter :: n_domains = (/2, 2/)
!-------------------------------------------------------------------------------
! priting flag (0,1)
#define I_PRINT 1
!
#if I_PRINT >= 1
   write(*,*)'customize_geometry: customizing geometry for SHELLS.'
#endif
   call add_layer2plane(n_plane,nr_bound,ns_bound,n_layers,thickness,n_domains)
#if I_PRINT >= 1
   write(*,*)'customize_geometry: done!'
#endif
!
end subroutine customize_geometry_shells
!-------------------------------------------------------------------------------
