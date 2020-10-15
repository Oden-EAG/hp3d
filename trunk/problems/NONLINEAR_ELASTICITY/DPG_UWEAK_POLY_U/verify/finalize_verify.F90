subroutine finalize_verify
  use data_structure3D_poly
  use geometry_polydpg
  !
  call dealloc_physics
  call deallocds
  call dealloc_geom_poly
  !
end subroutine finalize_verify