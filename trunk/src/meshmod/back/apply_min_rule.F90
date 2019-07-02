!> Purpose : uniform refinement on the mesh
subroutine apply_p_rule(Is_min_rule)
  use data_structure3D
  implicit none
  ! ** Locals
  integer :: iprint
  iprint=0

  if (iprint.eq.1) then
     write(*,*) 'apply_min_rule: Begin'
  endif

  ! collect elements
  nr_elements_to_refine = NRELES
  if (iprint.eq.1) then
     write(*,*) 'global_href_default: nr_elements_to_refine ', &
          nr_elements_to_refine
  endif
  !
  allocate(list(nr_elements_to_refine),stat=istat)
  if (istat.ne.SUCCESS) then
    call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
  endif
  !
  mdle=0
  do i=1,NRELES
     call nelcon(mdle, mdle)
     call elem_nodes(mdle, nodesl, norientl)

     do j=1, n_nodes
        norder_nodesl(j) = NODES(nodesl(j))%order
     enddo

     call set_mdle_order_to_nodes(type, nord, norientl, norder_mdle)
     call set_min_order(norder_mdle, norder_nodesl)
     ! if node is no visited, norder_mdle is noder_nodesl
     ! then visit again, apply min rule
     do j=1, n_nodes
        NODES(nodesl(j))%order = norder_nodesl(j)
     enddo
  enddo

end subroutine apply_min_rule
