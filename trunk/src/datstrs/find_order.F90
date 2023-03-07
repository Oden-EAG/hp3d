!> @brief find order of approximation for Mdle. It calls elem_nodes
!! @param[in]  Mdle   - middle node number
!! @param[out] Norder - order of approximation
!> @date Feb 2023
subroutine find_order(Mdle, Norder)
  use element_data
  use data_structure3D
  implicit none
  integer, intent(in)                 :: Mdle
  integer, intent(out), dimension(19) :: Norder
!
  integer, dimension(27) :: nodesl, norientl
  integer :: iprint
!
  iprint = 0
!
  Norder(1:19) = 0
  call elem_nodes(Mdle, nodesl,norientl)
  call find_order_from_list(NODES(Mdle)%Ntype,nodesl, Norder)
!
  if (iprint.eq.1) then
     write(*,7001) Mdle,S_Type(NODES(Mdle)%Ntype)
7001 format('find_order: Mdle = ',i8,' Mdle type = ',a5)
     write(*,7002) Norder(1:19)
7002 format('find_order: Norder = ',19i4)
  endif

end subroutine find_order
!
!--------------------------------------------------------------------------
!
!> @brief find order of approximation from the list of nodes
!! @param[in]  Ntype  - element type
!! @param[in]  Nodesl - middle node number
!! @param[out] Norder - order of approximation
!> @date Feb 2023
subroutine find_order_from_list(Ntype,Nodesl, Norder)
  use element_data
  use data_structure3D
  implicit none
  integer, intent(in)  :: Ntype
  integer, intent(in)  :: Nodesl(27)
  integer, intent(out) :: Norder(19)
  integer :: i, j, nrv, n
!
  nrv = NVERT(Ntype)
  n = NEDGE(Ntype) + NFACE(Ntype) + 1
!
  do i=1,n
     j = nrv+i
     Norder(i) = NODES(Nodesl(j))%order
  enddo
!
end subroutine find_order_from_list
!
!--------------------------------------------------------------------------
!
!> @brief find order of approximation in the element system of coords
!! @param[in]  Ntype        - element type
!! @param[in]  Norder       - order of approximation
!! @param[in]  Norient_face - orientation of element faces
!! @param[out] Norder_loc   - order of approx in element system of coords
!> @date Feb 2023
subroutine find_order_loc(Ntype,Norder,Norient_face, Norder_loc)
  use element_data, only: NFAXES
  use node_types
  implicit none
  integer, intent(in)  :: Ntype
  integer, intent(in)  :: Norder(19)
  integer, intent(in)  :: Norient_face(6)
  integer, intent(out) :: Norder_loc(19)
  integer :: j,nordh,nordv
!
  Norder_loc = Norder
!
  select case(Ntype)
  case(MDLB)
    do j=1,6
      if (NFAXES(3,Norient_face(j)).eq.1) then
        call decode(Norder(12+j), nordh,nordv)
        Norder_loc(12+j) = nordv*10+nordh
      endif
    enddo
  case(MDLP)
    do j=1,3
      if (NFAXES(3,Norient_face(2+j)).eq.1) then
        call decode(Norder(11+j), nordh,nordv)
        Norder_loc(11+j) = nordv*10+nordh
      endif
    enddo
  case(MDLD)
    if (NFAXES(3,Norient_face(1)).eq.1) then
      call decode(Norder(9), nordh,nordv)
      Norder_loc(9) = nordv*10+nordh
    endif
  end select
!
end subroutine find_order_loc
!
!--------------------------------------------------------------------------
!
!> @brief print order with format
!! @param[in]  Ntype  - element type
!! @param[in]  Norder - order of approximation
!> @date Feb 2023
subroutine print_order(Ntype,Norder)
  use element_data
  implicit none
  integer, intent(in) :: Ntype
  integer, intent(in) :: Norder(19)
  integer :: ii, nn
!
  ii=0; nn = nedge(Ntype)
  write(*,7001) Norder(ii+1:ii+nn)
7001 format('EDGE ORDERS = ',12i2)
!
  ii=ii+nn; nn = nface(Ntype)
  write(*,7002) Norder(ii+1:ii+nn)
7002 format('FACE ORDERS = ',6i3)
!
  ii=ii+nn+1
  write(*,7003) Norder(ii)
7003 format('MIDDLE NODE ORDER = ',i4)
!
end subroutine print_order
