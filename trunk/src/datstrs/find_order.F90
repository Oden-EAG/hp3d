!> Purpose : find order of approximation for Mdle. It calls elem_nodes
!! @param[in]  Mdle   - middle node number
!! @param[out] Norder - order of approximation
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
  call find_order_from_list(NODES(Mdle)%Type,nodesl, Norder)
!
  if (iprint.eq.1) then
     write(*,7001) Mdle,NODES(Mdle)%type
7001 format('find_order: Mdle = ',i8,' Mdle TYPE = ',a5)
     write(*,7002) Norder(1:19)
7002 format('find_order: Norder = ',19i4)
  endif

end subroutine find_order
!
!--------------------------------------------------------------------------
!
!> Purpose : find order of approximation from the list of nodes
!! @param[in]  Nodesl - middle node number
!! @param[out] Norder - order of approximation
subroutine find_order_from_list(Type,Nodesl, Norder)
  use element_data
  use data_structure3D
  implicit none
  character(len=4), intent(in)                 :: Type
  integer,          intent(in),  dimension(27) :: Nodesl
  integer,          intent(out), dimension(19) :: Norder
  integer :: i, j, nrv, n
!
  nrv = NVERT(Type)
  n = NEDGE(Type) + NFACE(Type) + 1
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
!> Purpose : find order of approximation in the element system of coords
!! @param[in]  Type         - element type
!! @param[in]  Norder       - order of approximation
!! @param[in]  Norient_face - orientation of element faces
!! @param[out] Norder_loc   - order of approx in element system of coords
subroutine find_order_loc(Type,Norder,Norient_face, Norder_loc)
  use element_data, only: NFAXES
  implicit none
  character(len=4), intent(in)                 :: Type
  integer,          intent(in),  dimension(19) :: Norder
  integer,          intent(in),  dimension(6)  :: Norient_face
  integer,          intent(out), dimension(19) :: Norder_loc
  integer :: j,nordh,nordv
!
  Norder_loc = Norder
!
  select case(Type)
  case('mdlb')
    do j=1,6
      if (NFAXES(3,Norient_face(j)).eq.1) then
        call decode(Norder(12+j), nordh,nordv)
        Norder_loc(12+j) = nordv*10+nordh
      endif
    enddo
  case('mdlp')
    do j=1,3
      if (NFAXES(3,Norient_face(2+j)).eq.1) then
        call decode(Norder(11+j), nordh,nordv)
        Norder_loc(11+j) = nordv*10+nordh
      endif
    enddo
  case('mdld')
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
!> Purpose : print order with format
!! @param[in]  Type   - element type
!! @param[in]  Norder - order of approximation
subroutine print_order(Type, Norder)
  use element_data
  implicit none
  character(len=4), intent(in) :: Type
  integer, intent(in), dimension(19) :: Norder
  integer :: in, nn
!
  in=0; nn = nedge(Type)
  write(*,7001) Norder(in+1:in+nn)
7001 format('EDGE ORDERS = ',12i2)
!
  in=in+nn; nn = nface(Type)
  write(*,7002) Norder(in+1:in+nn)
7002 format('FACE ORDERS = ',6i3)
!
  in=in+nn+1
  write(*,7003) Norder(in)
7003 format('MIDDLE NODE ORDER = ',i4)
!
end subroutine print_order
