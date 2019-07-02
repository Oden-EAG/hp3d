!> Purpose : print nodal connectivity
!! @param[in] Mdle - middle node
subroutine elem_show_var1(Mdle)
  use data_structure3D
  implicit none
  integer, intent(in)    :: Mdle
  integer, dimension(27) :: nodesl, norientl

  nodesl = 0; norientl = 0
  if (is_root(Mdle)) then
     call elem_dump(Mdle, nodesl, norientl)
  else
     call elem_nodes(Mdle, nodesl, norientl)
  endif

  write(*,*) 'elem_show: Mdle = ', Mdle, NODES(Mdle)%type
  select case(NODES(Mdle)%type)
  case('tetr','mdln')
     write(*,7103) nodesl(1:15)
     write(*,7103) norientl(1:15)
  case('pris','mdlp')
     write(*,7104) nodesl(1:21)
     write(*,7104) norientl(1:21)
  case('pyra','mdld')
     write(*,7105) nodesl(1:19)
     write(*,7105) norientl(1:19)
  case('bric','mdlb')
     write(*,7106) nodesl(1:27)
     write(*,7106) norientl(1:27)
  case default
     print *, "elem_show: unknown node type =", NODES(mdle)%type
     stop 1
  end select

7103 format(4i6,2x,6i6,2x,4i6,2x,i6)
7104 format(6i6,2x,9i6,2x,2i6,2x,3i6,2x,i6)
7105 format(5i6,2x,8i6,2x,i6,2x,4i6,2x,i6)
7106 format(8i6,2x,12i6,2x,6i6,2x,i6)
end subroutine elem_show_var1

!> Purpose : print nodal connectivity
!! @param[in] Mdle     - middle node
!! @param[in] Type     - middle type
!! @param[in] Nodesl   - nodes
!! @param[in] Norientl - orientation
subroutine elem_show_var2(Mdle, Type, Nodesl, Norientl)
  implicit none
  character(len=4),       intent(in) :: Type
  integer,                intent(in) :: Mdle
  integer, dimension(27), intent(in) :: Nodesl, Norientl

  write(*,*) 'elem_show: Mdle = ', Mdle, Type
  select case(Type)
  case('tetr','mdln')
     write(*,7103) nodesl(1:15)
     write(*,7103) norientl(1:15)
  case('pris','mdlp')
     write(*,7104) nodesl(1:21)
     write(*,7104) norientl(1:21)
  case('pyra','mdld')
     write(*,7105) nodesl(1:19)
     write(*,7105) norientl(1:19)
  case('bric','mdlb')
     write(*,7106) nodesl(1:27)
     write(*,7106) norientl(1:27)
  case default
     print *, "elem_show: unknown node type =", type
     stop 1
  end select

7103 format(4i6,2x,6i6,2x,4i6,2x,i6)
7104 format(6i6,2x,9i6,2x,2i6,2x,3i6,2x,i6)
7105 format(5i6,2x,8i6,2x,i6,2x,4i6,2x,i6)
7106 format(8i6,2x,12i6,2x,6i6,2x,i6)
end subroutine elem_show_var2

