!> Purpose : find out center coordinate of middle node
!! @param[in]  Mdle - middle node
!! @param[out] X    - center coordinate
subroutine find_element_size(Mdle, H)
  use data_structure3D
  implicit none
  integer, intent(in)  :: Mdle
  real(8), intent(out) :: H
  !
  integer :: nodesl(27),norientl(27)
  integer :: i, j, nv
  real(8) :: dist, x(3), y(3)
  !
  call elem_nodes(Mdle, nodesl,norientl)
  nv = nvert(NODES(Mdle)%ntype)
  !
  H = 0.d0
  do i=1,nv
     x = NODES(nodesl(i))%dof%coord(1:3,1)
     do j=1,nv
        y = NODES(nodesl(j))%dof%coord(1:3,1)
        call norm(x-y, dist)
        if (dist.gt.H) then
           H = dist
        end if
     end do
  end do
end subroutine find_element_size

subroutine find_element_size_minmax(Idom, Nelts, Hmin, Hmax)
  use data_structure3D
  implicit none
  integer, intent(in)    :: Idom
  integer, intent(out)   :: Nelts
  real(8), intent(inout) :: Hmin, Hmax
  !
  integer :: iel, mdle, ndom
  real(8) :: h
  !
  Nelts = 0
  mdle = 0; Hmin = 1e10; Hmax = 0.d0;
  do iel=1,NRELES
     call nelcon(mdle, mdle)
     call find_domain(mdle, ndom)
     if ((idom.eq.0).or.(idom.eq.ndom)) then
        call find_element_size(mdle, h)
        Nelts = Nelts + 1
        if (h.gt.Hmax) then
           Hmax = h
        end if
        if (h.lt.Hmin) then
           Hmin = h
        end if
     end if
  end do
end subroutine find_element_size_minmax
