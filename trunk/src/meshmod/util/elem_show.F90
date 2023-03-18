#if DEBUG_MODE

!> @brief     print nodal connectivity
!!
!! @param[in] Mdle - middle node
!!
!> @date      Feb 2023
subroutine elem_show_var1(Mdle)
   use data_structure3D
   implicit none
   integer, intent(in) :: Mdle
   integer :: nodesl(27),norientl(27)
   integer :: ntype
!
   nodesl = 0; norientl = 0
   if (is_root(Mdle)) then
      call elem_dump(Mdle, nodesl,norientl)
   else
      call elem_nodes(Mdle, nodesl,norientl)
   endif
!
   ntype = NODES(Mdle)%ntype
   select case(ntype)
   case(MDLB)
      write(*,7106) nodesl(1:27)
      write(*,7106) norientl(1:27)
   case(MDLN)
      write(*,7103) nodesl(1:15)
      write(*,7103) norientl(1:15)
   case(MDLP)
      write(*,7104) nodesl(1:21)
      write(*,7104) norientl(1:21)
   case(MDLD)
      write(*,7105) nodesl(1:19)
      write(*,7105) norientl(1:19)
   case default
      write(*,*) 'Error! elem_show: Mdle,type = ', Mdle, S_Type(ntype)
      stop 1
   end select
!
 7103 format(4i6,2x,6i6,2x,4i6,2x,i6)
 7104 format(6i6,2x,9i6,2x,2i6,2x,3i6,2x,i6)
 7105 format(5i6,2x,8i6,2x,i6,2x,4i6,2x,i6)
 7106 format(8i6,2x,12i6,2x,6i6,2x,i6)
!
end subroutine elem_show_var1
!
!----------------------------------------------------------------------
!
!> @brief     print nodal connectivity
!!
!! @param[in] Mdle     - middle node
!! @param[in] Ntype    - middle type
!! @param[in] Nodesl   - nodes
!! @param[in] Norientl - orientation
!!
!> @date      Feb 2023
subroutine elem_show_var2(Mdle,Ntype,Nodesl,Norientl)
   use node_types
   implicit none
   integer, intent(in) :: Mdle,Ntype
   integer, intent(in) :: Nodesl(27),Norientl(27)
!
   select case(Ntype)
   case(BRIC,MDLB)
      write(*,7106) nodesl(1:27)
      write(*,7106) norientl(1:27)
   case(TETR,MDLN)
      write(*,7103) nodesl(1:15)
      write(*,7103) norientl(1:15)
   case(PRIS,MDLP)
      write(*,7104) nodesl(1:21)
      write(*,7104) norientl(1:21)
   case(PYRA,MDLD)
      write(*,7105) nodesl(1:19)
      write(*,7105) norientl(1:19)
   case default
      write(*,*) 'Error! elem_show: Mdle,type = ', Mdle, S_Type(ntype)
      stop 1
   end select
!
 7103 format(4i6,2x,6i6,2x,4i6,2x,i6)
 7104 format(6i6,2x,9i6,2x,2i6,2x,3i6,2x,i6)
 7105 format(5i6,2x,8i6,2x,i6,2x,4i6,2x,i6)
 7106 format(8i6,2x,12i6,2x,6i6,2x,i6)
!
end subroutine elem_show_var2

#endif
