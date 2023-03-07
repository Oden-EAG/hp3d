!> @brief  displays active elements
!> @date   Feb 2023
subroutine display_act_elem
!
   use data_structure3D , only : NRELES,NODES
   use node_types
!
   implicit none
   integer :: mdle,iel
!
   write(*,*)''
   write(*,*)'Active elements:'
!
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
!
      write(*,7000) mdle,S_Type(NODES(mdle)%ntype),NODES(mdle)%order
 7000 format(' mdle = ',i4,' ; type = ',a4,' ; order = ',i3)
   enddo
   write(*,*) ''
!
end subroutine display_act_elem
