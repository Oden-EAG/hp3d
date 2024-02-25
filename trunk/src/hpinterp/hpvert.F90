!----------------------------------------------------------------------
!
!   routine name       - hpvert
!
!----------------------------------------------------------------------
!> @brief      routine determines geometry dof for a vertex node
!!
!> @param[in]  Iflag  - a flag specifying which of the objects the
!!                      vertex node is generated on
!!                       = 5  for a prism
!!                       = 6  for a hexahedron
!!                       = 7  for a tetrahedron
!!                       = 8  for a pyramid
!> @param[in]  No     - number of a specific object
!!                     (point,curve,triangle,rectangle...)
!> @param[in]  Xi     - reference coordinates of the point
!> @param[out] Xnod   - geometry dof
!!
!> @date       Aug 2019
!----------------------------------------------------------------------
!
subroutine hpvert(Iflag,No,Xi, Xnod)
!
   implicit none
!
   integer, intent(in)  :: Iflag
   integer, intent(in)  :: No
   real(8), intent(in)  :: Xi(3)
   real(8), intent(out) :: Xnod(3)
!
   real(8) :: void(3,3)
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
#endif
!
!----------------------------------------------------------------------
!
!..determine coordinates of the point
   select case(Iflag)
      case(5); call prism(No,Xi, Xnod,void)
      case(6); call hexa( No,Xi, Xnod,void)
      case(7); call tetra(No,Xi, Xnod,void)
      case(8); call pyram(No,Xi, Xnod,void)
      case default
         write(*,*) 'hpvert: unknown node type =', Iflag
         stop 1
   end select
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,7005) Iflag,No,Xi,Xnod
 7005 format('hpvert: Iflag,No = ',2i8,' Xi = ',3f8.3, &
             '  Xnod = ',3f8.3)
      call pause
   endif
#endif
!
end subroutine hpvert
