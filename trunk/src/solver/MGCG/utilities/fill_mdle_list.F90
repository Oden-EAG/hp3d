!----------------------------------------------------------------------
!
!   routine name       - fill_mdle_list
!
!----------------------------------------------------------------------
!
!   latest revision    - Apr 18
!
!   purpose            - file element mdle nod list for a grid
!
!   arguments
!       in
!              igrid   - grid index
!
!----------------------------------------------------------------------
!
   subroutine fill_mdle_list(igrid)
!
   use data_structure3D,  only: NRELES
   use mg_data_structure, only: GRID, NODES_MG
!
   implicit none
!
   integer :: igrid, mdle, iel
!
!-------------------------------------------------
!
   GRID(igrid)%nreles = NRELES
   allocate(GRID(igrid)%mdlel(NRELES))
   mdle = 0
   do iel = 1, GRID(igrid)%nreles
      call nelcon(mdle,mdle)
      GRID(igrid)%mdlel(iel) = mdle
      NODES_MG(mdle)%iel(igrid) = iel
   enddo

   end subroutine fill_mdle_list
