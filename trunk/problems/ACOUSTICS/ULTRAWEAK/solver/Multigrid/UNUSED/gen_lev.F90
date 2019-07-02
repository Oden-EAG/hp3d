!-----------------------------------------------------------------------------
!> Purpose : function determines the generation level for a node; the concept
!            is meaningful ONLY for isotropic refinements
!!
!! @param[in]  Nod - a node
!! @param[out] Gen_lev - the generation level for the node
!!
!! @revision Nov 17
!-----------------------------------------------------------------------------
!
      integer function Gen_lev(Nod)
!
      use data_structure3D
      implicit none
      integer, intent(in)  :: Nod
!
      integer :: igen, nfath
!-----------------------------------------------------------------------------
!
!  ...initiate the generation level
      igen=0
      nfath = NODES(Nod)%father
      do while(nfath.gt.0)
        nfath = NODES(nfath)%father
        igen = igen+1
      enddo
!
      Gen_lev = igen
!
!
      end function Gen_lev
