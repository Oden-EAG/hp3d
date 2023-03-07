!> @brief   show the mesh data structure
!> @date    Feb 2023
subroutine neig_show
   use data_structure3D
   use refinements
   implicit none
   integer :: mdle,ntype,i
   integer :: neig(4,6)
!
   mdle = 0
   do i=1,NRELES
      call nelcon(mdle, mdle)
      call find_neig(mdle, neig)
      ntype = NODES(mdle)%ntype
      write(*,*) 'neig_show: mdle, type, neig= ', &
          mdle, S_Type(ntype), neig(1,1:nface(ntype))
   enddo
!
end subroutine neig_show
