!--------------------------------------------------------------------
!
!     routine name      - random_refineDPG
!
!--------------------------------------------------------------------
!
!     latest revision:  - Feb 18
!
!     purpose:          - performs random refinements using Pseudo-random
!                         number generator functions of Fortran
!
!---------------------------------------------------------------------
!
   subroutine random_refineDPG
!
   use data_structure3D, ONLY: NRELES, NODES
!
   IMPLICIT NONE
!
   integer, allocatable :: mdle_list(:), mdle_ref_list(:)
   integer :: mdle, iel, kref
   real*8  :: x
   integer :: i, j, n, loc, nr_elem_to_refine
   integer, allocatable :: seed(:)
!
!---------------------------------------------------------------------
!
!..create the list of current elements
   allocate(mdle_list(NRELES))
!
   mdle = 0
   do iel = 1, NRELES
      call nelcon(mdle,mdle)
      mdle_list(iel) = mdle
   enddo
!
   nr_elem_to_refine = NRELES/4
   ! nr_elem_to_refine = 1

   allocate(mdle_ref_list(nr_elem_to_refine)) ; mdle_ref_list = 0

   i = 0
   do while (i .lt. nr_elem_to_refine)
!
      call random_number(x)
      j = 1 + int(x*NRELES)
      call locate(mdle_list(j), mdle_ref_list, i, loc)
      if (loc .eq. 0) then 
         i = i+1
         mdle_ref_list(i) = mdle_list(j)
      endif   
   enddo
!
!..refine the elements from the list
   do iel=1,nr_elem_to_refine
      mdle = mdle_ref_list(iel)
      if (NODES(mdle)%type .ne. 'mdlb') then
         write(*,*) 'random_refineDPG: type = ', NODES(mdle)%type
         stop 1
      endif
      kref = 111
!
      call refine(mdle,kref)
   enddo
!
!..close the mesh
   call close
   call enforce_max_rule
   call update_gdof
   call update_ddof

   deallocate(mdle_ref_list)
   deallocate(mdle_list)

   end subroutine random_refineDPG


   subroutine manual_pref

   use data_structure3D, ONLY: NRELES, NODES

!..manually raise the order of a face and one edge


   call nodmod(11,4)
   call nodmod(25,43)




   end subroutine manual_pref
