!--------------------------------------------------------------------
!> Purpose : routine performs a global isotropic h-refinement
!
!> @date Aug 2019
!--------------------------------------------------------------------
subroutine global_href
!
   use error
   use data_structure3D
   use environment , only : QUIET_MODE
   use mpi_param   , only : RANK,ROOT
!
   implicit none
!
   integer :: mdle_list(NRELES)
   integer :: nr_elements_to_refine,mdle,i,kref,ierr
!
!..collect elements
   nr_elements_to_refine = NRELES
!
   mdle=0
   do i=1,NRELES
      call nelcon(mdle, mdle)
      mdle_list(i) = mdle
   enddo
!
!..break the elements
   do i=1,nr_elements_to_refine
      mdle = mdle_list(i)
      if (is_leaf(mdle)) then
         call get_isoref(mdle, kref)
         call break(mdle,kref)
      endif
   enddo
!
   call refresh
!
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,100) ' # of elements broken         = ', nr_elements_to_refine
      write(*,200) ' # of current elements, nodes = ', NRELES, NRNODS
   endif
 100 format(A,I8)
 200 format(A,I8,', ',I9)
!
end subroutine global_href
!
!
!--------------------------------------------------------------------
!> Purpose : routine performs a global anisotropic h-refinement
!            note: works only for meshes with hexas/prisms,
!                  and breaking along the same axis
!                  (no mesh irregularity check are done here)
!
!> Arguments:
!     in:    Kref -   1 (break in z)  - hexa/prism
!                    10 (break in y)  - hexa only mesh
!                   100 (break in x)  - hexa only mesh
!                    11 (break in yz) - hexa only mesh
!                   101 (break in xz) - hexa only mesh
!                   110 (break in xy) - hexa/prism
!
!> @date Aug 2019
!--------------------------------------------------------------------
subroutine global_href_aniso(Kref)
!
   use error
   use data_structure3D
   use environment , only : QUIET_MODE
   use mpi_param   , only : RANK,ROOT
!
   implicit none
!
   integer, intent(in) :: Kref
!
   integer :: mdle_list(NRELES)
   integer :: nr_elements_to_refine,mdle,i,ierr
   integer :: kref_mdlb,kref_mdlp
!
!..check if valid refinement
   select case(Kref)
      case(1,110)
         kref_mdlb = Kref
         kref_mdlp = MOD(Kref,100)
      case(10,100,11,101)
         kref_mdlb = Kref
         kref_mdlp = 0
      case default
         write(*,*) 'global_href_aniso: invalid Kref. returning...'
         return
   end select
!
!..collect elements
   nr_elements_to_refine = NRELES
!
   mdle=0
   do i=1,NRELES
      call nelcon(mdle, mdle)
      mdle_list(i) = mdle
   enddo
!
!..break the elements
   do i=1,nr_elements_to_refine
      mdle = mdle_list(i)
      if (is_leaf(mdle)) then
         select case(NODES(mdle)%type)
            case('mdlb')
               call break(mdle,kref_mdlb)
            case('mdlp')
               if (kref_mdlp > 0) then
                  call break(mdle,kref_mdlp)
               else
                  write(*,*) 'global_href_aniso: invalid prism Kref. stop.'
                  stop
               endif
            case default
               write(*,*) 'global_href_aniso: unexpected node type. stop.'
               stop
         end select
      endif
   enddo
!
   call refresh
!
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,100) ' # of elements broken         = ', nr_elements_to_refine
      write(*,200) ' # of current elements, nodes = ', NRELES, NRNODS
   endif
 100 format(A,I8)
 200 format(A,I8,', ',I9)
!
end subroutine global_href_aniso
