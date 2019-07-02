!> Purpose : uniform refinement on the mesh
subroutine global_href
  use error
  use data_structure3D
  use environment , only : QUIET_MODE
  implicit none
  ! ** Locals
  integer, allocatable, dimension(:)   :: list
  integer :: iprint, nr_elements_to_refine, mdle, i, kref, istat
  
  iprint=0

  if (iprint.eq.1) then
     write(*,*) 'global_href: Begin'
  endif

  ! collect elements
  nr_elements_to_refine = NRELES
  if (iprint.eq.1) then
     write(*,*) 'global_href_default: nr_elements_to_refine ', &
          nr_elements_to_refine
  endif
  !
  allocate(list(nr_elements_to_refine),stat=istat)
  if (istat.ne.SUCCESS) then
    call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
  endif
  !
  mdle=0
  do i=1,NRELES
     call nelcon(mdle, mdle)
     list(i) = mdle
  enddo
  !
  ! break the elements
  do i=1,nr_elements_to_refine
     mdle = list(i)
     if (is_leaf(mdle)) then
        call get_isoref(mdle, kref)
        if (iprint.eq.1) then
           write(*,*) 'global_href_default: mdle, kref ', mdle, kref
        endif
        call break(mdle,kref)
     endif
  enddo
  !
  deallocate(list, stat=istat)
  !
  if (istat.ne.SUCCESS) then
    call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
  endif
  ! 
  call refresh
IF (.NOT. QUIET_MODE) write(*,*) '# of elements broken         = ', nr_elements_to_refine
IF (.NOT. QUIET_MODE) write(*,*) '# of current elements, nodes = ', NRELES, NRNODS
!
end subroutine global_href
