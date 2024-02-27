#if HP3D_DEBUG

!-----------------------------------------------------------------------
!> Purpose : routine employs the greedy strategy to perform adaptive
!!           h-refinements
!!
!> @param[in] Idom  - target domain, if 0, examine all elements
!> @param[in] H     - control element size
!> rev@Dec 13
!-----------------------------------------------------------------------
subroutine adapt_mesh_size(Idom,Hmax)
  !
  use data_structure3D, only : NRELES
  !
  implicit none
  !
  integer, intent(in) :: Idom
  real(8), intent(in) :: Hmax
  !
  integer, allocatable :: nlist(:)
  integer :: mdle, iel, ic, istat, kref, ndom, nelts_prev
  real(8) :: h
  !
  ic = 1;
  do while (ic > 0)
     allocate(nlist(NRELES), STAT=istat)
     if (istat.ne.0) then
        write(*,*) 'adapt_mesh_size: allocation error', istat
        stop 1
     end if
     !
     mdle=0; ic=0
     do iel=1,NRELES
        call nelcon(mdle, mdle)
        call find_domain(mdle, ndom)
        if ((Idom.eq.0).or.(Idom.eq.ndom)) then
           call find_element_size(mdle, h)
           if (h.gt.Hmax) then
              ic = ic + 1
              nlist(ic) = mdle
           end if
        end if
     enddo
     !
     nelts_prev = NRELES
     do iel=1,ic
        mdle = nlist(iel)
        call get_isoref(mdle, kref)
        call refine(mdle, kref)
     end do
     !
     call close_mesh
     call update_gdof
     call update_ddof
     !
     write(*,*) '# of elements (before, after, inc) ', nelts_prev, NRELES, (NRELES-nelts_prev)
     !
     deallocate(nlist, STAT=istat)
     if (istat.ne.0) then
        write(*,*) 'adapt_mesh_size: deallocation error', istat
        stop 1
     end if
  end do
  !
end subroutine adapt_mesh_size

#endif
