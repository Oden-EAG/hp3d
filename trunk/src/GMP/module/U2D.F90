!==============================================================================
module U2D                                                                    !
!==============================================================================
  use GMP                                                                     !
!------------------------------------------------------------------------------
! SHARED VARIABLES                                                            |
!------------------------------------------------------------------------------
  IMPLICIT NONE                                                               !
  SAVE                                                                        !
! surfaces a point is lying on                                                |
  integer, dimension(:,:), allocatable  :: POINT_TYPE                         !
!------------------------------------------------------------------------------
! SUBROUTINES                                                                 |
!------------------------------------------------------------------------------
    CONTAINS                                                                  !
!     allocate_EXTRA_PLANES                                                   |
!     deallocate_EXTRA_PLANES                                                 |
!     generate_POINT_TYPE                                                     |
!     deallocate_POINT_TYPE                                                   |
!==============================================================================
!
!
!
!
!
!------------------------------------------------------------------------------
subroutine allocate_EXTRA_PLANES
!------------------------------------------------------------------------------
  IMPLICIT NONE
  integer :: status
!------------------------------------------------------------------------------
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'allocate_EXTRA_PLANES: allocating planes...'
#endif
! ..1st plane
    NRSURFS = NRSURFS + 1
    allocate(SURFACES(NRSURFS)%Rdata(6), STAT = status)
    if (status .ne. 0) then
      write(*,*)'allocated_EXTRA_PLANES: Rdata not allocated for 1st plane.'
      stop
    endif
    SURFACES(NRSURFS)%Type = 'VecPt'
!
! ..2nd plane
    NRSURFS = NRSURFS + 1
    allocate(SURFACES(NRSURFS)%Rdata(6), STAT = status)
    if (status .ne. 0) then
      write(*,*)'allocated_EXTRA_PLANES: Rdata not allocated for 2nd plane.'
      stop
    endif
    SURFACES(NRSURFS)%Type = 'VecPt'
#if I_PRINT >= 1
    write(*,*)'allocate_EXTRA_PLANES: done!'
#endif
!
end subroutine allocate_EXTRA_PLANES
!-----------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------
subroutine deallocate_EXTRA_PLANES
!--------------------------------------------------------------------------------
  IMPLICIT NONE
  integer :: status
!--------------------------------------------------------------------------------
! ..1st plane
    deallocate(SURFACES(NRSURFS)%Rdata, STAT = status)
    if (status .ne. 0) then
      write(*,*)'deallocated_EXTRA_PLANES: Rdata not deallocated for 1st plane.'
      stop
    endif
    SURFACES(NRSURFS)%Type = 'Void'
    NRSURFS = NRSURFS - 1
!
! ..2nd plane
    deallocate(SURFACES(NRSURFS)%Rdata, STAT = status)
    if (status .ne. 0) then
      write(*,*)'deallocated_EXTRA_PLANES: Rdata not deallocated for 2nd plane.'
      stop
    endif
    SURFACES(NRSURFS)%Type = 'Void'
    NRSURFS = NRSURFS - 1
!
end subroutine deallocate_EXTRA_PLANES
!-----------------------------------------------------------------------------
!
!
!==============================================================================
!  subroutine - generate_POINT_TYPE
!==============================================================================
!  LATEST REVISION: Jul 09
!
!  PURPOSE: routine generates array POINT_TYPE;
!    POINT_TYPE(              1, np) = number of surfaces point np lies on
!    POINT_TYPE(2 : N_SURFS + 1, np) = list of surfaces
!==============================================================================
subroutine generate_POINT_TYPE
!------------------------------------------------------------------------------
  IMPLICIT NONE
!------------------------------------------------------------------------------
! LOCAL VARIABLES
  integer            :: nt,ns,iv,np,idec,i
  integer            :: status
!------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'generate_POINT_TYPE: generating POINT_TYPE...'
#endif
! ..allocate POINT_TYPE and initialize to 0
    allocate(POINT_TYPE(MAXSU + 1,NRPOINT), STAT = status)
    if (status .ne. 0) then
      write(*,*)'generate_POINT_TYPE: POINT_TYPE not allocated.'
      stop
    endif
    POINT_TYPE = 0
! ..loop through triangles
    do nt = 1, NRTRIAN
! ....cycle if triangle is not on an algebraic surface
      if (TRIANGLES(nt)%Type .ne. 'PTITri') cycle
#if I_PRINT >= 2
      write(*,*)'******************************************************************************'
      write(*,*)'generate_POINT_TYPE: visiting triangle nt = ',nt
#endif
! ....get surface number
      ns = TRIANGLES(nt)%Idata(1)
#if I_PRINT >= 2
      write(*,*)'generate_POINT_TYPE: ns = ',ns
#endif
! ....loop through triangle vertices
      do iv = 1, 3
        np = TRIANGLES(nt)%VertNo(iv)
! ......if vertex has not been visited
        if (POINT_TYPE(1,np) .eq. 0) then
#if I_PRINT >= 2
          write(*,*)'generate_POINT_TYPE: 1st visit to point np = ',np
#endif
          POINT_TYPE(1,np) = 1
          POINT_TYPE(2,np) = ns
#if I_PRINT >= 2
          write(*,*)'generate_POINT_TYPE: stored ns = ',POINT_TYPE(2,np)
#endif
! ......else vertex has been visited
        else
#if I_PRINT >= 2
          write(*,*)'generate_POINT_TYPE: revisiting point np = ',np
#endif
! ........determine if surface needs to be added
          idec = 0
          do i = 1, POINT_TYPE(1,np)
            if (POINT_TYPE(i + 1,np) .eq. ns) then
              idec = 1
              exit
            endif
          enddo
! ........if surface needs to be added
          if (idec .eq. 0) then
#if I_PRINT >= 2
            write(*,*)'generate_POINT_TYPE: new surface found.'
#endif
! ..........update number of surfaces
            POINT_TYPE(1,np) = POINT_TYPE(1,np) + 1
! ..........store surface number
            POINT_TYPE(POINT_TYPE(1,np) + 1,np) = ns
#if I_PRINT >= 2
            write(*,*)'generate_POINT_TYPE: surfaces = ',POINT_TYPE(2:(POINT_TYPE(1,np) + 1),np)
#endif
          endif
        endif
! ....end of loop through triangle vertices
      enddo
! ..end of loop through triangles
    enddo
#if I_PRINT >= 1
    write(*,*)'generate_POINT_TYPE: done!'
#endif
!
end subroutine generate_POINT_TYPE
!----------------------------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------------
subroutine deallocate_POINT_TYPE
!-----------------------------------------------------------------------------------
  IMPLICIT NONE
  integer :: status
!-----------------------------------------------------------------------------------
! printing flag (0,1)
#define I_PRINT 0
!
#if I_PRINT >= 1
     write(*,*)'deallocate_POINT_TYPE: deallocating POINT_TYPE.'
#endif
     deallocate(POINT_TYPE, STAT = status)
     if (status .ne. 0) then
       write(*,*)'deallocate_POINT_TYPE: POINT_TYPE not deallocated.'
       stop
     endif
#if I_PRINT >= 2
     write(*,*)'deallocate_POINT_TYPE: POINT_TYPE deallocated successfully.'
#endif
!
end subroutine deallocate_POINT_TYPE
!------------------------------------------------------------------------------------
!
end module U2D
