!------------------------------------------------------------------------
subroutine split_interface
!------------------------------------------------------------------------
! LATEST REVISION: Jul 09
!
! PURPOSE: routine generates twin planes.
!
! REMARKS: original plane must NOT be removed!
!------------------------------------------------------------------------
! MODULES
  use SPLIT_SURF
  use GMP
!------------------------------------------------------------------------
IMPLICIT NONE
!------------------------------------------------------------------------
! VARIABLES
  real*8  :: s
  integer :: i,ns
  integer :: status
!------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'split_interface: splitting interface...'
#endif
! ..add one subdomain
    NRDOMAIN = NRDOMAIN + 1
#if I_PRINT >= 2
    write(*,1) NRDOMAIN
1   format(' split_interface: adding subdomain, NRDOMAINS = ',i2)
#endif
! ..check if split surface is a plane
    if (SURFACES(Nsplit_MOD)%Type .ne. 'VecPt') then
      write(*,*) 'split_interface: CANNOT SPLIT ns = ',ns
      stop
    endif
! ..determine normal of split plane
    call norm(SURFACES(Nsplit_MOD)%Rdata(4:6), s)
    NORMAL(1:3) = SURFACES(Nsplit_MOD)%Rdata(4:6)/s
#if I_PRINT >= 2
    write(*,*)'split_interface: NORMAL = ',NORMAL(1:3)
#endif
! ..add 2 planes to GMP data structure
    do i = 1, 2
      NRSURFS = NRSURFS + 1
      if (NRSURFS .gt. MAXSU) then
        write(*,*) 'split_interface: increase MAXSU.'
        stop
      endif
      SURFACES(NRSURFS)%Type = 'VecPt'
      allocate(SURFACES(NRSURFS)%Rdata(6), STAT = status)
      if (status .ne. 0 ) then
        write(*,*)'split_interface: Rdata not allocated!'
        stop
      endif
      select case(i)
! ......plane on NEG side of split plane
        case(1)
          SURFACES(NRSURFS)%Rdata(1:3) = SURFACES(Nsplit_MOD)%Rdata(1:3) - NORMAL(1:3)*Dh1_MOD
          SURFACES(NRSURFS)%Rdata(4:6) = NORMAL(1:3)
! ......plane on POS side of split plane
        case(2)
          SURFACES(NRSURFS)%Rdata(1:3) = SURFACES(Nsplit_MOD)%Rdata(1:3) + NORMAL(1:3)*Dh2_MOD
          SURFACES(NRSURFS)%Rdata(4:6) = NORMAL(1:3)
      end select
    enddo
#if I_PRINT >= 2
    write(*,*)'************************************************************************'
    write(*,*)'split_interface: 2 planes added to GMP data structure.'
    write(*,*)'plane 1 = ', (NRSURFS - 1)
    write(*,*)'plane 2 = ', NRSURFS
    call print_GMP
#endif
#if I_PRINT >= 1
    write(*,*)'split_interface: done!'
#endif
!
end subroutine split_interface
