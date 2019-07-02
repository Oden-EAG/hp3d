!--------------------------------------------------------------------------------
subroutine create_TW_POINT_AND_CONN_CURVE(Np)
!--------------------------------------------------------------------------------
! LATEST REVISION: Jul 09
!
! PURPOSE: given a point, routine:
!   1. updates new_point(np) = twin_np;
!   2. if necessary:
!      a. generates twin point;
!      b. moves original point;
!      c. generates curve connecting twin points.
!--------------------------------------------------------------------------------
! MODULES
  use GMP
  use SPLIT_SURF
  use U2D
!--------------------------------------------------------------------------------
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! DUMMY ARGUMENTS
  integer, intent(in) :: Np
!-------------------------------------------------------------------------------- 
! VARIABLES
  real*8, dimension(3)    :: xp
  integer                   :: nrsrf  
  integer, dimension(MAXSU) :: SURFS
  integer                   :: status
!--------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'create_TW_POINT_AND_CONN_CURVE: Np = ',Np
#endif  
! ..if point does not need to be duplicated, twin point is point itself
    if (new_point(np) .eq. 0) then
      new_point(np) = np      
#if I_PRINT >= 1
      write(*,*)'create_TW_POINT_AND_CONN_CURVE: no need to duplicate point!'
#endif
      return
    endif  
! *************  MOVE ORIGINAL POINT  *******************************************
    xp(1:3) = POINTS(Np)%Rdata(1:3)
! ..determine surfaces point must conform to
    call give_surf(xp,Nr_confm_MOD,Ns_confm_MOD, nrsrf,SURFS)
! ..add plane shifted on NEG side to conforming surfaces
    nrsrf = nrsrf + 1
    SURFS(nrsrf) = NRSURFS - 3
    call reproject_point(Np,nrsrf,SURFS)
!  
! *************  GENERATE TWIN POINT  *******************************************
    NRPOINT = NRPOINT + 1
    if (NRPOINT .gt. MAXNP) then
      write(*,*) 'create_TW_POINT_AND_CONN_CURVE: increase MAXNP.'
      stop
    endif
    POINTS(NRPOINT)%Type  = 'Regular'
    allocate(POINTS(NRPOINT)%Rdata(1:3), STAT = status)
    if (status .ne. 0) then
      write(*,*)'create_TW_POINT_AND_CONN_CURVE: Rdata not allocated for Np = ',Np
      stop
    endif
    POINTS(NRPOINT)%Rdata(1:3) = xp(1:3)
! ..add plane shifted on POS side to conforming surfaces
    SURFS(nrsrf) = NRSURFS - 2
    call reproject_point(NRPOINT,nrsrf,SURFS)
#if I_PRINT >= 2   
    write(*,*)'**************  Np = ',Np
    write(*,*)'**  conf. surfaces = ',SURFS(1:nrsrf - 1)
    write(*,*)'*****  orig. point = ',xp(1:3)
    write(*,*)'**  point NEG side = ',POINTS(Np)%Rdata(1:3)
    write(*,*)'**  point POS side = ',POINTS(NRPOINT)%Rdata(1:3)
#endif 
! ..store new point to its old point connectivities
    new_point(Np) = NRPOINT
#if I_PRINT >= 2
    write(*,*)'************  twin point generated.'
#endif 
!
! *************  GENERATE CONNECTING CURVE  *************************************
    NRCURVE = NRCURVE + 1
    if (NRCURVE .gt. MAXNC) then
      write(*,*)'create_TW_POINT_AND_CONN_CURVE: increase MAXNC = ',MAXNC
      stop
    endif
    CURVES(NRCURVE)%Type       = 'Seglin'
    CURVES(NRCURVE)%EndPoNo(1) = Np
    CURVES(NRCURVE)%EndPoNo(2) = NRPOINT
#if I_PRINT >= 2
    write(*,*)'************  connecting curve generated.'
#endif 
#if I_PRINT >= 1
    write(*,*)'create_TW_POINT_AND_CONN_CURVE: done!'
#endif    
!    
end subroutine create_TW_POINT_AND_CONN_CURVE
