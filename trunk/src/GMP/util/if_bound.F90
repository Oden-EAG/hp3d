!--------------------------------------------------------------------------
!   function name       - if_bound
!--------------------------------------------------------------------------
!   latest revision    - Aug 08
!
!   purpose            - Routine checks whether a point lies within
!                        a domain enclosed by specified bounding 
!                        surfaces
!
!   arguments         
!     in:
!            Xp        - coordinates of a point
!            Nr_bound  - number of bounding surfaces
!            Ns_bound  - list of the bounding surfaces
!     out:
!            If_bound  = 1 if the point is within the area
!                      = 0 if the point is on the boundary
!                      = -1 if the point is outside
!--------------------------------------------------------------------------
integer function if_bound(Xp,Nr_bound,Ns_bound)
!--------------------------------------------------------------------------
! MODULES      
  use control
!--------------------------------------------------------------------------
! DUMMY ARGUMENTS
  real*8, dimension(3), intent(in)         :: Xp
  integer, intent(in)                      :: Nr_bound
  integer, dimension(Nr_bound), intent(in) :: Ns_bound
!--------------------------------------------------------------------------
! LOCAL VARIABLES
  integer              :: is
  real*8               :: fval
  real*8, dimension(3) :: dfdx  
!--------------------------------------------------------------------------
! 
    if_bound = 1
      do is = 1, Nr_bound
        call surf(Ns_bound(is),Xp, fval,dfdx)
        if (abs(fval) .le. GEOM_TOL)  if_bound = min(if_bound,0)
        if (fval .gt. GEOM_TOL)  if_bound = min(if_bound,-1)
      enddo
!
end function if_bound
