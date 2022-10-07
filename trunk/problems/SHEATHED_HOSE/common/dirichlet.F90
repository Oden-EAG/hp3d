!--------------------------------------------------------------------
!> Purpose : calculate dirichlet boundary condition
!! @param[in]  X      - physical coordinates of a point
!! @param[in]  Icase  - node case
!!
!! @param[out] ValH   - value of the H1 solution
!! @param[out] DvalH  - H1 corresponding first derivatives
!! @param[out] ValE   - value of the H(curl) solution
!! @param[out] DvalE  - H(curl) corresponding first derivatives
!! @param[out] ValV   - value of the H(div) solution
!! @param[out] DvalV  - H(div) corresponding first derivatives
!--------------------------------------------------------------------
subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
  use sheathed_isotropic_materials, only : R_inside,R_middle,R_outside,P_inner,P_outer,X_1,X_2
  use control    , only : NEXACT
  use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO
  implicit none
  real*8, dimension(3),          intent(in)  :: X
  integer,                       intent(in)  :: Icase,Mdle
! exact solution
  real*8,dimension(  MAXEQNH    ),intent(out) ::   ValH
  real*8,dimension(  MAXEQNH,3  ),intent(out) ::  DvalH
  real*8,dimension(  MAXEQNH,3,3)             :: d2valH
  real*8,dimension(3,MAXEQNE    ),intent(out) ::   ValE
  real*8,dimension(3,MAXEQNE,3  ),intent(out) ::  DvalE
  real*8,dimension(3,MAXEQNE,3,3)             :: d2valE
  real*8,dimension(3,MAXEQNV    ),intent(out) ::   ValV
  real*8,dimension(3,MAXEQNV,3  ),intent(out) ::  DvalV
  real*8,dimension(3,MAXEQNV,3,3)             :: d2valV
  real*8,dimension(  MAXEQNQ    )             ::   valQ
  real*8,dimension(  MAXEQNQ,3  )             ::  dvalQ
  real*8,dimension(  MAXEQNQ,3,3)             :: d2valQ
! printing flag
  integer :: iprint
! for exact solution
!--------------------------------------------------------------------
  real*8, dimension(3)   :: u
  real*8, dimension(3,3) :: gradU,epsilon,sigma
!--------------------------------------------------------------------
! miscellaneous
  integer :: ndom
  real(8) :: r,theta,tol
!--------------------------------------------------------------------
!
! printing flag : 0 - silent ; 1 - verbose
  iprint = 0
!
! initialize
  ValH = ZERO; DvalH = ZERO
  ValE = ZERO; DvalE = ZERO
  ValV = ZERO; DvalV = ZERO
  select case(NEXACT)
!==============================================================================
!  UNKNOWN EXACT SOLUTION                                                     |
!==============================================================================
    case(0)

      ValV(1,1)=1.0d0
      ValV(2,2)=1.0d0
      ValV(3,3)=1.0d0

      tol=1.0d-10
      r = dsqrt(X(2)**2+X(3)**2)
      theta = datan2(X(3),X(2))

      !  INNER FACE
      ! if (r.lt.R_middle) then
      if (r.lt.0.75d0) then
        ValV = P_inner*dcos(1.d0*theta)**2*ValV
      !  OUTER FACE
      else
        ValV = P_outer*ValV
      endif

      ! ENDS
      if ((abs(X(1)-X_1).le.tol).or.(abs(X(1)-X_2).le.tol)) then
        ValV=0.d0
      endif

!==============================================================================
!  KNOWN EXACT SOLUTION                                                       |
!==============================================================================
    case(1)
!     use the exact solution to determine Dirichlet data
      call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                          ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
    case(2)
      call find_domain(Mdle, ndom)

      !  add value of U (for RBM points)
      call exact_solution(X,ndom, u,gradU,epsilon,sigma)
      ValH = u
      ! ValH = u + 0.01d0 ! a constant shift in displacement should only affect convergence of displacement
      DvalH = gradU
      ValV = sigma

    case default
     write(*,*)'dirichlet: UNKNOWN EXACT SOLUTION FLAG', NEXACT
    stop 1
  end select

  if (iprint == 1) then
    write(*,1001)X(1:3),ValH(1:MAXEQNH)
 1001 format(' dirchlet: X,zvalH = ',3(e12.5,2x),2x,10(e12.5,2x))
  endif

end subroutine dirichlet
