!------------------------------------------------------------------------------
!> Purpose : source term
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[in]  X     - physical coordinates
!! @param[out] Fval - rhs
!------------------------------------------------------------------------------
!
subroutine getf(Mdle,X, Fval)
  use control   , only: NEXACT
  use assembly  , only: NR_RHS
  use parameters, only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ
  use isotropic_elast_material
  implicit none
!------------------------------------------------------------------------------
  integer,                      intent(in)  :: Mdle
  real*8, dimension(3),         intent(in)  :: X
  real*8, dimension(1:3,NR_RHS),intent(out) :: Fval
!------------------------------------------------------------------------------
! exact solution
  real*8, dimension(  MAXEQNH    ) ::   valH
  real*8, dimension(  MAXEQNH,3  ) ::  dvalH
  real*8, dimension(  MAXEQNH,3,3) :: d2valH
  real*8, dimension(3,MAXEQNE    ) ::   valE
  real*8, dimension(3,MAXEQNE,3  ) ::  dvalE
  real*8, dimension(3,MAXEQNE,3,3) :: d2valE
  real*8, dimension(3,MAXEQNV    ) ::   valV
  real*8, dimension(3,MAXEQNV,3  ) ::  dvalV
  real*8, dimension(3,MAXEQNV,3,3) :: d2valV
  real*8, dimension(  MAXEQNQ    ) ::   valQ
  real*8, dimension(  MAXEQNQ,3  ) ::  dvalQ
  real*8, dimension(  MAXEQNQ,3,3) :: d2valQ
! Elasticity Tensor
  real*8, dimension(3,3,3,3) :: C
! node case number
  integer :: icase
! counters
  integer :: icomp,j,k,l,iload
! printing flag
  integer :: iprint
!------------------------------------------------------------------------------
!
  iprint=0
!
  if (iprint == 1) then
    write(*,9999) Mdle,X
 9999 format(' get_f: Mdle,X = ',i8,2x,3(e12.5,2x))
  endif
! initialize source terms
  Fval(1:3,1:NR_RHS) = 0.d0
!
  select case(NEXACT)
!==============================================================================
!  UNKNOWN EXACT SOLUTION                                                     |
!==============================================================================
  case(0)
!
!==============================================================================
!  KNOWN EXACT SOLUTION                                                       |
!==============================================================================
  case(1)
!   compute exact soution
    call exact(X,icase, valH,dvalH,d2valH,valE,dvalE,d2valE, &
                        valV,dvalV,d2valV,valQ,dvalQ,d2valQ)
!   get elasticity tensor
    call getmat(X, C)
!
!   loop over rhs's
    do iload=1,NR_RHS
!
!     Assuming C is constant,
!       f_i = C_ijkl du_k/dx_lj
      do icomp=1,3
        do j=1,3; do k=1,3; do l=1,3
          Fval(icomp,iload) = Fval(icomp,iload) - C(icomp,j,k,l)*d2valH(k,l,j)
        enddo; enddo; enddo
        Fval(icomp,iload) = Fval(icomp,iload) - RHO*OMEGA**2*valH(icomp)
      enddo
    enddo
!==============================================================================
!  KNOWN EXACT SOLUTION , HOMOGENEOUS RHS                                     |
!==============================================================================
  case(2)
!
  endselect
!
endsubroutine getf
!
!
!
!------------------------------------------------------------------------------
!> Purpose : Neumann load
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[in]  Ibc   - boundary conditions
!! @param[in]  X     - physical coordinates
!! @param[in]  Rn    - outward unit vector
!! @param[out] Gval - Neumann loads
!------------------------------------------------------------------------------
!
subroutine getg(Mdle,Ibc,X,Rn, Gval)
  use control    , only: NEXACT
  use assembly   , only: NR_RHS
  use parameters , only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ
  use isotropic_elast_material
  implicit none
!------------------------------------------------------------------------------
  integer,                   intent(in)  :: Mdle
  real*8,dimension(3),       intent(in)  :: X,Rn
  integer,                   intent(in)  :: Ibc
  real*8,dimension(3,NR_RHS),intent(out) :: Gval
!------------------------------------------------------------------------------
! exact solution
  real*8,dimension(  MAXEQNH    ) ::   valH
  real*8,dimension(  MAXEQNH,3  ) ::  dvalH
  real*8,dimension(  MAXEQNH,3,3) :: d2valH
  real*8,dimension(3,MAXEQNE    ) ::   valE
  real*8,dimension(3,MAXEQNE,3  ) ::  dvalE
  real*8,dimension(3,MAXEQNE,3,3) :: d2valE
  real*8,dimension(3,MAXEQNV    ) ::   valV
  real*8,dimension(3,MAXEQNV,3  ) ::  dvalV
  real*8,dimension(3,MAXEQNV,3,3) :: d2valV
  real*8,dimension(  MAXEQNQ    ) ::   valQ
  real*8,dimension(  MAXEQNQ,3  ) ::  dvalQ
  real*8,dimension(  MAXEQNQ,3,3) :: d2valQ
!
! miscellaneous
  integer :: icomp,j,k,l,iload
  real*8, dimension(3,3,3,3) :: C
!------------------------------------------------------------------------------
!
! initialize source terms
  Gval(1:3,1:NR_RHS) = 0.d0
  !
  select case(NEXACT)
!==============================================================================
!  UNKNOWN EXACT SOLUTION                                                      |
!==============================================================================
  case(0)
!
!==============================================================================
!  KNOWN EXACT SOLUTION                                                       |
!==============================================================================
  case(1,2)
!
!   compute exact solution
    call exact(X,Mdle, valH,dvalH,d2valH, valE,dvalE,d2valE, &
                       valV,dvalV,d2valV, valQ,dvalQ,d2valQ)
!
!   compute the elastictiy tensor
    call getmat(X, C)
!
    select case(Ibc)
!   NEUMANN BC
    case(2)
!
!     loop over rhs's
      do iload=1,NR_RHS
!
!       g_i = C_ijkl du_k/dx_l n_j
        do icomp=1,3
          do k=1,3; do l=1,3; do j=1,3
            Gval(icomp,iload) = Gval(icomp,iload) + C(icomp,j,k,l)*dvalH(k,l)*Rn(j)
          enddo; enddo; enddo
        enddo
!     loop over rhs's
      enddo
!   INTERFACE BC
    case(3)
      ! nothing here yet
    endselect
  endselect
!
endsubroutine getg
