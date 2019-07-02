!------------------------------------------------------------------------------
!> Purpose : source term
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[in]  Xi    - master element coordinates
!! @param[in]  X     - physical coordinates
!! @param[out] Zfval - rhs
!------------------------------------------------------------------------------
!
subroutine getf(Mdle,Xi,X, Zfval)
  use control    , only: NEXACT
  use parameters , only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ
  use lapl       , only: NR_RHS_PROB 

  implicit none
  !------------------------------------------------------------------------------
  integer,                      intent(in)  :: Mdle
  real*8,dimension(3),          intent(in)  :: Xi,X
  real*8,dimension(NR_RHS_PROB),intent(out) :: Zfval
  !------------------------------------------------------------------------------
  ! exact solution
  real*8,dimension(  MAXEQNH    ) ::   zvalH
  real*8,dimension(  MAXEQNH,3  ) ::  zdvalH
  real*8,dimension(  MAXEQNH,3,3) :: zd2valH
  real*8,dimension(3,MAXEQNE    ) ::   zvalE
  real*8,dimension(3,MAXEQNE,3  ) ::  zdvalE
  real*8,dimension(3,MAXEQNE,3,3) :: zd2valE
  real*8,dimension(3,MAXEQNV    ) ::   zvalV
  real*8,dimension(3,MAXEQNV,3  ) ::  zdvalV
  real*8,dimension(3,MAXEQNV,3,3) :: zd2valV
  real*8,dimension(  MAXEQNQ    ) ::   zvalQ
  real*8,dimension(  MAXEQNQ,3  ) ::  zdvalQ
  real*8,dimension(  MAXEQNQ,3,3) :: zd2valQ
  !
  ! node case number
  integer :: icase

  Zfval(1:NR_RHS_PROB) = 0.d0 
  select case(NEXACT)
  case(0)          
     !******************************************************************************    
     ! Y O U R   F A V O R I T E   R H S   H E R E                                 !
     Zfval(1:NR_RHS_PROB) = 1.d0
     !******************************************************************************    
  case(1)
     ! compute exact soution      
     icase = 1
     call exact( &
          X,icase, zvalH,zdvalH,zd2valH,zvalE,zdvalE,zd2valE, &
          zvalV,zdvalV,zd2valV,zvalQ,zdvalQ,zd2valQ)
     ! evaluate rhs's                    
     Zfval(1:NR_RHS_PROB) = &
          -(zd2valH(1,1,1)**2 +   &
          zd2valH(1,2,2)**2 +   &
          zd2valH(1,3,3)**2)                  
  endselect
  !
endsubroutine getf
!------------------------------------------------------------------------------
!> Purpose : Neumann load
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[in]  Ibc   - boundary conditions 
!! @param[in]  Xi    - master element coordinates
!! @param[in]  X     - physical coordinates
!! @param[in]  Rn    - outward unit vector
!! @param[out] Zgval - Neumann loads
!------------------------------------------------------------------------------
subroutine getg(Mdle,Ibc,Xi,X,Rn, Zgval)
  use control    , only: NEXACT
  use parameters , only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ
  use lapl       , only: NR_RHS_PROB 
  implicit none
  !------------------------------------------------------------------------------
  integer,                      intent(in)  :: Mdle
  real*8,dimension(3),          intent(in)  :: Xi,X,Rn
  integer,dimension(6),         intent(in)  :: Ibc
  real*8,dimension(NR_RHS_PROB),intent(out) :: Zgval
  !------------------------------------------------------------------------------
  ! exact solution
  real*8,dimension(  MAXEQNH    ) ::   zvalH
  real*8,dimension(  MAXEQNH,3  ) ::  zdvalH
  real*8,dimension(  MAXEQNH,3,3) :: zd2valH
  real*8,dimension(3,MAXEQNE    ) ::   zvalE
  real*8,dimension(3,MAXEQNE,3  ) ::  zdvalE
  real*8,dimension(3,MAXEQNE,3,3) :: zd2valE
  real*8,dimension(3,MAXEQNV    ) ::   zvalV
  real*8,dimension(3,MAXEQNV,3  ) ::  zdvalV
  real*8,dimension(3,MAXEQNV,3,3) :: zd2valV
  real*8,dimension(  MAXEQNQ    ) ::   zvalQ
  real*8,dimension(  MAXEQNQ,3  ) ::  zdvalQ
  real*8,dimension(  MAXEQNQ,3,3) :: zd2valQ
  !  
  ! node case number
  integer :: icase
  !------------------------------------------------------------------------------
  ! initialize source terms
  Zgval(1:NR_RHS_PROB) = 0.d0 
  !
  select case(NEXACT)
  case(0)          
     !******************************************************************************    
     !   Y O U R   F A V O R I T E   N E U M A N N   L O A D   H E R E             !
     Zgval(1:NR_RHS_PROB) = 0.d0
     !******************************************************************************    
  case(1)
     ! compute exact soution      
     icase = 1
     call exact(X,icase, zvalH,zdvalH,zd2valH,zvalE,zdvalE,zd2valE, &
          zvalV,zdvalV,zd2valV,zvalQ,zdvalQ,zd2valQ)
     ! evaluate rhs's                    
     Zgval(1:NR_RHS_PROB) = ZdvalH(1,1)*Rn(1) +  &
          ZdvalH(1,2)*Rn(2) +  &
          ZdvalH(1,3)*Rn(3)                  
  endselect
  !
endsubroutine getg
