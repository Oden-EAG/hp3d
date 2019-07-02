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
  use assembly   , only: NR_RHS 
  use parameters , only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ
  implicit none
!------------------------------------------------------------------------------
  integer,                     intent(in)  :: Mdle
  real*8,dimension(3),         intent(in)  :: Xi,X
  real*8,dimension(1:3,NR_RHS),intent(out) :: Zfval
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
!
! initialize source terms
  Zfval(1:3,1:NR_RHS) = 0.d0 
  !
  select case(NEXACT)
  ! unknown exact solution  
  case(0)          
  !
  !
!******************************************************************************    
! Y O U R   F A V O R I T E   R H S   H E R E                                 !
  Zfval(3,1:NR_RHS) = -0.1d0
!******************************************************************************    
  !
  ! 
!!!  ! known exact solution
!!!  case(1)
!!!    ! compute exact soution      
!!!    icase = 1
!!!    call exact(X,icase, zvalH,zdvalH,zd2valH,zvalE,zdvalE,zd2valE, &
!!!                        zvalV,zdvalV,zd2valV,zvalQ,zdvalQ,zd2valQ)
!!!    ! evaluate rhs's                    
!!!    Zfval(1:NR_RHS) = -(zd2valH(1,1,1)**2 +   &
!!!                        zd2valH(1,2,2)**2 +   &
!!!                        zd2valH(1,3,3)**2)                  
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
!! @param[in]  Xi    - master element coordinates
!! @param[in]  X     - physical coordinates
!! @param[in]  Rn    - outward unit vector
!! @param[out] Zgval - Neumann loads
!------------------------------------------------------------------------------
!
subroutine getg(Mdle,Ibc,Xi,X,Rn, Zgval)
  use control    , only: NEXACT
  use assembly   , only: NR_RHS 
  use parameters , only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ
  implicit none
!------------------------------------------------------------------------------
  integer,                   intent(in)  :: Mdle
  real*8,dimension(3),       intent(in)  :: Xi,X,Rn
  integer,dimension(6),      intent(in)  :: Ibc
  real*8,dimension(3,NR_RHS),intent(out) :: Zgval
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
!
! initialize source terms
  Zgval(1:3,1:NR_RHS) = 0.d0 
  !
  select case(NEXACT)
  ! unknown exact solution  
  case(0)          
  !
  !
!******************************************************************************    
!   Y O U R   F A V O R I T E   N E U M A N N   L O A D   H E R E             !
    Zgval(1:3,1:NR_RHS) = 0.d0
!******************************************************************************    
  !
  ! 
!!!  ! known exact solution
!!!  case(1)
!!!    ! compute exact soution      
!!!    icase = 1
!!!    call exact(X,icase, zvalH,zdvalH,zd2valH,zvalE,zdvalE,zd2valE, &
!!!                        zvalV,zdvalV,zd2valV,zvalQ,zdvalQ,zd2valQ)
!!!    ! evaluate rhs's                    
!!!    Zgval(1:NR_RHS) = ZdvalH(1,1)*Rn(1) +  &
!!!                      ZdvalH(1,2)*Rn(2) +  &
!!!                      ZdvalH(1,3)*Rn(3)                  
  endselect        
!
endsubroutine getg
!
!
!
!------------------------------------------------------------------------------
!> Purpose : exact (manufactured, no miracle!) solution
!!
!! @param[in]  X       - a point in physical space 
!! @param[in]  Icase   - node case (specifies what variables are supported)
!! @param[out] zvalH   - value of the H1 solution
!! @param[out] ZdvalH  - corresponding first derivatives
!! @param[out] Zd2valH - corresponding second derivatives
!! @param[out] ZvalE   - value of the H(curl) solution
!! @param[out] ZdvalE  - corresponding first derivatives
!! @param[out] Zd2valE - corresponding second derivatives
!! @param[out] ZvalV   - value of the H(div) solution
!! @param[out] ZdvalV  - corresponding first derivatives
!! @param[out] Zd2valV - corresponding second derivatives
!! @param[out] ZvalQ   - value of the H(div) solution
!! @param[out] ZdvalQ  - corresponding first derivatives
!! @param[out] Zd2valQ - corresponding second derivatives
!------------------------------------------------------------------------------
!
subroutine exact(X,Icase, ZvalH,ZdvalH,Zd2valH, &
                          ZvalE,ZdvalE,Zd2valE, &
                          ZvalV,ZdvalV,Zd2valV, &
                          ZvalQ,ZdvalQ,Zd2valQ)
  use data_structure3D
  implicit none
!------------------------------------------------------------------------------
  real*8,dimension(3),             intent(in)  :: X
  integer,                         intent(in)  :: Icase
  real*8,dimension(  MAXEQNH    ), intent(out) ::   ZvalH
  real*8,dimension(  MAXEQNH,3  ), intent(out) ::  ZdvalH
  real*8,dimension(  MAXEQNH,3,3), intent(out) :: Zd2valH
  real*8,dimension(3,MAXEQNE    ), intent(out) ::   ZvalE
  real*8,dimension(3,MAXEQNE,3  ), intent(out) ::  ZdvalE
  real*8,dimension(3,MAXEQNE,3,3), intent(out) :: Zd2valE
  real*8,dimension(3,MAXEQNV    ), intent(out) ::   ZvalV
  real*8,dimension(3,MAXEQNV,3  ), intent(out) ::  ZdvalV
  real*8,dimension(3,MAXEQNV,3,3), intent(out) :: Zd2valV
  real*8,dimension(  MAXEQNQ    ), intent(out) ::   ZvalQ
  real*8,dimension(  MAXEQNQ,3  ), intent(out) ::  ZdvalQ
  real*8,dimension(  MAXEQNQ,3,3), intent(out) :: Zd2valQ
!------------------------------------------------------------------------------
!
! initialize exact solution
  ZvalH = 0.d0 ; ZdvalH = 0.d0 ; Zd2valH = 0.d0  
  ZvalE = 0.d0 ; ZdvalE = 0.d0 ; Zd2valE = 0.d0 
  ZvalV = 0.d0 ; ZdvalV = 0.d0 ; Zd2valV = 0.d0 
  ZvalQ = 0.d0 ; ZdvalQ = 0.d0 ; Zd2valQ = 0.d0 
  !
  !
!******************************************************************************    
! Y O U R   F A V O R I T E   H1   S O L U T I O N   H E R E                  !
  ZvalH(1:3) = 0.d0
!******************************************************************************    
  !
  !
end subroutine exact
