!--------------------------------------------------------------------------------------------
!   routine name       - quadB
!--------------------------------------------------------------------------------------------
!   latest revision    - Jul 09
!
!   purpose            - evaluate quad bubble
!
!   arguments :
!     in:
!        No            - quad number
!        Xi            - coordinates of a point in the quad
!                        in a system of coordinates
!        Norient       - quad orientation wrt the system
!                        of coordinates
!     out:
!        X             - physical coordinates
!        dX_dXi        - derivatives of physical coordinates 
!--------------------------------------------------------------------------------------------
subroutine quadB(No,Xi,Norient, X,dX_dXi)
!--------------------------------------------------------------------------------------------
! MODULES
  use kinds
  use control
  use GMP    
  use element_data
!--------------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 2
!--------------------------------------------------------------------------------------------
  IMPLICIT NONE
!--------------------------------------------------------------------------------------------
! DUMMY ARGUMENTS
  integer,                  intent(in)    :: No
  real*8, dimension(2),   intent(in)    :: Xi
  integer,                  intent(in)    :: Norient
  real*8, dimension(3),   intent(out)   :: X
  real*8, dimension(3,2), intent(out)   :: dX_dXi
!--------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  real*8, dimension(3,4)                :: xv
  real*8                                :: zeta
  real*8, dimension(2)                  :: dzeta_deta
  real*8, dimension(2)                  :: eta
  real*8, dimension(3,2)                :: dX_deta
  real*8, dimension(3)                  :: alpha
  real*8, dimension(3)                  :: dalpha_dzeta
  real*8, dimension(2,2,0:7), parameter :: deta_dXi = reshape                             &
                                                        ((/ 1.d0,  0.d0,  0.d0,  1.d0,      &
                                                            0.d0, -1.d0,  1.d0,  0.d0,      &
                                                           -1.d0,  0.d0,  0.d0, -1.d0,      &
                                                            0.d0,  1.d0, -1.d0,  0.d0,      &
                                                            0.d0,  1.d0,  1.d0,  0.d0,      &
                                                           -1.d0,  0.d0,  0.d0,  1.d0,      &
                                                            0.d0, -1.d0, -1.d0,  0.d0,      &
                                                            1.d0,  0.d0,  0.d0, -1.d0/),    &
                                                                                  (/2,2,8/))
  real*8, dimension(4)                  :: blend
  real*8, dimension(4,2)                :: dblend  
  real*8, dimension(4)                  :: shap_v
  real*8, dimension(4,2)                :: dshap_v  
  integer                                 :: norientc
  integer                                 :: iv,nc,ie,ivar,np,i,j
#if I_PRINT >= 2  
  real*8, dimension(3)                  :: aux
  real*8, dimension(3,2)                :: void
  real*8                                :: smax
#endif  
!--------------------------------------------------------------------------------------------
!
! ..check quad type
    if (RECTANGLES(No)%Type .ne. 'TraQua') then
      write(*,*)'quadB: inconsistent type = ',RECTANGLES(No)%Type
      stop
    endif
! ..transform Xi --> eta
    select case(Norient)
      case (0)
        eta(1) = Xi(1);         eta(2) = Xi(2)
      case (1)
        eta(1) = Xi(2);         eta(2) = 1.d0 - Xi(1)       
      case (2)
        eta(1) = 1.d0 - Xi(1);  eta(2) = 1.d0 - Xi(2)
      case (3)
        eta(1) = 1.d0 - Xi(2);  eta(2) = Xi(1)
      case (4)
        eta(2) = Xi(1);         eta(1) = Xi(2)
      case (5)
        eta(2) = Xi(2);         eta(1) = 1.d0 - Xi(1)       
      case (6)
        eta(2) = 1.d0 - Xi(1);  eta(1) = 1.d0 - Xi(2)
      case (7)
        eta(2) = 1.d0 - Xi(2);  eta(1) = Xi(1)
      case default
        write(*,*)'quadB: unknown orientation.'
        stop
    end select 
! ..blending functions
    blend(1)      = 1.d0 - eta(2)
    blend(2)      = eta(1)
    blend(3)      = eta(2)
    blend(4)      = 1.d0 - eta(1)
    dblend(1,1:2) = (/ 0.d0, -1.d0/)
    dblend(2,1:2) = (/ 1.d0,  0.d0/)
    dblend(3,1:2) = (/ 0.d0,  1.d0/)
    dblend(4,1:2) = (/-1.d0,  0.d0/)
! ..vertex shape functions
    shap_v(1)      = (1.d0 - eta(1))*(1.d0 - eta(2))
    shap_v(2)      = eta(1)*(1.d0 - eta(2))
    shap_v(3)      = eta(1)*eta(2)
    shap_v(4)      = (1.d0 - eta(1))*eta(2)
    dshap_v(1,1:2) = (/eta(2) - 1.d0, eta(1) - 1.d0/)
    dshap_v(2,1:2) = (/1.d0 - eta(2),       -eta(1)/)
    dshap_v(3,1:2) = (/       eta(2),        eta(1)/)
    dshap_v(4,1:2) = (/      -eta(2), 1.d0 - eta(1)/)
#if I_PRINT >= 1      
    write(*,*) 'quadB: computing quad bubble'
    write(*,1) No,Xi(1:2),Norient,eta(1:2)
1   format(' *****  No,Xi,Norient,eta = ',i6,2e12.5,i3,2x,2e12.5)
#endif
! ..evaluate rectangle parametrization
    call recta(No,eta, X,dX_deta)
#if I_PRINT >= 2    
    write(*,*) 'quadB: ORIGINAL X,dX_deta'
    do ivar = 1, 3
      write(*,7011) X(ivar),dX_deta(ivar,1:2)
7011  format(e12.5,2x,2e12.5)
    enddo
#endif      
! ..get the vertex coordinates
    do iv = 1, 4
      np = RECTANGLES(No)%VertNo(iv)
      xv(1:3,iv) = POINTS(np)%Rdata(1:3)
    enddo
!
! *****************************  SUBTRACT LINEAR INTERPOLANT  ***************************** |
! ..loop through vertices      
    do iv = 1, 4
#if I_PRINT >= 2    
!============================================================================================     
! Check consistency of parametrizations btw vertices an 'rect' routine                     |      
!--------------------------------------------------------------------------------------------      
! ....call 'rect' routine at master triangle vertices                                      !
      call recta(No,QUADR_COORD(1:2,iv), aux,void)                                          !
      smax = 0.d0                                                                           !
! ....accumulate error                                                                      !
      do ivar = 1, 3                                                                        !
        smax = max(smax,abs(aux(ivar) - xv(ivar,iv)))                                       !
      enddo                                                                                 !
! ....check whether GEOM_TOL is exceeded                                                    ! 
      if (smax .gt. GEOM_TOL) then                                                          !
        write(*,7001) No,iv,smax                                                            !
 7001   format('trianB: No,iv,smax = ',i5,i2,e12.5)                                         !
        write(*,*) 'aux = ',aux                                                             !
        write(*,*) 'xv  = ',xv(1:3,iv)                                                      !
        call pause                                                                          !
      endif                                                                                 !
!============================================================================================
#endif
      X = X - xv(1:3,iv)*shap_v(iv)
      do i = 1, 3
        do j = 1, 2
          dX_deta(i,j) = dX_deta(i,j) - xv(i,iv)*dshap_v(iv,j)
        enddo
      enddo
! ..end of loop through vertices        
    enddo
#if I_PRINT >= 2      
    write(*,*) 'quadB: AFTER VERTICES  X,dX_deta = '
    do ivar = 1, 3
      write(*,7011) X(ivar),dX_deta(ivar,1:2)
    enddo
#endif
!
! *****************************  SUBTRACT EDGE CONTRIBUTIONS  ***************************** |
! ..loop through edges      
    do ie = 1, 4
! ....get the curve number
      nc = RECTANGLES(No)%EdgeNo(ie)
      norientc = 0
      if (nc .lt. 0) then
        nc = -nc;  norientc = 1
      endif
      if (CURVES(nc)%Type .eq.' Seglin') cycle
      call proj_quad2edge(eta,ie, zeta,dzeta_deta)
#if I_PRINT >= 2       
      write(*,9001) ie,zeta
9001  format(' quadB: ie, zeta = ',I2,' ; ',E12.5)
#endif
      if ((zeta .lt. GEOM_TOL) .or. (zeta .gt. (1.d0 - GEOM_TOL))) cycle
#if I_PRINT >= 2 
      write(*,7003) ie,nc,CURVES(nc)%Type
7003  format(' quadB: ie,nc,Type = ',i2,i5,2x,a5)
#endif
! ....compute bubble function     
      call curveB(nc,zeta,norientc, alpha,dalpha_dzeta)
! ....subtract edge bubble        
      X = X - alpha*blend(ie)
      do i = 1, 3
        do j = 1, 2
          dX_deta(i,j) = dX_deta(i,j) - dalpha_dzeta(i)*dzeta_deta(j)*blend(ie)          &
                                      - alpha(i)*dblend(ie,j)
        enddo
      enddo
#if I_PRINT >= 2        
      write(*,*) 'quadB: X,dX_deta AFTER EDGE = ',ie
      do ivar = 1, 3
        write(*,7011) X(ivar),dX_deta(ivar,1:2)
      enddo
      call pause
#endif
! ..end of loop through edges        
    enddo
#if I_PRINT >= 2    
    write(*,*) 'quadB: AFTER EDGES  X,dX_deta = '
    do ivar = 1, 3
      write(*,7011) X(ivar),dX_deta(ivar,1:2)
    enddo
    call pause
#endif
! ..account for orientations in derivatives
    dX_dXi(1:3,1:2) = matmul(dX_deta(1:3,1:2), deta_dXi(1:2,1:2,Norient))
#if I_PRINT >= 2    
    write(*,*) 'trianB: FINAL  X,dX_dXi = '
    do ivar = 1, 3
      write(*,7011) X(ivar),dX_dXi(ivar,1:2)
    enddo
    call pause
#endif
!
end subroutine quadB
!----------------------------------------------------------------------------------------------------  
