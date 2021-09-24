!
! REMARK: THIS ROUTINE MUST BE OMP THREAD-SAFE
!         DUE TO OMP PARALLELIZATION OF UPDATE_DDOF->DIRICHLET->EXACT
!------------------------------------------------------------------------------
!> Purpose : exact (manufactured) solution
!> last mod: May 2021
!
!> @param[in]  X      - a point in physical space
!> @param[in]  Mdle   - element (middle node) number
!> @param[out] ValH   - value of the H1 solution
!> @param[out] DvalH  - corresponding first derivatives
!> @param[out] D2valH - corresponding second derivatives
!> @param[out] ValE   - value of the H(curl) solution
!> @param[out] DvalE  - corresponding first derivatives
!> @param[out] D2valE - corresponding second derivatives
!> @param[out] ValV   - value of the H(div) solution
!> @param[out] DvalV  - corresponding first derivatives
!> @param[out] D2valV - corresponding second derivatives
!> @param[out] ValQ   - value of the H(div) solution
!> @param[out] DvalQ  - corresponding first derivatives
!> @param[out] D2valQ - corresponding second derivatives
!------------------------------------------------------------------------------
!
      subroutine exact(Xp,Mdle, ValH,DvalH,D2valH, ValE,DvalE,D2valE, &
                       ValV,DvalV,D2valV, ValQ,DvalQ,D2valQ)
!
      use data_structure3D
!
      implicit none
      common /cexact/ NP1,NP2,NP3,ICOMP
      integer :: NP1,NP2,NP3,ICOMP
!
      real(8), intent(in)  :: Xp(3)
      integer, intent(in)  :: Mdle
!
      real(8),dimension(  MAXEQNH    ), intent(out) ::   ValH
      real(8),dimension(  MAXEQNH,3  ), intent(out) ::  DvalH
      real(8),dimension(  MAXEQNH,3,3), intent(out) :: D2valH
      real(8),dimension(3,MAXEQNE    ), intent(out) ::   ValE
      real(8),dimension(3,MAXEQNE,3  ), intent(out) ::  DvalE
      real(8),dimension(3,MAXEQNE,3,3), intent(out) :: D2valE
      real(8),dimension(3,MAXEQNV    ), intent(out) ::   ValV
      real(8),dimension(3,MAXEQNV,3  ), intent(out) ::  DvalV
      real(8),dimension(3,MAXEQNV,3,3), intent(out) :: D2valV
      real(8),dimension(  MAXEQNQ    ), intent(out) ::   ValQ
      real(8),dimension(  MAXEQNQ,3  ), intent(out) ::  DvalQ
      real(8),dimension(  MAXEQNQ,3,3), intent(out) :: D2valQ
!
!  ...Locals
      integer :: nsol
      real(8) :: x,y,z, fx,dfx,d2fx, fy,dfy,d2fy, fz,dfz,d2fz
!
      integer :: iprint = 0
!
!------------------------------------------------------------------------------
!
!  ...initialize exact solution
      ValH=ZERO ; DvalH=ZERO ; D2valH=ZERO
      ValE=ZERO ; DvalE=ZERO ; D2valE=ZERO
      ValV=ZERO ; DvalV=ZERO ; D2valV=ZERO
      ValQ=ZERO ; DvalQ=ZERO ; D2valQ=ZERO
!
      x = Xp(1); y = Xp(2); z = Xp(3)
      nsol=10
      select case(nsol)
!
!  .....a quadratic H1 function
        case(1)
          ValH(1:3) = x**2 + 1.d0
          DvalH(1:3,1) = 2.d0*x
          D2valH(1:3,1,1) = 2.d0
!
!  .....a quadratic function with zero Dirichlet BCs
        case(2)
          ValH(1:3) = x*(1.d0-x)
          DvalH(1:3,1) = -2.d0*x + 1.d0
          D2valH(1:3,1,1) = -2.d0
!
!  .....monomial solutions to test shape functions routines
        case(10)
          fx = x**NP1
          select case(NP1)
          case(0); dfx = 0.d0; d2fx = 0.d0
          case(1); dfx = 1.d0; d2fx = 0.d0
          case default; dfx = NP1*x**(NP1-1); d2fx = NP1*(NP1-1)*x**(NP1-2)
          end select
          fy = y**NP2
          select case(NP2)
          case(0); dfy = 0.d0; d2fy = 0.d0
          case(1); dfy = 1.d0; d2fy = 0.d0
          case default; dfy = NP2*y**(NP2-1); d2fy = NP2*(NP2-1)*y**(NP2-2)
          end select
          fz = z**NP3
          select case(NP3)
          case(0); dfz = 0.d0; d2fz = 0.d0
          case(1); dfz = 1.d0; d2fz = 0.d0
          case default; dfz = NP3*z**(NP3-1); d2fz = NP3*(NP3-1)*z**(NP3-2)
          end select
          if (iprint.ge.2) then
            write(*,*) 'exact: NP1,NP2,NP3 = ',NP1,NP2,NP3
            write(*,*) '       fx,dfx,d2fx = ',fx,dfx,d2fx
            write(*,*) '       fy,dfy,d2fy = ',fy,dfy,d2fy
            write(*,*) '       fz,dfz,d2fz = ',fz,dfz,d2fz
          endif
!
          ValH(1)       =   fx*  fy*  fz
          DvalH(1,1)    =  dfx*  fy*  fz
          DvalH(1,2)    =   fx* dfy*  fz
          DvalH(1,3)    =   fx*  fy* dfz
          D2valH(1,1,1) = d2fx*  fy*  fz
          D2valH(1,1,2) =  dfx* dfy*  fz
          D2valH(1,1,3) =  dfx*  fy* dfz
          D2valH(1,2,2) =   fx*d2fy*  fz
          D2valH(1,2,3) =   fx* dfy* dfz
          D2valH(1,3,3) =   fx*  fy*d2fz
          D2valH(1,2,1) = D2valH(1,1,2)
          D2valH(1,3,1) = D2valH(1,1,3)
          D2valH(1,3,2) = D2valH(1,2,3)
!
          ValE(ICOMP,1)       =   fx*  fy*  fz
          DvalE(ICOMP,1,1)    =  dfx*  fy*  fz
          DvalE(ICOMP,1,2)    =   fx* dfy*  fz
          DvalE(ICOMP,1,3)    =   fx*  fy* dfz
          D2valE(ICOMP,1,1,1) = d2fx*  fy*  fz
          D2valE(ICOMP,1,1,2) =  dfx* dfy*  fz
          D2valE(ICOMP,1,1,3) =  dfx*  fy* dfz
          D2valE(ICOMP,1,2,2) =   fx*d2fy*  fz
          D2valE(ICOMP,1,2,3) =   fx* dfy* dfz
          D2valE(ICOMP,1,3,3) =   fx*  fy*d2fz
          D2valE(ICOMP,1,2,1) = D2valE(ICOMP,1,1,2)
          D2valE(ICOMP,1,3,1) = D2valE(ICOMP,1,1,3)
          D2valE(ICOMP,1,3,2) = D2valE(ICOMP,1,2,3)
!
          ValV(ICOMP,1)       =   fx*  fy*  fz
          DvalV(ICOMP,1,1)    =  dfx*  fy*  fz
          DvalV(ICOMP,1,2)    =   fx* dfy*  fz
          DvalV(ICOMP,1,3)    =   fx*  fy* dfz
          D2valV(ICOMP,1,1,1) = d2fx*  fy*  fz
          D2valV(ICOMP,1,1,2) =  dfx* dfy*  fz
          D2valV(ICOMP,1,1,3) =  dfx*  fy* dfz
          D2valV(ICOMP,1,2,2) =   fx*d2fy*  fz
          D2valV(ICOMP,1,2,3) =   fx* dfy* dfz
          D2valV(ICOMP,1,3,3) =   fx*  fy*d2fz
          D2valV(ICOMP,1,2,1) = D2valV(ICOMP,1,1,2)
          D2valV(ICOMP,1,3,1) = D2valV(ICOMP,1,1,3)
          D2valV(ICOMP,1,3,2) = D2valV(ICOMP,1,2,3)
!
          ValQ(1)       =   fx*  fy*  fz
          DvalQ(1,1)    =  dfx*  fy*  fz
          DvalQ(1,2)    =   fx* dfy*  fz
          DvalQ(1,3)    =   fx*  fy* dfz
          D2valQ(1,1,1) = d2fx*  fy*  fz
          D2valQ(1,1,2) =  dfx* dfy*  fz
          D2valQ(1,1,3) =  dfx*  fy* dfz
          D2valQ(1,2,2) =   fx*d2fy*  fz
          D2valQ(1,2,3) =   fx* dfy* dfz
          D2valQ(1,3,3) =   fx*  fy*d2fz
          D2valQ(1,2,1) = D2valQ(1,1,2)
          D2valQ(1,3,1) = D2valQ(1,1,3)
          D2valQ(1,3,2) = D2valQ(1,2,3)
!
        end select
!
      if (iprint.ge.1) then
        write(*,7100) Xp
 7100   format('exact: Xp = ',3f8.3)
        write(*,7010) ValH(1),DvalH(1,1:3),D2valH(1,1:3,1:3)
 7010   format('       ValH = ',e12.5,' DvalH = ',3e12.5,' D2valH = ',9e12.5)
        call pause
      endif
      
      
!
      end subroutine exact
