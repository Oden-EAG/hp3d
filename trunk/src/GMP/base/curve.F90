!-----------------------------------------------------------------------
!> Purpose : physical coordinates for a curve parametrization, and 
!!           their derivative wrt to a GIVEN coordinate
!!
!! @param[in ] No      - curve number
!! @param[in ] Norient - orientation of GLOBAL REFERENCE coordinate 
!!                       wrt a GIVEN coordinate system
!! @param[in ] T       - coordinate of a point on the curve
!!                       according to GIVEN system
!! @param[out] X       - physical coordinates
!! @param[out] Dxdt    - derivatives of physical coordinate
!!
!! @revision Nov 12
!-----------------------------------------------------------------------
!
subroutine curve_local(No,Norient,T, X,Dxdt)
!      
      implicit none
      integer,            intent(in ) :: No,Norient
      real*8             ,intent(in ) :: T
      real*8,dimension(3),intent(out) :: X
      real*8,dimension(3),intent(out) :: Dxdt
!
!  ...GLOBAL REFERENCE coordinate
      real*8              :: eta
!  ...derivative of GLOBAL REFERENCE coordinate wrt GIVEN coordinate      
      real*8              :: detadt
!  ...derivative of PHYSICAL coordinates wrt GLOBAL REFERENCE coordinate      
      real*8,dimension(3) :: dxdeta
!-----------------------------------------------------------------------
!
!  ...GIVEN -> GLOBAL REFERENCE
      select case(Norient)
      case(0) ; eta=T      ; detadt= 1.d0
      case(1) ; eta=1.d0-T ; detadt=-1.d0
      case default
        write(*,7001) Norient
 7001   format(' curve_local: Norient = ',i3)
        stop
      endselect
!
!  ...GLOBAL REFERENCE -> PHYSICAL
      call curve(No,eta, X,dxdeta)
!
!  ...derivative of PHYSICAL wrt GIVEN
      Dxdt(1:3)=dxdeta(1:3)*detadt
!
!
endsubroutine curve_local
!
!
!
!----------------------------------------------------------------------------
!> Purpose :  physical coordinates for a curve parametrization, and
!!            their derivative wrt to reference coordinate
!!
!! @param[in ] No     - curve number
!! @param[in ] Eta    - reference coordinate  (between 0 and 1)
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives of the physical coordinates
!!
!! @revision Nov 12      
!----------------------------------------------------------------------------
!
subroutine curve(No,Eta, X,Dxdeta)
!
      use GMP
      use bezier
!      
      implicit none
      integer            ,intent(in ) :: No
      real*8             ,intent(in ) :: Eta
      real*8,dimension(3),intent(out) :: X,Dxdeta
!
      real*8, dimension(3,2) :: xv
      integer                :: iprint,i,np,j,ierror
      character(len=10)      :: type
!----------------------------------------------------------------------------
!
      iprint=0
!
!  ...curve type
      type=CURVES(No)%Type  
      if (iprint.eq.1) then
        write(*,7000) No, type, Eta
 7000   format(' curve: No,type,Eta = ',i4,2x,a10,2x,f8.3)
      endif
!
      select case(type)
!
!  ...straight segment
      case('Seglin')
        do i=1,2
          np=CURVES(No)%EndPoNo(i) ; call pointr(np, xv(1:3,i))
        enddo
        do j=1,3
          X(j)=Eta*xv(j,2)+(1.d0-Eta)*xv(j,1) ; Dxdeta(j)=xv(j,2)-xv(j,1)
        enddo
!
!  ...quarter of a circle..................................................
      case('QuaCir')    ; call curve_QuaCir(No,Eta, X,Dxdeta)
!
!  ...segment of a circle..................................................
      case('SegCir')    ; call curve_SegCir(No,Eta, X,Dxdeta)
!
!  ...quarter of an ellipse (0 to pi/2)....................................
      case('QuaEl1')    ; call curve_QuaEl1(No,Eta, X,Dxdeta)
!
!  ...quarter of an ellipse (-pi/4 to pi/4)................................
      case('QuaEl2')    ; call curve_QuaEl2(No,Eta, X,Dxdeta)
!
!  ...quarter of an ellipse (0 to pi/2)....................................
      case('QuaSEl')    ; call curve_QuaSEl(No,Eta, X,Dxdeta)
!
!  ...implicit curve.......................................................
      case('ImpCir')    ; call curve_ImpCir(No,Eta, X,Dxdeta)   
!
!  ...deg 3 Hermite curve..................................................
      case('3HermCur')  ; call curve_3HermCur(No,Eta, X,Dxdeta)   
!
!  ...deg 5 Hermite curve..................................................
      case('HermCur')   ; call curve_HermCur(No,Eta, X,Dxdeta)   
!   
!  ...geodesic curve on algebraic surface (Cylinder, Cone, Sphere).........
      case('1SurfsCur') ; call curve_1SurfsCur(No,Eta, X,Dxdeta)
!   
!  ...intersection of 2 algebraic surfaces (computed numerically)..........
      case('2SurfsCur') ; call curve_2SurfsCur(No,1,2,Eta, X,Dxdeta,ierror)
        if (ierror.eq.1) stop
!        
!  ...intersection of 3 algebraic surfaces (computed numerically)..........
      case('3SurfsCur') ; call curve_3SurfsCur(No,Eta, X,Dxdeta)        
!              
!  ...Bezier curve.........................................................
      case('5Bezier','7Bezier') ; call curve_bezier(No,Eta, X=X,dXdEta=Dxdeta)
!   
!  ...curve on reconstructed surface (NOT DEVELOPED).......................
      case('1SurfrCur') ; call curve_1SurfrCur(No,Eta, X,Dxdeta)
!   
!  ...intersection of 2 reconstructed surfaces (NOT DEVELOPED).............
      case('2SurfrCur') ; call curve_2SurfrCur(No,Eta, X,Dxdeta)
!              
!  ...cylinder geodesic (Cynthia's LEGACY version).........................
      case('CylGeod')   ; call curve_CylGeod(No,Eta, X,Dxdeta)
!              
!  ...part of spherical arc (NOT WORKING!)
!!!   case('SegCir') ; call curve_SegCir(No,Eta, X,Dxdeta)
!
#ifdef _PYGMP
      case('FourSurCur')
         call curve_fourier_surf(no, eta, x, dxdeta)
#endif
!
!  ...image of a straight line segment through a global system of 
!     cylindrical coordinates
      case('CylCoord')
!      case('CylCur')
         call curve_CylCoord(no, eta, x, dxdeta)
      case default
        write(*,7001) type
 7001   format(' curve: UNKNOWN CURVE, type = ',a10)
        stop
!
      endselect
!
!
end subroutine curve
