!------------------------------------------------------------------------------------
!> Purpose : implicit surface parameterization
!!
!! @param[in]  No   - surface number
!! @param[in]  X    - physical coordinates of a point
!! @param[out] Fval - function value
!! @param[out] Dfdx - function gradient
!!
!! @revision Nov 12
!------------------------------------------------------------------------------------
subroutine surf(No,X, Fval,Dfdx)
!
      use GMP
!
      implicit none
      integer             , intent(in ) :: No
      real*8, dimension(3), intent(in ) :: X
      real*8              , intent(out) :: Fval
      real*8, dimension(3), intent(out) :: Dfdx
!------------------------------------------------------------------------------------
      real*8, dimension(3) :: x0
      character(len=10)    :: stype
      integer              :: iprint,no1,nsign,np
!------------------------------------------------------------------------------------
!
      iprint=0
!
!  ...determine surface orientations
      no1=iabs(No) ; nsign=No/no1
      if (No.eq.0) then
        write(*,*) 'following : ', No, no1, nsign
      endif
!
!  ...find the surface type
      stype = SURFACES(no1)%Type
      if (iprint.eq.1) then
        write(*,7000) No,stype
 7000   format(' surf: No,type = ',i5,2x,a10)
      endif
!
!  ...select the surface
      select case(stype)
!
!  ...PLANE............................................................
      case('VecPt')
        call plane1(   X,SURFACES(no1)%Rdata(1:3),         &
                         SURFACES(no1)%Rdata(4:6), Fval,Dfdx)
!
!  ...PLANE............................................................
      case('ThrPt')
        call plane2(   X,SURFACES(no1)%Rdata(1:3),         &
                         SURFACES(no1)%Rdata(4:6),         &
                         SURFACES(no1)%Rdata(7:9), Fval,Dfdx)
!
!  ...PLANE............................................................
      case('PtNoVec')
        np=SURFACES(no1)%Idata(1)
        call pointr(np, x0)
        call plane1(   X,x0                      ,         &
                         SURFACES(no1)%Rdata(1:3), Fval,Dfdx)
!
!  ...PLANE............................................................
      case('PtNo2Pt')
        np=SURFACES(no1)%Idata(1)
        call pointr(np, x0)
        call plane2(   X,x0                      ,         &
                         SURFACES(no1)%Rdata(1:3),         &
                         SURFACES(no1)%Rdata(4:6), Fval,Dfdx)
!
!  ...PLANE PARAMETRIZED WITH CYLINDRICAL COORDINATES..................
      case('PPwCC')
        write(*,*) 'surf: stype = ',stype
        stop 1
!
!  ...SPHERE...........................................................
      case('Sphere')
        call sphere1(  X,SURFACES(no1)%Rdata(1:3),         &
                         SURFACES(no1)%Rdata(  4), Fval,Dfdx)
!
!  ...CYLINDER.........................................................
      case('Cylinder')
        call cylinder( X,SURFACES(no1)%Rdata(1:3),         &
                         SURFACES(no1)%Rdata(4:6),         &
                         SURFACES(no1)%Rdata(  7), Fval,Dfdx)
!
!  ...ELLIPSOID........................................................
      case('Ellipsoid')
        call ellipsoid(X,SURFACES(no1)%Rdata(1:3),         &
                         SURFACES(no1)%Rdata(  4),         &
                         SURFACES(no1)%Rdata(  5),         &
                         SURFACES(no1)%Rdata(  6), Fval,Dfdx)
!
!  ...CONE.............................................................
      case('Cone')
        call cone(     X,SURFACES(no1)%Rdata(1:3),         &
                         SURFACES(no1)%Rdata(4:6),         &
                         SURFACES(no1)%Rdata(  7), Fval,Dfdx)
!
      case default
           write(*,*)' surf: WRONG SURFACE TYPE'
           stop
      endselect
!
!  ...account for surface orientation
      Fval=Fval*nsign ; Dfdx(1:3)=Dfdx(1:3)*nsign
!
!
endsubroutine surf
