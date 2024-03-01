!------------------------------------------------------------------------------------
!> Purpose : parametrization of the hexahedron
!!
!! @param[in ] No     - hexa number
!! @param[in ] Eta    - reference coordinates of a point in the reference hexahedron
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives of the physical coordinates
!!
!! @revision Nov 12
!------------------------------------------------------------------------------------
!
subroutine hexa(No,Eta, X,Dxdeta)
!
      use GMP          , only : HEXAS,POINTS,CURVES,RECTANGLES
      use node_types   , only : BRIC
!
      implicit none
      integer               ,intent(in ) :: No
      real(8),dimension(3  ),intent(in ) :: Eta
      real(8),dimension(3  ),intent(out) :: X
      real(8),dimension(3,3),intent(out) :: Dxdeta
!------------------------------------------------------------------------------------
      real(8),dimension(  8) :: vshape
      real(8),dimension(3,8) :: dvshape
!
      integer :: i,nc,nr,lab,np,j
      integer :: iprint
!------------------------------------------------------------------------------------
!
      iprint=0
!
      select case(HEXAS(No)%Type)
!
!  ...trilinear Hexa.................................................................
      case('Linear')
!
!  .....check compatibility
        do i=1,12
          nc=abs(HEXAS(No)%EdgeNo(i))
          if (CURVES(nc)%Type.ne.'Seglin') then
            write(*,7002) No,i,nc,CURVES(nc)%Type
 7002       format(' hexa: incompatible edge! No,i,nc,type = ',i7,2x,i2,2x,i7,2x,a10)
            stop
          endif
        enddo
        do i=1,6
          if (HEXAS(No)%FigNo(i).eq.0) cycle
          call decode(HEXAS(No)%FigNo(i), nr,lab)
          if (RECTANGLES(nr)%Type.ne.'BilQua') then
            write(*,7003) No,i,nr,RECTANGLES(nr)%Type
 7003       format(' hexa: incompatible face! No,i,nr,type = ',i7,2x,i2,2x,i7,2x,a10)
            stop
          endif
        enddo
!
!  .....vertex shape functions
        call vshape3(BRIC,Eta, vshape,dvshape)
!
!  .....accumulate
        X(1:3)=0.d0 ; Dxdeta(1:3,1:3)=0.d0
        do i=1,8
          np=HEXAS(No)%VertNo(i)
          X(1:3) = X(1:3) + POINTS(np)%Rdata(1:3)*vshape(i)
          do j=1,3
            Dxdeta(j,1:3) = Dxdeta(j,1:3) + POINTS(np)%Rdata(j)*dvshape(1:3,i)
          enddo
        enddo
!
!  ....transfinite interpolation Hexa................................................
       case('TraHex','TrInHex') ; call hexa_TraHex(No,Eta, X,Dxdeta)
!
!  ....cylindrical coordinates Hexa
       case('CylHex') ; call hexa_CylHex(No,Eta, X,Dxdeta)
!
       case default
         write(*,7004) HEXAS(No)%Type
 7004    format(' hexa: unknown type! Type = ',a10)
         stop
       endselect
!
!
end subroutine hexa
!
!
!
!------------------------------------------------------------------------------------
!> Purpose : parametrization of the transfinite interpolation hexahedron
!!
!! @param[in ] No     - hexa number
!! @param[in ] Eta    - reference coordinates of a point in the reference hexahedron
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives of the physical coordinates
!!
!! @revision Nov 12
!------------------------------------------------------------------------------------
!
subroutine hexa_TraHex(No,Eta, X,Dxdeta)
!
      use GMP          , only : HEXAS,POINTS
      use element_data , only : BRICK_EDGE_TO_VERT,BRICK_FACE_TO_VERT
      use node_types   , only : BRIC
      implicit none
      integer               ,intent(in ) :: No
      real(8),dimension(3  ),intent(in ) :: Eta
      real(8),dimension(3  ),intent(out) :: X
      real(8),dimension(3,3),intent(out) :: Dxdeta
!------------------------------------------------------------------------------------
!  ...vertex shape functions
      real(8),dimension(  8)  :: vshape
      real(8),dimension(3,8)  :: dvshape
!  ...edge blending functions
      real(8),dimension(  12) :: blende
      real(8),dimension(3,12) :: dblende
!  ...face blending functions
      real(8),dimension(   6) :: blendf
      real(8),dimension(3, 6) :: dblendf
!  ...local edge coordinate
      real(8)                 :: te
      real(8),dimension(3)    :: dtedeta
!  ...curve
      real(8),dimension(3)    :: xe
      real(8),dimension(3)    :: dxedte
!  ...local face coodinates
      real(8),dimension(2)    :: tf
      real(8),dimension(2,3)  :: dtfdeta
!  ...quad
      real(8),dimension(3)    :: xf
      real(8),dimension(3,2)  :: dxfdtf
!
      integer                 :: i,j,np,nc,norient,nr,iv1,iv2,iv4
      integer                 :: iprint
!------------------------------------------------------------------------------------
!
      iprint=0
!
!  ...vertex shape functions
      call vshape3(BRIC,Eta, vshape,dvshape)
!
!  ...edge blending functions
      blende( 1)=(1.d0-Eta(2))*(1.d0-Eta(3))
      blende( 2)=      Eta(1) *(1.d0-Eta(3))
      blende( 3)=      Eta(2) *(1.d0-Eta(3))
      blende( 4)=(1.d0-Eta(1))*(1.d0-Eta(3))
      blende( 5)=(1.d0-Eta(2))*      Eta(3)
      blende( 6)=      Eta(1) *      Eta(3)
      blende( 7)=      Eta(2) *      Eta(3)
      blende( 8)=(1.d0-Eta(1))*      Eta(3)
      blende( 9)=(1.d0-Eta(1))*(1.d0-Eta(2))
      blende(10)=      Eta(1) *(1.d0-Eta(2))
      blende(11)=      Eta(1) *      Eta(2)
      blende(12)=(1.d0-Eta(1))*      Eta(2)
!
      dblende(1, 1)=  0.d0         ; dblende(2, 1)=-(1.d0-Eta(3)) ; dblende(3, 1)=-(1.d0-Eta(2))
      dblende(1, 2)= (1.d0-Eta(3)) ; dblende(2, 2)=  0.d0         ; dblende(3, 2)= -Eta(1)
      dblende(1, 3)=  0.d0         ; dblende(2, 3)= (1.d0-Eta(3)) ; dblende(3, 3)= -Eta(2)
      dblende(1, 4)=-(1.d0-Eta(3)) ; dblende(2, 4)=  0.d0         ; dblende(3, 4)=-(1.d0-Eta(1))
      dblende(1, 5)=  0.d0         ; dblende(2, 5)= -Eta(3)       ; dblende(3, 5)= (1.d0-Eta(2))
      dblende(1, 6)=  Eta(3)       ; dblende(2, 6)=  0.d0         ; dblende(3, 6)=  Eta(1)
      dblende(1, 7)=  0.d0         ; dblende(2, 7)=  Eta(3)       ; dblende(3, 7)=  Eta(2)
      dblende(1, 8)= -Eta(3)       ; dblende(2, 8)=  0.d0         ; dblende(3, 8)= (1.d0-Eta(1))
      dblende(1, 9)=-(1.d0-Eta(2)) ; dblende(2, 9)=-(1.d0-Eta(1)) ; dblende(3, 9)=  0.d0
      dblende(1,10)= (1.d0-Eta(2)) ; dblende(2,10)= -Eta(1)       ; dblende(3,10)=  0.d0
      dblende(1,11)=  Eta(2)       ; dblende(2,11)=  Eta(1)       ; dblende(3,11)=  0.d0
      dblende(1,12)= -Eta(2)       ; dblende(2,12)= (1.d0-Eta(1)) ; dblende(3,12)=  0.d0
!
!  ...face blending functions
      blendf(1)=1.d0-Eta(3) ; dblendf(1:3,1)=(/ 0.d0, 0.d0,-1.d0/)
      blendf(2)=     Eta(3) ; dblendf(1:3,2)=(/ 0.d0, 0.d0, 1.d0/)
      blendf(3)=1.d0-Eta(2) ; dblendf(1:3,3)=(/ 0.d0,-1.d0, 0.d0/)
      blendf(5)=     Eta(2) ; dblendf(1:3,5)=(/ 0.d0, 1.d0, 0.d0/)
      blendf(6)=1.d0-Eta(1) ; dblendf(1:3,6)=(/-1.d0, 0.d0, 0.d0/)
      blendf(4)=     Eta(1) ; dblendf(1:3,4)=(/ 1.d0, 0.d0, 0.d0/)
!
!------------------------------------------------------------------------------------
!     L I N E A R    I N T E R P O L A T I O N
!------------------------------------------------------------------------------------
      X(1:3)=0.d0 ; Dxdeta(1:3,1:3)=0.d0
      do i=1,8
        np=HEXAS(No)%VertNo(i)
        X(1:3) = X(1:3) + POINTS(np)%Rdata(1:3)*vshape(i)
        do j=1,3
          Dxdeta(j,1:3) = Dxdeta(j,1:3) +  &
                          POINTS(np)%Rdata(j)*dvshape(1:3,i)
        enddo
      enddo
!
!  ...printing
      if (iprint.eq.1) then
        write(*,*)'after vertices:'
        write(*,2000) X(1:3)
2000    format(' X = ',3(e12.5,2x))
        do j=1,3
          write(*,1000) j,Dxdeta(j,1:3)
1000      format(' j,Dxdeta(j,:) = ',i1,2x,3(e12.5,2x))
        enddo
      endif
!
!------------------------------------------------------------------------------------
!     E D G E S    C O N T R I B U T I O N S
!------------------------------------------------------------------------------------
      do i=1,12
        nc=HEXAS(No)%EdgeNo(i) ; norient=0
        if (nc.lt.0) then
          nc=-nc ; norient=1
        endif
!
!  .....get the edge vertices specifying the local edge orientation
        iv1=BRICK_EDGE_TO_VERT(1,i) ; iv2=BRICK_EDGE_TO_VERT(2,i)
!
!  .....project Eta onto the edge
        call proj_b2e(iv1,iv2,vshape,dvshape, te,dtedeta)
!
!  .....evaluate edge function
        call curve_local(nc,norient,te, xe,dxedte)
!
!  .....add (negative) edge contribution
        X(1:3) = X(1:3) - xe(1:3)*blende(i)
        do j=1,3
          Dxdeta(1:3,j) = Dxdeta(1:3,j)                    -  &
                          dxedte(1:3)*dtedeta(j)*blende(i) -  &
                          xe(1:3)*dblende(j,i)
        enddo
      enddo
!
!  ...printing
      if (iprint.eq.1) then
        write(*,*)'after edges:'
        write(*,2000) X(1:3)
        do j=1,3 ; write(*,1000) j,Dxdeta(j,1:3) ; enddo
      endif
!
!------------------------------------------------------------------------------------
!     F A C E S    C O N T R I B U T I O N S
!------------------------------------------------------------------------------------
      do i=1,6
        call decode(HEXAS(No)%FigNo(i), nr,norient)
!
!  .....get the face vertices specifying the local face orientation
        iv1=BRICK_FACE_TO_VERT(1,i)
        iv2=BRICK_FACE_TO_VERT(2,i)
        iv4=BRICK_FACE_TO_VERT(4,i)

!  .....project Eta onto face
        call proj_b2f(iv1,iv2,iv4,vshape,dvshape, tf,dtfdeta)
!
!  .....evaluate face function
        call recta_local(nr,norient,tf, xf,dxfdtf)
!
!  .....add face contribution
        X(1:3) = X(1:3) + xf(1:3)*blendf(i)
        do j=1,3
          Dxdeta(1:3,j) = Dxdeta(1:3,j)                           +  &
                          (dxfdtf(1:3,1)*dtfdeta(1,j)+               &
                           dxfdtf(1:3,2)*dtfdeta(2,j) )*blendf(i) +  &
                          xf(1:3)*dblendf(j,i)
        enddo
      enddo
!
!  ...printing
      if (iprint.eq.1) then
        write(*,*)'after faces:'
        write(*,2000) X(1:3)
        do j=1,3 ; write(*,1000) j,Dxdeta(j,1:3) ; enddo
      endif
!
!
end subroutine hexa_TraHex

!
!-----------------------------------------------------------------------
!
!   routine name       - hexa_CylHex
!
!-----------------------------------------------------------------------
!
!   latest revision    - Jan 17
!
!   purpose            - routine evaluates physical coordinates
!                        and its derivatives wrt to reference
!                        coordinates for a point on the image
!                        of a linear hexa through a global system
!                        of coordinates: x,y=rcos(\theta),z=rsin(\theta)
!                        and their  derivative wrt to reference
!                        coordinates

!
!   arguments :
!     in:
!               No     - a GMP hexahedron block number
!               Eta    - reference coordinates of a point
!                        in the reference hexa
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        to the reference coordinates
!
!-----------------------------------------------------------------------
!
   subroutine hexa_CylHex(No,Eta, X,Dxdeta)
!
      use control
      use GMP          , only : HEXAS,POINTS,NDIM
      use node_types   , only : BRIC
      implicit none
!----------------------------------------------------------------------
      integer,                 intent(in)  :: No
      real(8), dimension(3),   intent(in)  :: Eta(3)
      real(8), dimension(3),   intent(out) :: X(3)
      real(8), dimension(3,3), intent(out) :: Dxdeta(3,3)
!----------------------------------------------------------------------
!  ...vertex shape functions
      real(8), dimension(8)   :: vshape
      real(8), dimension(3,8) :: dvshape
!  ...cylindrical coordinates
      real(8)               :: r,theta,rp,thetap,thetaTmp
      real(8), dimension(3) :: drdeta,dthetadeta
!----------------------------------------------------------------------
!     misc.
      integer :: iprint,iv,np,i
      real(8) :: pi,twopi,costhet,sinthet
!----------------------------------------------------------------------
!
      iprint=0
!
      if ((HEXAS(No)%Type.ne.'CylHex'.or.(NDIM.ne.3))) then
        write(*,7001) HEXAS(No)%Type
 7001   format('Hexa_CylHex: WRONG HEXA TYPE = ',a10)
        stop 1
      endif
!
      if (iprint.eq.1) then
        write(*,7002) No,Eta
 7002   format('Hexa_CylHex: No,Eta = ',i4,2x,3f8.3)
      endif
!
      pi = acos(-1.d0)
      twopi = pi*2.d0
!
!  ...initiate
      X(1:3) = 0.d0; Dxdeta(1:3,1:3) = 0.d0
      r = 0.d0; drdeta(1:3) = 0.d0
      theta = 0.d0; dthetadeta(1:3) = 0.d0
!
!  ...vertex shape functions
      call vshape3(BRIC,Eta, vshape,dvshape)
!
!  ...interpolate in x,r,theta
      X(1)  = 0.d0; Dxdeta(1,1:3)     = 0.d0
      r     = 0.d0; Drdeta(1:3)       = 0.d0
      theta = 0.d0; Dthetadeta(1:3)   = 0.d0
      do iv=1,8
        np=HEXAS(No)%VertNo(iv)
        !  x-coord
        X(1) = X(1) + POINTS(np)%Rdata(1)*vshape(iv)
        Dxdeta(1,1:3) = Dxdeta(1,1:3)  &
                      + POINTS(np)%Rdata(1)*dvshape(1:3,iv)
        !  y,z-coords
        call coord_cart2polar(POINTS(np)%Rdata(2:3), rp,thetaTmp)
!  .....adjust the angle, if necessary
        select case(iv)
        case(1)
          thetap = thetaTmp
        case default
          if (thetaTmp-thetap.gt.pi) then
            thetap = thetaTmp - twopi
          elseif (thetaTmp-thetap.lt.-pi) then
            thetap = thetaTmp + twopi
          else
            thetap = thetaTmp
          endif
        end select
        ! write(*,*) 'thetap = ', thetap
        r     = r     + rp*vshape(iv)
        theta = theta + thetap*vshape(iv)
        Drdeta(1:3)     = Drdeta(1:3)     + rp*dvshape(1:3,iv)
        Dthetadeta(1:3) = Dthetadeta(1:3) + thetap*dvshape(1:3,iv)
      enddo

      costhet = cos(theta); sinthet = sin(theta)
      X(2) = r*costhet
      X(3) = r*sinthet
      Dxdeta(2,1:3) = Drdeta(1:3)*costhet - r*sinthet*Dthetadeta(1:3)
      Dxdeta(3,1:3) = Drdeta(1:3)*sinthet + r*costhet*Dthetadeta(1:3)
      if (iprint.eq.1) then
        write(*,*) 'theta = ', theta
        write(*,*) 'r = ', r
        write(*,7003) X
        write(*,7004) (Dxdeta(i,1:3),i=1,3)
 7003   format('Hexa_CylHex: ',/,'X      = ',3f8.3)
 7004   format('Dxdeta = ',3f8.3,/,  &
               '         ',3f8.3,/,  &
               '         ',3f8.3)
        call pause
      endif
!
!
   end subroutine hexa_CylHex
