!----------------------------------------------------------------------
!> @brief Curve kernel parametrization
!!
!> @param[in ] No      - curve number
!> @param[in ] T       - local curve coordinate
!> @param[in ] Norient - curve orientation
!> @param[out] X       - physical coordinates of the point
!> @param[out] Dxdt    - derivatives of the physical coordinates wrt
!!                       to the local curve coordinate
!!
!> @date Mar 2023
!----------------------------------------------------------------------
subroutine curveK(No,T,Norient, X,Dxdt)
!
      use control , only : GEOM_TOL
!
      implicit none
!
      integer,             intent(in ) :: No,Norient
      real(8),             intent(in ) :: T
      real(8),dimension(3),intent(out) :: X,Dxdt
!
      real(8) :: blend,dblend
      integer :: iprint
!----------------------------------------------------------------------
!
      iprint=0
!
!  ...check input
      if ((T.lt.GEOM_TOL).or.(T.gt.1.d0-GEOM_TOL)) then
        write(*,7001) No,T
 7001   format(' curveK: No,T = ',i5,2x,e12.5)
        call pause
      endif
!
!  ...evaluate the bubble function
      call curveB(No,T,Norient, X,Dxdt)
!
!  ...evaluate the kernel function
      blend=(1.d0-T)*T ; dblend=1.d0-2.d0*T
      Dxdt(1:3) = Dxdt(1:3)/blend - X(1:3)*dblend/blend**2
      X(1:3) = X(1:3)/blend
!
!
end subroutine curveK
!
!
!
!----------------------------------------------------------------------
!> @brief curve bubble parametrization
!!
!> @param[in ] No      - curve number
!> @param[in ] T       - local curve coordinate
!> @param[in ] Norient - curve orientation
!> @param[out] X       - physical coordinates of the point
!> @param[out] Dxdt    - derivatives of the physical coordinates wrt
!!                       to the local curve coordinate
!!
!> @date Nov 12
!----------------------------------------------------------------------
subroutine curveB(No,T,Norient, X,Dxdt)
!
      use element_data , only : EDGE_L2G
      use GMP          , only : POINTS , CURVES
      use control      , only : GEOM_TOL
!
      implicit none
      integer,             intent(in ) :: No,Norient
      real(8),             intent(in ) :: T
      real(8),dimension(3),intent(out) :: X,Dxdt
!
      real(8),dimension(3,2) :: xv
      real(8) :: t_aux,smax
      integer :: i,iv,np,ivar
      integer :: icheck
      integer :: iprint
!----------------------------------------------------------------------
!
      icheck=1
      iprint=0
!
!  ...collect curve endpoints wrt LOCAL coordinate T
      do i=1,2
        iv=EDGE_L2G(i,Norient) ; np=CURVES(No)%EndPoNo(iv)
        xv(1:3,i)=POINTS(np)%Rdata(1:3)
!
!  @@@@ CHECK CONSISTENCY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (icheck.ne.0) then
!  .......compute parameterization at endpoint
          t_aux=(i-1)*1.d0 ; call curve_local(No,Norient,t_aux, x,dxdt)
!  .......compare to vertex coordinates
          smax=0.d0
          do ivar=1,3 ; smax=max(smax,abs(x(ivar)-xv(ivar,i))) ; enddo
          if (smax.gt.GEOM_TOL) then
            write(*,7001) No,i,smax
 7001       format(' curveB: No,i,smax = ',i5,i2,e12.5)
            write(*,7002) x(1:3)
 7002       format(' x  = ',3(e12.5,2x))
            write(*,7003) xv(1:3,i)
 7003       format(' xv = ',3(e12.5,2x))
          endif
        endif
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
      enddo
!
!  ...curve parametrization
      call curve_local(No,Norient,T, X,Dxdt)
!
!  ...subtract the linear interpolant
      X(   1:3) = X(   1:3) - ( xv(1:3,1)*(1.d0-T) + xv(1:3,2)*T )
      Dxdt(1:3) = Dxdt(1:3) - ( xv(1:3,2) - xv(1:3,1))
!
!
end subroutine curveB
!
!
!
!----------------------------------------------------------------------
!> @brief triangle kernel parameterization
!!
!> @param[in ] No      - triangle number
!> @param[in ] T       - local triangle coordinates
!> @param[in ] Norient - triangle orientation
!> @param[out] X       - physical coordinates of the point
!> @param[out] Dxdt    - derivatives of the physical coordinates wrt
!!                       to the local curve coordinate
!!
!> @date Nov 12
!----------------------------------------------------------------------
subroutine trianK(No,T,Norient, X,Dxdt)
!
      implicit none
      integer,               intent(in ) :: No,Norient
      real(8),dimension(2  ),intent(in ) :: T
      real(8),dimension(3  ),intent(out) :: X
      real(8),dimension(3,2),intent(out) :: Dxdt
!
      real(8) :: blend
      real(8) :: dblend(2)
      integer :: i
      integer :: iprint
!
!----------------------------------------------------------------------
!
      iprint=0
!
!  ...evaluate blending shape functions
      blend = (1.d0 - T(1) - T(2))*T(1)*T(2)
      dblend(1) = -T(1)*T(2) + (1.d0 - T(1) - T(2))*T(2)
      dblend(2) = -T(1)*T(2) + (1.d0 - T(1) - T(2))*T(1)
!
!  ...evaluate the bubble
      call trianB(No,T,Norient, X,Dxdt)
!
!  ...divide by the blending function
      do i=1,2
        Dxdt(1:3,i) = Dxdt(1:3,i)/blend - X(1:3)*dblend(i)/blend**2
      enddo
      X(1:3) = X(1:3)/blend
!
end subroutine trianK
!
!
!
!----------------------------------------------------------------------
!> @brief triangle bubble parameterization
!!
!> @param[in ] No      - triangle number
!> @param[in ] T       - local triangle coordinates
!> @param[in ] Norient - triangle orientation
!> @param[out] X       - physical coordinates of the point
!> @param[out] Dxdt    - derivatives of the physical coordinates wrt
!!                       to the local curve coordinate
!!
!> @date Nov 12
!----------------------------------------------------------------------
subroutine trianB(No,T,Norient, X,Dxdt)
!
      use control , only : GEOM_TOL
      use GMP
      use element_data
      use node_types, only : TRIA
!
      implicit none
      integer,               intent(in ) :: No,Norient
      real(8),dimension(2  ),intent(in ) :: T
      real(8),dimension(3  ),intent(out) :: X
      real(8),dimension(3,2),intent(out) :: Dxdt
!
!  ...vertex shape functions
      real(8),dimension(  4) :: shapH
      real(8),dimension(2,4) :: dshapH
!
      real(8),dimension(2  ) :: eta,eta_aux,dtedeta,detadte,dblend
      real(8),dimension(3  ) :: x_aux,xc,dxcdte,dxdte
      real(8),dimension(2,2) :: detadt
      real(8),dimension(3,2) :: dxdeta,dxcdeta,dxdeta_aux
!
      real(8),dimension(3,3) :: xv
      real(8) :: smax,dmax,te,blend
      integer :: i,j,ivar,np,nc,norientc,iv1,iv2
!
!     consistency of parameterizations:
!     0 - no checking
!     1 - checking at vertices
!     2 - checking at vertices and edges
      integer :: icheck
      integer :: iprint
!-----------------------------------------------------------------------
!
      iprint=0
      icheck=2
!
 10   continue
!
!  ...if bubble is not needed, return
      if ( (TRIANGLES(No)%Type.eq.'TransTri').or.         &
           (TRIANGLES(No)%Type.eq.'PlaneTri')     ) then
        X(1:3)=0.d0 ; Dxdt(1:3,1:2)=0.d0
!        write(*,*)'trianB: warning! Bubble not needed!'
        write(*,6000) No,TRIANGLES(No)%Type
 6000   format(' No,Type = ',i6,2x,a10)
!        call pause
        return
      endif
!
!  ...master element, baricentric, and physical coordinates
      call local2global(TRIA,T,Norient, eta,detadt)
      call vshape2(TRIA,eta, shapH,dshapH)
      call trian(No,eta, X,dxdeta)
!
!----------------------------------------------------------------------
!  STEP 1 : subtract linear interpolant                               |
!----------------------------------------------------------------------
!  ...loop through vertices
      do i=1,3
!
!  .....get triangle vertex
        np=TRIANGLES(No)%VertNo(i) ; xv(1:3,i)=POINTS(np)%Rdata(1:3)
!
!  @@@@ CHECK CONSISTENCY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (icheck.ge.1) then
!  .......compute parameterization at vertex
          eta_aux(1:2)=TRIAN_COORD(1:2,i)
          call trian(No,eta_aux, x_aux,dxdeta_aux)
!  .......compare to vertex coordinates
          smax=0.d0
          do ivar=1,3
            smax=max(smax,abs(x_aux(ivar) - xv(ivar,i)))
          enddo
          if (smax.gt.GEOM_TOL) then
            write(*,7001) No,i,smax
 7001       format(' trianB: No,i,smax = ',i5,i2,e12.5)
            write(*,7002) x_aux(1:3)
 7002       format(' x  = ',3(e12.5,2x))
            write(*,7003) xv(1:3,i)
 7003       format(' xv = ',3(e12.5,2x))
            call pause
            iprint=1 ; goto 10
          endif
        endif
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  .....subtract linear interpolant
        X(1:3)=X(1:3) - xv(1:3,i)*shapH(i)
        do j=1,2
          dxdeta(1:3,j)=dxdeta(1:3,j) - xv(1:3,i)*dshapH(j,i)
        enddo

!  ...end of loop through vertices
      enddo
!
!  ...printing statement
      if (iprint.eq.1) then
        write(*,*) 'trianB: AFTER VERTICES'
        do ivar=1,3
          write(*,8000) ivar,X(ivar),dxdeta(ivar,1:2)
 8000     format(' i = ',i1,' ; X(i),dxdeta(i,1:2) = ',3(e12.5,2x))
        enddo
      endif
!
!----------------------------------------------------------------------
!  STEP 2 : subtract edge contributions                               |
!----------------------------------------------------------------------
!  ...loop through edges
      do i=1,3
        nc=TRIANGLES(No)%EdgeNo(i) ; norientc=0
        if (nc.lt.0) then ; nc=-nc ; norientc=1 ; endif
!
!  .....no contribution needed for straight segments
        if (CURVES(nc)%Type.eq.'Seglin') cycle
!
!  .....get edge coordinate
        iv1=TRIAN_EDGE_TO_VERT(1,i) ; iv2=TRIAN_EDGE_TO_VERT(2,i)
        call proj_t2e(iv1,iv2,shapH(1:3),dshapH(1:2,1:3), te,dtedeta)
!
!  .....no contribution needed at endpoints
        if ((te.lt.GEOM_TOL).or.(te.gt.(1.d0-GEOM_TOL))) cycle
!
!  @@@@ CHECK CONSISTENCY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (icheck.ge.2) then
!  .......(1) curve bubble from curve parameterization
          call curveB(nc,te,norientc, xc,dxcdte)
!
!  .......(2) curve bubble from triangle parameterization
          call edge_param(TRIA,i,te, eta_aux,detadte)
          call trian(No,eta_aux, x_aux,dxdeta_aux)
          dxdte(1:3)=dxdeta_aux(1:3,1)*detadte(1) + &
                     dxdeta_aux(1:3,2)*detadte(2)
          x_aux(1:3)=x_aux(1:3) - ( xv(1:3,iv1)*(1.d0-te) +   &
                                    xv(1:3,iv2)*te          )
          dxdte(1:3)=dxdte(1:3) - (xv(1:3,iv2) -  xv(1:3,iv1))
!
!  .......compare (1) & (2)
          smax=0.d0 ; dmax=0.d0
          do ivar=1,3
            smax=max(smax,abs(x_aux(ivar) - xc(ivar)))
            dmax=max(dmax,abs(dxdte(ivar) - dxcdte(ivar)))
          enddo
          if ((smax.gt.GEOM_TOL).or.(dmax.gt.GEOM_TOL)) then
            write(*,7004) No,i,smax,dmax
 7004       format(' trianB: No,i,smax,dmax = ',i5,i2,2e12.5)
            iprint=1 ; goto 10
          endif
        endif
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  .....construct blending function
        blend      = shapH(iv1)*shapH(iv2)
        dblend(1:2)=dshapH(1:2,iv1)*shapH(iv2) + shapH(iv1)*dshapH(1:2,iv2)

!  .....compute curve kernel
        call curveK(nc,te,norientc, xc,dxcdte)
        dxcdeta(1:3,1)=dxcdte(1:3)*dtedeta(1)
        dxcdeta(1:3,2)=dxcdte(1:3)*dtedeta(2)

!  .....subtract edge bubble
        X(1:3) = X(1:3) - xc(1:3)*blend
        do j=1,2
          dxdeta(1:3,j) = dxdeta(1:3,j) - dxcdeta(1:3,j)* blend    &
                                        -      xc(1:3  )*dblend(j)
        enddo
!
!  .....printing
        if (iprint.eq.1) then
          write(*,8001) i
 8001     format(' trianB: AFTER EDGE ',i1)
          do ivar=1,3
            write(*,8000) ivar,X(ivar),dxdeta(ivar,1:2)
          enddo
          call pause
        endif
!  ...end of loop through edges
      enddo
!
!  ...account for the orientation
      do j=1,2
        Dxdt(1:3,j) = dxdeta(1:3,1)*detadt(1,j) + &
                      dxdeta(1:3,2)*detadt(2,j)
      enddo
!
!  ...printing
      if (iprint.eq.1) then
        write(*,*) 'trianB: FINAL'
        do ivar=1,3
          write(*,8002) ivar,X(ivar),Dxdt(ivar,1:2)
 8002     format(' i = ',i1,' ; X(i),dxdt(i,1:2) = ',3(e12.5,2x))
        enddo
        call pause
      endif
!
!
end subroutine trianB
!
!
!
!----------------------------------------------------------------------
!> @brief rectangle bubble parameterization
!!
!> @param[in ] No      - rectangle number
!> @param[in ] T       - local rectangle coordinates
!> @param[in ] Norient - rectangle orientation
!> @param[out] X       - physical coordinates of the point
!> @param[out] Dxdt    - derivatives of the physical coordinates wrt
!!                       to the local curve coordinate
!!
!> @date Nov 12
!----------------------------------------------------------------------
subroutine rectaB(No,T,Norient, X,Dxdt)
!
      use control , only : GEOM_TOL
      use GMP
      use element_data
      use node_types, only : QUAD
!
      implicit none
      integer,               intent(in ) :: No,Norient
      real(8),dimension(2  ),intent(in ) :: T
      real(8),dimension(3  ),intent(out) :: X
      real(8),dimension(3,2),intent(out) :: Dxdt
!
!  ...vertex shape functions
      real(8),dimension(  4) :: shapH
      real(8),dimension(2,4) :: dshapH
!
      real(8),dimension(2  ) :: eta,eta_aux,dtedeta,detadte,dblend
      real(8),dimension(3  ) :: x_aux,xc,dxcdte,dxdte
      real(8),dimension(2,2) :: detadt
      real(8),dimension(3,2) :: dxdeta,dxcdeta,dxdeta_aux
!
      real(8),dimension(3,4) :: xv
      real(8) :: smax,dmax,te,blend
      integer :: i,j,ivar,np,nc,norientc,iv1,iv2
!
!     consistency of parameterizations:
!     0 - no checking
!     1 - checking at vertices
!     2 - checking at vertices and edges
      integer :: icheck
      integer :: iprint
!
!-----------------------------------------------------------------------
!
      iprint=0
      icheck=2
!
 10   continue
!
!  ...if bubble is not needed, return
      if (RECTANGLES(No)%Type.ne.'PTIRec') then
        X(1:3)=0.d0 ; Dxdt(1:3,1:2)=0.d0
!        write(*,*)'rectaB: warning! Bubble is not needed!'
        write(*,6000) No,RECTANGLES(No)%Type
 6000   format(' No,Type = ',i6,2x,a10)
!        call pause
        return
      endif
!
!  ...master element, baricentric, and physical coordinates
      call local2global(QUAD,T,Norient, eta,detadt)
      call vshape2(QUAD,eta, shapH,dshapH)
      call recta(No,eta, X,dxdeta)
!
!----------------------------------------------------------------------
!  STEP 1 : subtract bilinear interpolant                             |
!----------------------------------------------------------------------
!  ...loop over vertices
      do i=1,4
!
!  .....get rectangle vertex
        np=RECTANGLES(No)%VertNo(i) ; xv(1:3,i)=POINTS(np)%Rdata(1:3)
!
!  @@@@ CHECK CONSISTENCY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (icheck.ge.1) then
!  .......compute parameterization at vertex
          eta_aux(1:2)=QUADR_COORD(1:2,i)
          call recta(No,eta_aux, x_aux,dxdeta_aux)
!  .......compare to vertex coordinates
          smax=0.d0
          do ivar=1,3
            smax=max(smax,abs(x_aux(ivar) - xv(ivar,i)))
          enddo
          if (smax.gt.GEOM_TOL) then
            write(*,7001) No,i,smax
 7001       format(' rectaB: No,i,smax = ',i5,i2,e12.5)
            write(*,7002) x_aux(1:3)
 7002       format(' x  = ',3(e12.5,2x))
            write(*,7003) xv(1:3,i)
 7003       format(' xv = ',3(e12.5,2x))
            call pause
            iprint=1 ; goto 10
          endif
        endif
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  .....subtract bilinear interpolant
        X(1:3)=X(1:3) - xv(1:3,i)*shapH(i)
        do j=1,2
          dxdeta(1:3,j)=dxdeta(1:3,j) - xv(1:3,i)*dshapH(j,i)
        enddo
!
!  ...end of loop through vertices
      enddo
!
!  ...printing statement
      if (iprint.eq.1) then
        write(*,*) 'rectaB: AFTER VERTICES'
        do ivar=1,3
          write(*,8000) ivar,X(ivar),dxdeta(ivar,1:2)
 8000     format(' i = ',i1,' ; X(i),dxdeta(i,1:2) = ',3(e12.5,2x))
        enddo
      endif
!
!----------------------------------------------------------------------
!  STEP 2 : subtract edge contributions                               |
!----------------------------------------------------------------------
!  ...loop over edges
      do i=1,4
        nc=RECTANGLES(No)%EdgeNo(i) ; norientc=0
        if (nc.lt.0) then ;  nc=-nc ; norientc=1 ; endif
!
!  .....no contribution needed for straight segmentes
        if (CURVES(nc)%Type.eq.'Seglin') cycle
!
!  .....get edge coordinate
        iv1=QUADR_EDGE_TO_VERT(1,i) ; iv2=QUADR_EDGE_TO_VERT(2,i)
        call proj_r2e(iv1,iv2,shapH,dshapH, te,dtedeta)
!
!  .....no contribution needed at endpoints
        if ((te.lt.GEOM_TOL).or.(te.gt.(1.d0-GEOM_TOL))) cycle
!
!  @@@@ CHECK CONSISTENCY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (icheck.ge.2) then
!  .......(1) curve bubble from curve parameterization
          call curveB(nc,te,norientc, xc,dxcdte)
!
!  .......(2) curve bubble from triangle parameterization
          call edge_param(QUAD,i,te, eta_aux,detadte)
          call recta(No,eta_aux, x_aux,dxdeta_aux)
          dxdte(1:3)=dxdeta_aux(1:3,1)*detadte(1) + &
                     dxdeta_aux(1:3,2)*detadte(2)
          x_aux(1:3)=x_aux(1:3) - ( xv(1:3,iv1)*(1.d0-te) +   &
                                    xv(1:3,iv2)*te          )
          dxdte(1:3)=dxdte(1:3) - (xv(1:3,iv2) -  xv(1:3,iv1))
!
!  .......compare (1) & (2)
          smax=0.d0 ; dmax=0.d0
          do ivar=1,3
            smax=max(smax,abs(x_aux(ivar) - xc(ivar)))
            dmax=max(dmax,abs(dxdte(ivar) - dxcdte(ivar)))
          enddo
          if ((smax.gt.GEOM_TOL).or.(dmax.gt.GEOM_TOL)) then
            write(*,7004) No,i,smax,dmax
 7004       format(' rectaB: No,i,smax,dmax = ',i5,i2,2e12.5)
            iprint=1 ; goto 10
          endif
        endif
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  .....construct blending function
        select case(i)
        case(1) ; blend=1.d0 - eta(2) ; dblend(1:2)=(/ 0.d0, -1.d0/)
        case(2) ; blend=eta(1)        ; dblend(1:2)=(/ 1.d0,  0.d0/)
        case(3) ; blend=eta(2)        ; dblend(1:2)=(/ 0.d0,  1.d0/)
        case(4) ; blend=1.d0 - eta(1) ; dblend(1:2)=(/-1.d0,  0.d0/)
        endselect
!
!  .....compute curve bubble
        call curveB(nc,te,norientc, xc,dxcdte)
        dxcdeta(1:3,1)=dxcdte(1:3)*dtedeta(1)
        dxcdeta(1:3,2)=dxcdte(1:3)*dtedeta(2)
!
!  .....subtract edge bubble
        X(1:3) = X(1:3) - xc(1:3)*blend
        do j=1,2
          dxdeta(1:3,j) = dxdeta(1:3,j) - dxcdeta(1:3,j)* blend    &
                                        -      xc(1:3  )*dblend(j)
        enddo
!
!  .....printing
        if (iprint.eq.1) then
          write(*,8001) i
 8001     format(' rectaB: AFTER EDGE ',i1)
          do ivar=1,3
            write(*,8000) ivar,X(ivar),dxdeta(ivar,1:2)
          enddo
          call pause
        endif
!  ...loop over edges
      enddo
!
!  ...account for the orientation
      do j=1,2
        Dxdt(1:3,j) = dxdeta(1:3,1)*detadt(1,j) + &
                      dxdeta(1:3,2)*detadt(2,j)
      enddo
!
!  ...printing
      if (iprint.eq.1) then
        write(*,*) 'rectaB: FINAL'
        do ivar=1,3
          write(*,8002) ivar,X(ivar),Dxdt(ivar,1:2)
 8002     format(' i = ',i1,' ; X(i),dxdt(i,1:2) = ',3(e12.5,2x))
        enddo
        call pause
      endif
!
!
end subroutine rectaB
!
!----------------------------------------------------------------------
!> @date Mar 2023
!----------------------------------------------------------------------
subroutine rectaB_back(No,T,Norient, X,Dxdt)
!
      use control
      use GMP
      use element_data
!
      implicit none
!
      integer :: No,Norient
      real(8) :: T(2),X(3),Dxdt(3,2)
!
!  ...derivatives wrt rectangle coordinates
      real(8) :: dxds(3,2)
!
!  ...rectangle coordinates
      real(8) :: s(2),dsdt(2,2,0:7)
!
!  ...vertex shape functions
      real(8) :: shapH(4),dshapH(2,4)
!
!  ...vertex coordinates
      real(8) :: sv(2,4), xv(3,4)
!
!  ...edge projection, blending function
      real(8) :: dseds(2), dblend(2)
!
!  ...work space
      real(8) :: sw(2),dswdse(2), &
                 xw(3),dxwds(3,2),dxwdse(3),xc(3),dxcdse(3)
!
      real(8) :: blend,dmax,se,smax
      integer :: ie,is,iv,iv1,iv2,ivar,nc,norientc,nss,nv
!
      integer :: iprint
!
      data dsdt / 1.d0,  0.d0,  0.d0,  1.d0, &
                  0.d0, -1.d0,  1.d0,  0.d0, &
                 -1.d0,  0.d0,  0.d0, -1.d0, &
                  0.d0,  1.d0, -1.d0,  0.d0, &
                  0.d0,  1.d0,  1.d0,  0.d0, &
                 -1.d0,  0.d0,  0.d0,  1.d0, &
                  0.d0, -1.d0, -1.d0,  0.d0, &
                  1.d0,  0.d0,  0.d0, -1.d0 /
!
!  ...coordinates of master triangle vertices
      data sv /0.d0,0.d0, 1.d0,0.d0, 1.d0,1.d0, 0.d0,1.d0/
!
!-----------------------------------------------------------------------
!
      iprint=0
!
 10   continue
!
!  ...rectangle bubbles are needed only for 'PTIRec'
      if (RECTANGLES(No)%Type.ne.'PTIRec') then
        X(1:3) = 0.d0 ; Dxdt(1:3,1:2) = 0.d0
        write(*,6000) No,RECTANGLES(No)%Type
 6000   format(' rectaB: No,Type = ',i7,2x,a10)
!        write(*,*)'Warning: bubble is not needed!'
!        call pause
        return
      endif
!
!  ...account for orientation
      select case(Norient)
      case (0); s(1) = T(1)        ; s(2) = T(2)
      case (1); s(1) = T(2)        ; s(2) = 1.d0 - T(1)
      case (2); s(1) = 1.d0 - T(1) ; s(2) = 1.d0 - T(2)
      case (3); s(1) = 1.d0 - T(2) ; s(2) = T(1)
      case (4); s(2) = T(1)        ; s(1) = T(2)
      case (5); s(2) = T(2)        ; s(1) = 1.d0 - T(1)
      case (6); s(2) = 1.d0 - T(1) ; s(1) = 1.d0 - T(2)
      case (7); s(2) = 1.d0 - T(2) ; s(1) = T(1)
      case default
        write(*,7000)Norient
 7000   format(' rectaB: UNKNOWN ORIENTATION = ',i1)
        stop
      endselect
!
!  ...printing statement
      if (iprint .eq. 1) then
        write(*,*) '---------------------------------------------------'
        nss = RECTANGLES(No)%Idata(1)
        write(*,*)'surface = ',nss
        write(*,7004) No,T(1:2),Norient,s(1:2)
 7004   format(' rectaB: No,T,Norient,s = ',i6,2e12.5,i3,2x,2e12.5)
      endif
!
!  ...evaluate bilinear shape functions...
      shapH(1) = (1.d0 - s(1)) * (1.d0 - s(2))
      shapH(2) = s(1)          * (1.d0 - s(2))
      shapH(3) = s(1)          *  s(2)
      shapH(4) = (1.d0 - s(1)) *  s(2)
!  ...and their derivatives
      dshapH(1,1) = - (1.d0 - s(2))
      dshapH(2,1) = - (1.d0 - s(1))
      dshapH(1,2) =   (1.d0 - s(2))
      dshapH(2,2) = - s(1)
      dshapH(1,3) =   s(2)
      dshapH(2,3) =   s(1)
      dshapH(1,4) = - s(2)
      dshapH(2,4) =   (1.d0 - s(1))
!
!  ...start by evaluating rectangle parametrization
      call recta(No,s(1:2), X(1:3),dxds(1:3,1:2))
      if (iprint .eq. 1) then
        write(*,*)'rectaB: ORIGINAL X,dxds'
        do ivar = 1, 3
          write(*,7011) X(ivar),dxds(ivar,1:2)
 7011     format(e12.5,2x,2e12.5)
        enddo
      endif
!
!  ...get the vertex coordinates
      do iv=1,4
        nv = RECTANGLES(No)%VertNo(iv)
        xv(1:3,iv) = POINTS(nv)%Rdata(1:3)
      enddo
!
!----------------------------------------------------------------------
!     S U B T R A C T    V E R T E X    I N T E R P O L A T I O N
!----------------------------------------------------------------------
!
!  ...loop over vertices
      do iv=1,4
!
!======================================================================
!  check consistency of parametrizations btw vertices an 'recta'      |
!  routine.                                                           |
!----------------------------------------------------------------------
!  .....call 'recta' routine at master triangle vertices              !
        call recta(No,sv(1:2,iv), xw,dxwds)                           !
        smax = 0.d0                                                   !
!  .....accumulate error                                              !
        do ivar=1,3                                                   !
          smax = max(smax,abs(xw(ivar) - xv(ivar,iv)))                !
        enddo                                                         !
!  .....check whether GEOM_TOL is exceeded                            !
        if (smax .gt. GEOM_TOL) then                                  !
          write(*,7001) No,iv,smax                                    !
 7001     format(' rectaB: No,iv,smax = ',i5,i2,e12.5)                 !
          write(*,*) 'xw = ',xw                                       !
          write(*,*) 'xv = ',xv(1:3,iv)                               !
          call pause                                                  !
          iprint=1                                                    !
          goto 10                                                     !
        endif                                                         !
!======================================================================
!
!  .....subtract bilinear interpolant
        X(1:3) = X(1:3) - xv(1:3,iv)*shapH(iv)
        do is=1,2
          dxds(1:3,is) = dxds(1:3,is) - xv(1:3,iv)*dshapH(is,iv)
        enddo
!
!  ...loop over vertices
      enddo
!
!  ...printing statement
      if (iprint.eq.1) then
        write(*,*) 'rectaB: AFTER VERTICES  X,dxds = '
        do ivar=1,3
          write(*,7011) X(ivar),dxds(ivar,1:2)
        enddo
      endif
!
!----------------------------------------------------------------------
!     S U B T R A C T    E D G E    B U B B L E S
!----------------------------------------------------------------------
!
!  ...loop over edges
      do ie=1,4
!  .....get the curve number
        nc = RECTANGLES(No)%EdgeNo(ie) ; norientc=0
        if (nc.lt.0) then
          nc = -nc ;  norientc = 1
        endif
!  .....skip if straight segment
        if (CURVES(nc)%Type.eq.'Seglin') cycle
!  .....get the edge vertices specifying the local edge orientation
        iv1 = QUADR_EDGE_TO_VERT(1,ie) ;  iv2 = QUADR_EDGE_TO_VERT(2,ie)
!  .....project s onto the edge
        call proj_r2e(iv1,iv2,shapH,dshapH, se,dseds)
        if (iprint.eq.1) then
          write(*,9001) ie,se
 9001     format(' rectaB: ie, se = ',I2,' ; ',E12.5)
        endif
!  .....edge contributions (bubbles) are null at the end points
        if ((se.lt.GEOM_TOL).or.(se.gt.(1.d0-GEOM_TOL))) cycle
        if (iprint.eq.1) then
          write(*,7003) ie,nc,CURVES(nc)%Type
 7003     format(' rectaB: ie,nc,Type = ',i2,2x,i8,2x,a5)
          call print_GMP
        endif
!
!======================================================================
!  check consistency of curve and triangle parametrizations by        |
!  comparing edge bubble with triangle bubble.                        |
!----------------------------------------------------------------------
!       1st TERM OF COMPARISON: xc, dxcdse                            !
!  .....evaluate curve bubble at projected point se                   !
        call curveB(nc,se,norientc, xc,dxcdse)                        !
!----------------------------------------------------------------------
!       2nd TERM OF COMPARISON: xw, dxwdse                            !
!  .....evaluate local edge parameterization of master triangle       !
        call edge_param(QUAD,ie,se, sw,dswdse)                        !
!  .....evaluate rectangle parameterization                           !
        call recta(No,sw, xw,dxwds)                                   !
!  .....compute derivative wrt edge parameter                         !
        dxwdse(1:3) = dxwds(1:3,1)*dswdse(1) + dxwds(1:3,2)*dswdse(2) !
!  .....evaluate the edge bubble                                      !
        xw(1:3) = xw(1:3)                           &                 !
                - (xv(1:3,iv1)*(1.d0 - se) + xv(1:3,iv2)*se)          !
        dxwdse(1:3) = dxwdse(1:3)                            &        !
                - (xv(1:3,iv2) -  xv(1:3,iv1))                        !
!----------------------------------------------------------------------
!       COMPARE                                                       !
        smax = 0.d0;  dmax = 0.d0                                     !
        do ivar = 1, 3                                                !
          smax = max(smax,abs(xw(ivar) - xc(ivar)))                   !
          dmax = max(dmax,abs(dxcdse(ivar) - dxwdse(ivar)))           !
        enddo                                                         !
        if ((smax.gt.GEOM_TOL).or.(dmax.gt.GEOM_TOL)) then            !
          write(*,7002) No,ie,smax,dmax                               !
 7002     format(' rectaB: No,ie,smax,dmax = ',i5,i2,2e12.5)          !
!!          call pause                                                !
!!!          if (iprint .eq. 1) call my_tests                         !
!!!          iprint = 1                                               !
!!!          write(*,9000) No                                         !
 9000     format(' trianB: rerunning rountine for triangle No = ', &  !
                           I6,' with printing flag on.')              !
!!!          call pause                                               !
!!!          go to 10                                                 !
        endif                                                         !
!======================================================================
!
!  .....construct blending function
        select case(ie)
        case(1); blend = 1.d0 - s(2) ; dblend(1:2) = (/ 0.d0, -1.d0/)
        case(2); blend = s(1)        ; dblend(1:2) = (/ 1.d0,  0.d0/)
        case(3); blend = s(2)        ; dblend(1:2) = (/ 0.d0,  1.d0/)
        case(4); blend = 1.d0 - s(1) ; dblend(1:2) = (/-1.d0,  0.d0/)
        endselect
!
!  .....compute curve bubble
        call curveB(nc,se,norientc, xc,dxcdse)
!
!  .....subtract edge bubble
        X(1:3) = X(1:3) - xc(1:3)*blend
        do is = 1, 2
          dxds(1:3,is) = dxds(1:3,is)  &
                       - dxcdse(1:3)*dseds(is)*blend &
                       - xc(1:3)*dblend(is)
        enddo
!
!  .....printing
        if (iprint .eq. 1) then
          write(*,*)'rectaB: X,dxds AFTER EDGE = ',ie
          do ivar = 1, 3
            write(*,7011) X(ivar),dxds(ivar,1:2)
          enddo
          call pause
        endif
!
!  ...loop over edges
      enddo
!
!  ...printing statement
      if (iprint.eq.1) then
        write(*,*) 'rectaB: AFTER EDGES  X,dxds = '
        do ivar=1,3
          write(*,7011) X(ivar),dxds(ivar,1:2)
        enddo
        call pause
      endif
!
!---------------------------------------------------------------------
!  ...account for the orientation
      do is=1,2
        Dxdt(1:3,is) = dxds(1:3,1)*dsdt(1,is,Norient) &
                     + dxds(1:3,2)*dsdt(2,is,Norient)
      enddo
      if (iprint .eq. 1) then
        write(*,*)'rectaB: FINAL  X,Dxdt = '
        do ivar = 1, 3
          write(*,7011) X(ivar),Dxdt(ivar,1:2)
        enddo
        call pause
      endif
!
!
end subroutine rectaB_back
