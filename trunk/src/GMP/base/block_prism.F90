!------------------------------------------------------------------------------------
!> @brief Controls parametrization of prisms
!!
!! @param[in ] No     - prism number
!! @param[in ] Eta    - reference coordinates of a point in the reference prism
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives of the physical coordinates
!!
!> @date Mar 2023
!------------------------------------------------------------------------------------
subroutine prism(No,Eta, X,Dxdeta)
!
      use GMP
      use node_types, only: PRIS
!
      implicit none
      integer,               intent(in ) :: No
      real(8),dimension(3  ),intent(in ) :: Eta
      real(8),dimension(3  ),intent(out) :: X
      real(8),dimension(3,3),intent(out) :: Dxdeta
!
!  ...vertex coordinates
      real(8),dimension(3,6) :: xv
!  ...vertex shape functions
      real(8),dimension(  8) :: vshape
      real(8),dimension(3,8) :: dvshape
      integer :: iv,np,k,i
      integer :: iprint
!------------------------------------------------------------------------------------
!
      iprint=0
!
! ....printing statement
      if (iprint.eq.1) then
        write(*,7001) No,PRISMS(No)%Type
 7001   format('prism: No,PRISMS(No)%Type = ',i8,2x,a10)
        call pause
      endif
!
! ....select prism type
      select case(PRISMS(No)%Type)
!
!------------------------------------------------------------------------------------
!     L I N E A R
!------------------------------------------------------------------------------------
      case('Linear')
!  .....get the vertex coordinates
        do iv=1,6
          np=PRISMS(No)%VertNo(iv) ; call pointr(np, xv(1:3,iv))
        enddo
!  .....evaluate vertex shape functions
        call vshape3(PRIS,Eta, vshape,dvshape)
!  .....accumulate
        X(1:3)=0.d0; Dxdeta(1:3,1:3) = 0.d0
        do k=1,6
          X(1:3) = X(1:3) + xv(1:3,k)*vshape(k)
          do i=1,3
            Dxdeta(1:3,i) = Dxdeta(1:3,i) + xv(1:3,k)*dvshape(i,k)
          enddo
        enddo
!
!------------------------------------------------------------------------------------
!     T R A N S F I N I T E    I N T E R P O L A T I O N
!------------------------------------------------------------------------------------
      case('TIprism')
        call prism_TI(No,Eta, X,Dxdeta)
!
      case default
        write(*,7002) PRISMS(No)%Type
 7002   format('prism: unknown prism type = ',a10)
        stop
      endselect
!
!
   end subroutine prism
!
!
!
!------------------------------------------------------------------------------------
!> @brief Transfinite interpolation prism
!!
!! @param[in] No     - prism number
!! @param[in] Eta    - reference coordinates of a point in the reference prism
!! @param[in] X      - physical coordinates of the point
!! @param[in] Dxdeta - derivatives of the physical coordinates
!!
!! @date Mar 2023
!------------------------------------------------------------------------------------
!
subroutine prism_TI(No,Eta, X,Dxdeta)
!
      use control
      use GMP
      use element_data
      use node_types, only: PRIS
!
      implicit none
!
      integer :: No
      real(8) :: Eta(3),X(3),Dxdeta(3,3)
!
!  ...vertex shape functions
      real(8) :: shapH(6),dshapH(3,6)
!  ...2D and 1D barycentric coordinates
      real(8) :: vshapt(3),dvshapt(2,3),vshap(2),dvshap(2)
!
!  ...derivatives of edge coordinate
      real(8) :: dtedeta(3)
!
!  ...face coordinates
      real(8) :: tf(2),dtfdeta(2,3)
!
!  ...edge kernels
      real(8) :: xe(3),dxedt(3)
!
!  ...face kernels
      real(8) :: xf(3),dxfdtf(3,2)
!
!  ...blending function
      real(8) :: blend,dblend(3),fact,dfact(2)
!
      integer :: i,j,k,ie,ifig,iv,iv1,iv2,ivar
      integer :: nc,norient,np,nr,nt
      real(8) :: te
!
      integer :: iprint
      iprint=0
!
!------------------------------------------------------------------------
!
 10   continue
      if (iprint.eq.1) then
        write(*,7001) No,Eta(1:3)
 7001   format(' prism_TI: No,Eta = ',i5,2x,3e12.5)
      endif
!
!  ...affine coordinates for the triangular faces
      vshapt(1) = 1.d0 - Eta(1) - Eta(2)
      dvshapt(1,1) = -1.d0 ; dvshapt(2,1) = -1.d0
      vshapt(2) = Eta(1)
      dvshapt(1,2) =  1.d0 ; dvshapt(2,2) =  0.d0
      vshapt(3) = Eta(2)
      dvshapt(1,3) =  0.d0 ; dvshapt(2,3) =  1.d0
!
!  ...1D shape functions in the Eta_3 direction
      vshap(1) = 1.d0 - Eta(3) ; dvshap(1) = -1.d0
      vshap(2) = Eta(3)        ; dvshap(2) =  1.d0
!
!  ...vertex shape functions (products of affine coordinates)
      k=0
      do j=1,2
        do i=1,3
          k=k+1
           shapH(    k) =  vshapt(    i)* vshap(j)
          dshapH(1:2,k) = dvshapt(1:2,i)* vshap(j)
          dshapH(  3,k) =  vshapt(    i)*dvshap(j)
        enddo
      enddo
!
!------------------------------------------------------------------------
!     L I N E A R    I N T E R P O L A T I O N
!------------------------------------------------------------------------
      X(1:3)=0.d0 ; Dxdeta(1:3,1:3)=0.d0
      do iv=1,6
        np = PRISMS(No)%VertNo(iv)
        X(1:3) = X(1:3) + POINTS(np)%Rdata(1:3)*shapH(iv)
        do ivar=1,3
          Dxdeta(ivar,1:3) = Dxdeta(ivar,1:3) +  &
                             POINTS(np)%Rdata(ivar)*dshapH(1:3,iv)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'prism_TI: VERTEX INTERPOLANT = '
        do ivar=1,3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
 7011     format(e12.5,3x,3e12.5)
        enddo
      endif
!
!------------------------------------------------------------------------
!     H O R I Z O N T A L    E D G E    B U B B L E S
!------------------------------------------------------------------------
!  ...loop over horizontal edges
      ie=0
      do j=1,2
        do i=1,3
          ie=ie+1
!
!  .......get curve
          nc=PRISMS(No)%EdgeNo(ie) ; norient=0
          if (nc.lt.0) then
            nc = -nc; norient=1
          endif
          if (CURVES(nc)%Type.eq.'Seglin') cycle
!
!  .......get the edge vertices specifying the local edge orientation
          iv1=TRIAN_EDGE_TO_VERT(1,i) ; iv2=TRIAN_EDGE_TO_VERT(2,i)
!
!  .......project Eta onto the edge
          call proj_t2e(iv1,iv2,vshapt,dvshapt, te,dtedeta(1:2))
          dtedeta(3)=0
!
          if ((abs(te).lt.GEOM_TOL).or.(abs(1.d0-te).lt.GEOM_TOL)) cycle
!
          if (iprint.eq.1) then
            write(*,7012) ie,nc,CURVES(nc)%Type
 7012       format('prism_TI: ie,nc,Type = ',i2,i8,2x,a5)
          endif
!
!  .......evaluate edge kernel function
          call curveK(nc,te,norient, xe,dxedt)
          if (iprint.eq.1) then
            write(*,*) 'xe = ',xe
            write(*,*) 'dxedt = ',dxedt
          endif
!
!  .......blending function
          blend       =  vshapt(    iv1)* vshapt(    iv2)* vshap(j)
          dblend(1:2) = dvshapt(1:2,iv1)* vshapt(    iv2)* vshap(j) +  &
                         vshapt(    iv1)*dvshapt(1:2,iv2)* vshap(j)
          dblend(3)   =  vshapt(    iv1)* vshapt(    iv2)*dvshap(j)
!
!  .......add edge contribution
          X(1:3) = X(1:3) + xe(1:3)*blend
          do ivar=1,3
            Dxdeta(1:3,ivar) = Dxdeta(1:3,ivar)                  &
                             + dxedt(1:3)*dtedeta(ivar)*blend    &
                             + xe(1:3)*dblend(ivar)
          enddo
          if (iprint.eq.1) then
            write(*,*) 'prism_TI: AFTER EDGE ie = ',ie
            do ivar=1,3
              write(*,7011) X(ivar),Dxdeta(ivar,1:3)
            enddo
            call pause
          endif
        enddo
!
!  ...loop over horizontal edges
      enddo
!
!------------------------------------------------------------------------
!     V E R T I C A L    E D G E    B U B B L E S
!------------------------------------------------------------------------
!  ...loop over vertical edges
      do i=1,3
        ie=ie+1
!
!  .....get curve
        nc=PRISMS(No)%EdgeNo(ie) ; norient=0
        if (nc.lt.0) then
          nc = -nc; norient=1
        endif
        if (CURVES(nc)%Type.eq.'Seglin') cycle
!
! ......simple projection
        te=Eta(3)
        dtedeta(1:2)=0.d0 ; dtedeta(3)=1.d0
!
!  .....if at endpoint cycle
        if ((abs(te).lt.GEOM_TOL).or.(abs(1.d0-te).lt.GEOM_TOL)) cycle
!
        if (iprint.eq.1) then
          write(*,7012) ie,nc,CURVES(nc)%Type
        endif
!  .....evaluate edge bubble function
        call curveB(nc,te,norient, xe,dxedt)
        if (iprint.eq.1) then
          write(*,*) 'xe = ',xe
          write(*,*) 'dxedt = ',dxedt
        endif
!  .....blending function
        blend=vshapt(i) ; dblend(1:2)=dvshapt(1:2,i) ; dblend(3)=0.d0
!  .....add edge contribution
        X(1:3) = X(1:3) + xe(1:3)*blend
        do ivar=1,3
          Dxdeta(1:3,ivar) = Dxdeta(1:3,ivar)                  &
                           + dxedt(1:3)*dtedeta(ivar)*blend    &
                           + xe(1:3)*dblend(ivar)
        enddo
!  .....printing
        if (iprint.eq.1) then
          write(*,*) 'prism_TI: AFTER EDGE ie = ',ie
          do ivar=1,3
            write(*,7011) X(ivar),Dxdeta(ivar,1:3)
          enddo
          call pause
        endif
!
!  ...loop over vertical edges
      enddo
!
!------------------------------------------------------------------------
!     H O R I Z O N T A L    F A C E    B U B B L E S
!------------------------------------------------------------------------
      ifig=0
!  ...loop over horizontal faces
      do i=1,2
!
        ifig=ifig+1
        call decode(PRISMS(No)%FigNo(ifig), nt,norient)
!  .....skip if triangle type does not apply
        if ((TRIANGLES(nt)%Type.eq.'TransTri') .or.              &
            (TRIANGLES(nt)%Type.eq.'PlaneTri')      ) cycle
!
!  .....printing
        if (iprint.eq.1) then
          write(*,7013) ifig,nt,TRIANGLES(nt)%Type
 7013     format(' prism_TI: ifig,nt,Type = ',i2,i8,2x,a10)
        endif
!  .....project Eta onto the face
        tf(1:2)=Eta(1:2)
        dtfdeta(1:2,1)= (/1.d0,0.d0/)
        dtfdeta(1:2,2)= (/0.d0,1.d0/)
        dtfdeta(1:2,3)= 0.d0
!  .....if point is on the edge, then the bubble contribution is zero
        if ((abs(tf(2))           .lt.GEOM_TOL) .or.             &
            (abs(tf(1))           .lt.GEOM_TOL) .or.             &
            (abs(1.d0-tf(1)-tf(2)).lt.GEOM_TOL)      ) cycle
!  .....compute the face bubble
        call trianB(nt,tf,norient, xf,dxfdtf)
!  .....blending function
        blend=vshap(i) ; dblend(1:2)=0.d0 ; dblend(3)=dvshap(i)
!  .....add face contribution
        X(1:3) = X(1:3) + xf(1:3)*blend
        do ivar=1,3
          Dxdeta(1:3,ivar) = Dxdeta(1:3,ivar)                        &
                           + (dxfdtf(1:3,1)*dtfdeta(1,ivar)          &
                             +dxfdtf(1:3,2)*dtfdeta(2,ivar))*blend   &
                           + xf(1:3)*dblend(ivar)
        enddo
!  .....printing
        if (iprint.eq.1) then
          write(*,*)'prism_TI: AFTER FACE ifig = ',ifig
          do ivar=1,3
            write(*,7011) X(ivar),Dxdeta(ivar,1:3)
          enddo
          call pause
        endif
!
!  ...loop over horizontal faces
      enddo
!
!------------------------------------------------------------------------
!     V E R T I C A L    F A C E    B U B B L E S
!------------------------------------------------------------------------
!  ...loop over vertical faces
      do i=1,3
!
        ifig=ifig+1
        call decode(PRISMS(No)%FigNo(ifig), nr,norient)
!  .....skip if rectangle type does not apply
        if ((RECTANGLES(nr)%Type.eq.'TraQua') .or.              &
            (RECTANGLES(nr)%Type.eq.'BilQua')      ) cycle
!
!  .....printing
        if (iprint.eq.1) then
          write(*,7014) ifig,nr,RECTANGLES(nr)%Type
 7014     format(' prism_TI: ifig,nr,Type = ',i2,i8,2x,a10)
        endif
!  .....get the edge vertices specifying the local horizontal edge orientation
        iv1=TRIAN_EDGE_TO_VERT(1,i) ; iv2=TRIAN_EDGE_TO_VERT(2,i)
!  .....project Eta onto face
        call proj_t2e(iv1,iv2,vshapt,dvshapt, tf(1),dtfdeta(1,1:2))
        dtfdeta(1,3) = 0.d0
        tf(2)=Eta(3) ; dtfdeta(2,1:3) = (/0.d0,0.d0,1.d0/)
!  .....if point is on the edge, then the bubble contribution is zero
        if ((abs(tf(2)     ).lt.GEOM_TOL) .or.             &
            (abs(tf(1)     ).lt.GEOM_TOL) .or.             &
            (abs(1.d0-tf(2)).lt.GEOM_TOL) .or.             &
            (abs(1.d0-tf(1)).lt.GEOM_TOL)      ) cycle
!  .....compute the face bubble
        call rectaB(nr,tf,norient, xf,dxfdtf)

!!!        write(*,*)'ifig,xf = ',ifig,xf

!  .....compute semikernel
        fact=(1.d0-tf(1))*tf(1) ; dfact(1)=1.d0-2.d0*tf(1) ; dfact(2) = 0.d0
        do ivar=1,2
          dxfdtf(1:3,ivar) = (dxfdtf(1:3,ivar)*fact - xf(1:3)*dfact(ivar))/fact**2
        enddo
        xf(1:3) = xf(1:3)/fact
!  .....blending function
         blend      =  vshapt(    iv1)* vshapt(    iv2)
        dblend(1:2) = dvshapt(1:2,iv1)* vshapt(    iv2) +     &
                       vshapt(    iv1)*dvshapt(1:2,iv2)
        dblend(3)   = 0.d0
!  .....add face contribution

!!!        write(*,*)'X before = ',X

        X(1:3) = X(1:3) + xf(1:3)*blend

!!!        write(*,*)'X after = ',X

        do ivar=1,3
          Dxdeta(1:3,ivar) = Dxdeta(1:3,ivar)                        &
                           + (dxfdtf(1:3,1)*dtfdeta(1,ivar)          &
                             +dxfdtf(1:3,2)*dtfdeta(2,ivar))*blend   &
                           + xf(1:3)*dblend(ivar)
        enddo
!  .....printing
        if (iprint.eq.1) then
          write(*,*)'prism_TI: AFTER FACE ifig = ',ifig
          do ivar=1,3
            write(*,7011) X(ivar),Dxdeta(ivar,1:3)
          enddo
          call pause
        endif
!
!  ...loop over vertical faces
      enddo
!
!
end subroutine prism_TI
