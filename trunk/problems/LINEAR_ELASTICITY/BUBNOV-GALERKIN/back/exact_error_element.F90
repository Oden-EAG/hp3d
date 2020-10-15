!-----------------------------------------------------------------
!> Purpose : compute element contributions to global H1 error
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[out] Derr  - element error (squared)
!! @param[out] Dnorm - element norm of the exact solution (squared)
!-----------------------------------------------------------------
subroutine exact_error_element(Mdle, Derr,Dnorm)
  use data_structure3D
  use element_data, only: nface, nsign_param, face_param, face_type, face_order
  use parameters, only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ
  use common_prob_data
  use isotropic_elast_material
  !---------------------------------------------------------------
  implicit none
  integer, intent(in)  :: Mdle
  real*8,  intent(out) :: Derr, Dnorm
  !---------------------------------------------------------------

  ! local error
  real*8 :: derr_l2, derr_h1(3), derr_en_intr, derr_en_bndry

  ! local energy
  real*8 :: jsol

  !  element type
  character(len=4) :: etype

  ! face type
  character(len=4) :: ftype

  !  element order of approximation
  integer :: norder(19), nordf(5), nedge_orient(12), nface_orient(6)

  ! face node (middle node)
  integer :: mdlf

  ! face boundary condition
  integer :: ibc(6), fbc

  !  geometry
  real*8,dimension(3,MAXbrickH) :: xnod
  real*8,dimension(3)           :: xi,x
  real*8,dimension(2)           :: t
  real*8,dimension(3,3)         :: dxdxi,dxidx
  real*8,dimension(3,2)         :: dxdt,dxidt

  !  approximate solution dof's
  real*8,dimension(MAXEQNH,MAXbrickH) :: zdofH
  real*8,dimension(MAXEQNE,MAXbrickE) :: zdofE
  real*8,dimension(MAXEQNV,MAXbrickV) :: zdofV
  real*8,dimension(MAXEQNQ,MAXbrickQ) :: zdofQ

  !  approximate solution and its derivatives
  real*8,dimension(MAXEQNH  ) :: zsolH
  real*8,dimension(MAXEQNH,3) :: zdsolH

  !  shape functions and their derivatives
  real*8,dimension(  MAXbrickH) :: shapH
  real*8,dimension(3,MAXbrickH) :: dshapH,dshapHx

  !  3D quadrature data
  real*8,dimension(3,MAX_NINT3) :: xiloc
  real*8,dimension(MAX_NINT3)   :: wxi

  !  2D quadrature data for boundary terms
  double precision :: tloc(2,MAXquadH),wt(MAXquadH)

  !  external unit vector
  double precision :: rn(3)

  !  exact solution
  real*8,    dimension(  MAXEQNH    ) ::   zvalH
  real*8,    dimension(  MAXEQNH,3  ) ::  zdvalH
  real*8,    dimension(  MAXEQNH,3,3) :: zd2valH
  real*8,    dimension(3,MAXEQNE    ) ::   zvalE
  real*8,    dimension(3,MAXEQNE,3  ) ::  zdvalE
  real*8,    dimension(3,MAXEQNE,3,3) :: zd2valE
  real*8,    dimension(3,MAXEQNV    ) ::   zvalV
  real*8,    dimension(3,MAXEQNV,3  ) ::  zdvalV
  real*8,    dimension(3,MAXEQNV,3,3) :: zd2valV
  real*8,    dimension(  MAXEQNQ    ) ::   zvalQ
  real*8,    dimension(  MAXEQNQ,3  ) ::  zdvalQ
  real*8,    dimension(  MAXEQNQ,3,3) :: zd2valQ

  !  number of H1, Hcurl, Hdiv, L2 dof's
  integer :: nrdofH,nrdofE,nrdofV,nrdofQ
  !
  !  miscellaneous
  integer :: i,j,k,l,nint,icase,icomp,iflag,icompH,iload,nf,nsign
  real*8  :: wa,rjac,weight
  !
  integer :: iprint
  !---------------------------------------------------------------
  !
  iprint = 0
  if (iprint.eq.1) then
     write(*,7001) Mdle,NODES(Mdle)%type
7001 format('exact_error_element: Mdle,type = ',i10,2x,a5)
  endif
  !
  !  order of approximation, orientations, geometry dof's
  call find_order (Mdle, norder)
  call find_orient(Mdle, nedge_orient,nface_orient)
  call nodcor     (Mdle,xnod)
  call solelm     (Mdle, zdofH,zdofE,zdofV,zdofQ)
  !
  !  initialize
  Derr    = 0.d0; Dnorm   = 0.d0
  derr_l2 = 0.d0; derr_h1 = 0.d0; derr_en_intr = 0.d0; derr_en_bndry = 0.d0


  !  set up the element quadrature
  etype = NODES(Mdle)%type
  icase = NODES(Mdle)%case
  call set_3Dint(etype,norder, nint,xiloc,wxi)

  !  loop through integration points
  do l=1,nint
    xi(1:3) = xiloc(1:3,l); wa = wxi(l)

    !  evaluate appropriate shape functions at the point
    call shape3H(etype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,dshapH)

    !  geometry map
    x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0
    do k=1,nrdofH
      x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
      do i=1,3
        dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
      enddo
    enddo

    !  evaluate the inverse derivatives and jacobian
    call geom(dxdxi, dxidx,rjac,iflag)
    if (iflag.ne.0) then
      write(*,*) 'elem_elem_error: NEGATIVE JACOBIAN FOR Mdle = ', &
           Mdle, '        rjac = ',rjac
      stop 1
    endif

    !  compute derivatives wrt physical coordinates
    do k=1,nrdofH
      dshapHx(1:3,k) = 0.d0
      do icomp=1,3
        dshapHx(1:3,k) = dshapHx(1:3,k) + dshapH(icomp,k)*dxidx(icomp,1:3)
      enddo
    enddo

    !  total weight
    weight = wa*rjac

    !  evaluate the approximate solution
    zsolH(1:MAXEQNH) = ZERO; zdsolH(1:MAXEQNH,1:3) = ZERO
    do k=1,nrdofH
      zsolH(1:MAXEQNH) = zsolH(1:MAXEQNH) + zdofH(1:MAXEQNH,k)*shapH(k)
      do icomp=1,3
        zdsolH(1:MAXEQNH,icomp) = zdsolH(1:MAXEQNH,icomp) +  zdofH(1:MAXEQNH,k)*dshapHx(icomp,k)
      enddo
    enddo

    select case(IERROR_PROB)
    case(IERROR_L2,IERROR_H1)
      !  compute the exact solution
      icase = 1
      call exact(x,icase, &
          zvalH,zdvalH,zd2valH, &
          zvalE,zdvalE,zd2valE, &
          zvalV,zdvalV,zd2valV, &
          zvalQ,zdvalQ,zd2valQ)

      if (iprint.eq.1) then
        write(*,7003) l,x(1:3),zsolH(1:3),zvalH(1:3)
      7003    format('exact_error_element: l,x = ',i4,2x,3f8.3,2x,  &
             /,' APPROX SOLUTION = ',3e12.5,&
             /,' EXACT  SOLUTION = ',3e12.5)
      endif

      ! Loop through components
      do icompH=1,MAXEQNH
        ! L2 Error
        Dnorm   = Dnorm    + abs(zvalH(icompH)               )**2*weight
        derr_l2 = derr_l2  + abs(zvalH(icompH) - zsolH(icompH))**2*weight

        ! H1 Error
        if (IERROR_PROB.eq.IERROR_H1) then
          do icomp=1,3
             Dnorm          = Dnorm          + abs(zdvalH(icompH,icomp))**2*weight
             derr_h1(icomp) = derr_h1(icomp) + abs(zdvalH(icompH,icomp) - zdsolH(icompH,icomp))**2*weight
          enddo
        endif
      enddo

    ! Energy Error
    !     1/2 ||u_h - u||^2  =  J(u_h) - J(u)
    ! where
    !     J(v)  =  1/2 b(v,v) - l(v)
    case(IERROR_ENERGY)

      !  compute interior integral contributions to J(u_h)
      iload=1
      call energy_integrand_interior(Mdle,x,zsolH,zdsolH,iload, jsol)

      Dnorm   = Dnorm  +1.d0*weight
      derr_en_intr = derr_en_intr+jsol*weight

    case default
      write(*,*) 'exact_error_element : IERROR_PROB = ', IERROR_PROB
      return
    endselect
  enddo

  select case(IERROR_PROB)

  case(IERROR_L2,IERROR_H1)
    !  boundary integrals do not contibute to these norms
    Derr = derr_l2 + sum(derr_h1(1:3))

  case(IERROR_ENERGY)
    nint = 0
    !  find boundary conditions on element
    call find_bc(mdle, ibc)
    !  compute boundary integrals
    do nf=1,nface(etype)
      ! !  get middle node for face
      ! call elem_face(Mdle,nf, mdlf)
      ! !  get boundary condition at face
      ! fbc = NODES(mdlf)%bcond
      !  get boundary condition at face
      fbc = ibc(nf)
      !  skip if this face is not Neumann in any component
      if ((fbc.ne.2).and.(fbc.ne.8)) cycle
      !  determine the type of the face
      ftype = face_type(etype,nf)
      !  determine order for the face
      call face_order(etype,nf,norder, nordf)
      !  begin quadrature on the face
      call set_2Dint(ftype,nordf, nint,tloc,wt)
      if (iprint.eq.1) then
        call pause
        write(*,4001) nf,fbc,nint,nordf
4001    format('exact_error_element: nf,fbc,nint,nordf = ', &
               i2,2x,i2,2x,i3,2x,4i2,i3)
        write(*,*)
      endif
      nsign = Nsign_param(etype,nf)
      !  loop through face integration points
      do l=1,nint
        t(1:2)=tloc(1:2,l) ; wa=wt(l)
        !  determine the master element coordinates
        call face_param(etype,nf,t, xi,dxidt)
        !  derivatives and values of the shape functions
        call shape3H(etype,xi,norder,nedge_orient,nface_orient, &
                     nrdofH,shapH,dshapH)
        !  geometry map
        x(1:3)=0.d0 ; dxdxi(1:3,1:3)=0.d0
        do k=1,nrdofH
          x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
          do i=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
          enddo
        enddo
        if (iprint.eq.1) then
          write(*,4002) (x(i),i=1,3)
4002      format('exact_error_element: x = ',3(/,f8.3))
          write(*,4003) (dxdxi(i,1:3),i=1,3)
4003      format('exact_error_element: dxdxi = ',3(/,3f8.3))
          write(*,4004) (dxidt(i,1:2),i=1,3)
4004      format('exact_error_element: dxidt = ',3(/,2f8.3))
        endif
        !  determine normal unit vector and Jacobian
        dxdt(1:3,1:2) = 0.d0
        do i=1,2
          do j=1,3
            dxdt(1:3,i) = dxdt(1:3,i) + dxdxi(1:3,j)*dxidt(j,i)
          enddo
        enddo
        call cross_product(dxdt(1:3,1),dxdt(1:3,2), rn)
        call norm(rn, rjac)
        rn(1:3) = rn(1:3)*nsign/rjac
        if (iprint.eq.1) then
          write(*,4005) nf, nsign, rn
4005      format('exact_error_element: nf,nsign,rn = ',i2,2x,i2,2x,3f8.3)
        endif
        !  evaluate the approximate solution
        zsolH(1:MAXEQNH) = ZERO; zdsolH(1:MAXEQNH,1:3) = ZERO
        do k=1,nrdofH
          !  only the solution values are necessary
          zsolH(1:MAXEQNH) = zsolH(1:MAXEQNH) + zdofH(1:MAXEQNH,k)*shapH(k)
        enddo
        !  total weight
        weight = wa*rjac
        !  compute integrand
        iload = 1
        call energy_integrand_boundary(Mdle,fbc,x,zsolH,rn,iload, jsol)
        derr_en_bndry = derr_en_bndry+jsol*weight

        ! write(*,*)
        ! write(*,*) 'zsolH  = ', zsolH
        ! write(*,*) 'jsol   = ', jsol
        ! write(*,*) 'weight = ', weight

      !  end of loop through integration points
      enddo
    !  end loop through faces
    enddo
    if (iprint.eq.1) then
      write(*,*)
      write(*,4006) derr_en_intr, derr_en_bndry
4006  format('exact_error_element: derr_en_intr, derr_en_bndry = ',f8.6,2x,f10.8)
    endif

    Derr = derr_en_intr + derr_en_bndry

  endselect

  NODES(Mdle)%error(0,  0) = derr_l2
  NODES(Mdle)%error(1:3,0) = derr_h1(1:3)
  ! NODES(Mdle)%error(4,  0) = derr_en



end subroutine exact_error_element

!-----------------------------------------------------------------
!> Purpose : compute interior integrands in energy functional for PVW linear elasticity
!!
!! @param[in]  Mdle   - element (middle node) number
!! @param[in]  X      - physical coordinates
!! @param[in]  ZsolH  - value of the H1 solution
!! @param[in]  ZdsolH - corresponding first derivatives
!! @param[in]  iload  - load number
!! @param[out] Jsol   - energy
!-----------------------------------------------------------------
subroutine energy_integrand_interior(Mdle,X,ZsolH,ZdsolH,Iload, Jsol)
  use isotropic_elast_material
  use assembly  , only: NR_RHS
  use parameters, only: MAXEQNH
!-----------------------------------------------------------------
  implicit none
  integer,                     intent(in)  :: Mdle
  real*8,dimension(3),         intent(in)  :: X
  real*8,dimension(MAXEQNH  ), intent(in)  ::  ZsolH
  real*8,dimension(MAXEQNH,3), intent(in)  :: ZdsolH
  integer,                     intent(in)  :: Iload
  real*8,                      intent(out) :: Jsol
!-----------------------------------------------------------------
  integer :: icomp,j,k,l
! Elasticity Tensor
  real*8,dimension(3,3,3,3) :: EE
! Body Force
  real*8,dimension(3,NR_RHS) :: zfval
!-----------------------------------------------------------------

! NOTE THAT THIS ONLY COMPUTES THE INTERIOR INTEGRALS
! WE STILL HAVE TO COMPUTE THE BOUNDARY INTEGRALS!!!!

! get elasticity tensor
  call getC(X, EE)

! get l
  call getf(Mdle,X, zfval)

  Jsol=0.d0
! b(u,v) = \int_\Omega EE_ijkl du_k/dx_l dv_i/dx_j
! l(u)   = \int_\Omega F_i u_i + \int_\Gamma_{Neumann} EE_ijkl du_k/dx_l n_j v_i
! J(u)   = 1/2 b(u,u)-l(u)
  do icomp=1,3
    do j=1,3; do k=1,3; do l=1,3
      Jsol = Jsol + 0.5d0*EE(icomp,j,k,l)*ZdsolH(k,l)*ZdsolH(icomp,j)
    enddo; enddo; enddo
    Jsol = Jsol - zfval(icomp,Iload)*ZsolH(icomp)
  enddo

end subroutine energy_integrand_interior


!-----------------------------------------------------------------
!> Purpose : compute boundary integrands in energy functional for PVW linear elasticity
!!
!! @param[in]  Mdle   - element (middle node) number
!! @param[in]  Fbc    - boundary condition on given face
!! @param[in]  X      - physical coordinates
!! @param[in]  ZsolH  - value of the H1 solution
!! @param[in]  iload  - load number
!! @param[out] Jsol   - energy
!-----------------------------------------------------------------
subroutine energy_integrand_boundary(Mdle,Fbc,X,ZsolH,Rn,Iload, Jsol)
  use isotropic_elast_material
  use assembly  , only: NR_RHS
  use parameters, only: MAXEQNH
!-----------------------------------------------------------------
  implicit none
  integer,                   intent(in)  :: Mdle
  integer,                   intent(in)  :: fbc
  real*8,dimension(3),       intent(in)  :: X
  real*8,dimension(MAXEQNH), intent(in)  :: ZsolH
  real*8,dimension(3),       intent(in)  :: Rn
  integer,                   intent(in)  :: Iload
  real*8,                    intent(out) :: Jsol
!-----------------------------------------------------------------
  integer :: icomp
! Neumann data
  real*8,dimension(3,NR_RHS) :: gval
!-----------------------------------------------------------------

! NOTE THAT THIS ONLY COMPUTES THE BOUNDARY INTEGRALS

! get the Neumann data
  call getg(Mdle,Fbc,X,Rn, gval)

  Jsol=0.d0
! b(u,v) = \int_\Omega EE_ijkl du_k/dx_l dv_i/dx_j
! l(u)   = \int_\Omega F_i u_i + \int_\Gamma_{Neumann} g_i u_i
! J(u)   = 1/2 b(u,u)-l(u)

  if (Fbc.eq.2) then
    do icomp=1,3
      Jsol = Jsol - gval(icomp,Iload)*ZsolH(icomp)
    enddo
  elseif (Fbc.eq.8) then
    !  Dirichlet on 3rd component
    do icomp=1,2
      Jsol = Jsol - gval(icomp,Iload)*ZsolH(icomp)
    enddo
  endif

end subroutine energy_integrand_boundary
