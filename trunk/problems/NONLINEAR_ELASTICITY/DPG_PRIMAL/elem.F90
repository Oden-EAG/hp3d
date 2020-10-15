!--------------------------------------------------------------------------
!> Purpose : return stiffness matrix and load vector for element
!!
!! @param[in]  Mdle   - an element (middle node) number
!! @param[out] Nrdof  - number of dof for a single component
!! @param[out] Itest  - index for assembly
!! @param[out] Itrial - index for assembly
!--------------------------------------------------------------------------
!
subroutine elem(Mdle, Itest,Itrial)
  use parameters, only : ZERO
  use physics   , only : NR_PHYSA
  use assembly  , only : ALOC,BLOC,NR_RHS
  use data_structure3D
 !--------------------------------------------------------------------------
  implicit none
  integer,                    intent(in)  :: Mdle
  integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
  ! integer :: clock1,clock2,clockrate,clockmax
!--------------------------------------------------------------------------
  ! call system_clock(clock1,clockrate,clockmax)
  Itest (1:NR_PHYSA) = 0
  Itrial(1:NR_PHYSA) = 0

  ALOC(1,1)%array = ZERO; ALOC(1,2)%array = ZERO; BLOC(1)%array = ZERO
  ALOC(2,1)%array = ZERO; ALOC(2,2)%array = ZERO; BLOC(2)%array = ZERO

  select case(NODES(Mdle)%case)
  !  we wish to support both the elasticity field variables and the elasticity flux variables at each point
  !  (all physical attributes are supported when NODES(Mdle)%case == 2**NR_PHYSA-1)
  case(3)

    Itest(1:2)=1; Itrial(1:2)=1
    ! write(*,*) ''
    ! write(*,*) 'before elem...'
    call elem_DPG_PRIMAL(Mdle,  &
                BLOC(1)%nrow,BLOC(2)%nrow,  &
                ALOC(1,1)%array,ALOC(1,2)%array,BLOC(1)%array,  &
                ALOC(2,1)%array,ALOC(2,2)%array,BLOC(2)%array)
    ! write (*,*) 'after elem...'
  case default
    write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',  &
               Mdle,NODES(Mdle)%case
    call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
  end select
! !
!    call system_clock(clock2,clockrate,clockmax)
!    write(*,*),'Elapsed time in elem() subroutine: ', &
!    real(clock2-clock1, kind=8)/real( clockrate, kind=8)
end subroutine elem

!------------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for primal DPG elasticity problem
!! @param[in]  Mdle      - middle node number
!! @param[in]  Nrow_B1   - number of rows of 1-component of load vector
!! @param[in]  Nrow_B2   - number of rows of 2-component of load vector
!!
!! @param[out] Aloc11,Aloc12,Aloc11,Aloc12    - elem stiffness matrix
!! @param[out] Bloc1,Bloc2                    - elem load vector(s)
!------------------------------------------------------------------------------------------
!
!   3 equations.
!
!                    v \in H(grad)^3
!
!    \int_\Omega [  C grad(u) : grad(v)  ]
!  -             <   v , \hat \sigma_n   >
!  = \int_\Omega [     f \cdot v         ]
!
!------------------------------------------------------------------------------------------
!
subroutine elem_DPG_PRIMAL(Mdle,Nrow_Bloc1,Nrow_Bloc2,  &
                           Aloc11,Aloc12,Bloc1,Aloc21,Aloc22,Bloc2)
      use control, only: INTEGRATION
      ! use primal_module
      use parameters
      use parametersDPG
      use data_structure3D
      use element_data
      ! use isotropic_elast_material
      use assembly, only: NR_RHS
      use common_prob_data, only: SYMMETRY_TOL, TEST_NORM
      use hyperelasticity
      use nl_solver_module, only: LINESEARCH_FACTOR
      ! use m_assembly, only:  mdle_list,norderList,nedge_orientList,nface_orientList,xnodList

!------------------------------------------------------------------------------------------
      implicit none
      integer,                                  intent(in)  :: Mdle
      integer,                                  intent(in)  :: Nrow_Bloc1
      integer,                                  intent(in)  :: Nrow_Bloc2
      real*8, dimension(Nrow_Bloc1,Nrow_Bloc1), intent(out) :: Aloc11
      real*8, dimension(Nrow_Bloc1,Nrow_Bloc2), intent(out) :: Aloc12
      real*8, dimension(Nrow_Bloc2,Nrow_Bloc1), intent(out) :: Aloc21
      real*8, dimension(Nrow_Bloc2,Nrow_Bloc2), intent(out) :: Aloc22
      real*8, dimension(Nrow_Bloc1,NR_RHS),     intent(out) :: Bloc1
      real*8, dimension(Nrow_Bloc2,NR_RHS),     intent(out) :: Bloc2
!------------------------------------------------------------------------------------------
!
!  ...element and face type
      character(len=4) :: etype,ftype
!
!  ...element and face order, enriched order
      integer, dimension(19) :: norder
      integer, dimension(5)  :: nordf
      integer                :: nordP
!
!  ...node edge and face orientations
      integer, dimension(12) :: nedge_orient
      integer, dimension(6)  :: nface_orient
!
!  ...SHAPE FUNCTIONS
!     H1
      real*8, dimension(  MAXbrickH)  :: shapH
      real*8, dimension(3,MAXbrickH)  :: gradH
      integer                         :: nrdofH
!     discontinuous H1
      real*8, dimension(  MAXbrickHH) :: shapHH
      real*8, dimension(3,MAXbrickHH) :: gradHH
      integer                         :: nrdofHH
!     H(div)
      real*8, dimension(3,MAXbrickV)  :: shapV
      real*8, dimension(  MAXbrickV)  :: dshapV ! never used
      integer                         :: nrdofV
!  ...flux
      real*8, dimension(  MAXbrickV)  :: shapV_n
!
!
!  ...approximate solution
      real*8, dimension(MAXEQNH,MAXbrickH) :: dofH
      real*8, dimension(MAXEQNE,MAXbrickE) :: dofE
      real*8, dimension(MAXEQNV,MAXbrickV) :: dofV
      real*8, dimension(MAXEQNQ,MAXbrickQ) :: dofQ
      real*8, dimension(  MAXEQNH  )       :: solH
      real*8, dimension(  MAXEQNH,3)       :: dsolH
      real*8, dimension(3,MAXEQNE  )       :: solE
      real*8, dimension(3,MAXEQNE  )       :: curlE
      real*8, dimension(3,MAXEQNV  )       :: solV
      real*8, dimension(  MAXEQNV  )       :: divV
      real*8, dimension(  MAXEQNQ  )       :: solQ
      real*8, dimension(  NRHVAR  )       :: solH_up
      real*8, dimension(  NRHVAR,3)       :: dsolH_up
      ! real*8, dimension(3,NREVAR  )       :: solE_up
      ! real*8, dimension(3,NREVAR  )       :: curlE_up
      real*8, dimension(3,NRVVAR  )       :: solV_up
      ! real*8, dimension(  NRVVAR  )       :: divV_up
      ! real*8, dimension(  NRQVAR  )       :: solQ_up
!     flux
      real*8, dimension(3)                  :: sigma_n
! 
!  ...load vector for the enriched space
      real*8, dimension(3*MAXbrickHH,MAXNRHS) :: EnrLoad
!     stiffnes matrices (and copies) for the enriched test space
      real*8, dimension(3*MAXbrickHH,3*MAXbrickH) :: EnrField,EnrFieldc
      real*8, dimension(3*MAXbrickHH,3*(MAXbrickV-MAXmdlbV)) :: EnrTrace,EnrTracec
!     Gram matrix for the local Riesz matrix in LAPACK format
      real*8, dimension(3*MAXbrickHH*(3*MAXbrickHH+1)/2) :: Gram
!
!  ...geometry
      real*8, dimension(3,MAXbrickH) :: xnod
      real*8, dimension(3)           :: xi,x,rn
      real*8, dimension(3,3)         :: dxdxi,dxidx
      real*8, dimension(2)           :: t
      real*8, dimension(3,2)         :: dxidt,dxdt,rt
!
!  ...stiffness tensors in master coordinates and physical coordinates
      ! real*8, dimension(3,3,3,3) :: C,Symm,CC
      ! real*8, dimension(3,3)     :: kwave
!
! ....hyperelasticity 
      real*8 :: W,dWdF(3,3),d2WdF(3,3,3,3), Ftensor(3,3)
      integer:: imat
!  ...source term (don't need Neumann term)
      real*8, dimension(3,MAXNRHS) :: fval
!
!  ...3D quadrature data
      real*8, dimension(3,MAXNINT3ADD) :: xiloc
      real*8, dimension(MAXNINT3ADD)   :: wxi
!
!  ...2D quadrature data for boundary terms
      real*8, dimension(2,MAXNINT2ADD) :: tloc
      real*8, dimension(MAXNINT2ADD)   :: wt
!
!  ...miscellaneous
      integer :: i,j,k,l,m,n,mm,nn,k1,k2,k3,m1,m2,m3,ipt,icomp,jcomp,  &
                 nint,ifc,nsign,iprint,iflag,info,info1,info2,info3,   &
                 kH,kV,lH,lV,iel,                                      &
                 nrHbub,nrEbub,nrVbub,nrQbub
      real*8  :: weight,wa,rjac,brjac,tmp,diffmax,dmax
! 
      integer :: clock1,clock2,clockrate,clockmax
!
!  ...LAPACK stuff
      character uplo
! NOTE: nk is a "statement function"
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
      iprint=0
!
!  ...element type
      etype = NODES(Mdle)%type
!
!  ...order of approximation, orientations, geometry dof's (don't need bc flags)
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
      ! do i=1,NRELES
      !   if (mdle_list(i) == Mdle) iel = i
      ! enddo
      ! norder = norderList(1:19,iel)
      ! nedge_orient = nedge_orientList(1:12,iel)
      ! nface_orient = nface_orientList(1:6,iel)
      ! xnod = xnodList(1:3,1:MAXbrickH,iel)
!
!  ...set the enriched order of appoximation
      select case(etype)
      case('mdlb')        ; nordP = NODES(Mdle)%order + NORD_ADD*111
      case('mdln','mdld') ; nordP = NODES(Mdle)%order + NORD_ADD*1
      case('mdlp')        ; nordP = NODES(Mdle)%order + NORD_ADD*11
      end select
! 
      ! call system_clock(clock1,clockrate,clockmax)
      ! do iter=1,20
!
!  ...initialize the local element matrices and load vectors
      Aloc11 = ZERO; Aloc12 = ZERO; Bloc1 = ZERO
      Aloc21 = ZERO; Aloc22 = ZERO; Bloc2 = ZERO
!  ...initialize the enriched local element stiffness matrices and load vectors
      EnrField=ZERO; EnrTrace=ZERO; EnrLoad=ZERO
!  ...initialize the Gram matrix
      Gram=ZERO
!
      ! call getSymm(Symm)
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!
!
!  ...determine solution dof
      call solelm(Mdle, dofH,dofE,dofV,dofQ)      



!  ...set up the element quadrature
      INTEGRATION = NORD_ADD
      call set_3dint_DPG(etype,norder, nint,xiloc,wxi)
      INTEGRATION = 0

      ! do l=1,nint
      !      write(*,*) 'l1, xiloc (1,l) = ', l, xiloc (1,l)
      !      write(*,*) 'l2, xiloc (2,l) = ', l, xiloc (2,l)
      !      write(*,*) 'l3, xiloc (3,l) = ', l, xiloc (3,l)
      !      write(*,*) 'l , wxi   (l)   = ', l, wxi   (l)     
      ! enddo
! 
!
!  ...loop through integration points
      do ipt=1,nint
        xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt)
!
!  .....Compute shape functions needed for test/trial field variables
!       H1 (trial/geometry)
        call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                     nrdofH,shapH,gradH)
!       H1 (test)
        call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!  .....geometry map
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
                    x,dxdxi,dxidx,rjac,iflag)
! 
!                   Determine D1-6, keeping notation by Jason Kurtz
                    ! D1=rjac*(dxidx(1,1)**2+dxidx(1,2)**2+dxidx(1,3)**2)
                    ! D2=rjac*(dxidx(1,1)*dxidx(2,1)+ &
                    !     dxidx(1,2)*dxidx(2,2)+dxidx(1,3)*dxidx(2,3))
                    ! D3=rjac*(dxidx(1,1)*dxidx(3,1)+ &
                    !     dxidx(1,2)*dxidx(3,2)+dxidx(1,3)*dxidx(3,3))
                    ! D4=rjac*(dxidx(2,1)**2+dxidx(2,2)**2+dxidx(2,3)**2)
                    ! D5=rjac*(dxidx(2,1)*dxidx(3,1)+ &
                    !     dxidx(2,2)*dxidx(3,2)+dxidx(2,3)*dxidx(3,3))
                    ! D6=rjac*(dxidx(3,1)**2+dxidx(3,2)**2+dxidx(3,3)**2)
                    ! E =rjac;
! 
                    ! write(*,*) 'ipt = ', ipt
                    ! write(*,*) 'D1 = ', D1
                    ! write(*,*) 'D2 = ', D2
                    ! write(*,*) 'D3 = ', D3
                    ! write(*,*) 'D4 = ', D4
                    ! write(*,*) 'D5 = ', D5
                    ! write(*,*) 'D6 = ', D6
                    ! write(*,*) 'E  = ', E
                    ! call pause
! 
        if (iflag.ne.0) then
          write(*,1000) Mdle,rjac
 1000     format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
          stop
        endif
!
!
!  .....compute the approximate solution on the PHYSICAL element
        call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                     x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)


        ! solH_up (1:NRHVAR  ) = solH (NRHVAR+1:2*NRHVAR  ) + LINESEARCH_FACTOR*solH (1:NRHVAR  )
        dsolH_up(1:NRHVAR,:) = dsolH(NRHVAR+1:2*NRHVAR,:) + LINESEARCH_FACTOR*dsolH(1:NRHVAR,:)

        Ftensor = DEL + dsolH_up

        call find_material(Mdle,imat)

        if (MATERIALS(Imat)%FLAG_INCOM.and.MATERIALS(Imat)%CONSTIT.eq.0) then
          write(*,*) 'elem: primal formulation does not support an incompressible linear elastic material!'
          stop
        endif

        call eval_strain_energy_W_F(imat,x,Ftensor,W,dWdF,d2WdF)
! 
! 
!  .....integration weight
        weight = wa*rjac
!
!  .....Change coordinates so the shape functions are on the physical element
!       H1 (trial)
        do k=1,nrdofH
          gradH(1:3,k) = gradH(1,k)*dxidx(1,1:3)  &
                       + gradH(2,k)*dxidx(2,1:3)  &
                       + gradH(3,k)*dxidx(3,1:3)
        enddo
!       H1 (test)
        do k=1,nrdofHH
          gradHH(1:3,k) = gradHH(1,k)*dxidx(1,1:3)  &
                        + gradHH(2,k)*dxidx(2,1:3)  &
                        + gradHH(3,k)*dxidx(3,1:3)
        enddo
!
! !  .....compute the stiffness tensor
!         call getC(x, C)
! !
! !  .....compute the stiffness tensor squared (for gram matrix)
!         if (TEST_NORM.eq.3) then
!           call getCC(x, CC)
!         else if (TEST_NORM.eq.5) then
! !       H1 (projection)
!           call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!         endif
!
!  .....get the source term
        call getf(Mdle,x, fval)
!
!  .....OUTER loop through enriched dofs
        do k1=1,nrdofHH
          do icomp=1,3
            m1 = (k1-1)*3+icomp
! !
! !           E N R I C H E D   L O A D   V E C T O R
! !
! !   f \cdot v
! !
            EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS) +                    &
                                   ( fval(icomp,1:NR_RHS)*shapHH(k1)         &
                                    -dWdF(icomp,1)*gradHH(1,k1)              &
                                    -dWdF(icomp,2)*gradHH(2,k1)              &
                                    -dWdF(icomp,3)*gradHH(3,k1)      )*weight
!           ! if (m1.eq.49) then
!           !   write(*,*) 'ipt = ', ipt
!           !     ! write(*,*) 'weight = ', weight
!           !     ! call pause
!           !   write(*,*) 'm1, icomp, fval(icomp,1) = ',m1, icomp, fval(icomp,1)
!           !   write(*,*) 'k1,  shapHH(k1) = ', k1,  shapHH(k1)
!           !   write(*,*) 'weight = ', weight
!           ! endif
!
!           G R A M   M A T R I X
!
            select case(TEST_NORM)
! GALERKIN
!   (C:grad(v_2),grad(v)) + (v_2,v)
            case(1)
              do m2=m1,3*nrdofHH
                jcomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-1)/3)+1
                  tmp = ZERO
                  do m=1,3; do n=1,3
                  tmp = tmp  &
                      + d2WdF(jcomp,n,icomp,m)*gradHH(n,k2)*gradHH(m,k1)
                  enddo; enddo
                  Gram(k) = Gram(k)  &
                          + ( tmp + shapHH(k1)*shapHH(k2) )*weight
                endif
              enddo
! 
! MATHEMATICIANS
!   (grad(v_2),grad(v)) + (v_2,v)
            case(2)
              do m2=m1,3*nrdofHH
                jcomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( shapHH(k1)*shapHH(k2)  &
                            + gradHH(1,k1)*gradHH(1,k2)  &
                            + gradHH(2,k1)*gradHH(2,k2)  &
                            + gradHH(3,k1)*gradHH(3,k2) )*weight
                endif
              enddo
! PSEUDO-STRAIN
!   (C:eps(v_2),C:eps(v)) + (v_2,v)
              case(3)
                ! do m2=m1,3*nrdofHH
                !   jcomp = mod(m2-1,3)+1
                !   if (icomp.eq.jcomp) then
                !     k = nk(m1,m2)
                !     k2 = int((m2-1)/3)+1
                !     do m=1,3; do n=1,3
                !     tmp = tmp  &
                !         + CC(jcomp,n,icomp,m)*gradHH(n,k2)*gradHH(m,k1)
                !     enddo; enddo
                !     Gram(k) = Gram(k)  &
                !             + ( tmp + shapHH(k1)*shapHH(k2) )*weight
                !   endif
                ! enddo
! PSEUDO-STRESS
!   (eps(v_2),eps(v)) + (v_2,v)
            case(4)
                ! do m2=m1,3*nrdofHH
                !   jcomp = mod(m2-1,3)+1
                !   if (icomp.eq.jcomp) then
                !     k = nk(m1,m2)
                !     k2 = int((m2-1)/3)+1
                !     do m=1,3; do n=1,3
                !     tmp = tmp  &
                !         + Symm(jcomp,n,icomp,m)*gradHH(n,k2)*gradHH(m,k1)
                !     enddo; enddo
                !     Gram(k) = Gram(k)  &
                !             + ( tmp + shapHH(k1)*shapHH(k2) )*weight
                !   endif
                ! enddo
! STRAIN
!   (P_N(C:eps(v_2)),C:eps(v)) + (v_2,v)
            case(5)
! STRESS
!   (P_M(eps(v_2)),eps(v)) + (v_2,v)
            case(6)
! 
            end select
!
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!
!   C grad(u) : grad(v)
!
!  .........INNER loop through trial dofs
            do k3=1,nrdofH
              do jcomp=1,3
                m3 = (k3-1)*3+jcomp
                tmp = ZERO
                do m=1,3; do n=1,3
                  tmp = tmp  &
                      + d2WdF(jcomp,n,icomp,m)*gradH(n,k3)*gradHH(m,k1)
                enddo; enddo
                EnrField(m1,m3) = EnrField(m1,m3) + tmp*weight
!  .........INNER loop
              enddo
            enddo
!  .....OUTER loop
          enddo
        enddo
! !
!  ...end of loop through integration points
      enddo
! 
! enddo
! 
   !   call system_clock(clock2,clockrate,clockmax)
! 
   ! write(*,*) 'Elapsed time in Gram matrix integration: ', &
   ! real(clock2-clock1, kind=8)/real( clockrate, kind=8)/20.d0
!
!-----------------------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L                                           |
!-----------------------------------------------------------------------------------
!
    !write(*,*) 'elem: before boundary integral...'
!
!  ...get number of Hdiv bubbles
    call find_ndof(mdle,nrHbub,nrEbub,nrVbub,nrQbub)
!
!  ...loop through element faces
      do ifc=1,nface(etype)
!
!  .....sign factor to determine the OUTWARD normal unit vector
        nsign = nsign_param(etype,ifc)
!
!  .....face type
        ftype = face_type(etype,ifc)
!
!  .....face order of approximation
        call face_order(etype,ifc,norder, nordf)
!
!  .....set up the face quadrature
        INTEGRATION = NORD_ADD
        call set_2dint_DPG(ftype,nordf, nint,tloc,wt)
        INTEGRATION = 0
!
!  .....loop through face integration points
        do ipt=1,nint
!
!  .......face coordinates
          t(1:2) = tloc(1:2,ipt); wa = wt(ipt)
!
!  .......master element coordinates using face parameterization
          call face_param(etype,ifc,t, xi,dxidt)
!
!  .......Compute shape functions needed for test/trial field variables and geometry
!         H1 (geometry)
          call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                       nrdofH,shapH,gradH)
!         H(div) (trial)
          call shape3DV(etype,xi,norder,nface_orient,  &
                       nrdofV,shapV,dshapV)
!         H1 (test)
          call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!  .......geometry map
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,  &
                       x,dxdxi,dxidx,rjac,dxdt,rn,brjac)
          weight = wa*brjac

!
!  .......compute the approximate solution on the PHYSICAL element
          call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                       x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)

          solV_up(:,1:NRVVAR) = solV(:,NRVVAR+1:2*NRVVAR) + LINESEARCH_FACTOR*solV(:,1:NRHVAR)

!
!  .......compute approximate flux vector
          sigma_n(1:3) = solV_up(1,1:3)*rn(1)  &
                       + solV_up(2,1:3)*rn(2)  &
                       + solV_up(3,1:3)*rn(3)


!
!  .......Change coordinates so the shape functions are on the physical element.
!         Also construct trial fluxes
!         H(div) (trial)
          do k=1,nrdofV-nrVbub
            shapV(1:3,k) = dxdxi(1:3,1)*shapV(1,k)  &
                         + dxdxi(1:3,2)*shapV(2,k)  &
                         + dxdxi(1:3,3)*shapV(3,k)
            shapV_n(k) = shapV(1,k)*rn(1)  &
                       + shapV(2,k)*rn(2)  &
                       + shapV(3,k)*rn(3)
          enddo
          shapV_n(1:nrdofV-nrVbub) = shapV_n(1:nrdofV-nrVbub)/rjac
!
!  .......OUTER loop through enriched H1 test functions
          do k1=1,nrdofHH
            do icomp=1,3
              m1 = (k1-1)*3+icomp
!
!
!           E N R I C H E D   L O A D  -  T R A C E   C O N T R I B U T I O N
! 
              EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS) &
                                   + sigma_n(icomp)*shapHH(k1)*weight
! 
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!
!   - <  v , \hat \sigma_n  >
!
!  ...........INNER loop through enriched dofs
              do k2=1,nrdofV-nrVbub
                !  contribute only when components match (populate only the diagonal)
                m2 = (k2-1)*3+icomp
                EnrTrace(m1,m2) = EnrTrace(m1,m2)  &
                                - shapV_n(k2)*shapHH(k1)*weight
              enddo
!
!  .......OUTER loop
            enddo
          enddo
!
!  .....end of loop over integration points
        enddo
!
!  ...end of loop over faces
      enddo

      !write(*,*) 'elem: after boundary integral...'

!
      if (iprint.eq.1) then
! 
        write(*,*) 'EnrLoad = '
        do i=1,3*nrdofHH+1
          write(*,6000) i,EnrLoad(i,1)
 6000     format('i = ',i4,'  ',e12.5)
        enddo
! 
        call pause
! 
        write(*,*) 'Gram = '
        do i=1,25
          write(*,6001) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 6001     format('i = ',i3,'  ',25e12.5)
        enddo
! 
        call pause
! 
        write(*,*) 'EnrField = '
        do i=1,3*nrdofHH+1
          write(*,6002) i,EnrField(i,1:3*nrdofH)
 6002     format('i = ',i4,'  ',15(/,10e12.5))
        enddo
! 
        call pause
! 
        write(*,*) 'EnrTrace = '
        do i=1,3*nrdofHH+1
          write(*,6002) i,EnrTrace(i,1:3*nrdofV)
        enddo
! 
      endif
! 
! 
!
      !write(*,*) 'elem: before local Riesz problem...'
!-----------------------------------------------------------------------------------
!     D P G    L O C A L    A S S E M B L Y                                        |
!-----------------------------------------------------------------------------------
!
!  ...factor the Gram matrix
      uplo = 'U'
      call DPPTRF(uplo, 3*nrdofHH, Gram, info)
      if (info.ne.0) then
        write(*,*) 'elem_DPG_PRIMAL: info = ',info ; stop
      endif
!
!  ...save copies of enriched stiffness matrices
      EnrFieldc=EnrField ; EnrTracec=EnrTrace
!
!  ...G^-1 * Load
      call DPPTRS(uplo,3*nrdofHH,NR_RHS,Gram,EnrLoad,3*MAXbrickHH,info1)
      if (info1.ne.0) then
        write(*,*) 'elem_DPG_PRIMAL: info1 = ',info1 ; stop
      endif
!
!  ...G^-1 * EnrField
      call DPPTRS(uplo,3*nrdofHH,3*nrdofH,Gram,EnrFieldc,3*MAXbrickHH,info2)
      if (info2.ne.0) then
        write(*,*) 'elem_DPG_PRIMAL: info2 = ',info2 ; stop
      endif
!
!  ...G^-1 * EnrTrace
      call DPPTRS(uplo,3*nrdofHH,3*(nrdofV-nrVbub),Gram,EnrTracec,3*MAXbrickHH,info3)
      if (info3.ne.0) then
        write(*,*) 'elem_DPG_PRIMAL: info3 = ',info3 ; stop
      endif

      !write(*,*) 'elem: after local Riesz problem...'
!
!
!  ...ultimate DPG load vector and stiffness matrices
      do kH=1,3*nrdofH
!
!       EnrField^T * (G^-1 * EnrLoad)
        do m=1,NR_RHS
          do k=1,3*nrdofHH
            Bloc1(kH,m) = Bloc1(kH,m)  &
                        + EnrField(k,kH)*EnrLoad(k,m)
          enddo
        enddo
!
!       EnrField^T * (G^-1 * EnrField)
        do lH=1,3*nrdofH
          do k=1,3*nrdofHH
            Aloc11(kH,lH) = Aloc11(kH,lH)  &
                          + EnrField(k,kH)*EnrFieldc(k,lH)
          enddo
        enddo
!
!       EnrField^T * (G^-1 * EnrTrace)
        do lV=1,3*(nrdofV-nrVbub)
          do k=1,3*nrdofHH
            Aloc12(kH,lV) = Aloc12(kH,lV)  &
                          + EnrField(k,kH)*EnrTracec(k,lV)
          enddo
        enddo
      enddo
!
      !write(*,*) 'elem: before Hdiv part...'

      do kV=1,3*(nrdofV-nrVbub)
!
!       EnrTrace^T * (G^-1 * EnrLoad)
        do m=1,NR_RHS
          do k=1,3*nrdofHH
            Bloc2(kV,m) = Bloc2(kV,m)  &
                        + EnrTrace(k,kV)*EnrLoad(k,m)
          enddo
        enddo
!
!       EnrTrace * (G^-1 * EnrField)^T
        do lH=1,3*nrdofH
          do k=1,3*nrdofHH
            Aloc21(kV,lH) = Aloc21(kV,lH)  &
                          + EnrTrace(k,kV)*EnrFieldc(k,lH)
          enddo
        enddo
!
!       EnrTrace * (G^-1 * EnrTrace)^T
        do lV=1,3*(nrdofV-nrVbub)
          do k=1,3*nrdofHH
            Aloc22(kV,lV) = Aloc22(kV,lV)  &
                          + EnrTrace(k,kV)*EnrTracec(k,lV)
          enddo
        enddo
      enddo

      !write(*,*) 'elem: before test part...'
!
!-----------------------------------------------------------------------------------
!     T E S T S    A N D    P R I N T    S T A T E M E N T S                       |
!-----------------------------------------------------------------------------------
!
      if (iprint.ge.2) then
  !  ...check symmetry
        if (iprint.ge.2) then
          diffmax = ZERO; dmax = ZERO
          do k1=1,3*nrdofH
            do k2=k1,3*nrdofH
              diffmax = max(diffmax,abs(Aloc11(k1,k2)-Aloc11(k2,k1)))
              dmax = max(dmax,abs(Aloc11(k1,k2)))
            enddo
          enddo
          if (diffmax/dmax.gt.SYMMETRY_TOL) then
            write(*,7021) diffmax, dmax
     7021   format('elem_DPG_PRIMAL: diffmax,dmax FOR Aloc11 = ',2e12.5)
            call pause
          endif
          diffmax = ZERO; dmax = ZERO
          do k1=1,nrdofV
            do k2=k1,nrdofV
              diffmax = max(diffmax,abs(Aloc22(k1,k2)-Aloc22(k2,k1)))
              dmax = max(dmax,abs(Aloc22(k1,k2)))
            enddo
          enddo
          if (diffmax/dmax.gt.SYMMETRY_TOL) then
            write(*,7022) diffmax, dmax
     7022   format('elem_DPG_PRIMAL: diffmax,dmax FOR Aloc22 = ',2e12.5)
            call pause
          endif
          diffmax = ZERO; dmax = ZERO
          do k1=1,nrdofH
            do k2=1,nrdofv
              diffmax = max(diffmax,abs(Aloc12(k1,k2)-Aloc21(k2,k1)))
              dmax = max(dmax,abs(Aloc12(k1,k2)))
            enddo
          enddo
          if (diffmax/dmax.gt.SYMMETRY_TOL) then
            write(*,7023) diffmax, dmax
     7023   format('elem_DPG_PRIMAL: diffmax,dmax FOR Aloc12 = ',2e12.5)
            call pause
          endif
        endif
      endif
!
!  ...print statments
      if (iprint.ge.1) then
        write(*,7010)0
        
 7010   format('elem_DPG_PRIMAL: Bloc1,Bloc2 = ')
        write(*,7011) Bloc1(1:3*NrdofH,1)
        write(*,7011) Bloc2(1:3*NrdofV,1)
 7011   format(10e12.5)
        write(*,7012)
 7012   format('elem_DPG_PRIMAL: Aloc11 = ')
        do i=1,3*NrdofH
          write(*,7013) i,Aloc11(i,1:3*NrdofH)
 7013     format('i = ',i3,10(/,10e12.5))
        enddo
        write(*,7014)
 7014   format('elem_DPG_PRIMAL: Aloc12 = ')
        do i=1,3*NrdofH
          write(*,7013) i,Aloc12(i,1:3*NrdofV)
        enddo
        write(*,7015)
 7015   format('elem_DPG_PRIMAL: Aloc22 = ')
        do i=1,3*NrdofV
          write(*,7013) i,Aloc22(i,1:3*NrdofV)
        enddo
      endif
!
!
end subroutine elem_DPG_PRIMAL
! 