!--------------------------------------------------------------------------
!> Purpose : return stiffness matrix and Residual vector for element
!!
!! @param[in]  Mdle      - an element (middle node) number
!! @param[out] Resid     - element residual (squared)
!! @param[out] Nref_flag - suggested h-refinement flag
!------------------------------------------------------------------------------------------
!       L I N E A R I Z E D    E L A S T O S T A T I C S    E Q U A T I O N S    
!                                   R E S I D U A L
! 
!     DISPLACEMENT GRADIENT   | LINEAR MOMENTUM BALANCE | PK1 CONSTITUTIVE LAW 
!
!  |   \tau \in H(div)^3      |       v \in (H1)^3      | \chi \in (L^2)^3x3
! 
!  = 
!                             + (f , v)                 + ( K(I+D),\chi)     
!  - (u,div(\tau))            + <\hat t ,v>             - ( P,\chi )
!  + < \hat u,(\tau n)>       - (P,grad(v))
!  - (D,\tau)
!------------------------------------------------------------------------------------------
!
      subroutine elem_residual(Mdle, Resid,Nref_flag)
!
      use control, only: INTEGRATION
      ! use primal_module
      use parametersDPG
      use element_data
      use data_structure3D
      ! use isotropic_elast_material
      use common_prob_data, only: TEST_NORM, IP
      use assembly  , only: NR_RHS
      use hyperelasticity
      use nl_solver_module, only: LINESEARCH_FACTOR     
!------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in)  :: Mdle
      real*8,  intent(out) :: Resid
      integer, intent(out) :: Nref_flag
!------------------------------------------------------------------------------------------
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
!     H1  (geometry and trial)
      real*8, dimension(  MAXbrickH)  :: shapH
      real*8, dimension(3,MAXbrickH)  :: gradH
      integer                         :: nrdofH
!     H(div)  (trial)
      real*8, dimension(3,MAXbrickV)  :: shapV
      real*8, dimension(  MAXbrickV)  :: divV
      real*8, dimension(  MAXbrickV)  :: shapV_n
      integer                         :: nrdofV
!     L2  (trial)
      real*8, dimension(  MAXbrickQ)  :: shapQ
      integer                         :: nrdofQ
!     H1   (test)
      real*8, dimension(  MAXbrickHH) :: shapHH
      real*8, dimension(3,MAXbrickHH) :: gradHH
      integer                         :: nrdofHH
!     H(div)  (test)
      real*8, dimension(3,MAXbrickVV) :: shapVV
      real*8, dimension(  MAXbrickVV) :: divVV
      real*8, dimension(  MAXbrickVV) :: shapVV_n
      integer                         :: nrdofVV
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
      real*8, dimension(3,MAXEQNE  )       :: dsolE
      real*8, dimension(3,MAXEQNV  )       :: solV
      real*8, dimension(  MAXEQNV  )       :: dsolV
      real*8, dimension(  MAXEQNQ  )       :: solQ
      real*8, dimension(  NRHVAR  )       :: solH_up
      ! real*8, dimension(  NRHVAR,3)       :: dsolH_up
      ! real*8, dimension(3,NREVAR  )       :: solE_up
      ! real*8, dimension(3,NREVAR  )       :: curlE_up
      real*8, dimension(3,NRVVAR  )       :: solV_up
      ! real*8, dimension(  NRVVAR  )       :: divV_up
      real*8, dimension(  NRQVAR  )       :: solQ_up
!     flux
      ! real*8, dimension(3)                  :: sigma_n
! 
!  ...Necessary constant since there is no dynamic allocation
      integer, parameter :: MAXNRHS_MOD = 1
!
!  ...residual for the enriched test space 
      real*8, dimension(3*MAXbrickVV+3*MAXbrickHH) :: EnrResid,EnrResidc
!     Gram matrix for the local Riesz matrix in LAPACK format
      real*8, dimension((3*MAXbrickVV+3*MAXbrickHH)*(3*MAXbrickVV+3*MAXbrickHH+1)/2) :: Gram
!     residual for Least Squares equations
      real*8      :: FOSLS_Resid
!
!  ...geometry
      real*8, dimension(3,MAXbrickH) :: xnod
      real*8, dimension(3)           :: xi,x,rn
      real*8, dimension(3,3)         :: dxdxi,dxidx
      real*8, dimension(2)           :: t
      real*8, dimension(3,2)         :: dxidt,dxdt,rt
!
! ....hyperelasticity 
      real*8 :: W,dWdF(3,3),d2WdF(3,3,3,3), gradu(3,3) , PK1(3,3), u(3), uhat(3), trac(3)
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
      integer :: i,j,k,l,m,n,k1,k2,k3,k4,k5,m1,m2,m3,m4,m5,n1,n2,n3,nsign,ipt,ifc,  &
                 icomp,jcomp,kcomp,lcomp,nint,iprint,iflag,info,info1,info2,info3,  &
                 kH,kV,kQ,lH,lV,lQ,kmin,kmax,enrdof,                                &
                 nrHbub,nrEbub,nrVbub,nrQbub,       gdump,bdump6,bdump7
      real*8  :: weight,wa,rjac,bjac,tmp,diff
!
!  ...LAPACK stuff
      character uplo
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
!  ...Load number for which to calculate residual from
      ! iload = 1
!
      ! select case(Mdle)
      ! case(1)
      !   iprint=0
      ! case default
      !   iprint=0
      ! end select
      iprint=0
!
!  ...element type
      etype = NODES(Mdle)%type
!
!  ...order of approximation, orientations, geometry dof's
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
!
!  ...set the enriched order of appoximation
      ! nordtmp = NORD_ADD
      ! nordtmp = 4 - IP !max(NORD_ADD,2)
      select case(etype)
      case('mdlb')        ; nordP = NODES(Mdle)%order + NORD_ADD*111
      case('mdln','mdld') ; nordP = NODES(Mdle)%order + NORD_ADD*1
      case('mdlp')        ; nordP = NODES(Mdle)%order + NORD_ADD*11
      end select
!
!  ...initialize the enriched local element residual vector
      EnrResid = ZERO; FOSLS_Resid = ZERO
!  ...initialize the Gram matrix
      Gram = ZERO

      ! call getSymm(Symm)
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!
!
!  ...determine solution dof
      call solelm(Mdle, dofH,dofE,dofV,dofQ)      
! 
!  ...set up the element quadrature
      INTEGRATION = NORD_ADD
      call set_3dint_DPG(etype,norder, nint,xiloc,wxi)
      INTEGRATION = 0
!
!  ...loop through integration points
      do ipt=1,nint
        xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt)
!
!  .....Compute shape functions needed for test/trial field variables and geometry
!       H1 (geometry)
        call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                     nrdofH,shapH,gradH)
!       L2 (field trial)
        call shape3DQ(etype,xi,norder, nrdofQ,shapQ)
!       H1 (test)
        call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!       H(div) (test)
        call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!
!  .....geometry map
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
                    x,dxdxi,dxidx,rjac,iflag)
        if (iflag.ne.0) then
          write(*,1000) Mdle,rjac
 1000     format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
          stop
        endif
!
!  .....Change coordinates so the shape functions are on the physical element
!       L2 (trial)
        shapQ(1:nrdofQ) = shapQ(1:nrdofQ)/rjac
!       H1 (test)
        do k=1,nrdofHH
          gradHH(1:3,k) = gradHH(1,k)*dxidx(1,1:3)  &
                        + gradHH(2,k)*dxidx(2,1:3)  &
                        + gradHH(3,k)*dxidx(3,1:3)
        enddo
!       H(div) (test)
        do k=1,nrdofVV
          shapVV(1:3,k) = dxdxi(1:3,1)*shapVV(1,k)  &
                        + dxdxi(1:3,2)*shapVV(2,k)  &
                        + dxdxi(1:3,3)*shapVV(3,k)
        enddo
        shapVV(1:3,1:nrdofVV) = shapVV(1:3,1:nrdofVV)/rjac
        divVV(1:nrdofVV) = divVV(1:nrdofVV)/rjac
!
!
!  .....compute the approximate solution on the PHYSICAL element
        call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                     x,dxdxi,solH,dsolH,solE,dsolE,solV,dsolV,solQ)


        solQ_up = solQ(NRQVAR+1:2*NRQVAR) + LINESEARCH_FACTOR*solQ(1:NRQVAR)
        u = solQ_up(1:3)
        PK1 = reshape(solQ_up(4 :12),(/3,3/))
        gradu = reshape(solQ_up(13:21),(/3,3/))


        call find_material(Mdle,imat)

        ! if (MATERIALS(Imat)%FLAG_INCOM.and.MATERIALS(Imat)%CONSTIT.eq.0) then
        !   write(*,*) 'elem: primal formulation does not support an incompressible linear elastic material!'
        !   stop
        ! endif

        call eval_strain_energy_W_F(imat,x,gradu+DEL,W,dWdF,d2WdF)
! 
! 
!  .....integration weight
        weight = wa*rjac
!
!  .....compute the stiffness tensor
        ! call getC(x, C)
!
!  .....get the source term
        call getf(Mdle,x, fval)
! 
! 
!    P A R T  1 : go through \tau\in H(div)^3 test space (this fills the first set of rows)
!
!
!  .....FIRST OUTER loop through enriched H(div) dofs
        do k1=1,nrdofVV
!  .......OUTER loop through components
          do jcomp=1,3
            m1 = (k1-1)*3+jcomp
!
!           E N R I C H E D   R E S I D U A L   V E C T O R
!
!    - (u,div(\tau)) - (D,\tau)
            EnrResid(m1)  = EnrResid(m1) + weight                 &
                          * (- u(jcomp)*divVV(k1)                          &
                             - dot_product(gradu(jcomp,:),shapvv(:,k1)) )
! 
! 
!
!           G R A M   M A T R I X
!
            select case(TEST_NORM)
!
!   MATHEMATICIAN'S
!   (\tau_2,\tau)+(div(\tau_2),div(tau))
!
            case(2)
              do m2=m1,3*nrdofVV
                icomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( divVV(k1)*divVV(k2)  &
                            + shapVV(1,k1)*shapVV(1,k2)  &
                            + shapVV(2,k1)*shapVV(2,k2)  &
                            + shapVV(3,k1)*shapVV(3,k2) )*weight
                endif
              enddo
            end select
!      ...end of test component loop
          enddo
!    ...end of first outer loop
        enddo
!
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows) 
!
!  .....SECOND OUTER loop through enriched H1 dofs
        do k1=1,nrdofHH
          do jcomp=1,3
            m1 = 3*nrdofVV + (k1-1)*3+jcomp
!
!           E N R I C H E D   R E S I D U A L   V E C T O R
!
!           ( f , v ) - ( P, Grad v )
!
            EnrResid(m1) = EnrResid(m1)                              &
                         + ( fval(jcomp,1)*shapHH(k1)                &
                            -PK1 (jcomp,1)*gradHH(1,k1)              &
                            -PK1 (jcomp,2)*gradHH(2,k1)              &
                            -PK1 (jcomp,3)*gradHH(3,k1)     )*weight
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
! 
! MATHEMATICIANS
!   (grad(v_2),grad(v)) + (v_2,v)
            case(2)
              do m2=m1,3*nrdofHH + 3*nrdofVV
                icomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2 - 3*nrdofVV-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( shapHH(k1)*shapHH(k2)  &
                            + gradHH(1,k1)*gradHH(1,k2)  &
                            + gradHH(2,k1)*gradHH(2,k2)  &
                            + gradHH(3,k1)*gradHH(3,k2) )*weight
                endif
              enddo
! 
            end select
!  .....end of second OUTER loop
          enddo
        enddo
! 
!
!    P A R T  3 : constitutive law as a FOSLS implementation
!
!  FOSLS_Resid = || K(I+D) - P ||^2
!
        FOSLS_Resid = FOSLS_Resid                               &
                    + double_dot_prod(dWdF-PK1,dWdF-PK1)*weight
!  ...end of loop through integration points
      enddo
! 
! 
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
                       nrdofV,shapV,divV)
!         H1 (test)
          call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!         H(div) (test)
          call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!
!  .......geometry map
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,  &
                       x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
          weight = wa*bjac

!
!  .......compute the approximate solution on the PHYSICAL element
          call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                       x,dxdxi,solH,dsolH,solE,dsolE,solV,dsolV,solQ)

          uhat (1:NRHVAR ) = solH (NRHVAR+1:2*NRHVAR ) + LINESEARCH_FACTOR*solH (1:NRHVAR )

          solV_up(:,1:NRVVAR) = solV(:,NRVVAR+1:2*NRVVAR) + LINESEARCH_FACTOR*solV(:,1:NRHVAR)
!
!  .......compute approximate flux vector (solV_up is the transposed stress)
          trac(1:3) = solV_up(1,1:3)*rn(1)  &
                    + solV_up(2,1:3)*rn(2)  &
                    + solV_up(3,1:3)*rn(3)
!
!  .......Change coordinates so the shape functions are on the physical element.
!         H(div) (test)
          do k=1,nrdofVV
            shapVV(1:3,k) = dxdxi(1:3,1)*shapVV(1,k)  &
                          + dxdxi(1:3,2)*shapVV(2,k)  &
                          + dxdxi(1:3,3)*shapVV(3,k)
            shapVV_n(k) = shapVV(1,k)*rn(1)  &
                        + shapVV(2,k)*rn(2)  &
                        + shapVV(3,k)*rn(3)
          enddo
          shapVV_n(1:nrdofVV) = shapVV_n(1:nrdofVV)/rjac
!
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
!
!
!    P A R T  1 : go through \tau\in H(div)^3 test space (this fills the first set of rows)
!
!
!  .......OUTER loop through enriched H(div) test functions
          do k1=1,nrdofVV
            do jcomp=1,3
              m1 = (k1-1)*3+jcomp
! 
! 
!           E N R I C H E D   R E S I D U A L  -  T R A C E   C O N T R I B U T I O N
! 
!           < \hat u,(\tau n)>
! 
              EnrResid(m1) = EnrResid(m1)                   &
                           + uhat(jcomp)*shapVV_n(k1)*weight
!
!        ...end of outer Hdiv loop
            enddo
          enddo
!
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows) 
!
!  .......SECOND OUTER loop through enriched H1 test function
          do k1=1,nrdofHH
!  .........OUTER loop through components
            do jcomp=1,3
              ! counter of row
              m1 = 3*nrdofVV+(k1-1)*3+jcomp
!
!
!           E N R I C H E D   R E S I D U A L  -  T R A C T I O N   C O N T R I B U T I O N
! 
!             <\hat t , v>
!
              EnrResid(m1) = EnrResid(m1)                 &
                           + trac(jcomp)*shapHH(k1)*weight
! 
!
!  .......OUTER loop
            enddo
          enddo
!
!  .....end of loop over integration points
        enddo
        ! if (iprint.eq.1) call pause
      enddo
      if (iprint.ge.1) then
        write(*,7015) EnrResid(1:nrdofHH)
 7015   format('elem_residual: FINAL EnrResid = ',10(/,10(e12.5,2x)))
        call pause
      endif
!
!-----------------------------------------------------------------------------------
!     D P G    L O C A L    A S S E M B L Y                                        |
!-----------------------------------------------------------------------------------
!
!  ...Compact enriched number of rows (total enriched test dof)
      enrdof = 3*nrdofVV+3*nrdofHH
! 
 !            gdump=75
 !        open(unit=gdump,file='output/gram', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(gdump,5999)Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 ! 5999     format(189(e12.5,","))
 !        enddo
 !        close(gdump)
!
!  ...factor the Gram matrix
      uplo = 'U'
      call DPPTRF(uplo, enrdof, Gram, info)
      if (info.ne.0) then
        write(*,*) 'elem_DPG_UWEAK: info = ',info ; stop
      endif
!
!  ...save copies of the RHS to compute later the residual
      EnrResidc = EnrResid
!
!  ...compute the product of inverted test Gram matrix with RHS,
!     EnrResid is overwritten with the solution
      call DPPTRS(uplo, enrdof, 1, Gram, EnrResid, 3*MAXbrickVV+3*MAXbrickHH, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem_residual: info1 = ',info1
        stop 1
      endif
!
!  ...compute the residual
      Resid = 0.d0
!  ...contribution from enriched test functions
      do k=1,3*nrdofHH
        Resid = Resid  &
              + EnrResidc(k)*EnrResid(k)
      enddo
!  ...contribution from L2 test functions
      Resid = Resid + FOSLS_Resid
! !
! !-----------------------------------------------------------------------------------
! !     E R R O R    R E P R E S E N T A T I O N    F U N C T I O N                  |
! !-----------------------------------------------------------------------------------
! !
! !
! !  ...recompute the element residual through direct integration to
! !     establish anistropy flags
!       residd(0:3) = 0.d0
!       do ipt=1,nint
!         xi(1:3) = xiloc(1:3,ipt); wa = wxiloc(ipt)
! !
! !  .....H1 shape functions
!         call shape3H(etype,xi,norder,nedge_orient,nface_orient,  &
!                      nrdofH,shapH,gradH)
! !
! !  .....discontinuous H1 shape functions
!         call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
! !
! !  .....geometry map
!         call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
!                     x,dxdxi,dxidx,rjac,iflag)
!         if (iflag.ne.0) then
!           write(*,1000) Mdle,rjac
!  1000     format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
!           stop
!         endif
! !
! !  .....integration weight
!         weight = rjac*wa
! !
! !  .....Change coordinates so the shape functions are on the physical element
!         do k=1,nrdofHH
!           gradHH(1:3,k) = gradHH(1,k)*dxidx(1,1:3)  &
!                         + gradHH(2,k)*dxidx(2,1:3)  &
!                         + gradHH(3,k)*dxidx(3,1:3)
!         enddo
! !
! !  .....compute the error representation function (for only one load vector)
!         psi(1:3) = 0.d0; gradpsi(1:3,1:3) = 0.d0
!         do k1=1,nrdofHH
!           do icomp=1,3
!             m1 = (k1-1)*3+icomp
! !
!             psi(icomp)         = psi(icomp)         + EnrResid(m1)*shapHH(k1)
!             gradpsi(icomp,1:3) = gradpsi(icomp,1:3) + EnrResid(m1)*gradHH(1:3,k1)
!           enddo
!         enddo
! !
! !  .....compute the stiffness tensor
!         call getC(x, C)
! !
! !  .....compute the semi inner product tensor (depends upon chosen norm)
!         call getSemiIPTensor(x, SIP)

!         !  a weird way of partitioning the H1 norm of psi...
!         do k=1,3; do l=1,3
!           if (k.eq.l) then
!             do m=1,3; do n=1,3
!               residd(k) = residd(k)  &
!                         + SIP(m,k,n,k)*gradpsi(m,k)*gradpsi(n,k)*weight
!             enddo; enddo
!           else
!             do m=1,3; do n=1,3
!               residd(0) = residd(0)  &
!                         + SIP(m,k,n,l)*gradpsi(m,k)*gradpsi(n,l)*weight
!             enddo; enddo
!           endif
!         enddo; enddo
!         do icomp=1,3
!           residd(0) = residd(0) + psi(icomp)**2*weight
!         enddo
! !
! !  ...end of loop through integration points
!       enddo
! !
! !  ...a test
!       diff = residd(0)+residd(1)+residd(2)+residd(3) - Resid
!      if ((abs(diff).gt.1.d-8*abs(Resid)).and.(abs(diff).gt.1.d-14)) then
!        write(*,*) 'Resid = ', Resid
!        write(*,*) 'sum(residd) = ', residd(0)+residd(1)+residd(2)+residd(3)
!        write(*,*) 'residd = ', residd(0:3)
!        call pause
!      endif
! !
! !  ...determine the refinement flag (I see no reason to change this)
!       select case(etype)
!       case('mdlb')
!         if (residd(0).lt..1d0*Resid) then
!           nref(1:3) = 1
!           do i=1,3
!             if (residd(i).lt..1d0*Resid) nref(i)=0
!           enddo
!           Nref_flag = nref(1)*100+nref(2)*10+nref(3)
!         else
!           Nref_flag = 111
!         endif
!       case('mdln','mdld')
!         Nref_flag = 1
!       case('mdlp')
!         if (residd(0).lt..1d0*Resid) then
!           nref(1:2) = 1
!           if (residd(1)+residd(2).lt..2d0*Resid) nref(1)=0
!           if (residd(3).lt..1d0*Resid) nref(2)=0
!           Nref_flag = nref(1)*10+nref(2)
!         else
!           Nref_flag = 111
!         endif
!       end select
!
!-----------------------------------------------------------------------
!
      if (iprint.ge.1) then
        write(*,7010) Mdle, Resid
 7010   format('elem_residual: Mdle, Resid = ',i5,3x,e12.5)
        call pause
      endif
!
!-----------------------------------------------------------------------------------
!             R E F I N E M E N T  F L A G S                                       |
!-----------------------------------------------------------------------------------
!  ...determines what refinement flag to use
!  ...if isotropic h refinements
      call get_isoref(Mdle, Nref_flag)
      ! Nref_flag = 110
!  ...if anisotropic h refinements -> missing
!
!
      end subroutine elem_residual

!--------------------------------------------------------------------------
!> Purpose : returns global residual
!!
!--------------------------------------------------------------------------
!
      subroutine compute_residual
!
      use data_structure3D
      use environment, only : QUIET_MODE
!------------------------------------------------------------------------------------------
      implicit none
!
!  ...middle node
      integer :: mdle
!
!  ...residual
      real*8 :: resid,residual
!
!  ...rate
      real*8 :: rate
!
!  ...number of dof for higher order node
      integer :: ndofH,ndofE,ndofV,ndofQ
      integer :: nrdof_total
!
!  ...refinement flag
      integer :: nref_flag
!
!  ...field variables flag
      integer, dimension(NR_PHYSA) :: nflag
!
!  ...visitation flag
      integer, save :: ivis = 0
!
!  ...total number of dof for the old mesh
      integer, save :: nrdof_total_old
!
!  ...residual for the old mesh
      real*8,  save :: residual_old
!
!  ...residuals and rates to display
      real*8 , dimension(2,10), save :: rwork
!
!  ...number of dof to display
      integer , dimension(10), save :: iwork
!
!  ...miscellaneous
      integer :: i,iel,iphys,iprint
!------------------------------------------------------------------------------------------
!
      iprint=0
!
!  ...field variables flag
      nflag(1:NR_PHYSA)=(/1,0/)
!  ...compute total residual and number of dof
      nrdof_total = NRDOFSH + NRDOFSE + NRDOFSV + NRDOFSQ
      residual = 0.d0
      mdle = 0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_residual(mdle, resid,nref_flag)
        if (iprint.eq.1) then
          write(*,7010) iel, mdle, resid
 7010     format('compute_residual: iel, mdle = ',2i5,  &
                 ' element residual = ',e12.5)
        endif
        residual = residual + resid
        call find_ndof(mdle, ndofH,ndofE,ndofV,ndofQ)
!  .....subtract the bubbles from trace variables
        do iphys=1,NR_PHYSA
!      ...skip if field variable, otherwise remove middle dof and
!         leave trace dof only
          if (nflag(iphys).eq.1) cycle
          select case(DTYPE(iphys))
          case('contin'); nrdof_total=nrdof_total-ndofH*NR_COMP(iphys)
          case('tangen'); nrdof_total=nrdof_total-ndofE*NR_COMP(iphys)
          case('normal'); nrdof_total=nrdof_total-ndofV*NR_COMP(iphys)
          case('discon'); nrdof_total=nrdof_total-ndofQ*NR_COMP(iphys)
          end select
        enddo
      enddo
      residual = sqrt(residual)
!
!  ...compute rate
      rate = 0.d0
      if (ivis.ne.0) then
        if (nrdof_total.gt.nrdof_total_old) then
          rate = log(residual_old/residual)  &
                /log(float(nrdof_total_old)/float(nrdof_total))
        endif
      endif
!
!  ...save current data
      ivis = ivis+1
      nrdof_total_old = nrdof_total
      residual_old = residual
!
!  ...store data to display
      iwork(ivis) = nrdof_total
      rwork(1,ivis) = residual
      rwork(2,ivis) = rate
!
!  ...display the convergence history
      if (.NOT. QUIET_MODE) then
        write(*,*)''
        write(*,*)'         -- Error Report --'
        write(*,7100)
 7100   format(' Mesh  ','  Nrdof  ', ' Residual   ','     Rate ')
        do i=1,ivis
          write(*,7110)i,iwork(i),rwork(1:2,i)
 7110     format(i3,4x,i7,2x,e12.5,2x,f8.3)
        enddo
        write(*,*)''
      endif
!
!
      end subroutine compute_residual
