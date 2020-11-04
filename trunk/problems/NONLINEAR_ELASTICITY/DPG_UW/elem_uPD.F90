!------------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for UW DPG elasticity problem
!! @param[in]  Mdle      - middle node number
! 
!------------------------------------------------------------------------------------------
!       L I N E A R I Z E D    E L A S T O S T A T I C S    E Q U A T I O N S    
!
!     DISPLACEMENT GRADIENT   | LINEAR MOMENTUM BALANCE | PK1 CONSTITUTIVE LAW 
!
!  |   \tau \in H(div)^3      |       v \in (H1)^3      | \chi \in (L^2)^3x3
!
!  - <\delta \hat u,(\tau n)> +            0
!  +               0          -   < \delta \hat t ,v >                       
!  + ( \delta u , div(\tau) ) +            0           
!  +               0          + ( \delta P , grad(v) )  + (    \delta P , \chi  )
!  + ( \delta D , \tau )      +            0            - (L : \delta D , \chi)
!  = 
!                             + (f , v)                 + ( K(I+D),\chi)     
!  - (u,div(\tau))            + <\hat t ,v>             - ( P,\chi )
!  + < \hat u,(\tau n)>       - (P,grad(v))
!  - (D,\tau)
!------------------------------------------------------------------------------------------
!
subroutine elem_DPG_UW_uPD(Mdle)
      use control, only: INTEGRATION
      use parameters
      use parametersDPG
      use data_structure3D
      use element_data
      ! use isotropic_elast_material
      use physics   , only : ADRES,NR_PHYSA      
      use assembly, only: ALOC,BLOC,NR_RHS
      use common_prob_data, only: SYMMETRY_TOL, TEST_NORM
      use hyperelasticity
      use nl_solver_module, only: LINESEARCH_FACTOR
      ! use m_assembly, only:  mdle_list,norderList,nedge_orientList,nface_orientList,xnodList

!------------------------------------------------------------------------------------------
      implicit none
      integer,                                  intent(in)  :: Mdle
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
!  ...MATRICES
!     (TRANSPOSED) stiffness and load matrices for the enriched test space 
      real*8, dimension(3*(MAXbrickH-MAXmdlbH),3*MAXbrickVV+3*MAXbrickHH) :: EnrTraceDispl!,EnrTraceDisplc
      real*8, dimension(3*(MAXbrickV-MAXmdlbV),3*MAXbrickVV+3*MAXbrickHH) :: EnrTraceStress!,EnrTraceStressc
      real*8, dimension(3*MAXbrickQ,3*MAXbrickVV+3*MAXbrickHH) :: EnrFieldDispl!,EnrFieldDisplc
      real*8, dimension(9*MAXbrickQ,3*MAXbrickVV+3*MAXbrickHH) :: EnrFieldStress!,EnrFieldStressc
      real*8, dimension(9*MAXbrickQ,3*MAXbrickVV+3*MAXbrickHH) :: EnrFieldGrad!,EnrFieldStressc
      real*8, dimension(3*MAXbrickVV+3*MAXbrickHH,MAXNRHS_MOD) :: EnrLoad
!     ENRICHED STIFFNESS matrix (excludes least squares part)
      ! real*8, dimension(3*MAXbrickVV+3*MAXbrickHH,3*(MAXbrickH-MAXmdlbH)+3*(MAXbrickV-MAXmdlbV)+21*MAXbrickQ) :: EnrStiffness
      real*8, allocatable :: EnrStiffness(:,:),EnrEverything(:,:)
!     full enriched [stiffness | load] matrix
      ! real*8, dimension(3*MAXbrickVV+3*MAXbrickHH,3*(MAXbrickH-MAXmdlbH)+3*(MAXbrickV-MAXmdlbV)+21*MAXbrickQ+MAXNRHS_MOD) :: EnrEverything
!     Gram matrix for the local Riesz matrix in LAPACK format
      real*8, dimension((3*MAXbrickVV+3*MAXbrickHH)*(3*MAXbrickVV+3*MAXbrickHH+1)/2) :: Gram
!     (TRANSPOSED) stiffness and load matrices for Least Squares equations
      real*8, dimension(9*MAXbrickQ,18*MAXbrickQ) :: FOSLS_Stress
      real*8, dimension(9*MAXbrickQ,18*MAXbrickQ) :: FOSLS_Grad
      real*8, dimension(18*MAXbrickQ,MAXNRHS_MOD) :: FOSLS_Load
!     Assembled local optimal stiffness matrix
      real*8, allocatable :: FullDPG(:,:)
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
                 nrHbub,nrEbub,nrVbub,nrQbub,       gdump,bdump6,bdump7 , trialdof
      integer, dimension(NR_PHYSA) :: ndofphysics
      real*8  :: weight,wa,rjac,brjac,tmp,diffmax,dmax
! 
      ! integer :: clock1,clock2,clockrate,clockmax
!
!  ...LAPACK stuff
      character uplo,transa,transb
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
!
!  ...initialize the enriched local element stiffness matrices and load vectors
      EnrTraceDispl=ZERO; EnrTraceStress=ZERO; EnrFieldDispl=ZERO; EnrFieldStress=ZERO
      EnrFieldGrad=ZERO;FOSLS_Stress=ZERO; FOSLS_Grad=ZERO
      EnrLoad=ZERO; FOSLS_Load=ZERO
! !  ...initialize the Gram matrix
      Gram=ZERO
!
      ! call getSymm(Symm)
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!
!
 !  write(*,1999) NODES(mdle)%dof%zdofQ(1:42,1)
 ! 1999 format(42(es8.1))
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

        call eval_strain_energy_W_F(imat,x,DEL+gradu,W,dWdF,d2WdF)
! 
! 
!  .....integration weight
        weight = wa*rjac
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
!           E N R I C H E D   L O A D   V E C T O R
!
!    - (u,div(\tau)) - (D,\tau)
            EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS) + weight                  &
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
!
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                    
!           T E R M   1 :    ( \delta u , div(\tau) )
! 
!  .........INNER loop through trial dofs for displacement
            do k4=1,nrdofQ
              do icomp=1,3
                m4 = (k4-1)*3+icomp
                if (icomp.eq.jcomp) then
                  EnrFieldDispl(m4,m1) = EnrFieldDispl(m4,m1)       &
                                       + shapQ(k4)*divVV(k1)*weight
                endif
              enddo
            enddo
!         
! 
!           T E R M   2 :    ( \delta D , \tau )
! 
! 
!        ...INNER loop through trial dofs for displacement gradient
            do k5 = 1, nrdofQ
!             icomp ranges 1 through 9, corresponds to number of copy
              icomp = 0
!             j is the column of \delta D, i is the row
              do j=1,3; do i=1,3
                icomp = icomp + 1
!               dof number
                m5 = (k5-1)*9 + icomp
!               contributions come from matching rows of \tau and \delta D (jcomp and i)
                if (jcomp.ne.i) cycle
                EnrFieldGrad(m5,m1) = EnrFieldGrad(m5,m1)             &
                                       + shapQ(k5)*shapVV(j,k1)*weight
              enddo; enddo
!        ...end of displacement gradient inner loop
            enddo
!      ...end of test component loop
          enddo
!    ...end of first outer loop
        enddo
!
!
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows) 
!
!  .....SECOND OUTER loop through enriched H1 dofs
        do k1=1,nrdofHH
          do jcomp=1,3
            m1 = 3*nrdofVV + (k1-1)*3+jcomp
!
!           E N R I C H E D   L O A D   V E C T O R
!
!           ( f , v ) - ( P, Grad v )
!
            EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS) +                    &
                                   ( fval(jcomp,1:NR_RHS)*shapHH(k1)         &
                                    -PK1 (jcomp,1)*gradHH(1,k1)              &
                                    -PK1 (jcomp,2)*gradHH(2,k1)              &
                                    -PK1 (jcomp,3)*gradHH(3,k1)      )*weight
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
                  k2 = int((m2 - 3*nrdofVV -1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( shapHH(k1)*shapHH(k2)  &
                            + gradHH(1,k1)*gradHH(1,k2)  &
                            + gradHH(2,k1)*gradHH(2,k2)  &
                            + gradHH(3,k1)*gradHH(3,k2) )*weight
                endif
              enddo
! 
            end select
!
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!
!           ( \delta P , grad(v) )
!
!  .........INNER loop through trial dofs
            do k3=1,nrdofQ
!             icomp ranges 1 through 9, corresponds to number of copy
              icomp = 0
!             j is the column of \delta P, i is the row
              do j=1,3; do i=1,3
                icomp = icomp + 1
!               dof number
                m3 = (k3-1)*9 + icomp
!               contributions come from matching row of \delta P and comp of v (i and jcomp)
                if (jcomp.ne.i) cycle
                EnrFieldStress(m3,m1) = EnrFieldStress(m3,m1)             &
                                      + shapQ(k3)*gradHH(j,k1)*weight
              enddo; enddo
!  .........INNER loop
            enddo
!  .....OUTER loop
          enddo
        enddo
! 
!
!    P A R T  3 : linearized constitutive law as a FOSLS implementation
!
!  + (     \delta P , \chi  )
!  - ( L : \delta D , \chi  )
! =
!    ( K(I+D) - P , \chi )
!     
!    by FOSLS, \chi = \delta P` - L : \delta D`
!
!  ...loop over optimal test function \delta P' 
        do k1=1,nrdofQ
          jcomp = 0
          do l=1,3; do k=1,3
            jcomp = jcomp + 1
            m1 = (k1-1)*9 + jcomp
!
!
!     F O S L S   L O A D   V E C T O R
!
!     T E R M   1 :     ( K(I+D) - P , \delta P` )
!
            FOSLS_Load(m1,1) = FOSLS_Load(m1,1)                        &
                             + ( dWdF(k,l)-PK1(k,l) )*shapQ(k1)*weight
!
!
!     F O S L S   S T I F F N E S S   M A T R I X
!
!     T E R M   1:    (\delta P, \delta P`)
!
            do k3=1,nrdofQ
              icomp = 0
              do j=1,3; do i=1,3
                icomp = icomp + 1
                m3 = (k3-1)*9 + icomp
                if (icomp.ne.jcomp) cycle
!
                FOSLS_Stress(m3,m1) = FOSLS_Stress(m3,m1)              &
                                    + (shapQ(k1)*shapq(k3))*weight
!
              enddo; enddo
            enddo
!
!     T E R M   2:   - ( L : \delta D , \delta P`)
!
            do k5=1,nrdofQ
              icomp = 0
              do j=1,3; do i=1,3
                icomp = icomp + 1
                m5 = (k5 - 1)*9 + icomp
!
                FOSLS_Grad(m5,m1) = FOSLS_Grad(m5,m1)                     &
                                  - d2WdF(i,j,k,l)*shapQ(k1)*shapQ(k5)*weight
!
              enddo; enddo
            enddo
!    ...end of \delta P` loop
          enddo;enddo
        enddo
!
!
!  ...loop over optimal test function \delta D` 
        do k1=1,nrdofQ
          jcomp = 0
          do l=1,3; do k=1,3
            jcomp = jcomp + 1
            m1 = (nrdofQ + k1-1)*9 + jcomp
!
!
!     L O A D   V E C T O R
!
!     T E R M   2 :    -( K(I+D) - P , L : \delta D` )
!
            FOSLS_Load(m1,1) = FOSLS_Load(m1,1)                         &
                             - double_dot_prod(dWdF-PK1,d2WdF(:,:,k,l)) &
                             * shapQ(k1)*weight
!
!
!     F O S L S   S T I F F N E S S   M A T R I X
!
!     T E R M   4:    ( L : \delta D , L : \delta D`)
!
            do k3=1,nrdofQ
              icomp = 0
              do j=1,3; do i=1,3
                icomp = icomp + 1
                m3 = (k3-1)*9 + icomp
!
                FOSLS_Grad(m3,m1) = FOSLS_Grad(m3,m1)                              &
                                  + shapQ(k1)*shapQ(k3)                            &
                                  * double_dot_prod(d2WdF(:,:,i,j),d2WdF(:,:,k,l)) &
                                  * weight
!
              enddo; enddo
            enddo
!      ...end of outer loop over \delta F`
          enddo;enddo
        enddo
!  ...end of quadrature point loop
      enddo
!
!     F O S L S   S T I F F N E S S   M A T R I X
!
!     T E R M   3 :   we take the transpose of term 2
!
      FOSLS_Stress(1:9*nrdofQ,9*nrdofQ+1:18*nrdofQ) = transpose(FOSLS_Grad(1:9*nrdofQ,1:9*nrdofQ))
!
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
                       x,dxdxi,dxidx,rjac,dxdt,rn,brjac)
          weight = wa*brjac

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
!           E N R I C H E D   L O A D  -  T R A C E   C O N T R I B U T I O N
! 
!           < \hat u,(\tau n)>
! 
              EnrLoad(m1,1) = EnrLoad(m1,1)                   &
                            + uhat(jcomp)*shapVV_n(k1)*weight
!
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!
!         - <\delta \hat u,(\tau n)>
! 
!  ...........INNER loop through H1 bdry trial dofs
              do k2=1,nrdofH-nrHbub
                do icomp=1,3
                  m2 = (k2-1)*3+icomp
                  if (icomp.eq.jcomp) then
                    EnrTraceDispl(m2,m1) = EnrTraceDispl(m2,m1)          &
                                         - shapH(k2)*shapVV_n(k1)*weight
                  endif
                enddo
              enddo
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
!           E N R I C H E D   L O A D  -  T R A C T I O N   C O N T R I B U T I O N
! 
!             <\hat t , v>
!
              EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS)         &
                                   + trac(jcomp)*shapHH(k1)*weight
! 
!
!           E N R I C H E D   T R A C T I O N   S T I F F N E S S   M A T R I X
!
!   - <  \delta \hat t , v >
!
!  ...........INNER loop through enriched dofs
              do k2=1,nrdofV-nrVbub
            !  contribute only when components match (populate only the diagonal)
                m2 = (k2-1)*3+jcomp
                EnrTraceStress(m2,m1) = EnrTraceStress(m2,m1)         &
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
      ! iprint=2
      if (iprint.eq.2) then
        write(*,*) 'Gram = '
        do i=1,25
          write(*,6000) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 6000     format('i = ',i3,'  ',25e12.5)
        enddo

        call pause

        write(*,*) 'EnrFieldDispl = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrFieldDispl(1:3*nrdofQ,i)
 6001     format('i = ',i4,'  ',19(/,10e12.5))
        enddo

        call pause

        write(*,*) 'EnrFieldStress = '
        do i=1,3*nrdofVV
          write(*,6001) i,EnrFieldStress(1:9*nrdofQ,i)
        enddo

        call pause

        do i=1+3*nrdofVV,3*nrdofVV+3*nrdofHH
          write(*,6001) i,EnrFieldStress(1:9*nrdofQ,i)
        enddo

        call pause

        write(*,*) 'EnrFieldOmega = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrFieldGrad(1:9*nrdofQ,i)
        enddo

        call pause

        write(*,*) 'EnrTraceDispl = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrTraceDispl(1:3*(nrdofH-nrHbub),i)
        enddo

        call pause

        write(*,*) 'EnrTraceStress = '
        do i=3*nrdofVV,3*nrdofVV+3*nrdofHH
          write(*,6001) i,EnrTraceStress(1:3*(nrdofV-nrVbub),i)
        enddo

        call pause

        write(*,*) 'EnrLoad = '
        do i=1,189
          write(*,6001) i,EnrLoad(i,1)
        enddo

        call pause
        write(*,*) 'FOSLS_Load = '
        do i=1,18
          write(*,6001) i,FOSLS_Load(i,1)
        enddo

        call pause

      endif
! 
! 
!
      !write(*,*) 'elem: before local Riesz problem...'
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
 !        close(gdump); call pause
!
!  ...factor the Gram matrix
      uplo = 'U'
      call DPPTRF(uplo, enrdof, Gram, info)
      if (info.ne.0) then
        write(*,*) 'elem_DPG_UWEAK: info = ',info ; stop
      endif
! 
 !            gdump=76
 !        open(unit=gdump,file='output/gram_fact', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(gdump,5999)Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 ! ! 5999     format(189(e12.5,","))
 !        enddo
 !        close(gdump)
!
!  ...Create vector of indices with dof of each physics variable (in the same order of physics file)
      ndofphysics = (/3*(nrdofH-nrHbub),3*(nrdofV-nrVbub),3*nrdofQ,9*nrdofQ,9*nrdofQ/)

      trialdof = sum(ndofphysics)
!
!  ...Construct EnrEverything by packing all Enriched Matrices in the order H1,Hcurl,Hdiv,L2 
!     (in the same order of physics file) and load
      allocate(EnrStiffness(enrdof,trialdof))
      EnrStiffness=ZERO
      allocate(EnrEverything(enrdof,trialdof+NR_RHS))
      EnrEverything=ZERO
!  ...EnrTraceDispl (H1)
      kmin=0
      kmax=kmin+ndofphysics(1)
      EnrStiffness(1:enrdof,kmin+1:kmax)=transpose(EnrTraceDispl(1:(kmax-kmin),1:enrdof))
!  ...EnrTraceStress (H(div))
      kmin=kmax
      kmax=kmin+ndofphysics(2)
      EnrStiffness(1:enrdof,kmin+1:kmax)=transpose(EnrTraceStress(1:(kmax-kmin),1:enrdof))
!  ...EnrFieldDispl (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(3)
      EnrStiffness(1:enrdof,kmin+1:kmax)=transpose(EnrFieldDispl(1:(kmax-kmin),1:enrdof))
!  ...EnrFieldStress (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(4)
      EnrStiffness(1:enrdof,kmin+1:kmax)=transpose(EnrFieldStress(1:(kmax-kmin),1:enrdof))
!  ...EnrFieldGrad (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(5)
      EnrStiffness(1:enrdof,kmin+1:kmax)=transpose(EnrFieldGrad(1:(kmax-kmin),1:enrdof))
!  ...store in right position in EnrEverything
      EnrEverything(1:enrdof,1:kmax) = EnrStiffness(1:enrdof,1:kmax)
!  ...EnrLoad
      kmin=kmax
      kmax=kmin+NR_RHS
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrLoad(1:enrdof,1:(kmax-kmin))
!
      if (iprint.eq.1) then
        write(*,*) 'EnrEverything = '
        do i=1,enrdof
          write(*,6001) i,EnrEverything(i,1:kmax)
        enddo
      endif
 !            bdump7=82
 !      open(unit=bdump7,file='output/EnrEverything', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(bdump7,5995) EnrEverything(i,1:kmax)
 ! 5995     format(69(e12.5,","))        
 !        enddo
 !      close(bdump7)

!
!
!  ...G^-1 * EnrEverything
      ! call DPPTRS(uplo,enrdof,kmax,Gram,EnrEverything(:,1:kmax),3*MAXbrickVV+3*MAXbrickHH,info1)
      call DPPTRS(uplo,enrdof,kmax,Gram,EnrEverything(:,1:kmax),enrdof,info1)
      if (info1.ne.0) then
        write(*,*) 'elem_DPG_UWEAK: info1 = ',info1 ; stop
      endif
! 
 !            bdump7=88
 !      open(unit=bdump7,file='output/EnrEverything_post', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(bdump7,5995) EnrEverything(i,1:kmax)
 ! ! 5995     format(66(e12.5,","))        
 !        enddo
 !      close(bdump7)
!
!  ...Build full DPG matrix (stiffness + load) in one go
      allocate(FullDPG(kmin,kmax))
      FullDPG = ZERO
      transa = 'T'
      transb = 'N'
      m = kmin
      n = kmax
      k = enrdof
      l = enrdof! 3*MAXbrickVV+3*MAXbrickHH
      ! call DGEMM(transa,transb,m,n,k,1.d0,EnrStiffness,l,EnrEverything,l,0.d0,FullDPG,m)!3*(MAXbrickH-MAXmdlbH)+3*(MAXbrickV-MAXmdlbV)+12*MAXbrickQ)
      call DGEMM(transa,transb,m,n,k,1.d0,EnrStiffness,l,EnrEverything,l,0.d0,FullDPG,m)
!
      !write(*,*) 'elem: after local Riesz problem...'
 !            bdump7=88
 !      open(unit=bdump7,file='output/FullDPG', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,kmin
 !          write(bdump7,5995) FullDPG(i,1:kmax)
 ! ! 5995     format(66(e12.5,","))        
 !        enddo
 !      close(bdump7)

      k3 = ndofphysics(1)+ndofphysics(2)+ndofphysics(3)
      k4 = k3 + ndofphysics(4)
      k5 = k4 + ndofphysics(5)

      FullDPG(k3+1:k4,k3+1:k4) =  FullDPG(k3+1:k4,k3+1:k4)                        &
                               +  FOSLS_Stress(1:ndofphysics(4),1:ndofphysics(4))

      FullDPG(k3+1:k4,k4+1:k5) =  FullDPG(k3+1:k4,k4+1:k5)                        &
                               +  FOSLS_Stress(1:ndofphysics(4),ndofphysics(4)+1:ndofphysics(4)+ndofphysics(5))

      FullDPG(k4+1:k5,k3+1:k4) =  FullDPG(k4+1:k5,k3+1:k4)                        &
                               +  FOSLS_Grad(1:ndofphysics(5),1:ndofphysics(4))

      FullDPG(k4+1:k5,k4+1:k5) =  FullDPG(k4+1:k5,k4+1:k5)                        &
                               +  FOSLS_Grad(1:ndofphysics(5),ndofphysics(4)+1:ndofphysics(4)+ndofphysics(5))

      FullDPG( k3+1: k5, k5+1) =  FullDPG( k3+1: k5, k5+1)            &
                               +  FOSLS_Load(1:k5-k3,1) 
      ! write(*,*) 'elem: after local Riesz problem...'
 !            bdump7=84
 !      open(unit=bdump7,file='output/FOSLS', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,kmin
 !          write(bdump7,5995) FullDPG(i,1:kmax)
 ! 5995     format(66(e12.5,","))        
 !        enddo
 !      close(bdump7)
 !      call pause
!
!
!     ULTIMATE DPG LOAD VECTORS AND STIFFNESS MATRICES THROUGH STATIC CONDENSATION
!
!  ...Populate the assembly matrices accordingly
      n1=0
      do i=1,NR_PHYSA
        n=ndofphysics(i)
        n2=n1+n
        m1=0
        do j=1,NR_PHYSA
          m=ndofphysics(j)
          m2=m1+m
!  .......First initialize
          ALOC(i,j)%array = ZERO
!  .......Then populate
          ALOC(i,j)%array(1:n,1:m) = FullDPG(n1+1:n2,m1+1:m2)
          m1=m2  
        enddo
        m2=m1+NR_RHS
!  .....First initialize
        BLOC(i)%array = ZERO
!  .....Then populate
        BLOC(i)%array(1:n,1:NR_RHS) = FullDPG(n1+1:n2,m1+1:m2)
        n1=n2
      enddo
      !write(*,*) 'elem: before test part...'
!
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'EnrStiffness = '
        do i=1,enrdof
          write(*,6001) i,EnrStiffness(i,1:kmin)
        enddo
        write(*,*) 'Gram^-1*EnrEverything = '
        do i=1,enrdof
          write(*,6001) i,EnrEverything(i,1:kmax)
        enddo
        write(*,*) 'FullDPG = '
        do i=1,kmin
          write(*,6001) i,FullDPG(i,1:kmax)
        enddo
      endif
      ! call pause
! 
      deallocate(FullDPG)
      deallocate(EnrStiffness,EnrEverything)
!
end subroutine elem_DPG_UW_uPD
! 
subroutine elem_DPG_UW_uPD_test(Mdle)
  integer :: Mdle
  end subroutine