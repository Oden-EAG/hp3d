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

  use physics   , only : NR_PHYSA
  use data_structure3D_poly
!--------------------------------------------------------------------------
  implicit none
  integer,                    intent(in)  :: Mdle
  integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!--------------------------------------------------------------------------

  Itest (1:NR_PHYSA) = 1
  Itrial(1:NR_PHYSA) = 1

! $OMP critical
  ! write(*,*) 'mdle = ', mdle
  call elem_DPG_UWEAK_poly(Mdle)
! $OMP end critical

!
end subroutine elem

!------------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for PolyDPG - UW elasticity problem
!! @param[in]  Mdle      - middle node number
!------------------------------------------------------------------------------------------
!
!                |   \tau \in H(div)^3   |       v \in (H1)^3     |
!
!                   - <\hat u,(\tau n)>  +            0
!  +                        0            + - <\hat (\sigma n),v>
!  + \int_\Omega [   u \cdot div(\tau)   +            0           ]
!  + \int_\Omega [    A \sigma : \tau    +    \sigma : grad(v)    ]
!  + \int_\Omega [     \omega : \tau     +            0           ]
!  = \int_\Omega [          0            +        f \cdot v       ]
!
!------------------------------------------------------------------------------------------
!
subroutine elem_DPG_UWEAK_poly(Mdle)

      use uweak_module_poly
      use connectivity_poly
      use control, only: INTEGRATION
      use parameters
      use parametersDPG
      use data_structure3D_poly
      use element_data_poly
      use isotropic_elast_material_compo
      use physics   , only : NR_PHYSA
      use assembly_poly  , only : ALOC,BLOC,NR_RHS
      use common_prob_data, only: SYMMETRY_TOL, TEST_NORM
      use geometry_polydpg,only: ALLC_DOM
!------------------------------------------------------------------------------------------
      implicit none
!
      integer, intent(in) :: Mdle

!  ...MATRICES
!     stiffnes and load matrices for the enriched test space
  real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),(3*MAXtriaH*MAX_NRFC)) :: EnrTraceDispl!,EnrTraceDisplc
  real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),3*MAXtriaQ*MAX_NRFC) :: EnrTraceStress!,EnrTraceStressc
  real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),3*MAXtetraQ) :: EnrFieldDispl!,EnrFieldDisplc
  real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),6*MAXtetraQ) :: EnrFieldStress!,EnrFieldStressc
  real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),3*MAXtetraQ) :: EnrFieldOmega!,EnrFieldOmegac
  real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),MAXNRHS_MOD) :: EnrLoad
  real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),                                               &
                    3*(MAXtriaQ+MAXtriaH)*MAX_NRFC+12*MAXtetraQ) :: EnrStiffness
  real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),                                               &
                    3*(MAXtriaQ+MAXtriaH)*MAX_NRFC+12*MAXtetraQ+MAXNRHS_MOD) :: EnrEverything
  ! real*8, dimension((3*(MAXtriaQ+MAXtriaH)*MAX_NRFC+12*MAXtetraQ),                             &
  !                   3*(MAXtriaQ+MAXtriaH)*MAX_NRFC+12*MAXtetraQ+MAXNRHS_MOD) :: FullDPG
!     Gram matrix for the local Riesz matrix in LAPACK format
  real*8, dimension((3*MAXtetraVV+3*MAXtetraHH)*(3*MAXtetraVV+3*MAXtetraHH+1)/2) :: Gram
!  ...element and face type
      character(len=4) :: etype,ftype
!
!  ...number of topological entities (vertices,edges,faces)
      integer :: nrv,nre,nrf
!
!  ...element and face order, enriched order
      integer, dimension(19) :: norder
      integer, dimension(5)  :: nordf
      integer                :: nordP
!
!  ...node edge and face orientations
      integer, dimension(12) :: nedge_orient
      integer, dimension(6)  :: nface_orient
      integer, dimension(4)  :: norie
!
!  ...SHAPE FUNCTIONS
!     H1  (geometry and trial)
      real*8, dimension(  MAXtetraH)  :: shapH
      real*8, dimension(3,MAXtetraH)  :: gradH
      integer                         :: nrdofH
!     H(div)  (trial)
      real*8, dimension(3,MAXtetraV)  :: shapV
      real*8, dimension(  MAXtetraV)  :: divV
      real*8, dimension(  MAXtetraV)  :: shapV_n
      integer                         :: nrdofV
!     Discont Face H1 (trial)      
      real*8, dimension(MAXtriaH   )  :: shapH_f
      real*8, dimension(2,MAXtriaH )  :: gradH_f
!     Face L2 (trial)      
      real*8, dimension(MAXtriaQ   )  :: shapQ_f
!     L2  (trial)
      real*8, dimension(  MAXtetraQ)  :: shapQ
      integer                         :: nrdofQ
!     H1   (test)
      real*8, dimension(  MAXtetraHH) :: shapHH
      real*8, dimension(3,MAXtetraHH) :: gradHH
      integer                         :: nrdofHH
!     H(div)  (test)
      real*8, dimension(3,MAXtetraVV) :: shapVV
      real*8, dimension(  MAXtetraVV) :: divVV
      real*8, dimension(  MAXtetraVV) :: shapVV_n
      integer                         :: nrdofVV
!
!
!  ...geometry
      real*8, dimension(3,MAXtetraH) :: xnod
      real*8, dimension(3)           :: xi,x,rn,fn,x_aux,b0,b1,b2,a0,a1,a2,a3,fn_a
      real*8, dimension(3,3)         :: dxdxi_e,dxidx_e,Qrot_e,Qrot_f,maptet,maptetinv,b_aff
      real*8, dimension(2,2)         :: dxdxi_f,dxidx_f,dxdxi_bf,maptri
      real*8, dimension(2)           :: t
      real*8, dimension(3,2)         :: dxidt,dxdt,rt,dxi_eta
      real*8, dimension(3,4)         :: a_aff
      integer                        :: nsign
!
!  ...tensors in physical coordinates
      real*8, dimension(3,3,3,3) :: A, AA
!
!  ...source term (don't need Neumann term)
      real*8, dimension(3,MAXNRHS) :: fval
!
!  ...3D quadrature data
      real*8, allocatable :: xiloc(:,:),xloc(:,:)
      real*8, allocatable :: wxi(:)
!
!  ...2D quadrature data for boundary terms
      real*8, allocatable :: tloc(:,:)
      real*8, allocatable :: wt(:)
!
!  ...miscellaneous
      integer :: i,j,k,l,m,n,k1,k2,k3,k4,k5,m1,m2,m3,m4,m5,n1,n2,n3,ipt,ifc,  &
                 icomp,jcomp,kcomp,lcomp,nint,iprint,iflag,info,info1,info2,info3,  &
                 kH,kV,kQ,lH,lV,lQ,kmin,kmax,enrdof,       &
                 Nrv_f,mdlf,jv,jf,loc,nrdofQ_f,nrdofH_f, &
                 gdump,bdump1,bdump2,bdump3,bdump4,bdump5,bdump6,bdump7,  &
                 nord_add_local,nord_add_ini,  imat , idom

      integer, dimension(NR_PHYSA) :: ndofphysics
      real*8  :: weight,wa,rjac,brjac,tmp,diffmax,dmax,rjac_e,rjac_f,rjac_bf,r_e, &
                 maptetdet,maptridet,Area_f,r_f, dir,b_area,a_vol , alpha

      integer, allocatable :: nfaces(:),Norientf(:),Nedges(:),Nverts(:),nverts_f(:)
      real*8, allocatable :: xverts(:,:),xverts_f(:,:),                             FullDPG(:,:)


!
!  ...LAPACK stuff
      character :: uplo,transa,transb

      integer :: liwork,lwork,ldz
      integer,allocatable :: iwork(:)
      real*8,allocatable ::work(:),eval(:),evec(:,:)
! NOTE: nk is a "statement function"
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
!       write(*,*) 'mdle = ', mdle
! goto 90
      iprint=0
      ! pause
!
!  ...Coefficient for L^2 term in adjoint-graph norm
      alpha=1.d0
!     
      etype = 'tetr'
      ftype = 'tria'
!     retrieve element node data
      call elem_nodes_poly(Mdle,Nfaces,Norientf,Nrf,Nedges,Nre,Nverts,Nrv)
! !     allocate array for vertices 
      allocate(xverts(3,Nrv))
!     get vertices coordinates
      do jv = 1, nrv
        xverts(1:3,jv) = NODES(nverts(jv))%coord(1:3,1)
      enddo
!
!     get element affine coordinates points a0,a1,a2,a3
      call poly_affine_3d(xverts,nrv,a0,a1,a2,a3)
!     construct array a_aff with a0,a1,a2,a3 as columns
      a_aff(:,1) = a0; a_aff(:,2) = a1; a_aff(:,3) = a2; a_aff(:,4) = a3
! ! computational of jacobian determinant
! ! get normal of face a1a2a3
!       call face_normal(a_aff(:,2:4),fn_a)
! ! get tetrahedron volume
!       call tetra_centr_vol(a0,a1,a2,a3,fn_a,x_aux,a_vol)
! ! correct
!       rjac_e = 6.d0*a_vol
! construct linear transformation associated to affine map. Interpreted as classic jacobian
      dxdxi_e(:,1) = a1 - a0
      dxdxi_e(:,2) = a2 - a0
      dxdxi_e(:,3) = a3 - a0
! get inverse and determinant of linear transformation
      call geom(dxdxi_e,dxidx_e,rjac_e,iflag)
!
      if (iflag.ne.0) then
        write(*,8999) Mdle,rjac
 8999   format('Negative Jacobian!Mdle,rjac=',i8,2x,e12.5)
        stop
      endif
!     get order, in the norder structure of a uniform tetrahedral element
      norder = 0
      norder(1:15) = NODES(Mdle)%order      
!
!
!     set nord_add_local
      nord_add_ini = NORD_ADD
      if (VARIABLE_DP.eq.1) then
        call correct_dp(Mdle,nord_add_ini,NODES(Mdle)%order,nrf,nord_add_local)
      ! write(*,*) 'mdle, nord_add_local = ',Mdle,nord_add_local
      else
        nord_add_local = NORD_ADD
      endif
!
!  ...set the enriched order of appoximation
      nordP = NODES(Mdle)%order + nord_add_local
! 
!  ...set up the element quadrature
      INTEGRATION = max(0, nord_add_local-1)
      allocate(xiloc(3,(Nrf-3)*(MAX_NRVF-2)*MAXNINT3ADD),wxi((Nrf-3)*(MAX_NRVF-2)*MAXNINT3ADD))
      call set_quadrature_affine_elem(mdle,Nrf,NODES(Mdle)%order,a_aff,nint,xiloc,wxi)
      INTEGRATION = 0
      allocate(xloc(3,nint))
      xloc = 0.d0
      call poly_affine_3d_to_x(xiloc(:,1:nint),nint,a_aff,xloc)
!  ...initialize the enriched local element stiffness matrices and load vectors
      EnrTraceDispl=ZERO; EnrTraceStress=ZERO; EnrFieldDispl=ZERO; EnrFieldStress=ZERO; EnrFieldOmega=ZERO
      EnrLoad=ZERO
! !  ...initialize the Gram matrix
      Gram=ZERO
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!      

!        ...loop through integration points
      do ipt=1,nint
!              
        xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt); x(1:3) = xloc(1:3,ipt)
!       
!        .....Compute shape functions needed for test/trial field variables and geometry
!             L2 (field trial)
        call shape3Q(etype,xi,norder, nrdofQ,shapQ)              
!             H1 (test)
        call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)              
!             H(div) (test)
        call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!
!       .....Change coordinates so the shape functions are on the physical element
!             L2 (trial)
        shapQ(1:nrdofQ) = shapQ(1:nrdofQ)/rjac_e
!             H1 (test)
        do k=1,nrdofHH
          gradHH(1:3,k) = gradHH(1,k)*dxidx_e(1,1:3)  &
                        + gradHH(2,k)*dxidx_e(2,1:3)  &
                        + gradHH(3,k)*dxidx_e(3,1:3)
        enddo
!             H(div) (test)
        do k=1,nrdofVV
          shapVV(1:3,k) = dxdxi_e(1:3,1)*shapVV(1,k)  &
                        + dxdxi_e(1:3,2)*shapVV(2,k)  &
                        + dxdxi_e(1:3,3)*shapVV(3,k)
        enddo
        shapVV(1:3,1:nrdofVV) = shapVV(1:3,1:nrdofVV)/rjac_e
        divVV(1:nrdofVV) = divVV(1:nrdofVV)/rjac_e
!
!        .....integration weight
        weight = wa !*rjac_e
!
        ! call find_material(x,(/0.d0,0.d0,0.d0/),imat)
        idom = ALLC_DOM(mdle)
        select case(idom)
        case(1)
          imat = 1
        case default
          imat = 2
        end select
!        .....compute the compliance tensor
        call getA(imat , x, A)
!
!  .....need this for the adjoint graph norm
        call getAA(imat , X, AA)
!
!        .....get the source term
        call getf(Mdle,x, fval)
!
!
!    P A R T  1 : go through \tau\in H(div)^3 test space (this fills the first set of rows)
!
!
!        .....FIRST OUTER loop through enriched H(div) dofs
        do k1=1,nrdofVV
!        .......OUTER loop through components
          do jcomp=1,3
            m1 = (k1-1)*3+jcomp
!     
!                 E N R I C H E D   L O A D   V E C T O R
!
!   0
!
!           G R A M   M A T R I X
!
            select case(TEST_NORM)
!
!   ADJOINT GRAPH
!   (A:\tau_2,A:\tau)+(div(\tau_2),div(tau))+(\tau_2,\tau)+alpha*(\tau_2,\tau)
!
            case(1)
              do m2=m1,3*nrdofVV
                icomp = mod(m2-1,3)+1
                k = nk(m1,m2)
                k2 = int((m2-1)/3)+1
                do kcomp=1,3
                  do lcomp=1,3
                    if ((icomp.eq.jcomp.and.kcomp.eq.lcomp).or. &
                        (icomp.eq.kcomp.and.jcomp.eq.lcomp).or. &
                        (icomp.eq.lcomp.and.jcomp.eq.kcomp))    &
                      Gram(k) = Gram(k)  &
                              + ( AA(icomp,kcomp,jcomp,lcomp)                 &
                                 *shapVV(lcomp,k1)*shapVV(kcomp,k2) )*weight
                  enddo
                enddo
                if (icomp.eq.jcomp) then
                  Gram(k) = Gram(k)  + divVV(k1)*divVV(k2)*weight
                  Gram(k) = Gram(k)  &
                          + (1.d0+alpha)*(               &
                            + shapVV(1,k1)*shapVV(1,k2)  &
                            + shapVV(2,k1)*shapVV(2,k2)  &
                            + shapVV(3,k1)*shapVV(3,k2) )*weight
                endif
              enddo
!
!   (A:\tau_2, \grad v)
!
              do m2=3*nrdofVV+1,3*nrdofVV+3*nrdofHH
                icomp = mod(m2-1,3)+1
                k = nk(m1,m2)
                k2 = int((m2-3*nrdofVV-1)/3)+1
                do kcomp=1,3
                  do lcomp=1,3
                      if ((icomp.eq.jcomp.and.kcomp.eq.lcomp).or. &
                        (icomp.eq.kcomp.and.jcomp.eq.lcomp).or. &
                        (icomp.eq.lcomp.and.jcomp.eq.kcomp))    &
                        Gram(k) = Gram(k)  &
                                + ( A(icomp,kcomp,jcomp,lcomp)                 &
                                   *shapVV(lcomp,k1)*gradHH(kcomp,k2) )*weight
                  enddo
                enddo
              enddo
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
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   C A U C H Y    S T R E S S   )
!
!   A \sigma : \tau
!
!          ( sigma1   sigma4  sigma5 )
!  sigma = ( sigma4   sigma2  sigma6 )
!          ( sigma5   sigma6  sigma3 )
!
!  .........INNER loop through trial dofs for Cauchy stress
            do k3=1,nrdofQ
              do icomp=1,6
                m3 = (k3-1)*6+icomp
                if (icomp.eq.1) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(1,1,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                elseif (icomp.eq.2) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(2,2,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                elseif (icomp.eq.3) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(3,3,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                elseif (icomp.eq.4) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(1,2,jcomp,n)*shapQ(k3)*shapVV(n,k1) &
                        + A(2,1,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                elseif (icomp.eq.5) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(1,3,jcomp,n)*shapQ(k3)*shapVV(n,k1) &
                        + A(3,1,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                elseif (icomp.eq.6) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(2,3,jcomp,n)*shapQ(k3)*shapVV(n,k1) &
                        + A(3,2,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                endif
                EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                      + tmp*weight
              enddo
            enddo
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   u \cdot div(\tau)
!
!  .........INNER loop through trial dofs for displacement
            do k4=1,nrdofQ
              do icomp=1,3
                m4 = (k4-1)*3+icomp
                if (icomp.eq.jcomp) then
                  EnrFieldDispl(m1,m4) = EnrFieldDispl(m1,m4)  &
                                       + shapQ(k4)*divVV(k1)*weight
                endif
              enddo
            enddo
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                   (   L A G R A N G E    M U L T I P L I E R   )
!
!   + omega : \tau
!
!          (   0     omega3  -omega2 )    1
!  omega = (-omega3    0      omega1 ) = --- ( grad(u)-grad(u)^T)
!          ( omega2  -omega1    0    )    2
!
!  .........INNER loop through trial dofs for omega (antisymmetric part of displacement gradient)
            do k5=1,nrdofQ
              do icomp=1,3
                m5 = (k5-1)*3+icomp
                if (icomp.eq.1) then
                  if (jcomp.eq.2) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         + shapQ(k5)*shapVV(3,k1)*weight
                  elseif (jcomp.eq.3) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         - shapQ(k5)*shapVV(2,k1)*weight
                  endif
                elseif (icomp.eq.2) then
                  if (jcomp.eq.1) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         - shapQ(k5)*shapVV(3,k1)*weight
                  elseif (jcomp.eq.3) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         + shapQ(k5)*shapVV(1,k1)*weight
                  endif
                elseif (icomp.eq.3) then
                  if (jcomp.eq.1) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         + shapQ(k5)*shapVV(2,k1)*weight
                  elseif (jcomp.eq.2) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         - shapQ(k5)*shapVV(1,k1)*weight
                  endif
                endif
              enddo
            enddo
!  .....END OUTER LOOP through test stresses
          enddo
        enddo
!
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows)
!
!        .....SECOND OUTER loop through enriched H1 dofs
        do k1=1,nrdofHH
!        .......OUTER loop through components
          do jcomp=1,3
      ! counter of row
            m1 = 3*nrdofVV+(k1-1)*3+jcomp
!
!
!            E N R I C H E D   L O A D   V E C T O R
!
!   f \cdot v
!
            EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS) &
                           + fval(jcomp,1:NR_RHS)*shapHH(k1)*weight
!
!            G R A M   M A T R I X
!
            select case(TEST_NORM)
!
!   ADJOINT GRAPH
!   (\grad(v_1),A:\tau_2+\grad(v_2))+alpha*(v_1,v_2)
!
            case(1)

              do m2=m1,3*nrdofVV+3*nrdofHH
                icomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-3*nrdofVV-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( alpha*shapHH(k1)*shapHH(k2)  &
                            + gradHH(1,k1)*gradHH(1,k2)  &
                            + gradHH(2,k1)*gradHH(2,k2)  &
                            + gradHH(3,k1)*gradHH(3,k2) )*weight
                endif
              enddo
!
!   MATHEMATICIAN'S
!   (v_2,v) + (grad(v_2),grad(v))
!
            case(2)
              do m2=m1,3*nrdofVV+3*nrdofHH
                icomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-3*nrdofVV-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( shapHH(k1)*shapHH(k2)  &
                            + gradHH(1,k1)*gradHH(1,k2)  &
                            + gradHH(2,k1)*gradHH(2,k2)  &
                            + gradHH(3,k1)*gradHH(3,k2) )*weight
                endif
              enddo
            end select
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   C A U C H Y    S T R E S S   )
!
!   + \sigma : grad(v)
!
!          ( sigma1   sigma4  sigma5 )
!  sigma = ( sigma4   sigma2  sigma6 )
!          ( sigma5   sigma6  sigma3 )
!
!  .........INNER loop through trial dofs for Cauchy stress
            do k3=1,nrdofQ
              do icomp=1,6
                m3 = (k3-1)*6+icomp
                if (icomp.eq.1) then
                  if (jcomp.eq.1) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(1,k1)*weight
                  endif
                elseif (icomp.eq.2) then
                  if (jcomp.eq.2) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(2,k1)*weight
                  endif
                elseif (icomp.eq.3) then
                  if (jcomp.eq.3) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(3,k1)*weight
                  endif
                elseif (icomp.eq.4) then
                  if (jcomp.eq.1) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(2,k1)*weight
                  elseif (jcomp.eq.2) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(1,k1)*weight
                  endif
                elseif (icomp.eq.5) then
                  if (jcomp.eq.1) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(3,k1)*weight
                  elseif (jcomp.eq.3) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(1,k1)*weight
                  endif
                elseif (icomp.eq.6) then
                  if (jcomp.eq.2) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(3,k1)*weight
                  elseif (jcomp.eq.3) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(2,k1)*weight
                  endif
                endif
              enddo
            enddo
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   0
!
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                   (   L A G R A N G E    M U L T I P L I E R   )
!
!   0
!
!
!        .......END OUTER LOOP
          enddo
        enddo
!           end of integration point loop
      enddo
!
!-----------------------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L                                           |
!-----------------------------------------------------------------------------------
!

!     loop over faces
      do jf = 1,nrf
        mdlf = Nfaces(jf)
        call face_vert_list(mdlf,nrv_f,nverts_f)
        allocate(xverts_f(3,nrv_f))
        do jv=1,Nrv_f
          xverts_f(1:3,jv) = NODES(nverts_f(jv))%coord(1:3,1)
        enddo
!
! get face affine coordinates points b0,b1,b2
        call poly_affine_2d(xverts_f,nrv_f,b0,b1,b2)
        b_aff(:,1) = b0
        b_aff(:,2) = b1
        b_aff(:,3) = b2
! 
!       set up face order in usual structure
        nordf = 0
        nordf(1:4) = NODES(Mdlf)%order
        norie = 0
!  .....set up the face quadrature
        allocate(tloc(2,(nrv_f-2)*MAXNINT2ADD),wt((nrv_f-2)*MAXNINT2ADD))
! get quadrature points for polygon using the affine coordinates
        INTEGRATION = max(0, nord_add_local - 1)
        call set_quadrature_affine_face(mdlf,NODES(Mdlf)%order,nrv_f,xverts_f,b_aff,nint,tloc,wt)
        INTEGRATION = 0
! transform integration points' face affine coordinates to element affine coordinates
        xiloc = 0.d0
        dxi_eta = 0.d0
        call poly_affine_2d_to_3d(tloc,nint,b_aff,a_aff,xiloc,dxi_eta)
! get area of triangle b0b1b2
        call trian_centr_area(b0,b1,b2,x_aux,b_area)
        rjac_f = 2.d0*b_area
! 
        ! write(*,*) 'elem_poly:   mdlf,rjac_f=',mdlf,rjac_f
!       get face normal in physical coordinates
        call face_normal(xverts_f(:,1:3),fn)
! 
!       check if face normal locally goes outward (Norientf = 0). Also correct sign of integral (dir)
        dir = 1.d0
        if(Norientf(jf).ne.0) then
!       if normal fn goes inward, correct it
          ! fn = -1.d0*fn
          dir = -1.d0
        endif
 
        ! faceint = 0.d0
        ! write(*,*)'elem: mdle, jf,fn=',mdle,jf,fn

!    .....loop through face integration points
        do ipt=1,nint
!
!    .......quadrature point in face coordinates
          t(1:2) = tloc(1:2,ipt); wa = wt(ipt); xi(1:3) = xiloc(1:3,ipt)
!           evaluate trial (2d face) and test (3d element) functions
!         Discont Face H1
          call shape2DH(ftype,t,nordf, norie,nrdofH_f,shapH_f,gradH_f)
!         Face L2
          call shape2DQ(ftype,t,nordf, nrdofQ_f,ShapQ_f)
!           H1 (test)
          call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!           H(div) (test)
          call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!       .....Change coordinates so the shape functions are on the physical element
!           L2 (trial)
          shapQ_f(1:nrdofQ_f) = shapQ_f(1:nrdofQ_f)/rjac_f
!
!           H(div) (test)
          do k=1,nrdofVV
            shapVV(1:3,k) = dxdxi_e(1:3,1)*shapVV(1,k)  &
                          + dxdxi_e(1:3,2)*shapVV(2,k)  &
                          + dxdxi_e(1:3,3)*shapVV(3,k)
            shapVV_n(k) = shapVV(1,k)*fn(1)  &
                        + shapVV(2,k)*fn(2)  &
                        + shapVV(3,k)*fn(3)
          enddo

          shapVV_n(1:nrdofVV) = shapVV_n(1:nrdofVV)/rjac_e
!
          weight = wa !*rjac_f
          
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
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   - <\hat u,(\tau n)>
!
!  ...........INNER loop through \hat{u} bdry trial dofs
              do k2=1,nrdofH_f
                do icomp=1,3
                  m2 = (k2-1)*3+icomp + 3*nrdofH_f*(jf-1)

                  if (icomp.eq.jcomp) then
                    ! if (mdlf.eq.50) then
                    !   faceint((k2-1)*3+icomp) = faceint((k2-1)*3+icomp)  - shapQ_f(k2)*weight*dir
                    ! endif
                      EnrTraceDispl(m1,m2) = EnrTraceDispl(m1,m2)  &
                                         - shapH_f(k2)*shapVV_n(k1)*weight*dir
                    

                    ! write(*,*) 'elem    EnrTraceDispl integr:   m1,m2=',m1,m2
                    ! write(*,*) 'elem    EnrTraceDispl integr:   k1,k2=',k1,k2
                    ! write(*,*) 'elem    EnrTraceDispl integr:   jcomp,icomp=',jcomp,icomp
                    ! write(*,*) 'elem    EnrTraceDispl integr:   shapQ_f(k2),shapVV_n(k1)=',shapQ_f(k2),shapVV_n(k1)

                  endif
                enddo
              enddo
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   C A U C H Y    S T R E S S   )
!
!   0
!
!
!  .........END OUTER LOOP
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
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   0
!
!
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   - <\hat \sigma,v>
!
!  ...........INNER loop through H(div) bdry trial dofs
              do k3=1,nrdofQ_f
                do icomp=1,3
                  m3 = (k3-1)*3+icomp + 3*nrdofQ_f*(jf-1)
                  if (icomp.eq.jcomp) then
                    EnrTraceStress(m1,m3) = EnrTraceStress(m1,m3)  &
                                          - shapQ_f(k3)*shapHH(k1)*weight*dir
                  endif
                enddo
              enddo
!
!    .........END OUTER LOOP
            enddo
          enddo
!    .....end of loop over integration points
        enddo

 !        if (mdlf.eq.50) then
 !          write(*,*)'face 50 shape funs integral, mdle ',mdle
 !          write(*,7410) faceint(1:9)
 ! 7410     format(9e12.5)
 !        endif

        deallocate(xverts_f,nverts_f,tloc,wt)
!     end of face loop
      enddo
!
!
        enrdof = 3*nrdofVV+3*nrdofHH

        ! if (mdle.eq.0) then

 ! !        write(*,*) 'Gram       mdle = ',mdle
 ! !        do i=1,18
 ! !          write(*,6000) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 ! ! 6000     format('i = ',i3,'  ',25e12.5)
 ! !        enddo

 !        gdump=75
 !        open(unit=gdump,file='output/gram', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(gdump,5999)Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 ! 5999     format(3003(e22.15,","))
 !        enddo
 !        close(gdump)
 !        bdump1=76
 !        open(unit=bdump1,file='output/b_u', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(bdump1,5998) EnrFieldDispl(i,1:3*nrdofQ)
 ! 5998     format(3(e22.15,","))        
 !        enddo
 !        close(bdump1)
 !        bdump2=77
 !        open(unit=bdump2,file='output/b_sigma', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(bdump2,5997) EnrFieldStress(i,1:6*nrdofQ)
 ! 5997     format(6(e22.15,","))
 !        enddo
 !        close(bdump2)
 !        bdump3=78
 !        open(unit=bdump3,file='output/b_omega', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(bdump3,5998) EnrFieldOmega(i,1:3*nrdofQ) 
 !        enddo
 !        close(bdump3)
 !        bdump4=79
 !        open(unit=bdump4,file='output/b_uhat', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(bdump4,5996) EnrTraceDispl(i,1:3*nrdofH_f*Nrf)
 ! 5996     format(36(e22.15,","))        
 !        enddo
 !        close(bdump4)
 !        bdump5=80
 !        open(unit=bdump5,file='output/b_sigmahat', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(bdump5,5995) EnrTraceStress(i,1:3*nrdofQ_f*Nrf)
 ! 5995     format(12(e22.15,","))  
 !        enddo
 !        close(bdump5)
 !        bdump6=81
 !        open(unit=bdump6,file='output/f', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(bdump6,5994) EnrLoad(i,1)
 ! 5994     format(e22.15,",")  
 !        enddo
 !        close(bdump5)

 !      endif
        ! call pause

 !      if (iprint.eq.1) then
 !        write(*,*) 'Gram = '
 !        do i=1,25
 !          write(*,6000) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 ! ! 6000     format('i = ',i3,'  ',25e12.5)
 !        enddo

 !        call pause

 !        write(*,*) 'EnrFieldDispl = '
 !        do i=1,3*nrdofVV+1
 !          write(*,6001) i,EnrFieldDispl(i,1:3*nrdofQ)
 ! 6001     format('i = ',i4,'  ',15(/,10e12.5))
 !        enddo

 !        call pause

 !        write(*,*) 'EnrFieldStress = '
 !        do i=1,3*nrdofVV
 !          write(*,6001) i,EnrFieldStress(i,1:6*nrdofQ)
 !        enddo

 !        call pause

 !        do i=1+3*nrdofVV,3*nrdofVV+3*nrdofHH
 !          write(*,6001) i,EnrFieldStress(i,1:6*nrdofQ)
 !        enddo

 !        call pause

 !        write(*,*) 'EnrFieldOmega = '
 !        do i=1,3*nrdofVV+1
 !          write(*,6001) i,EnrFieldOmega(i,1:3*nrdofQ)
 !        enddo

 !        call pause

 !        write(*,*) 'EnrTraceDispl = '
 !        do i=1,3*nrdofVV+1
 !          write(*,6001) i,EnrTraceDispl(i,1:3*nrdofH_f*Nrf)
 !        enddo

 !        call pause

 !        write(*,*) 'EnrTraceStress = '
 !        do i=3*nrdofVV,3*nrdofVV+3*nrdofHH
 !          write(*,6001) i,EnrTraceStress(i,1:3*nrdofQ_f*Nrf)
 !        enddo

 !        call pause

 !      endif

!
!-----------------------------------------------------------------------------------
!     D P G    L O C A L    A S S E M B L Y                                        |
!-----------------------------------------------------------------------------------
!
!  ...Compact enriched number of rows (total enriched test dof)
      enrdof = 3*nrdofVV+3*nrdofHH
!
!  ...factor the Gram matrix
      uplo = 'U'


      ! allocate(Gramtmp(enrdof,enrdof))
 !      Gramtmp= ZERO
 !      call DTPTTR(uplo,enrdof,Gram,Gramtmp,enrdof,info)
 !      call DPOTRF(uplo,enrdof,Gramtmp,enrdof,info)

 !      gdump=75
 !        open(unit=gdump,file='output/gram_fact', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(gdump,5999)Gramtmp(i,1:75)
 ! ! 5990     format(75(e12.5,","))
 !        enddo
 !        close(gdump)

 !      call DTRTTP(uplo,enrdof,Gramtmp,enrdof,Gram,info)

! subroutine dspevd ( character   JOBZ,
! character   UPLO,
! integer   N,
! double precision, dimension( * )  AP,
! double precision, dimension( * )  W,
! double precision, dimension( ldz, * )   Z,
! integer   LDZ,
! double precision, dimension( * )  WORK,
! integer   LWORK,
! integer, dimension( * )   IWORK,
! integer   LIWORK,
! integer   INFO 
! ) 
      
      ! Checking eigenvalues
      ! liwork = 2*enrdof
      ! lwork = 2*enrdof
      ! allocate(work(lwork))
      ! allocate(iwork(liwork))
      ! ldz = enrdof
      ! allocate(eval(enrdof),evec(ldz,enrdof))
      ! call DSPEVD('N',uplo,enrdof,Gram(1:enrdof*(enrdof+1)/2), &
      !             Eval,Evec,ldz,Work,lwork,iwork,liwork,info)
      ! write (*,*) 'DSPEVD, info = ',info


      ! Cholesky Factorization
      call DPPTRF(uplo, enrdof, Gram, info)



 !      gdump=75
 !        open(unit=gdump,file='output/gram', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(gdump,5999)Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 ! ! 5999     format(75(e12.5,","))
 !        enddo
 !        close(gdump)
      if (info.ne.0) then
        write(*,*) 'elem_POLYDPG, mdle,enrdof=',mdle,enrdof
        write(*,*) 'elem_POLYDPG Gram Cholesky Factorization: info = ',info
        stop
      endif
!
!  ...Create vector of indices with dof of each physics variable (in the same order of physics file)
      ndofphysics = (/3*nrdofH_f*Nrf,3*nrdofQ_f*Nrf,3*nrdofQ,6*nrdofQ,3*nrdofQ/)
!
!  ...Construct EnrEverything by packing all Enriched Matrices in the order H1,Hcurl,Hdiv,L2
!     (in the same order of physics file) and load
      EnrEverything=ZERO
!  ...EnrTraceDispl (Face L2)
      kmin=0
      kmax=kmin+ndofphysics(1)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrTraceDispl(1:enrdof,1:(kmax-kmin))
!  ...EnrTraceStress (Face L2)
      kmin=kmax
      kmax=kmin+ndofphysics(2)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrTraceStress(1:enrdof,1:(kmax-kmin))
!  ...EnrFieldDispl (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(3)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrFieldDispl(1:enrdof,1:(kmax-kmin))
!  ...EnrFieldStress (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(4)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrFieldStress(1:enrdof,1:(kmax-kmin))
!  ...EnrFieldOmega (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(5)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrFieldOmega(1:enrdof,1:(kmax-kmin))
!  ...EnrLoad
      kmin=kmax
      kmax=kmin+NR_RHS
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrLoad(1:enrdof,1:(kmax-kmin))
!
!  ...Save copy of EnrStiffness which implies deleting the load part from EnrEverything
      EnrStiffness=ZERO
      EnrStiffness(1:enrdof,1:kmin)=EnrEverything(1:enrdof,1:kmin)
!
 !      bdump7=82
 !      open(unit=bdump7,file='output/EnrEverything', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,75
 !          write(bdump7,5995) EnrEverything(i,1:67)
 ! ! 5995     format(66(e12.5,","))        
 !        enddo
 !      close(bdump7)
! !  ...save copies of enriched stiffness matrices
!       EnrTraceDisplc=EnrTraceDispl; EnrTraceStressc=EnrTraceStress
!       EnrFieldDisplc=EnrFieldDispl; EnrFieldStressc=EnrFieldStress
!       EnrFieldOmegac=EnrFieldOmega
!
!  ...G^-1 * EnrEverything
      call DPPTRS(uplo,enrdof,kmax,Gram,EnrEverything(:,1:kmax),3*MAXtetraVV+3*MAXtetraHH,info1)
      if (info1.ne.0) then
        write(*,*) 'elem_DPG_UWEAK: info1 = ',info1 ; stop
      endif
!
!  ...Build full DPG matrix (stiffness + load) in one go
      allocate(FullDPG(kmin,kmax))
      FullDPG = ZERO
      transa = 'T'
      transb = 'N'
      m = kmin
      n = kmax
      k = enrdof
      l = 3*MAXtetraVV+3*MAXtetraHH
      call DGEMM(transa,transb,m,n,k,1.d0,EnrStiffness,l,EnrEverything,l,0.d0,FullDPG,m)!3*(MAXtriaH+MAXtriaQ)*MAX_NRFC+12*MAXtetraQ)
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
!

 !      if (mdle.eq.1) then
 !        bdump6=81      
 !      open(unit=bdump6,file='output/FullDPG_e1', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,210
 !          write(bdump6,5995) FullDPG(i,1:210)
 ! 5995     format(210(e22.15,","))        
 !        enddo
 !        close(bdump6)
 !      else
 !        bdump7=71      
 !      open(unit=bdump7,file='output/FullDPG_e2', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,174
 !          write(bdump7,5885) FullDPG(i,1:174)
 ! 5885     format(174(e22.15,","))        
 !        enddo
 !        close(bdump7)
 !      endif
      
      
      deallocate(Nfaces,Norientf,Nedges,Nverts)
      deallocate(xiloc,wxi,xloc)
      deallocate(FullDPG)
!
! 90 continue
end subroutine elem_DPG_UWEAK_poly