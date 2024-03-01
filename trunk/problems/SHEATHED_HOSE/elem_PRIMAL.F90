!------------------------------------------------------------------------------------------
!> @brief element stiffness matrix and load vector for primal DPG elasticity problem
!> @param[in]  Mdle      - middle node number
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
subroutine elem_DPG_PRIMAL(Mdle)
      use control, only: INTEGRATION
      use primal_module
      use parameters
      use parametersDPG
      use data_structure3D
      use element_data
      use sheathed_isotropic_materials
      use assembly, only: ALOC,BLOC,NR_RHS
      use common_prob_data, only: TEST_NORM

      ! use m_assembly, only:  mdle_list,norderList,nedge_orientList,nface_orientList,xnodList

!------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in) :: Mdle
!------------------------------------------------------------------------------------------
!
!  ...element and face type
      integer :: etype,ftype
!
!  ...element and face order, enriched order
      integer, dimension(19) :: norder
      integer, dimension(5)  :: nordf
      integer                :: nordP
!
!  ...domain number
      integer :: ndom
!
!  ...node edge and face orientations
      integer, dimension(12) :: nedge_orient
      integer, dimension(6)  :: nface_orient
!
!  ...SHAPE FUNCTIONS
!     H1
      real(8), dimension(  MAXbrickH)  :: shapH
      real(8), dimension(3,MAXbrickH)  :: gradH
      integer                         :: nrdofH
!     discontinuous H1
      real(8), dimension(  MAXbrickHH) :: shapHH
      real(8), dimension(3,MAXbrickHH) :: gradHH
      integer                         :: nrdofHH
!     H(div)
      real(8), dimension(3,MAXbrickV)  :: shapV
      real(8), dimension(  MAXbrickV)  :: divV ! never used
      integer                         :: nrdofV
!  ...flux
      real(8), dimension(  MAXbrickV)  :: shapV_n
! !
! !  ...load vector for the enriched space
!       real(8), dimension(3*MAXbrickHH,NRRHS) :: EnrLoad
!
!  ...geometry
      real(8), dimension(3,MAXbrickH) :: xnod
      real(8), dimension(3)           :: xi,x,rn
      real(8), dimension(3,3)         :: dxdxi,dxidx
      real(8), dimension(2)           :: t
      real(8), dimension(3,2)         :: dxidt,dxdt
!
!  ...stiffness tensors in master coordinates and physical coordinates
      real(8), dimension(3,3,3,3) :: C,Symm !,CC
!
!  ...source term (don't need Neumann term)
      real(8), dimension(3,NRRHS) :: fval
!
!  ...3D quadrature data
      real(8), dimension(3,MAXNINT3ADD) :: xiloc
      real(8), dimension(MAXNINT3ADD)   :: wxi
!
!  ...2D quadrature data for boundary terms
      real(8), dimension(2,MAXNINT2ADD) :: tloc
      real(8), dimension(MAXNINT2ADD)   :: wt
!
!  ...miscellaneous
      integer :: i,j,k,l,m,n,k1,k2,k3,m1,m2,m3,n1,n2,ipt,icomp,jcomp,  &
                 nint,ifc,nsign,iprint,iflag,info,info1,   &
                 enrdof,kmin,kmax
      integer :: ndofphysics(2)
      real(8)  :: weight,wa,rjac,brjac,tmp
!
!  ...LAPACK stuff
      character :: uplo,transa,transb
      integer, external :: ij_upper_to_packed
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
      iprint=0
      tmp = 0.d0
!
!  ...element type
      etype = NODES(Mdle)%ntype
!
!  ...order of approximation, orientations, geometry dof's, domain number
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
      call find_domain(Mdle, ndom)
!
!  ...set the enriched order of appoximation
      select case(etype)
      case(MDLB)      ; nordP = NODES(Mdle)%order + NORD_ADD*111
      case(MDLP)      ; nordP = NODES(Mdle)%order + NORD_ADD*11
      case(MDLN,MDLD) ; nordP = NODES(Mdle)%order + NORD_ADD*1
      end select
!
!  ...initialize the enriched local element stiffness matrices and load vectors
      EnrField=ZERO; EnrTrace=ZERO; EnrLoad=ZERO
!  ...initialize the Gram matrix
      Gram=ZERO
!
      call getSymm(Symm)
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!
!  ...set up the element quadrature
      INTEGRATION = NORD_ADD !+3
      call set_3Dint(etype,norder, nint,xiloc,wxi)
      INTEGRATION = 0
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
        if (iflag.ne.0) then
          write(*,1000) Mdle,rjac
 1000     format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
          stop
        endif
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
!  .....compute the stiffness tensor
        call getC(x,ndom, C)
!
!  .....compute the stiffness tensor squared (for gram matrix)
        if (TEST_NORM.eq.3) then
          ! call getCC(x, CC)
        else if (TEST_NORM.eq.5) then
!       H1 (projection)
          ! call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
        endif
!
!  .....get the source term
        call getf(Mdle,x, fval)
!
!  .....OUTER loop through enriched dofs
        do k1=1,nrdofHH
          do icomp=1,3
            m1 = (k1-1)*3+icomp
!
!           E N R I C H E D   L O A D   V E C T O R
!
!   f \cdot v
!
            EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS) &
                                 + fval(icomp,1:NR_RHS)*shapHH(k1)*weight
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
                  k = ij_upper_to_packed(m1,m2)
                  k2 = int((m2-1)/3)+1
                  do m=1,3; do n=1,3
                  tmp = tmp  &
                      + C(jcomp,n,icomp,m)*gradHH(n,k2)*gradHH(m,k1)
                  enddo; enddo
                  Gram(k) = Gram(k)  &
                          + ( tmp + shapHH(k1)*shapHH(k2) )*weight
                endif
              enddo

! MATHEMATICIANS
!   (grad(v_2),grad(v)) + (v_2,v)
            case(2)
              do m2=m1,3*nrdofHH
                jcomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = ij_upper_to_packed(m1,m2)
                  k2 = int((m2-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( shapHH(k1)*shapHH(k2)  &
                            + gradHH(1,k1)*gradHH(1,k2)  &
                            + gradHH(2,k1)*gradHH(2,k2)  &
                            + gradHH(3,k1)*gradHH(3,k2) )*weight
                endif
              enddo
! ! PSEUDO-STRAIN
! !   (C:eps(v_2),C:eps(v)) + (v_2,v)
!               case(3)
!                 do m2=m1,3*nrdofHH
!                   jcomp = mod(m2-1,3)+1
!                   if (icomp.eq.jcomp) then
!                     k = ij_upper_to_packed(m1,m2)
!                     k2 = int((m2-1)/3)+1
!                     do m=1,3; do n=1,3
!                     tmp = tmp  &
!                         + CC(jcomp,n,icomp,m)*gradHH(n,k2)*gradHH(m,k1)
!                     enddo; enddo
!                     Gram(k) = Gram(k)  &
!                             + ( tmp + shapHH(k1)*shapHH(k2) )*weight
!                   endif
!                 enddo
! PSEUDO-STRESS
!   (eps(v_2),eps(v)) + (v_2,v)
            case(4)
                do m2=m1,3*nrdofHH
                  jcomp = mod(m2-1,3)+1
                  if (icomp.eq.jcomp) then
                    k = ij_upper_to_packed(m1,m2)
                    k2 = int((m2-1)/3)+1
                    do m=1,3; do n=1,3
                    tmp = tmp  &
                        + Symm(jcomp,n,icomp,m)*gradHH(n,k2)*gradHH(m,k1)
                    enddo; enddo
                    Gram(k) = Gram(k)  &
                            + ( tmp + shapHH(k1)*shapHH(k2) )*weight
                  endif
                enddo
! STRAIN
!   (P_N(C:eps(v_2)),C:eps(v)) + (v_2,v)
            case(5)
! STRESS
!   (P_M(eps(v_2)),eps(v)) + (v_2,v)
            case(6)

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
                tmp = 0.d0
                do m=1,3; do n=1,3
                  tmp = tmp  &
                      + C(jcomp,n,icomp,m)*gradH(n,k3)*gradHH(m,k1)
                enddo; enddo
                EnrField(m1,m3) = EnrField(m1,m3) + tmp*weight
!  .........INNER loop
              enddo
            enddo
!  .....OUTER loop
          enddo
        enddo
!
!  ...end of loop through integration points
     enddo
!
!-----------------------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L                                           |
!-----------------------------------------------------------------------------------
!
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
        INTEGRATION = NORD_ADD !+3
        call set_2Dint(ftype,nordf, nint,tloc,wt)
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
!
!  .......geometry map
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,  &
                       x,dxdxi,dxidx,rjac,dxdt,rn,brjac)
          weight = wa*brjac
!
!  .......Change coordinates so the shape functions are on the physical element.
!         Also construct trial fluxes
!         H(div) (trial)
          do k=1,nrdofV
            shapV(1:3,k) = dxdxi(1:3,1)*shapV(1,k)  &
                         + dxdxi(1:3,2)*shapV(2,k)  &
                         + dxdxi(1:3,3)*shapV(3,k)
            shapV_n(k) = shapV(1,k)*rn(1)  &
                       + shapV(2,k)*rn(2)  &
                       + shapV(3,k)*rn(3)
          enddo
          shapV_n(1:nrdofV) = shapV_n(1:nrdofV)/rjac
!
!  .......OUTER loop through enriched H1 test functions
          do k1=1,nrdofHH
            do icomp=1,3
              m1 = (k1-1)*3+icomp
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!
!   - <  v , \hat \sigma_n  >
!
!  ...........INNER loop through enriched dofs
              do k2=1,nrdofV
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
!
      if (iprint.eq.1) then

        write(*,*) 'EnrLoad = '
        do i=1,3*nrdofHH+1
          write(*,6000) i,EnrLoad(i,1)
 6000     format('i = ',i4,'  ',e12.5)
        enddo

        call pause

        write(*,*) 'Gram = '
        do i=1,25
          write(*,6001) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 6001     format('i = ',i3,'  ',25e12.5)
        enddo

        call pause

        write(*,*) 'EnrField = '
        do i=1,3*nrdofHH+1
          write(*,6002) i,EnrField(i,1:3*nrdofH)
 6002     format('i = ',i4,'  ',15( / ,10e12.5))
        enddo

        call pause

        write(*,*) 'EnrTrace = '
        do i=1,3*nrdofHH+1
          write(*,6002) i,EnrTrace(i,1:3*nrdofV)
        enddo

      endif
!
!-----------------------------------------------------------------------------------
!     D P G    L O C A L    A S S E M B L Y                                        |
!-----------------------------------------------------------------------------------
!
!  ...Compact enriched number of rows (total enriched test dof)
      enrdof = 3*nrdofHH
!
!  ...factor the Gram matrix
      uplo = 'U'
      call DPPTRF(uplo, enrdof, Gram, info)
      if (info.ne.0) then
        write(*,*) 'elem_DPG_PRIMAL: info = ',info ; stop
      endif
!
!  ...Create vector of indices with dof of each physics variable (in the same order of physics file)
      ndofphysics = (/3*nrdofH,3*nrdofV/)
!
!  ...Construct EnrEverything by packing all Enriched Matrices in the order H1,Hcurl,Hdiv,L2
!     (in the same order of physics file) and load
      EnrEverything=ZERO
!  ...EnrTraceDispl (H1)
      kmin=0
      kmax=kmin+ndofphysics(1)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrField(1:enrdof,1:(kmax-kmin))
!  ...EnrTraceStress (H(div))
      kmin=kmax
      kmax=kmin+ndofphysics(2)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrTrace(1:enrdof,1:(kmax-kmin))
!  ...EnrLoad
      kmin=kmax
      kmax=kmin+NR_RHS
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrLoad(1:enrdof,1:(kmax-kmin))
!
!  ...Save copy of EnrStiffness which implies deleting the load part from EnrEverything
      EnrStiffness=ZERO
      EnrStiffness(1:enrdof,1:kmin)=EnrEverything(1:enrdof,1:kmin)
!
!  ...G^-1 * EnrEverything
      call DPPTRS(uplo,enrdof,kmax,Gram,EnrEverything,3*MAXbrickHH,info1)
      if (info1.ne.0) then
        write(*,*) 'elem_DPG_UWEAK: info1 = ',info1 ; stop
      endif     
!
!  ...Build full DPG matrix (stiffness + load) in one go
      FullDPG = ZERO
      transa = 'T'
      transb = 'N'
      m = kmin
      n = kmax
      k = enrdof
      l = 3*MAXbrickHH
      call DGEMM(transa,transb,m,n,k,1.d0,EnrStiffness,l,EnrEverything,l,0.d0,FullDPG,3*MAXbrickH+3*MAXbrickV)
!
!     ULTIMATE DPG LOAD VECTORS AND STIFFNESS MATRICES THROUGH STATIC CONDENSATION
!
!  ...Populate the assembly matrices accordingly
      n1=0
      do i=1,2
        n=ndofphysics(i)
        n2=n1+n
        m1=0
        do j=1,2
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
end subroutine elem_DPG_PRIMAL
