!------------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for primal DPG elasticity problem
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
subroutine elem_DPG_UWEAK(Mdle)
      use uweak_module
      use control, only: INTEGRATION
      use parameters
      use parametersDPG
      use data_structure3D
      use element_data
      use sheathed_isotropic_materials
      use physics   , only : NR_PHYSA
      use assembly  , only : ALOC,BLOC,NR_RHS
      use common_prob_data, only: SYMMETRY_TOL, TEST_NORM
!------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in) :: Mdle
!------------------------------------------------------------------------------------------
!
!  ...element and face type
      integer :: etype,ftype
!
!  ...number of topological entities (vertices,edges,faces)
      integer :: nrv,nre,nrf
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
!  ...geometry
      real*8, dimension(3,MAXbrickH) :: xnod
      real*8, dimension(3)           :: xi,x,rn
      real*8, dimension(3,3)         :: dxdxi,dxidx
      real*8, dimension(2)           :: t
      real*8, dimension(3,2)         :: dxidt,dxdt
      integer                        :: nsign
!
!  ...tensors in physical coordinates
      real*8, dimension(3,3,3,3) :: A,AA,Symm,Skew
!
!  ...source term (don't need Neumann term)
      real*8, dimension(3,NRRHS) :: fval
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
      integer :: i,j,k,l,m,n,k1,k2,k3,k4,k5,m1,m2,m3,m4,m5,n1,n2,ipt,ifc,  &
                 icomp,jcomp,nint,iprint,iflag,info,info1,  &
                 kH,lH,kmin,kmax,enrdof,weightedWeakSymmetryConstraint
      integer, dimension(NR_PHYSA) :: ndofphysics
      real*8  :: weight,wa,rjac,brjac,tmp,diffmax,dmax,omegaWeight,l2Weight,l2StressWeight
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
      weightedWeakSymmetryConstraint=1

      l2Weight = 1.d-3
      l2StressWeight = 0.d0
!
!  ...element type
      etype = NODES(Mdle)%ntype
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
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
      EnrTraceDispl=ZERO; EnrTraceStress=ZERO; EnrFieldDispl=ZERO; EnrFieldStress=ZERO; EnrFieldOmega=ZERO
      EnrLoad=ZERO
!  ...initialize the Gram matrix
      Gram=ZERO
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
!  ...Tensors for adjoint graph norm
      if (TEST_NORM.eq.1) then
        call getSymm(Symm)
        call getSkew(Skew)
      endif
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
!  .....integration weight
        weight = wa*rjac
!
!  .....compute the compliance tensor
        call getA(x,ndom, A)
!
!  .....One more tensor for adjoint graph norm
        if (TEST_NORM.eq.1) then
          call getAA(x, AA)
          ! call getWeightedSkew(X, Skew)
        endif
!
!  .....get the source term
        call getf(Mdle,x, fval)
!
!  .....Weight for omega term
        if (weightedWeakSymmetryConstraint.eq.1) then
          if (dsqrt(xi(2)**2+xi(3)**3).le.R_middle) then
            omegaWeight = 1.d0/E_rubber
          else
            omegaWeight = 1.d0/E_steel
          endif
        else
          omegaWeight = 1.d0
        endif
!
!  NOTE: Because nrddofVV \neq nrdofQ, there is no convenient way of organizing
!        the matrices testing all of the n-th basis functions in turn before
!        moving on to the (n+1)-th basis functions. Therefore, and which will lead
!        to a sparser system as well, we test all the basis functions representing
!        \tau and then all those representing v.
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
!   0
!
!           G R A M   M A T R I X
!
            select case(TEST_NORM)
!
!   (A:\tau_1,A:\tau_2+\nabla(v_2))_symm+(\tau_1,\tau_2)_skew+(div(\tau_1),div(tau_2))
!
            case(1)
              do m2=m1,3*nrdofVV
                icomp = mod(m2-1,3)+1
                k = nk(m1,m2)
                k2 = int((m2-1)/3)+1

                tmp=0.d0
                do m=1,3; do n=1,3
                  tmp = tmp  &
                      ! + Symm(icomp,m,jcomp,n)*shapVV(m,k2)*shapVV(n,k1)  &
                      + AA(icomp,m,jcomp,n)*shapVV(m,k2)*shapVV(n,k1)  &
                      + Skew(icomp,m,jcomp,n)*shapVV(m,k2)*shapVV(n,k1)*omegaWeight**2
                enddo; enddo

                !  div term
                if (icomp.eq.jcomp) then
                  tmp = tmp + divVV(k1)*divVV(k2)
                endif

                ! additional L2 term
                tmp = tmp + ( shapVV(1,k1)*shapVV(1,k2)  &
                            + shapVV(2,k1)*shapVV(2,k2)  &
                            + shapVV(3,k1)*shapVV(3,k2) )*l2StressWeight

                Gram(k) = Gram(k) + tmp*weight
              enddo

              do m2=3*nrdofVV+1,3*nrdofVV+3*nrdofHH
                k = nk(m1,m2)
                k2 = int((m2-3*nrdofVV-1)/3)+1
                icomp = mod(m2-1,3)+1

                tmp = 0.d0
                do m=1,3; do n=1,3
                  tmp = tmp + A(icomp,m,jcomp,n)*gradHH(m,k2)*shapVV(n,k1)
                enddo; enddo

                Gram(k) = Gram(k) + tmp*weight
              enddo
!
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
                                         + shapQ(k5)*shapVV(3,k1)*weight*omegaWeight
                  elseif (jcomp.eq.3) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         - shapQ(k5)*shapVV(2,k1)*weight*omegaWeight
                  endif
                elseif (icomp.eq.2) then
                  if (jcomp.eq.1) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         - shapQ(k5)*shapVV(3,k1)*weight*omegaWeight
                  elseif (jcomp.eq.3) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         + shapQ(k5)*shapVV(1,k1)*weight*omegaWeight
                  endif
                elseif (icomp.eq.3) then
                  if (jcomp.eq.1) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         + shapQ(k5)*shapVV(2,k1)*weight*omegaWeight
                  elseif (jcomp.eq.2) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         - shapQ(k5)*shapVV(1,k1)*weight*omegaWeight
                  endif
                endif
              enddo
            enddo
!  .....END OUTER LOOP through test stresses
          enddo
        enddo
!
!
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows)
!
!  .....SECOND OUTER loop through enriched H1 dofs
        do k1=1,nrdofHH
!  .......OUTER loop through components
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
!   (\nabla(v_1),\nabla(v_2))_symm + (v_1,v_2)
!
            case(1)
              do m2=m1,3*nrdofVV+3*nrdofHH
                icomp = mod(m2-1,3)+1
                k = nk(m1,m2)
                k2 = int((m2-3*nrdofVV-1)/3)+1

                tmp=0.d0
                do m=1,3; do n=1,3
                  tmp = tmp  &
                      + Symm(icomp,m,jcomp,n)*gradHH(m,k2)*gradHH(n,k1)
                enddo; enddo
                ! L2-term
                if (icomp.eq.jcomp) then
                  tmp = tmp + shapHH(k1)*shapHH(k2)*l2Weight
                endif

                Gram(k) = Gram(k) + tmp*weight
              enddo
!
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
!  .......END OUTER LOOP
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
!  ...loop through element faces
      do ifc=1,nrf
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
!         H1 (geometry and bdry trial)
          call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                        nrdofH,shapH,gradH)
!         H(div) (bdry trial)
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
!  .......Change coordinates (only of what is needed in this case) so the shape functions are on the physical element
!         H(div) (trial)
          do k=1,nrdofV
            shapV(1:3,k) = dxdxi(1:3,1)*shapV(1,k)  &
                         + dxdxi(1:3,2)*shapV(2,k)  &
                         + dxdxi(1:3,3)*shapV(3,k)
            shapV_n(k)  = shapV(1,k)*rn(1)  &
                        + shapV(2,k)*rn(2)  &
                        + shapV(3,k)*rn(3)
          enddo
          shapV_n(1:nrdofV) = shapV_n(1:nrdofV)/rjac
!
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
!  ...........INNER loop through H1 bdry trial dofs
              do k2=1,nrdofH
                do icomp=1,3
                  m2 = (k2-1)*3+icomp
                  if (icomp.eq.jcomp) then
                    EnrTraceDispl(m1,m2) = EnrTraceDispl(m1,m2)  &
                                         - shapH(k2)*shapVV_n(k1)*weight
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
              do k3=1,nrdofV
                do icomp=1,3
                  m3 = (k3-1)*3+icomp
                  if (icomp.eq.jcomp) then
                    EnrTraceStress(m1,m3) = EnrTraceStress(m1,m3)  &
                                          - shapV_n(k3)*shapHH(k1)*weight
                  endif
                enddo
              enddo
!
!  .........END OUTER LOOP
            enddo
          enddo
!  .....end of loop over integration points
        enddo
!
!  ...end of loop over faces
      enddo
!
      if (iprint.eq.1) then
        write(*,*) 'Gram = '
        do i=1,25
          write(*,6000) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 6000     format('i = ',i3,'  ',25e12.5)
        enddo

        call pause

        write(*,*) 'EnrFieldDispl = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrFieldDispl(i,1:3*nrdofQ)
 6001     format('i = ',i4,'  ',15(/,10e12.5))
        enddo

        call pause

        write(*,*) 'EnrFieldStress = '
        do i=1,3*nrdofVV
          write(*,6001) i,EnrFieldStress(i,1:6*nrdofQ)
        enddo

        call pause

        do i=1+3*nrdofVV,3*nrdofVV+3*nrdofHH
          write(*,6001) i,EnrFieldStress(i,1:6*nrdofQ)
        enddo

        call pause

        write(*,*) 'EnrFieldOmega = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrFieldOmega(i,1:3*nrdofQ)
        enddo

        call pause

        write(*,*) 'EnrTraceDispl = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrTraceDispl(i,1:3*nrdofH)
        enddo

        call pause

        write(*,*) 'EnrTraceStress = '
        do i=3*nrdofVV,3*nrdofVV+3*nrdofHH
          write(*,6001) i,EnrTraceStress(i,1:3*nrdofV)
        enddo

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
!  ...factor the Gram matrix
      uplo = 'U'
      call DPPTRF(uplo, enrdof, Gram, info)
      if (info.ne.0) then
        write(*,*) 'elem_DPG_UWEAK: info = ',info ; stop
      endif
!
!  ...Create vector of indices with dof of each physics variable (in the same order of physics file)
      ndofphysics = (/3*nrdofH,3*nrdofV,3*nrdofQ,6*nrdofQ,3*nrdofQ/)
!
!  ...Construct EnrEverything by packing all Enriched Matrices in the order H1,Hcurl,Hdiv,L2
!     (in the same order of physics file) and load
      EnrEverything=ZERO
!  ...EnrTraceDispl (H1)
      kmin=0
      kmax=kmin+ndofphysics(1)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrTraceDispl(1:enrdof,1:(kmax-kmin))
!  ...EnrTraceStress (H(div))
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
!  ...G^-1 * EnrEverything
      call DPPTRS(uplo,enrdof,kmax,Gram,EnrEverything,3*MAXbrickVV+3*MAXbrickHH,info1)
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
      l = 3*MAXbrickVV+3*MAXbrickHH
      call DGEMM(transa,transb,m,n,k,1.d0,EnrStiffness,l,EnrEverything,l,0.d0,FullDPG,3*MAXbrickH+3*MAXbrickV+12*MAXbrickQ)
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
!-----------------------------------------------------------------------------------
!     T E S T S    A N D    P R I N T    S T A T E M E N T S                       |
!-----------------------------------------------------------------------------------
!
!  ...check symmetry (not complete)
      diffmax = ZERO; dmax = ZERO
      do kH=1,3*nrdofH
        do lH=kH,3*nrdofH
          diffmax = max(diffmax,abs(ALOC(1,1)%array(kH,lH)-ALOC(1,1)%array(lH,kH)))
          dmax = max(dmax,abs(ALOC(1,1)%array(kH,lH)))
        enddo
      enddo
      if (diffmax/dmax.gt.SYMMETRY_TOL) then
        write(*,7021) diffmax, dmax
 7021   format('elem_DPG_UWEAK_SYMM: diffmax,dmax FOR ALOC(1,1)%array = ',2e12.5)
        call pause
      endif
!
end subroutine elem_DPG_UWEAK
