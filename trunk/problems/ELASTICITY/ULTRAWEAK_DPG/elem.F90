!--------------------------------------------------------------------------
!> Purpose : return stiffness matrix and load vector for element
!!
!! @param[in]  Mdle   - an element (middle node) number
!! @param[out] Nrdof  - number of dof for a single component
!! @param[out] Itest  - index for assembly
!! @param[out] Itrial - index for assembly
!--------------------------------------------------------------------------
   subroutine elem(Mdle, Itest,Itrial)
!
      use physics   , only : NR_PHYSA
      use assembly  , only : ALOC,BLOC,NR_RHS
      use data_structure3D
!
      implicit none
!
      integer, intent(in)  :: Mdle
      integer, intent(out) :: Itest(NR_PHYSA), Itrial(NR_PHYSA)
!
!--------------------------------------------------------------------------
!
      Itest (1:NR_PHYSA) = 0
      Itrial(1:NR_PHYSA) = 0
!
      select case(NODES(Mdle)%case)
      case(15)
         Itest(1:NR_PHYSA)=1; Itrial(1:NR_PHYSA)=1
         call elem_DPG_UWEAK_SYMM(Mdle)
      case default
         write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',  &
                     Mdle,NODES(Mdle)%case
         call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
      end select
!
   end subroutine elem




!--------------------------------------------------------------------------
!> Purpose :   Element stiffness matrix and load vector for UW DPG
!!             elasticity problem
!!
!! @param[in]  Mdle      - middle node number
!!
!> @date       July 2023
!--------------------------------------------------------------------------
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
!--------------------------------------------------------------------------
   subroutine elem_DPG_UWEAK_SYMM(Mdle)
!
      use control,                  only: INTEGRATION
      use parameters
      use parametersDPG
      use data_structure3D
      use element_data
      use isotropic_elast_material
      use physics,                  only: NR_PHYSA
      use assembly,                 only: ALOC, BLOC, NR_RHS
      use common_prob_data,         only: SYMMETRY_TOL, TEST_NORM, ALPHA
!
      implicit none
!
      integer, intent(in)  :: Mdle
!
!  ...element and face data
      integer :: ntype,ftype
      integer :: nrv,nre,nrf
      integer :: norder(19), nordf(5), nordP
      integer :: nedge_orient(12),nface_orient(6)
!
!  ...shape functions
      real(8) :: shapH  (  MAXbrickH), gradH(3,MAXbrickH)
      real(8) :: shapV  (3,MAXbrickV), divV (  MAXbrickV)
      real(8) :: shapV_n(  MAXbrickV)
      real(8) :: shapQ  (  MAXbrickQ)
      real(8) :: shapHH (  MAXbrickHH), gradHH(3,MAXbrickHH)
!
      integer :: nrdofH, nrdofV, nrdofQ, nrdofHH
!
!  ...geometry
      real(8) :: xnod(3,MAXbrickH)
      real(8) :: xi(3), x(3), rn(3)
      real(8) :: dxdxi(3,3), dxidx(3,3)
      real(8) :: t(2), rjac, bjac
      real(8) :: dxidt(3,2), dxdt(3,2), rt(3,2)
      integer :: nsign
!
!  ...tensors in physical coordinates
      real(8), dimension(3,3,3,3) :: A,AA,Symm
!
!  ...temporary variables
      real(8) :: tmp
      real(8) :: tmpDivTau1(3), tmpDivTau2(3), tmpTau_n(3)
!
!  ...source term (don't need Neumann term)
      real(8) :: fval(3,NR_RHS)
!
!  ...quadrature data
      real(8) :: xiloc(3,MAXNINT3ADD), wxi(MAXNINT3ADD)
      real(8) :: tloc(2,MAXNINT2ADD),  wt(MAXNINT2ADD)
      real(8) :: weight, wa
      integer :: nint
!
!  ...Matrices
!     stiffness and load matrices for the enriched test space
      real(8) :: EnrTraceDispl (9*MAXbrickHH,3*MAXbrickH)  !EnrTraceDisplc
      real(8) :: EnrTraceStress(9*MAXbrickHH,3*MAXbrickV)  !EnrTraceStressc
      real(8) :: EnrFieldDispl (9*MAXbrickHH,3*MAXbrickQ)  !EnrFieldDisplc
      real(8) :: EnrFieldStress(9*MAXbrickHH,6*MAXbrickQ)  !EnrFieldStressc
      real(8) :: EnrLoad       (9*MAXbrickHH,NR_RHS)
      real(8) :: Stiff_All(9*MAXbrickHH,                                &
                           3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ+NR_RHS)
      real(8) :: FullDPG(3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ,           &
                         3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ+NR_RHS)
!  ...Gram matrix
      real(8) :: Gram((9*MAXbrickHH)*(9*MAXbrickHH+1)/2)
!
!  ...miscellaneous
      integer :: i, j, k, l, m, n, mm, nn
      integer :: k1, k2, k3, k4, k5
      integer :: m1, m2, m3, m4, m5
      integer :: n1, n2, n3
      integer :: ic, jc, kc, lc, klc, ijc
      integer :: ipt, ifc, iflag, info
      integer :: kH, kV, kQ, lH, lV, lQ
      integer :: kmin, kmax, enrdof
      integer :: ndofphysics(NR_PHYSA)
      real(8) :: diffmax,dmax
!
      integer :: iprint = 0
!
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1
!
!--------------------------------------------------------------------------
!
!  ...element type
      ntype = NODES(Mdle)%ntype
      nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
!  ...order of approximation, orientations, geometry dof's (don't need bc flags)
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
!
!  ...set the enriched order of appoximation
      select case(ntype)
      case(MDLB);       nordP = NODES(Mdle)%order + NORD_ADD*111
      case(MDLP);       nordP = NODES(Mdle)%order + NORD_ADD*11
      case(MDLN,MDLD);  nordP = NODES(Mdle)%order + NORD_ADD*1
      end select
!
!  ...initialize the enriched local element stiffness matrices and load vectors
      EnrTraceDispl(:,:)  = ZERO
      EnrTraceStress(:,:) = ZERO
      EnrFieldDispl(:,:)  = ZERO
      EnrFieldStress(:,:) = ZERO
      EnrLoad(:) = ZERO
!
!  ...initialize the Gram matrix
      Gram(:) = ZERO
!
!--------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L
!--------------------------------------------------------------------------
!
!  ...set up the element quadrature
      INTEGRATION = NORD_ADD
      call set_3Dint(ntype,norder, nint,xiloc,wxi)
      INTEGRATION = 0
!
!  ...Auxiliary tensors
      call getSymm(Symm)
!
!  ...loop through integration points
      do ipt=1,nint
         xi(1:3) = xiloc(1:3,ipt)
         wa = wxi(ipt)
!
!     ...geometry shape functions
         call shape3DH(ntype,xi,norder,nedge_orient,nface_orient,  &
                      nrdofH,shapH,gradH)
!
!     ...L2 shape functions
         call shape3DQ(ntype,xi,norder, nrdofQ,shapQ)
!
!     ...enriched H1 shape functions
         call shape3HH(ntype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!     ...geometry map
         call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
                     x,dxdxi,dxidx,rjac,iflag)
!
         if (iflag.ne.0) then
            write(*,1000) Mdle,rjac
 1000       format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
            stop
         endif
!
!     ...Apply pullbacks
!        L2 (trial)
         shapQ(1:nrdofQ) = shapQ(1:nrdofQ)/rjac
!
! TODO: call Lapack routine
!     ...H1 (test)
         do k=1,nrdofHH
            gradHH(1:3,k) = gradHH(1,k) * dxidx(1,1:3)  &
                          + gradHH(2,k) * dxidx(2,1:3)  &
                          + gradHH(3,k) * dxidx(3,1:3)
         enddo
!
!     ...integration weight
         weight = wa*rjac
!
!     ...compute the compliance tensor
         call getA(x, A)
!
!     ...need this for the adjoint graph norm
         call getAA(X, AA)
!
!     ...get the source term
         call getf(Mdle,x, fval)
!
!    P A R T  1 : go through \tau\in H(div)^3 test space (this fills the first set of rows)
!
!     ...FIRST OUTER loop through enriched stress dofs
         do k1=1,nrdofHH
!        ...OUTER loop through components
            do ic=1,3
               do jc=ic,3
                  ijc = nk(ic,jc)
                  m1 = (k1-1)*6+ijc

                  tmpDivTau1=0.d0
                  tmpDivTau1(ic) = gradHH(jc,k1)
                  tmpDivTau1(jc) = gradHH(ic,k1)
!
!     --- Gram matrix ---
!
                  select case(TEST_NORM)
!
!   (A:\tau_1,A:\tau_2+\varepsilon(v_2))+(div(\tau_1),div(tau_2))
!
                  case(1)
!
!                    Tau loop
!
                     do m2=m1,6*nrdofHH
                        k = nk(m1,m2)
                        k2 = int((m2-1)/6)+1
                        klc = mod(m2-1,6)+1
!
                        select case(klc)
                        case(1); kc=1; lc=1;
                        case(2); kc=1; lc=2;
                        case(3); kc=2; lc=2;
                        case(4); kc=1; lc=3;
                        case(5); kc=2; lc=3;
                        case(6); kc=3; lc=3;
                        end select
!
                        tmpDivTau2=0.d0
                        tmpDivTau2(kc) = gradHH(lc,k2)
                        tmpDivTau2(lc) = gradHH(kc,k2)
!
                        tmp=0.d0
!
!                    ...L2-terms
                        if (ic.ne.jc) then
                           if (kc.ne.lc) then
                              tmp = AA(kc,lc,ic,jc) * shapHH(k2) * shapHH(k1)    &
                                  + AA(kc,lc,jc,ic) * shapHH(k2) * shapHH(k1)    &
                                  + AA(lc,kc,ic,jc) * shapHH(k2) * shapHH(k1)    &
                                  + AA(lc,kc,jc,ic) * shapHH(k2) * shapHH(k1)
                           else
                              tmp = AA(kc,lc,ic,jc) * shapHH(k2) * shapHH(k1)    &
                                  + AA(kc,lc,jc,ic) * shapHH(k2) * shapHH(k1)
                           endif
                        elseif (kc.ne.lc) then
                           tmp = AA(kc,lc,ic,jc) * shapHH(k2) * shapHH(k1)       &
                               + AA(lc,kc,ic,jc) * shapHH(k2) * shapHH(k1)
                        else
                           tmp = AA(kc,lc,ic,jc) * shapHH(k2) * shapHH(k1)
                        endif
!
!                    ...div term
                        do m=1,3
                           tmp = tmp + tmpDivTau1(m) * tmpDivTau2(m)
                        enddo
!
!                    ...additional L2-terms
                        if (ic.eq.jc) then
                           if ((kc.eq.ic).and.(lc.eq.jc)) then
                              tmp = tmp + ALPHA * shapHH(k1) * shapHH(k2)
                           endif
                        elseif ((kc.eq.ic).and.(lc.eq.jc)) then
                           tmp = tmp + 2 * ALPHA * shapHH(k1) * shapHH(k2)
                        endif
!
                        Gram(k) = Gram(k) + tmp*weight
                     enddo
!
!                    v-loop
!
                     do m2=6*nrdofHH+1,6*nrdofHH+3*nrdofHH
                        k = nk(m1,m2)
                        k2 = int((m2-6*nrdofHH-1)/3)+1
                        kc = mod(m2-6*nrdofHH-1,3)+1

                        tmp=0.d0
                        if (ic.eq.jc) then
                           do m=1,3
                              tmp = tmp  &
                                  + A(ic,jc,kc,m) * gradHH(m,k2) * shapHH(k1)
                           enddo
                        else
                           do m=1,3
                              tmp = tmp  &
                                  + A(ic,jc,kc,m) * gradHH(m,k2) * shapHH(k1)  &
                                  + A(jc,ic,kc,m) * gradHH(m,k2) * shapHH(k1)
                           enddo
                        endif
!
                        Gram(k) = Gram(k) + tmp*weight
                     enddo
!
!   (\tau_1,\tau_2)+(div(\tau_1),div(tau_2))
!
                  case(2)
                     do m2=m1,6*nrdofHH
                        k = nk(m1,m2)
                        k2 = int((m2-1)/6)+1
                        klc = mod(m2-1,6)+1
!
                        select case(klc)
                        case(1); kc=1; lc=1;
                        case(2); kc=1; lc=2;
                        case(3); kc=2; lc=2;
                        case(4); kc=1; lc=3;
                        case(5); kc=2; lc=3;
                        case(6); kc=3; lc=3;
                        end select
!
                        tmpDivTau2=0.d0
                        tmpDivTau2(kc) = gradHH(lc,k2)
                        tmpDivTau2(lc) = gradHH(kc,k2)
!
                        tmp=0.d0
!
!                    ...L2 terms
                        if (ic.eq.jc) then
                           if ((kc.eq.ic).and.(lc.eq.jc)) then
                              tmp = shapHH(k1) * shapHH(k2)
                           endif
                        elseif ((kc.eq.ic).and.(lc.eq.jc)) then
                           tmp = 2*shapHH(k1) * shapHH(k2)
                        endif
!
!                    ...Div term
                        do m=1,3
                           tmp = tmp + tmpDivTau1(m) * tmpDivTau2(m)
                        enddo
!
                        Gram(k) = Gram(k) + tmp*weight
                     enddo
                  end select
!
!     --- cauchy stress stiffness matrix ---
!
!   A \sigma : \tau
!
!          ( sigma1  sigma2  sigma4 )
!  sigma = ( sigma2  sigma3  sigma5 )
!          ( sigma4  sigma5  sigma6 )
!
!              ...INNER loop through trial dofs for Cauchy stress
                  do k3=1,nrdofQ
                     do klc=1,6
                        m3 = (k3-1)*6+klc
!
                        if (ic.ne.jc) then
                           if (klc.eq.1) then
                              tmp = A(1,1,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(1,1,jc,ic) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.2) then
                              tmp = A(1,2,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(1,2,jc,ic) * shapQ(k3) * shapHH(k1)  &
                                  + A(2,1,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(2,1,jc,ic) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.3) then
                              tmp = A(2,2,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(2,2,jc,ic) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.4) then
                              tmp = A(1,3,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(1,3,jc,ic) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,1,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,1,jc,ic) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.5) then
                              tmp = A(2,3,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(2,3,jc,ic) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,2,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,2,jc,ic) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.6) then
                              tmp = A(3,3,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,3,jc,ic) * shapQ(k3) * shapHH(k1)
                           endif
                        else
                           if (klc.eq.1) then
                              tmp = A(1,1,ic,jc) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.2) then
                              tmp = A(1,2,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(2,1,ic,jc) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.3) then
                              tmp = A(2,2,ic,jc) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.4) then
                              tmp = A(1,3,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,1,ic,jc) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.5) then
                              tmp = A(2,3,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,2,ic,jc) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.6) then
                              tmp = A(3,3,ic,jc) * shapQ(k3) * shapHH(k1)
                           endif
                        endif

                        EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                              + tmp * weight
                     enddo
                  enddo
!
!     --- displacement stiffness matrix --
!
!   u \cdot div(\tau)
!
!              ...INNER loop through trial dofs for displacement
                  do k4=1,nrdofQ
                     do kc=1,3
                        m4 = (k4-1)*3 + kc
!
                        EnrFieldDispl(m1,m4) = EnrFieldDispl(m1,m4)  &
                                             + shapQ(k4) * tmpDivTau1(kc) * weight
                     enddo
                  enddo
!           ...END OUTER LOOP through test stresses
               enddo
            enddo
         enddo
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows)
!
!     ...SECOND OUTER loop through enriched H1 dofs
         do k1=1,nrdofHH
!        ...OUTER loop through components
            do ic=1,3
!
!           ...counter of row
               m1 = 6*nrdofHH+(k1-1)*3+ic
!
!     --- load vector ---
!
!           ...f \cdot v
               EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS) &
                                    + fval(ic,1:NR_RHS) * shapHH(k1) * weight
!
!     --- Gram matrix ---
!
               select case(TEST_NORM)
!
!   (\varepsilon(v_1),A:\tau_2+\varepsilon(v_2))+(v_1,v_2)
!
               case(1)
                  do m2=m1,6*nrdofHH+3*nrdofHH
                     k = nk(m1,m2)
                     k2 = int((m2-6*nrdofHH-1)/3)+1
                     jc = mod(m2-1,3)+1

                     tmp=0.d0
!                 ...grad terms
                     do m=1,3
                        do n=1,3
                           tmp = tmp  &
                               + Symm(jc,m,ic,n) * gradHH(m,k2) * gradHH(n,k1)
                        enddo
                     enddo
!                 ...L2-term
                     if (ic.eq.jc) then
                        tmp = tmp + shapHH(k1)*shapHH(k2)
                     endif
!
                     Gram(k) = Gram(k) + tmp*weight
                  enddo
!
!   (v_2,v) + (grad(v_2),grad(v))
!
               case(2)
                  do m2=m1,6*nrdofHH+3*nrdofHH
                     jc = mod(m2-1,3)+1
                     if (ic.eq.jc) then
                        k = nk(m1,m2)
                        k2 = int((m2-6*nrdofHH-1)/3)+1
                        Gram(k) = Gram(k)  &
                                + ( shapHH(k1) * shapHH(k2)                &
                                  + gradHH(1,k1) * gradHH(1,k2)            &
                                  + gradHH(2,k1) * gradHH(2,k2)            &
                                  + gradHH(3,k1) * gradHH(3,k2) ) * weight
                     endif
                  enddo
               end select
!
!     --- cauchy stress stiffness matrix ---
!
!   + \sigma : \varepsilon(v)
!
!          ( sigma1  sigma2  sigma4 )
!  sigma = ( sigma2  sigma3  sigma5 )
!          ( sigma4  sigma5  sigma3 )
!
!           ...INNER loop through trial dofs for Cauchy stress
               do k3=1,nrdofQ
                  do klc=1,6
                     m3 = (k3-1)*6 + klc
!
                     select case(klc)
                     case(1); kc=1; lc=1;
                     case(2); kc=1; lc=2;
                     case(3); kc=2; lc=2;
                     case(4); kc=1; lc=3;
                     case(5); kc=2; lc=3;
                     case(6); kc=3; lc=3;
                     end select
!
                     tmp=0.d0
                     if (kc.eq.ic) then
                        tmp = gradHH(lc,k1) * shapQ(k3)
                     elseif (lc.eq.ic) then
                        tmp = gradHH(kc,k1) * shapQ(k3)
                     endif

                     EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                           + tmp * weight
                  enddo
               enddo
!
!        ...end outer loop
            enddo
         enddo
!
!  ...end of loop through integration points
      enddo
!
!--------------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L
!--------------------------------------------------------------------------
!
!  ...loop through element faces
      do ifc=1,nrf
!
!     ...sign factor to determine the OUTWARD normal unit vector
         nsign = nsign_param(ntype,ifc)
!
!     ...face type
         ftype = face_type(ntype,ifc)
!
!     ...face order of approximation
         call face_order(ntype,ifc,norder, nordf)
!
!     ...set up the face quadrature
         INTEGRATION = NORD_ADD
         call set_2Dint(ftype,nordf, nint,tloc,wt)
         INTEGRATION = 0
!
!     ...loop through face integration points
         do ipt=1,nint
!
!        ...face coordinates
            t(1:2) = tloc(1:2,ipt); wa = wt(ipt)
!
!        ...master element coordinates using face parameterization
            call face_param(ntype,ifc,t, xi,dxidt)
!
!        ...Compute shape functions needed for test/trial field variables and geometry
!           H1 (geometry and bdry trial)
            call shape3DH(ntype,xi,norder,nedge_orient,nface_orient,  &
                         nrdofH,shapH,gradH)
!        ...H(div) (bdry trial)
            call shape3DV(ntype,xi,norder,nface_orient,  &
                         nrdofV,shapV,divV)
!        ...H1 (test)
            call shape3HH(ntype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!        ...geometry map
            call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,  &
                         x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
            weight = wa*bjac
!
!        ...apply pullbacks
! TODO: use Lapack
!           H(div) (trial)
            do k=1,nrdofVi
               shapV(1:3,k) = dxdxi(1:3,1)*shapV(1,k)  &
                            + dxdxi(1:3,2)*shapV(2,k)  &
                            + dxdxi(1:3,3)*shapV(3,k)
               shapV_n(k) = shapV(1,k)*rn(1)  &
                          + shapV(2,k)*rn(2)  &
                          + shapV(3,k)*rn(3)
            enddo
            shapV_n(1:nrdofV) = shapV_n(1:nrdofV)/rjac
!
!
!    P A R T  1 : go through \tau test space
!
!
!        ...OUTER loop through enriched H(div) test functions
            do k1=1,nrdofHH
               do ic=1,3
                  do jc=ic,3
                     ijc = nk(ic,jc)
                     m1 = (k1-1)*6+ijc

                     tmpTau_n(1:3)=0.d0
                     tmpTau_n(ic) = shapHH(k1)*rn(jc)
                     tmpTau_n(jc) = shapHH(k1)*rn(ic)
!
!        --- displacement stiffness matrix ---
!
!   - <\hat u,(\tau n)>
!
!                 ...INNER loop through H1 bdry trial dofs
                     do k2=1,nrdofHi
                        do kc=1,3
                           m2 = (k2-1)*3+kc

                           EnrTraceDispl(m1,m2) = EnrTraceDispl(m1,m2)  &
                                                - shapH(k2)*tmpTau_n(kc)*weight
                        enddo
                     enddo
!
!              ...END OUTER LOOP
                  enddo
               enddo
            enddo
!
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows)
!
!        ...SECOND OUTER loop through enriched H1 test function
            do k1=1,nrdofHH
!           ...OUTER loop through components
               do jc=1,3
                  m1 = 6*nrdofHH+(k1-1)*3+jc
!
!        --- cauchy stress stiffness matrix ---
!
!   - <\hat \sigma,v>
!
!              ...INNER loop through H(div) bdry trial dofs
                  do k3=1,nrdofVi
                     do ic=1,3
                        m3 = (k3-1)*3+ic
                        if (ic.eq.jc) then
                           EnrTraceStress(m1,m3) = EnrTraceStress(m1,m3)  &
                                                 - shapV_n(k3)*shapHH(k1)*weight
                        endif
                     enddo
                  enddo
!
!           ...END OUTER LOOP
               enddo
            enddo
!     ...end of loop over integration points
         enddo
!
!  ...end of loop over faces
      enddo
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,*) 'Gram = '
         do i=1,25
            write(*,6000) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 6000       format('i = ',i3,'  ',25e12.5)
         enddo
!
         call pause
!
         write(*,*) 'EnrFieldDispl = '
         do i=1,6*nrdofHH+1
            write(*,6001) i,EnrFieldDispl(i,1:3*nrdofQ)
 6001       format('i = ',i4,'  ',15(/,10e12.5))
         enddo
!
         call pause
!
         write(*,*) 'EnrFieldStress = '
         do i=1,6*nrdofHH
           write(*,6001) i,EnrFieldStress(i,1:6*nrdofQ)
         enddo
!
         call pause
!
         do i=1+6*nrdofHH,6*nrdofHH+3*nrdofHH
            write(*,6001) i,EnrFieldStress(i,1:6*nrdofQ)
         enddo
!
         call pause
!
         write(*,*) 'EnrTraceDispl = '
         do i=1,6*nrdofHH+1
            write(*,6001) i,EnrTraceDispl(i,1:3*nrdofH)
         enddo
!
         call pause
!
         write(*,*) 'EnrTraceStress = '
         do i=6*nrdofHH,6*nrdofHH+3*nrdofHH
            write(*,6001) i,EnrTraceStress(i,1:3*nrdofV)
         enddo
!
         call pause
!
      endif
#endif
!
!--------------------------------------------------------------------------
!     D P G    F I N A L    A S S E M B L Y
!--------------------------------------------------------------------------
!
!  ...Compact enriched number of rows (total enriched test dof)
      enrdof = 9*nrdofHH
!
!  ...factor the Gram matrix
      call DPPTRF('U', enrdof, Gram, info)
      if (info.ne.0) then
        write(*,*) 'elem_DPG_UWEAK_SYMM: info = ',info ; stop
      endif
!
!  ...Create vector of indices with dof of each physics variable (in the same order of physics file)
      ndofphysics = (/3*nrdofH,3*nrdofV,3*nrdofQ,6*nrdofQ/)
!
!  ...Construct Stiff_ALL by packing all Enriched Matrices in the order H1,Hcurl,Hdiv,L2
!     (in the same order of physics file) and load
      Stiff_ALL(:,:) = ZERO
!
!  ...EnrTraceDispl (H1)
      kmin = 0
      kmax = kmin + ndofphysics(1)
      Stiff_ALL(1:enrdof,kmin+1:kmax) = EnrTraceDispl(1:enrdof,1:(kmax-kmin))
!  ...EnrTraceStress (H(div))
      kmin = kmax
      kmax = kmin + ndofphysics(2)
      Stiff_ALL(1:enrdof,kmin+1:kmax) = EnrTraceStress(1:enrdof,1:(kmax-kmin))
!  ...EnrFieldDispl (L2)
      kmin = kmax
      kmax = kmin + ndofphysics(3)
      Stiff_ALL(1:enrdof,kmin+1:kmax) = EnrFieldDispl(1:enrdof,1:(kmax-kmin))
!  ...EnrFieldStress (L2)
      kmin = kmax
      kmax = kmin + ndofphysics(4)
      Stiff_ALL(1:enrdof,kmin+1:kmax) = EnrFieldStress(1:enrdof,1:(kmax-kmin))
!  ...EnrLoad
      kmin = kmax
      kmax = kmin + NR_RHS
      Stiff_ALL(1:enrdof,kmin+1:kmax) = EnrLoad(1:enrdof,1:(kmax-kmin))
!
!  ...G^-1 * Stiff_ALL
      call DPPTRS('U',enrdof,kmax,Gram,Stiff_ALL,9*MAXbrickHH,info)
      if (info.ne.0) then
         write(*,*) 'elem_DPG_UWEAK_SYMM: info = ',info ; stop
      endif
!
!  ...Build full DPG matrix (stiffness + load) in one go
      FullDPG(:,:) = ZERO
      m = kmin
      n = kmax
      k = enrdof
      l = 9*MAXbrickHH
      call DGEMM('T','N',m,n,k,1.d0,EnrStiffness,l,Stiff_ALL,l,0.d0,FullDPG,3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ)
!
!     ULTIMATE DPG LOAD VECTORS AND STIFFNESS MATRICES THROUGH STATIC CONDENSATION
!
!  ...Populate the assembly matrices accordingly
      n1 = 0
      do i = 1,NR_PHYSA
         n = ndofphysics(i)
         n2 = n1+n
         m1 = 0
         do j = 1,NR_PHYSA
            m = ndofphysics(j)
            m2 = m1 + m
!        ...First initialize
            ALOC(i,j)%array(:,:) = ZERO
!        ...Then populate
            ALOC(i,j)%array(1:n,1:m) = FullDPG(n1+1:n2,m1+1:m2)
            m1 = m2
        enddo
        m2 = m1 + NR_RHS
!  .....First initialize
        BLOC(i)%array(:,:) = ZERO
!  .....Then populate
        BLOC(i)%array(1:n,1:NR_RHS) = FullDPG(n1+1:n2,m1+1:m2)
        n1 = n2
      enddo
!
   end subroutine elem_DPG_UWEAK_SYMM
