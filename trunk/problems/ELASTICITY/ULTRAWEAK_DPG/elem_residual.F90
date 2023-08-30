!--------------------------------------------------------------------------
!> Purpose : return stiffness matrix and Residual vector for element
!!
!! @param[in]  Mdle      - an element (middle node) number
!! @param[out] Resid     - element residual (squared)
!! @param[out] Nref_flag - suggested h-refinement flag
!--------------------------------------------------------------------------
   subroutine elem_residual(Mdle, Resid,Nref_flag)
!
      use control,                  only : INTEGRATION
      use parametersDPG
      use element_data
      use data_structure3D
      use isotropic_elast_material
      use common_prob_data,         only: SYMMETRY_TOL, TEST_NORM, ALPHA
      use assembly,                 only: NR_RHS
!
      implicit none
!
      integer, intent(in)  :: Mdle
      real(8), intent(out) :: Resid
      integer, intent(out) :: Nref_flag
!
!  ...element data
      integer :: ntype,ftype
      integer :: nrv, nre, nrf
      integer :: norder(19), nordf(5), nordP, norderP(19)
      integer :: nedge_orient(12), nface_orient(6)
!
!  ...shape functions
      real(8) :: shapH(MAXbrickH),   gradH(3,MAXbrickH)
      real(8) :: shapHH(MAXbrickHH), gradHH(3,MAXbrickHH)
      integer :: nrdofH, nrdofHH, nrdofEE, nrdofVV, nrdofQQ
!
!  ...geometry
      real(8) :: xnod(3,MAXbrickH)
      real(8) :: xi(3), x(3), rn(3)
      real(8) :: dxdxi(3,3), dxidx(3,3)
      real(8) :: t(2), rjac, brjac
      real(8) :: dxidt(3,2), dxdt(3,2),rt(3,2)
      integer :: nsign
!
!  ...Resid vector for the enriched space
      real(8), allocatable :: EnrResid(:), EnrResidc(:)
!
!  ...gram matrix
      real(8), allocatable :: Gram(:)
!
!  ...tensors in physical coordinates
      real(8), dimension(3,3,3,3) :: A,AA,Symm
!
!  ...temporary variables
      real(8) :: tmp, tmpDivTau1(3), tmpDivTau2(3), tmpTau_n(3)
!
!  ...source term (don't need Neumann term)
      real(8) :: fval(3,NR_RHS)
!
!  ...quadrature data
      real(8) :: xiloc(3,MAXNINT3ADD), wxi(MAXNINT3ADD)
!
!  ...2D quadrature data for boundary terms
      real(8) :: tloc(2,MAXNINT2ADD), wt(MAXNINT2ADD)
      real(8) :: weight, wa
!
!  ...approximate solution
      real(8) :: dofH(MAXEQNH,MAXbrickH)
      real(8) :: dofE(MAXEQNE,MAXbrickE)
      real(8) :: dofV(MAXEQNV,MAXbrickV)
      real(8) :: dofQ(MAXEQNQ,MAXbrickQ)
      real(8) :: solH (  MAXEQNH), dsolH(  MAXEQNH,3)
      real(8) :: solE (3,MAXEQNE), curlE(3,MAXEQNE  )
      real(8) :: solV (3,MAXEQNV), divV (  MAXEQNV  )
      real(8) :: solQ (  MAXEQNQ)
!
!  ...displacement
      real(8) :: u(3)
!
!  ...stress
      real(8) :: sigma(3,3)
      real(8) :: sigma_n(3)
!
!  ...misc
      integer :: i, j, k, m, n, k1, k2, m1, m2
      integer :: ic, jc, ijc, kc, lc, klc
      integer :: ipt, nint, ifc, iload, iflag, info, nordtmp
      integer :: nrTest
      real(8) :: diff, DDOT
!
      integer :: iprint = 0
!
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1
!
!--------------------------------------------------------------------------
!
!  ...Load number for which to calculate residual from
      iload = 1
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
!  ...set the enriched order of approximation
      nordtmp = NORD_ADD
      select case(ntype)
      case(MDLB);       nordP = NODES(Mdle)%order + nordtmp*111
      case(MDLP);       nordP = NODES(Mdle)%order + nordtmp*11
      case(MDLN,MDLD);  nordP = NODES(Mdle)%order + nordtmp*1
      end select
!
!  ...get enriched order vector
      call compute_enriched_order(ntype,nordP, norderP)
!
!  ...number of element test DOFs
      call celndof(ntype,norderP, nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!
!  ...determine solution dof
      call solelm(Mdle, dofH,dofE,dofV,dofQ)
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,7020) xnod(1,1:8),xnod(2,1:8),xnod(3,1:8)
 7020    format('elem_residual: xnod  = ',8(f8.3,2x),  &
            2(  /,'                       ',8(f8.3,2x)))
         write(*,7030) dofH(1,1:8),dofV(1,1:6),dofQ(1:2,1)
 7030    format('elem_residual: dofH(1,1:8) = ',8(e12.5,2x),  &
               /,'               dofV(1,1:6) = ',6(e12.5,2x),  &
               /,'               dofQ(1:2,1) = ',2(e12.5,2x))
      endif
#endif
!
!  ...initialize the enriched local element residual vector
      nrTest = 9*nrdofHH
      allocate(EnrResid(nrTest)); EnrResid(:) = ZERO
!
!  ...initialize the Gram matrix
      allocate(Gram(nrTest*(nrTest+1)/2))
      Gram(:) = ZERO
!
!--------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L
!--------------------------------------------------------------------------
!
!  ...set up the element quadrature
      INTEGRATION = nordtmp
      call set_3Dint(ntype,norder, nint,xiloc,wxi)
      INTEGRATION = 0
!
!  ...Auxiliary tensors
      call getSymm(Symm)
!
!  ...loop through integration points
      do ipt=1,nint
         xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt)
!
!     ...Compute shape functions needed for geometry
         call shape3DH(ntype,xi,norder,nedge_orient,nface_orient,  &
                      nrdofH,shapH,gradH)
!     ...H1 (test)
         call shape3HH(ntype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!     ...geometry map
         call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
                     x,dxdxi,dxidx,rjac,iflag)
         if (iflag.ne.0) then
            write(*,1000) Mdle,rjac
 1000       format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
            stop
         endif
!
!     ...pullbacks
         do k=1,nrdofHH
            gradHH(1:3,k) = gradHH(1,k)*dxidx(1,1:3)  &
                          + gradHH(2,k)*dxidx(2,1:3)  &
                          + gradHH(3,k)*dxidx(3,1:3)
         enddo
!
!     ...integration weight
         weight = wa * rjac
!
!     ...compute the approximate solution on the PHYSICAL element
         call soleval(mdle,xi,nedge_orient,nface_orient,             &
                      norder,xnod,dofH,dofE,dofV,dofQ,1,             &
                      x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)
!
         u(1:3)       = solQ(1:3)
         sigma(1,1:3) = (/ solQ(4), solQ(5), solQ(7) /)
         sigma(2,1:3) = (/ solQ(5), solQ(6), solQ(8) /)
         sigma(3,1:3) = (/ solQ(7), solQ(8), solQ(9) /)
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
!    P A R T  1 : go through \tau\in H(div)^3 test space
!
!     ...FIRST OUTER loop through enriched stress dofs
         do k1=1,nrdofHH
!        ...OUTER loop through components
            do ic=1,3
               do jc=ic,3
                  ijc = nk(ic,jc)
                  m1 = (k1-1)*6+ijc
!
                  tmpDivTau1=0.d0
                  tmpDivTau1(ic) = gradHH(jc,k1)
                  tmpDivTau1(jc) = gradHH(ic,k1)
!
!     --- residual ---
!
!   ( A:\sigma, \tau ) + ( u , div(\tau) )
!
                  tmp=0.d0
                  if (ic.eq.jc) then
                     do m=1,3
                        do n=1,3
                        tmp = tmp  &
                            + A(ic,jc,m,n) * sigma(m,n) * shapHH(k1)
                        enddo
                     enddo
                  else
                     do m=1,3
                        do n=1,3
                           tmp = tmp  &
                               + A(ic,jc,m,n) * sigma(m,n) * shapHH(k1)  &
                               + A(jc,ic,m,n) * sigma(m,n) * shapHH(k1)
                        enddo
                     enddo
                  endif
                  do m=1,3
                     tmp = tmp + u(m) * tmpDivTau1(m)
                  enddo
                  EnrResid(m1) = EnrResid(m1) + tmp*weight
!
!     --- Gram matrix ---
!
                  select case(TEST_NORM)
!
!   (A:\tau_1,A:\tau_2+\varepsilon(v_2))+(div(\tau_1),div(tau_2))
!
                  case(1)
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
                        tmpDivTau2 = ZERO
                        tmpDivTau2(kc) = gradHH(lc,k2)
                        tmpDivTau2(lc) = gradHH(kc,k2)
!
                        tmp = ZERO
!
!                    ...L2-terms
                        if (ic.ne.jc) then
                           if (kc.ne.lc) then
                              tmp = AA(kc,lc,ic,jc) * shapHH(k2) * shapHH(k1)  &
                                  + AA(kc,lc,jc,ic) * shapHH(k2) * shapHH(k1)  &
                                  + AA(lc,kc,ic,jc) * shapHH(k2) * shapHH(k1)  &
                                  + AA(lc,kc,jc,ic) * shapHH(k2) * shapHH(k1)
                           else
                              tmp = AA(kc,lc,ic,jc) * shapHH(k2) * shapHH(k1)  &
                                  + AA(kc,lc,jc,ic) * shapHH(k2) * shapHH(k1)
                           endif
                        elseif (kc.ne.lc) then
                           tmp = AA(kc,lc,ic,jc) * shapHH(k2) * shapHH(k1)  &
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
                        Gram(k) = Gram(k) + tmp * weight
                     enddo
!
                     do m2 = 6*nrdofHH+1,6*nrdofHH+3*nrdofHH
                        k = nk(m1,m2)
                        k2 = int((m2-6*nrdofHH-1)/3)+1
                        kc = mod(m2-6*nrdofHH-1,3)+1

                        tmp = ZERO
                        if (ic.eq.jc) then
                           do m = 1,3
                              tmp = tmp  &
                                  + A(ic,jc,kc,m) * gradHH(m,k2) * shapHH(k1)
                           enddo
                        else
                           do m = 1,3
                              tmp = tmp  &
                                  + A(ic,jc,kc,m) * gradHH(m,k2) * shapHH(k1)  &
                                  + A(jc,ic,kc,m) * gradHH(m,k2) * shapHH(k1)
                           enddo
                        endif
!
                        Gram(k) = Gram(k) + tmp * weight
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
                        tmpDivTau2 = ZERO
!
                        tmpDivTau2(kc) = gradHH(lc,k2)
                        tmpDivTau2(lc) = gradHH(kc,k2)
!
                        tmp = ZERO
!
!                    ...L2-terms
                        if (ic.eq.jc) then
                           if ((kc.eq.ic).and.(lc.eq.jc)) then
                              tmp = shapHH(k1) * shapHH(k2)
                           endif
                        elseif ((kc.eq.ic).and.(lc.eq.jc)) then
                           tmp = 2*shapHH(k1) * shapHH(k2)
                        endif
!
!                    ...div term
                        do m=1,3
                           tmp = tmp + tmpDivTau1(m) * tmpDivTau2(m)
                        enddo
!
                        Gram(k) = Gram(k) + tmp * weight
                     enddo
                  end select
!
!           ...END OUTER LOOP through test stresses
               enddo
            enddo
        enddo
!
!    P A R T  2 : go through v\in(H1)^3 test space
!
!     ...SECOND OUTER loop through enriched H1 dofs
         do k1=1,nrdofHH
!        ...OUTER loop through components
            do ic=1,3
               m1 = 6*nrdofHH+(k1-1)*3+ic
!
!     --- residual ---
!
!   + (\sigma, grad(v) ) - (f,v)
!
               tmp=0.d0
               do m=1,3
               tmp = tmp + sigma(ic,m) * gradHH(m,k1)
               enddo
               tmp = tmp - fval(ic,iload) * shapHH(k1)
               EnrResid(m1) = EnrResid(m1) + tmp*weight
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
!
                     tmp = ZERO
!
!                 ...grad terms
                     do m=1,3
                        do n=1,3
                           tmp = tmp  &
                               + Symm(jc,m,ic,n) * gradHH(m,k2) * gradHH(n,k1)
                        enddo
                     enddo
!
!                 ...L2-term
                     if (ic.eq.jc) then
                        tmp = tmp + shapHH(k1) * shapHH(k2)
                     endif
!
                     Gram(k) = Gram(k) + tmp * weight
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
                              + (  shapHH(k1) * shapHH(k2)                 &
                                 + gradHH(1,k1) * gradHH(1,k2)             &
                                 + gradHH(2,k1) * gradHH(2,k2)             &
                                 + gradHH(3,k1) * gradHH(3,k2) ) * weight
                     endif
                  enddo
               end select
!
!        ...END OUTER LOOP
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
         INTEGRATION = nordtmp
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
!        ...Compute shape functions needed for test variables and geometry
            call shape3DH(ntype,xi,norder,nedge_orient,nface_orient,     &
                        nrdofH,shapH,gradH)
            call shape3HH(ntype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!        ...geometry map
            call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,   &
                         x,dxdxi,dxidx,rjac,dxdt,rn,brjac)
            weight = wa*brjac
!
!        ...compute the approximate solution on the PHYSICAL element
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,      &
                         xnod,dofH,dofE,dofV,dofQ,1,                    &
                         x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)
!
            u(1:3)       = solH(1:3)
            sigma_n(1:3) = solV(1,1:3)*rn(1) + solV(2,1:3)*rn(2) + solV(3,1:3)*rn(3)
!
!    P A R T  1 : go through \tau\in H(div)^3 test space
!
!        ...OUTER loop through enriched test functions
            do k1=1,nrdofHH
               do ic=1,3
                  do jc=ic,3
                     ijc = nk(ic,jc)
                     m1 = (k1-1)*6+ijc

                     tmpTau_n(1:3)=0.d0
                     tmpTau_n(ic) = shapHH(k1) * rn(jc)
                     tmpTau_n(jc) = shapHH(k1) * rn(ic)
!
!     --- trace residual ---
!
!   - <\hat u,(\tau n)>
!
                     EnrResid(m1) = EnrResid(m1)                     &
                                  - ( u(1) * tmpTau_n(1)             &
                                    + u(2) * tmpTau_n(2)             &
                                    + u(3) * tmpTau_n(3) ) * weight
!
                  enddo
               enddo
!        ...END OUTER loop
            enddo
!
!    P A R T  2 : go through v\in(H1)^3 test space
!
!        ...SECOND OUTER loop through enriched H1 test function
            do k1=1,nrdofHH
               do kc=1,3
                  m1 = 6*nrdofHH+(k1-1)*3+kc
!
!     --- trace residual ---
!
!   - <sigma_n,v>
!
                  EnrResid(m1) = EnrResid(m1)                        &
                               - sigma_n(kc) * shapHH(k1) * weight
!
!           ...OUTER loop
               enddo
            enddo
!     ...end of loop over integration points
         enddo
!  ...end of loop over faces
      enddo
!
!--------------------------------------------------------------------------
!     D P G    F I N A L    A S S E M B L Y
!--------------------------------------------------------------------------
!
!  ...factor the Gram matrix
      call DPPTRF('U',9*nrdofHH,Gram, info)
      if (info.ne.0) then
        write(*,*) 'elem_residual: info = ',info ; stop
      endif
!
!  ...save copies of the RHS to compute later the residual
      allocate(EnrResidc(9*nrdofHH))
      EnrResidc(:) = EnrResid(:)
!
!  ...compute the product of inverted test Gram matrix with RHS
      call DPPTRS('U',9*nrdofHH,1,Gram,EnrResid,NrTest, info)
      if (info.ne.0) then
         write(*,*) 'elem_residual: info1 = ',info
         stop 1
      endif
!
!  ...compute the residual
      Resid = DDOT(9*nrdofHH,EnrResidc,1,EnrResid,1)
      deallocate(Gram,EnrResid,EnrResidc)
!
!  ...determines what refinement flag to use
!  ...if isotropic h refinements
      call get_isoref(Mdle, Nref_flag)
!
   end subroutine elem_residual





!--------------------------------------------------------------------------
!> @brief      Returns global residual
!!
!> @date       July 2023
!--------------------------------------------------------------------------
   subroutine compute_residual
!
      use data_structure3D
      use environment,        only: QUIET_MODE
      use par_mesh,           only: DISTRIBUTED
      use assembly_sc,        only: NRDOF_CON, NRDOF_TOT
      use MPI
      use mpi_param
!
      implicit none
!
!  ...middle node
      integer :: mdle
!
!  ...residual
      real(8) :: resid, residual, rate
!
!  ...number of dof for higher order node
      integer :: ndofH, ndofE, ndofV, ndofQ, nrdof_total
!
!  ...refinement flag
      integer :: nref_flag
!
!  ...old data
      integer, save :: ivis = 0
      real(8), save :: residual_old
      integer, save :: nrdof_total_old
      real(8), save :: rwork(2,10)
      integer, save :: iwork(10)
!
!  ...miscellaneous
      integer :: i, iel, iphys, ierr
!
      integer :: iprint = 0
!
!--------------------------------------------------------------------------
!
!  ...initialize
      residual = 0.d0
      if (.not.DISTRIBUTED) then
         NRELES_SUBD = NRELES
         ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      endif
!
!$OMP PARALLEL DO DEFAULT(SHARED)         &
!$OMP PRIVATE(iel,mdle,resid,nref_flag))  &
!$OMP REDUCTION(+:residual)
!  ...compute residual on subdomain
      do iel=1,NRELES_SUBD
         mdle = ELEM_SUBD(iel)
         call elem_residual(mdle, resid,nref_flag)
         residual = residual + resid
      enddo
!$OMP END PARALLEL DO
!
!  ...reduce if distributed
      if (DISTRIBUTED) then
         call MPI_ALLREDUCE(MPI_IN_PLACE,residual,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif
!
      if (RANK.le.ROOT) then
         write(*,7030) NRDOF_TOT, sqrt(residual)
      endif
 7030 format('compute_residual: TOTAL NUMBER OF DOF, RESIDUAL = ',i8,3x,1e12.5)
!
   end subroutine compute_residual
