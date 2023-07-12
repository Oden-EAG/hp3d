!----------------------------------------------------------------------
!> @brief       Wrapper for element stiffness matrix assembly
!!
!> @param[in]   mdle    - middle node number
!> @param[out]  Itest   - index for assembly
!!                         (trivial since not multiphysics application)
!> @param[out]  Itrial  - index for assembly
!!                         (trivial since not multiphysics application)
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine elem(Mdle, Itest,Itrial)
!
      use physics,       only : NR_PHYSA
      use data_structure3D
      use parametersDPG
      use control,       only : INTEGRATION
!
      implicit none
!
      integer, intent(in)  :: Mdle
      integer, intent(out) :: Itest(NR_PHYSA), Itrial(NR_PHYSA)
!
!  ...element order, orientation for edges and faces
      integer :: norder(19),norderP(19), nordP
      integer :: nordx,nordy,nordz
      integer :: nrdofH,nrdofE,nrdofV,nrdofQ
      integer :: nrdofHi,nrdofVi
      integer :: nrTest, nrTrial
      integer :: nrdofHH,nrdofEE,nrdofVV,nrdofQQ
      integer :: ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl
      integer :: nord_add_local
      integer :: nrv, nre, nrf
      integer :: ntype
!
!----------------------------------------------------------------------
!
!  ...get element type and number of vertices, edges, and faces
      ntype = NODES(Mdle)%ntype
      nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
      Itest (1:NR_PHYSA) = 0
      Itrial(1:NR_PHYSA) = 0
!
      select case(NODES(Mdle)%case)
      case(7)
         Itest(1:NR_PHYSA)=1; Itrial(1:NR_PHYSA)=1
!
!     ...get DPG enrichment order
         nord_add_local = NORD_ADD
!      
!     ...determine order of approximation
         call find_order(Mdle, norder)
!
 10      continue
!     ...set the enriched order of approximation
         select case(ntype)
            case(MDLB)      ; nordP = NODES(Mdle)%order + nord_add_local*111
            case(MDLN,MDLD) ; nordP = NODES(Mdle)%order + nord_add_local*1
            case(MDLP)      ; nordP = NODES(Mdle)%order + nord_add_local*11
         end select
!
         call compute_enriched_order(nordP, norderP)
!
!     ...number of element trial DOFs
         call celndof(ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!     ...number of element test DOFs
         call celndof(ntype,norderP,nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!     ...number of element bubble DOFs of each type
         call ndof_nod(ntype,norder(nre+nrf+1),ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
!
!     ...number of test functions
         nrTest = nrdofHH + nrdofVV
!
!     ...number of trial functions of each type (H,V,Q)
         nrdofHi = nrdofH - ndofHmdl
         nrdofVi = nrdofV - ndofVmdl
         nrTrial = nrdofHi + nrdofVi + 4*nrdofQ
!      
         if (nrTest .le. nrTrial-ndofHmdl-ndofVmdl) then
            nord_add_local = nord_add_local + 1
!$omp critical         
            write(*,*) 'elem: WARNING mdle =  ', mdle
            write(*,*) 'elem: nrTest, nrTrial = ', nrTest, nrTrial-ndofHmdl-ndofVmdl
            write(*,*) 'elem: nord_add_local  = ', nord_add_local
!$omp end critical                  
            go to 10
         endif
!
         call elem_DPG_UWEAK_ACOUSTICS(mdle,nord_add_local,nrTest,nrTrial,    &
                                       nrdofHH,nrdofVV,nrdofH,nrdofV,nrdofQ,  &
                                       ndofHmdl,ndofVmdl)
!
      case default
         write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',  &
         Mdle,NODES(Mdle)%case
         call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
      end select
!
   end subroutine elem





!----------------------------------------------------------------------
!> @brief       Assembles element stiffness matrix
!!
!> @param[in]   mdle             - middle node number
!> @param[in]   nord_add_local   - DPG enrichment order
!> @param[in]   nrTest           - number of element test functions
!> @param[in]   nrTrial          - number of element trial functions
!> @param[in]   nrdofHH          - number of enriched H^1 element DOFs
!> @param[in]   nrdofVV          - number of enriched H(div) element DOFs
!> @param[in]   nrdofH           - number of H^1 element DOFs
!> @param[in]   nrdofH           - number of H(div) element DOFs
!> @param[in]   nrdofH           - number of L^2 element DOFs
!> @param[in]   nrdofHmdl        - number of H(1) bubble DOFs
!> @param[in]   nrdofVmdl        - number of H(div) bubble DOFs
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine elem_DPG_UWEAK_ACOUSTICS(Mdle,nord_add_local,nrTest,nrTrial,    &
                                       nrdofHH,nrdofVV,nrdofH,nrdofV,nrdofQ,  &
                                       ndofHmdl,ndofVmdl)
!
      use data_structure3D
      use parametersDPG
      use control,       only: INTEGRATION
      use assembly,      only: ALOC,BLOC,NR_RHS
      use common_prob_data_UW
!
      implicit none
!
      integer, intent(in) :: mdle
      integer, intent(in) :: nord_add_local
      integer, intent(in) :: nrTest, nrTrial
      integer, intent(in) :: nrdofHH, nrdofVV
      integer, intent(in) :: nrdofH, nrdofV, nrdofQ
      integer, intent(in) :: ndofHmdl, ndofVmdl
!
!  ...element and face types
      integer :: ntype, ftype
!
!  ...element data
      integer :: norder(19), norderf(5)
      integer :: norient_edge(12), norient_face(6)
      integer :: ibc(6,NRINDEX)
      integer :: nordP, nre, nrv, nrf
!
!  ...geometry
      real(8) :: xnod(3,MAXbrickH)
      real(8) :: xi(3), x(3), t(2)
      real(8) :: dxdxi(3,3), dxidx(3,3), rjac
      real(8) :: dxidt(3,2), dxdt(3,2), bjac
      real(8) :: rt(3,2), rn(3)
!
!  ...quadrature
      real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
      real(8) :: tloc(2,MAXNINT2ADD),  wtloc(MAXNINT2ADD)
      real(8) :: wa, weight
!
!  ...trial and test variables
      real(8) :: u(3), v(3), div_u, div_v, un, vn, vec(3)
      real(8) :: q, p, dq(3), dp(3)
      integer :: nrdofHi, nrdofVi
!
!  ...load
      complex(8) :: zf(4), zg
!
!  ...shape function workspace
      real(8) :: shapH(MAXbrickH),     gradH(3,MAXbrickH)
      real(8) :: shapHH(MAXbrickHH),   gradHH(3,MAXbrickHH)
      real(8) :: shapV(3,MAXbrickV),   divV(MAXbrickV)
      real(8) :: shapVV(3,MAXbrickVV), divVV(MAXbrickVV)
      real(8) :: shapQ(MAXbrickQ)
!
!  ...load and intermediate matrices
      complex(8), allocatable :: BloadHV(:), Gram(:,:)
      complex(8), allocatable :: StiffHV_H(:,:), StiffHV_V(:,:), StiffHV_Q(:,:)
      complex(8), allocatable :: Stiff_ALL(:,:)
      complex(8), allocatable :: Zaloc(:,:)
!
!  ...lapack
      character(1) :: uplo
!
!  ...misc
      real(8) :: h_elem, omeg
      integer :: m, k1, k2, n1, n2, i, i1, j1, j2, j3, ifc, l
      integer :: nint, n, nrhs, nsign, iflag, info
!
      integer :: iprint = 0
!
!----------------------------------------------------------------------
!
!  ...element type
      ntype = NODES(Mdle)%ntype
      nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
!  ...determine order of approximation
      call find_order(Mdle, norder)
!
!  ...set the enriched order of approximation
      select case(ntype)
         case(MDLB);      nordP = NODES(Mdle)%order + nord_add_local*111
         case(MDLN,MDLD); nordP = NODES(Mdle)%order + nord_add_local*1
         case(MDLP);      nordP = NODES(Mdle)%order + nord_add_local*11
      end select
!
!  ...determine edge and face orientations
      call find_orient(Mdle, norient_edge,norient_face)
!                                                                     
!  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
!
!  ...determine element size and scale correction
      call find_hmin(Mdle, h_elem)
! 
!  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
!
      omeg = OMEGA
!
!  ...initialize matrices and load vector
      allocate(BloadHV(nrTest));                BloadHV(:)     = ZERO
      allocate(Gram(nrTest,nrTest));            Gram(:,:)      = ZERO
!
      allocate(StiffHV_H(nrdofH,nrTest));       StiffHV_H(:,:) = ZERO
      allocate(StiffHV_V(nrdofV,nrTest));       StiffHV_V(:,:) = ZERO
      allocate(StiffHV_Q(4*nrdofQ,nrTest));     StiffHV_Q(:,:) = ZERO
      allocate(Stiff_ALL(nrTest,nrTRIAL+1));    Stiff_ALL(:,:) = ZERO
      allocate(Zaloc(nrTRIAL+1,nrTRIAL+1));     Zaloc(:,:)     = ZERO
!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                             
!----------------------------------------------------------------------
!
!  ...use the enriched order to set the quadrature
      INTEGRATION = nord_add_local
      call set_3Dint_DPG(ntype,norder, nint,xiloc,waloc)
      INTEGRATION = 0
!
!  ...loop over integration points
      do l=1,nint
!      
         xi(1:3)=xiloc(1:3,l) ; wa=waloc(l)
!
!     ...H1 shape functions (for geometry)
         call shape3DH(ntype,xi,norder,norient_edge,norient_face,nrdofH,shapH,gradH)
!
!     ...L2 shape functions for the trial space
         call shape3DQ(ntype,xi,norder, nrdofQ,shapQ)
!
!     ...discontinuous H1 shape functions for the enriched test space
         call shape3HH(ntype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!     ...discontinuous H(div) shape functions for the enriched test space
         call shape3VV(ntype,xi,nordP, nrdofVV,shapVV,divVV)
!
!     ...geometry map
         call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,x,dxdxi,dxidx,rjac,iflag)
!
!     ...integration weight
         weight = rjac*wa
!
!     ...get the RHS
         call getf(Mdle,x, zf)
!
!   ...loop through enriched H1 test functions in the enriched space
         do k1=1,nrdofHH
!
!        ...Piola transformation
            q = shapHH(k1)
            dq(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                    + gradHH(2,k1)*dxidx(2,1:3) &
                    + gradHH(3,k1)*dxidx(3,1:3)
!
!
!     --- Load vector ---
!
            BloadHV(k1) = BloadHV(k1) + q*zf(1)*weight
!
!     --- Stiffness matrix ---
!
!        ...loop through L2 trial shape functions
            do k2=1,nrdofQ
!          
!           ...Piola transformation
               p = shapQ(k2)/rjac ; u(1:3) = p
!
               m = (k2-1)*4+1
               StiffHV_Q(m,k1) = StiffHV_Q(m,k1) + ZIMG*OMEGA*q*p*weight
!
               m = (k2-1)*4+2
               StiffHV_Q(m,k1) = StiffHV_Q(m,k1) - dq(1)*u(1)*weight
!
               m = (k2-1)*4+3
               StiffHV_Q(m,k1) = StiffHV_Q(m,k1) - dq(2)*u(2)*weight
!
               m = (k2-1)*4+4
               StiffHV_Q(m,k1) = StiffHV_Q(m,k1) - dq(3)*u(3)*weight
!
            enddo
!
!     --- Gram matrix ---
!
!        ...H1 test functions
            do k2=1,k1
!           ...Piola transformation
               p = shapHH(k2)
               dp(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                       + gradHH(2,k2)*dxidx(2,1:3) &
                       + gradHH(3,k2)*dxidx(3,1:3)
! 
!           ...accumulate for the gram matrix
               Gram(k2,k1) = Gram(k2,k1)                              &
                           + (dq(1)*dp(1)+dq(2)*dp(2)+dq(3)*dp(3)   &
                           + (omeg**2+ALPHA)*q*p)*weight
            enddo
!     ...end of loop through H1 test functions         
         enddo
!      
!     ...loop through enriched H(div) test functions in the enriched space
         do k1 = 1,nrdofVV
!
            n1 = nrdofHH+k1
!
!        ...Piola transformation
            v(1:3) = (dxdxi(1:3,1)*shapVV(1,k1)                 &
                   +  dxdxi(1:3,2)*shapVV(2,k1)                 &
                   +  dxdxi(1:3,3)*shapVV(3,k1))/rjac
            div_v  =  divVV(k1)/rjac
!
!     --- Load vector ---
!
            BloadHV(n1) = BloadHV(n1) + (v(1)*zf(2)+v(2)*zf(3)+v(3)*zf(4))*weight
!
!     --- Stiffness matrix ---
!
!        ...loop through L2 trial shape functions
            do k2=1,nrdofQ
!
!           ...Piola transformation
               p = shapQ(k2)/rjac; u(1:3)=p
!
               m = (k2-1)*4+1
               StiffHV_Q(m,n1) = StiffHV_Q(m,n1) - div_v*p*weight
!
               m = (k2-1)*4+2
               StiffHV_Q(m,n1) = StiffHV_Q(m,n1) + ZIMG*OMEGA*v(1)*u(1)*weight
!
               m = (k2-1)*4+3
               StiffHV_Q(m,n1) = StiffHV_Q(m,n1) + ZIMG*OMEGA*v(2)*u(2)*weight
!
               m = (k2-1)*4+4
               StiffHV_Q(m,n1) = StiffHV_Q(m,n1) + ZIMG*OMEGA*v(3)*u(3)*weight
!            
         enddo
!
!     --- Gram matrix ---
!
!        ...H1 test functions
            do k2=1,nrdofHH
!           ...Piola transformation
               p = shapHH(k2)
               dp(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                       + gradHH(2,k2)*dxidx(2,1:3) &
                       + gradHH(3,k2)*dxidx(3,1:3)
! 
!           ...accumulate for the gram matrix
               Gram(k2,n1) = Gram(k2,n1)                            &
                           + ZIMG*omeg*(dp(1)*v(1)+dp(2)*v(2)+dp(3)*v(3) - p*div_v)*weight
            enddo
            do k2 = 1,k1
               n2 = nrdofHH+k2
!
!           ...Piola transformation
               u(1:3) = (dxdxi(1:3,1)*shapVV(1,k2)             &
                      +  dxdxi(1:3,2)*shapVV(2,k2)             &
                      +  dxdxi(1:3,3)*shapVV(3,k2))/rjac
               div_u  =  divVV(k2)/rjac

               Gram(n2,n1) = Gram(n2,n1)        &
                           + ((omeg**2+ALPHA)*(v(1)*u(1)+v(2)*u(2)+v(3)*u(3))  &
                           + (div_v*div_u))*weight
            enddo
!     ...end of loop through H(div) test functions
         enddo
!  ...end of loop through integration points
      enddo
!
!----------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L S                           
!----------------------------------------------------------------------
!
      nrdofHi = nrdofH - ndofHmdl
      nrdofVi = nrdofV - ndofVmdl
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
         call face_order(ntype,ifc,norder, norderf)
!
!     ...set 2D quadrature
         INTEGRATION = nord_add_local
         call set_2Dint_DPG(ftype,norderf, nint,tloc,wtloc)
         INTEGRATION = 0
!
!     ...loop through integration points
         do l = 1,nint
!
!        ...face coordinates
            t(1:2) = tloc(1:2,l)
!
!        ...face parametrization
            call face_param(ntype,ifc,t, xi,dxidt)
! 
!        ...determine element H1 shape functions (for geometry)
            call shape3DH(ntype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
! 
!        ...determine element Hdiv shape functions (for fluxes)
            call shape3DV(ntype,xi,norder,norient_face, nrdofV,shapV,divV)
!
!        ...determine discontinuous H1 shape functions
            call shape3HH(ntype,xi,nordP, nrdofHH,shapHH,gradHH)
!         
!        ...determine discontinuous H(div) shape functions
            call shape3VV(ntype,xi,nordP, nrdofVV,shapVV,divVV)
!
!        ...geometry
            call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
            weight = bjac*wtloc(l)
!
!        ...check if on impedance boundary
            if (ibc(ifc,2) .eq. 3) then
!           ...get boundary data
               call getg(mdle,x,rn,ibc(ifc,2),zg)
!
!           ...loop through H1 test functions
               do k1 = 1,nrdofHH
!              ...value of the shape function at the point
                  q = shapHH(k1)
!              ...accumulate for the load vector
                  BloadHV(k1) = BloadHV(k1) + q*zg*weight
!              ...loop through H1 trial functions
                  do k2 = 1,nrdofHi
!                 ...value of the shape function at the point
                     p = shapH(k2)
!                 ...accumulate for the Stiffness matrix
                     StiffHV_H(k2,k1) = StiffHV_H(k2,k1) + q*p*weight
                  enddo
!
               enddo
!
!           ...regular boundary
            else
! 
!           ...loop through enriched H1 test functions
               do k1 = 1,nrdofHH
                  q  = shapHH(k1)
!
!              ...loop through H(div) trial functions
                  do k2=1,nrdofVi
! 
!                 ...normal component (Piola transformation)
                     vec(1:3) = dxdxi(1:3,1)*shapV(1,k2)   &
                              + dxdxi(1:3,2)*shapV(2,k2)   &
                              + dxdxi(1:3,3)*shapV(3,k2)
                     vec(1:3) = vec(1:3)/rjac
                     un = vec(1)*rn(1)+vec(2)*rn(2)+vec(3)*rn(3)
!
!                 ...accumulate for the Stiffness matrix
                     StiffHV_V(k2,k1) = StiffHV_V(k2,k1) + q*un*weight
!
!                 ...end of loop through H(div) trial functions
                  enddo
!              ...end of loop through H1 test functions
               enddo
            endif
!         
!        ...loop through H(div) enriched test functions
            do k1 = 1,nrdofVV
               n1 = nrdofHH+k1
!
!           ...normal component of the test function (Piola transformation at work!)
               vec(1:3) = dxdxi(1:3,1)*shapVV(1,k1)   &
                        + dxdxi(1:3,2)*shapVV(2,k1)   &
                        + dxdxi(1:3,3)*shapVV(3,k1)
               vec(1:3) = vec(1:3)/rjac
               vn = vec(1)*rn(1)+vec(2)*rn(2)+vec(3)*rn(3)
!
!           ...loop through H1 trial functions
               do k2=1,nrdofHi
!
!              ...value of the shape function at the point
                  p = shapH(k2)
!
!              ...accumulate for the rectangular Stiffness matrix
                  StiffHV_H(k2,n1) = StiffHV_H(k2,n1) + vn*p*weight
!
!              ...end of loop through H1 trial functions
               enddo
!           ...end of loop through H(div) test functions
            enddo
!        ...end of loop through integration points
         enddo
!     ...end of loop through faces
      enddo
!
!
!----------------------------------------------------------------------
!      Alternative construction of normal matrix
!----------------------------------------------------------------------
!
      i1 = nrTEST ; j1 = nrdofHi ; j2 = nrdofVi ; j3 = 4*nrdofQ
!
      Stiff_ALL(1:i1,1:j1)             = transpose(StiffHV_H(1:j1,1:i1))
      Stiff_ALL(1:i1,j1+1:j1+j2)       = transpose(StiffHV_V(1:j2,1:i1))
      Stiff_ALL(1:i1,j1+j2+1:j1+j2+j3) = transpose(StiffHV_Q(1:j3,1:i1))
      Stiff_ALL(1:i1,j1+j2+j3+1)       = BloadHV(1:i1)
!
      N     = nrTEST
      NRHS  = nrdofHi + nrdofVi + 4*nrdofQ + 1
! 
!  ...factorize the Gram matrix
      call ZPOTRF('U',N,Gram,N,info)
      if (info.ne.0) then
         write(*,*) 'elem: Gram ZPOTRF: Mdle,info = ',Mdle,info
         write(*,*) 'elem: nrTest, nrTrial = ',nrTest, nrTrial
         stop 1
      endif
!
      call ZTRSM('L','U','C','N',N,NRHS,ZONE,Gram,N,Stiff_ALL,N)
      call ZHERK('U','C',NRHS,N,ZONE,Stiff_ALL,N,ZERO,Zaloc,NRHS)
!
      do i=1,NRHS-1
         Zaloc(i+1:NRHS,i) = conjg(Zaloc(i,i+1:NRHS))
      enddo
!
      N = NRHS-1
!
      BLOC(1)%array(1:j1,1) = Zaloc(1:j1,j1+j2+j3+1)
      BLOC(2)%array(1:j2,1) = Zaloc(j1+1:j1+j2,j1+j2+j3+1)
      BLOC(3)%array(1:j3,1) = Zaloc(j1+j2+1:j1+j2+j3,j1+j2+j3+1)
!
      ALOC(1,1)%array(1:j1,1:j1) = Zaloc(1:j1,1:j1)
      ALOC(1,2)%array(1:j1,1:j2) = Zaloc(1:j1,j1+1:j1+j2)
      ALOC(1,3)%array(1:j1,1:j3) = Zaloc(1:j1,j1+j2+1:j1+j2+j3)
!
      ALOC(2,1)%array(1:j2,1:j1) = Zaloc(j1+1:j1+j2,1:j1)
      ALOC(2,2)%array(1:j2,1:j2) = Zaloc(j1+1:j1+j2,j1+1:j1+j2)
      ALOC(2,3)%array(1:j2,1:j3) = Zaloc(j1+1:j1+j2,j1+j2+1:j1+j2+j3)
!
      ALOC(3,1)%array(1:j3,1:j1) = Zaloc(j1+j2+1:j1+j2+j3,1:j1)
      ALOC(3,2)%array(1:j3,1:j2) = Zaloc(j1+j2+1:j1+j2+j3,j1+1:j1+j2)
      ALOC(3,3)%array(1:j3,1:j3) = Zaloc(j1+j2+1:j1+j2+j3,j1+j2+1:j1+j2+j3)
!
      deallocate(Zaloc,Stiff_ALL,Gram)
!   
   end subroutine elem_DPG_UWEAK_ACOUSTICS

 




!----------------------------------------------------------------------
!> @brief       Compute enriched order for all element nodes
!!
!> @param[in]   Nord    - enriched element order
!> @param[out]  Norder  - integer array containing order of elem nodes
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine compute_enriched_order(Nord,Norder)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: Nord
      integer, intent(out) :: Norder(19)
      integer :: temp(2)
      integer :: nordF(3)
      integer :: nordB(3),ndofH(3)
!
!----------------------------------------------------------------------
!
      call decod(Nord,MODORDER,2, temp)
      nordF(1) = temp(1) ; nordB(3) = temp(2)
      call decod(nordF(1),MODORDER,2, nordB(1:2))
      call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2))
      call encod(nordB(2:3),MODORDER,2, nordF(3))
!
!  ...edges
      norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
      norder(5:8)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
      norder(9:12)=nordB(3)
!
!  ...faces
      norder(13:14)=nordF(1)
      norder(15:18)=(/nordF(2),nordF(3),nordF(2),nordF(3)/)
!
!  ...middle
      norder(19)=Nord
!
   end subroutine compute_enriched_order





!----------------------------------------------------------------------
!> @brief       Decode order in x-,y-, and z-directions from element order
!!
!> @param[in]   ntype   - middle node type
!> @param[in]   nord    - element order
!> @param[out]  nordx   - order in x-direction
!> @param[out]  nordy   - order in y-direction
!> @param[out]  nordz   - order in z-direction
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine compute_1D_ord(ntype,nord,nordx,nordy,nordz)
!
      use parametersDPG,     only : MAXPP
      use control,           only : INTEGRATION
!
      implicit none
!
      integer, intent(in)  :: ntype
      integer, intent(in)  :: Nord
      integer, intent(out) :: nordx, nordy, nordz
      integer              :: nordh, nord1,nord2,nord3
!
      call decode(Nord, nordh,nord3)
      call decode(nordh, nord1,nord2)
!
      nordx=min(nord1+INTEGRATION,MAXPP)
      nordy=min(nord2+INTEGRATION,MAXPP)
      nordz=min(nord3+INTEGRATION,MAXPP)
!
   end subroutine compute_1D_ord
     
