!----------------------------------------------------------------------
!> @brief       Compute element residual
!!
!> @param[in]   mdle       - middle node number
!> @param[out]  Resid      - element residual (norm)
!> @param[out]  Nref_flag  - suggested element refinement
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine elem_residual(Mdle, Resid,Nref_flag)
!
      use data_structure3D
      use control
      use parametersDPG
      use common_prob_data_UW
!
      implicit none
!
!  ...element and face types
      integer :: ntype, ftype
!
!  ...element data
      integer :: norder(19), norderf(5)
      integer :: norient_edge(12), norient_face(6)
      integer :: ibc(6,NRINDEX)
      integer :: nre, nrv, nrf
!
!  ...geometry
      real(8) :: xnod(3,MAXbrickH)
      real(8) :: xi(3), x(3), t(2)
      real(8) :: dxdxi(3,3), dxidx(3,3), rjac
      real(8) :: dxidt(3,2), dxdt(3,2), bjac
      real(8) :: rt(3,2), rn(3)
!
!  ...solution DOFs
      complex(8) :: zdofH(MAXEQNH,MAXbrickH)
      complex(8) :: zdofE(MAXEQNE,MAXbrickE)
      complex(8) :: zdofV(MAXEQNV,MAXbrickV)
      complex(8) :: zdofQ(MAXEQNQ,MAXbrickQ)
!
!  ...solution values
      complex(8) :: zsolH(MAXEQNH), zdsolH(MAXEQNH,3)
      complex(8) :: zsolV(MAXEQNV,3), zdivV(MAXEQNV)
      complex(8) :: zsolQ(MAXEQNQ)
!
!  ...quadrature
      real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
      real(8) :: tloc(2,MAXNINT2ADD),  wtloc(MAXNINT2ADD)
      real(8) :: wa, weight
!
!  ...trial and test variables
      real(8)    :: u(3), v(3), div_v, vec(3), vn
      real(8)    :: p, dp(3), q, dq(3),
      complex(8) :: zu(3), zp, zvec(3), zun
      integer    :: nrdofHi, nrdofVi
!
      complex(8) :: zf(4)
!
!  ...shape function workspace
      real(8) :: shapH(MAXbrickH),     gradH(3,MAXbrickH)
      real(8) :: shapHH(MAXbrickHH),   gradHH(3,MAXbrickHH)
      real(8) :: shapV(3,MAXbrickV),   divV(MAXbrickV)
      real(8) :: shapVV(3,MAXbrickVV), divVV(MAXbrickVV)
      real(8) :: shapQ(MAXbrickQ)
!
!  ...load and intermediate matrices
      complex(8), allocatable :: BLOADHV(:), AP(:), temp(:)
!
!  ...lapack
      character(1) :: uplo
!
!  ...misc
      real(8) :: h_elem, omeg
      integer :: m, k, k1, k2, n1, n2, i, i1, j1, j2, j3, ifc, l
      integer :: nint, n, nrhs, nsign, iflag, info
!
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1
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
         case(MDLB);      nordP = NODES(Mdle)%order+NORD_ADD*111
         case(MDLN,MDLD); nordP = NODES(Mdle)%order+NORD_ADD
         case(MDLP);      nordP = NODES(Mdle)%order+NORD_ADD*11
      end select
!
      call compute_enriched_order(nordP, norderP)
!
!  ...determine edge and face orientations
      call find_orient( Mdle, norient_edge,norient_face)
!                                                                     
!  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
!
!  ...determine element size and scale correction
      call find_hmin(Mdle,h_elem)
!
      omeg = OMEGA
!
!  ...get number of DOFs, enriched DOFs, and bubble DOFs
      call celndof(NODES(Mdle)%ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
      call celndof(NODES(Mdle)%ntype,norderP,nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
      call ndof_nod(ntype,norder(nre+nrf+1), ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
!
!  ...number of H^1 and H(div) interface DOFs
      nrdofHi = nrdofH - ndofHmdl
      nrdofVi = nrdofV - ndofVmdl
!
!  ...number of text and trial dofs
      nrTest  = nrdofHH + nrdofVV
      nrTrial = nrdofHi + nrdofVi + 4*nrdofQ
!
!  ...memory for the matrices
      allocate(BloadHV(nrTest));          BloadHV(:)  = ZERO
      allocate(AP(nrTest*(nrTest+1)/2));  AP(:)       = ZERO
!
!  ...determine solution dof
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
!  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
!
!-----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                               |
!-----------------------------------------------------------------------
!
!  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD
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
!     ...evaluate the L2 solution
         zsolQ(1:4) = ZERO
         do k=1,NrdofQ
            zsolQ(1:4) = zsolQ(1:4) + zdofQ(1:4,k)*shapQ(k)
         enddo
!
!     ...Piola transformation
         zp = zsolQ(1)/rjac
         zu(1:3) = zsolQ(2:4)/rjac
!
!     ...loop through enriched H1 test functions in the enriched space
         do k1=1,nrdofHH
!
!        ...Piola transformation
            q = shapHH(k1)
            dq(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                    + gradHH(2,k1)*dxidx(2,1:3) &
                    + gradHH(3,k1)*dxidx(3,1:3)
!
!     --- Residual ---
!
            BloadHV(k1) = BloadHV(k1)                          &
                        + (ZIMG*OMEGA*zp*q                       &
                        -  zu(1)*dq(1)-zu(2)*dq(2)-zu(3)*dq(3) &
                        -  zf(1)*q)                            &
                        *  weight
!
!     --- Gram matrix ---
!
!        ...H1 test functions
            do k2=k1,nrdofHH
!           ...Piola transformation
               p = shapHH(k2)
               dp(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                       + gradHH(2,k2)*dxidx(2,1:3) &
                       + gradHH(3,k2)*dxidx(3,1:3)

!           ...determine index in triangular format
               k = nk(k1,k2)
!
!           ...accumulate for the gram matrix
               AP(k) = AP(k)                                     &
                     + (dq(1)*dp(1)+dq(2)*dp(2)+dq(3)*dp(3)      &
                     + (omeg**2+ALPHA)*q*p)*weight
            enddo
!        ...H(div) test functions
            do k2=1,nrdofVV
               n2 = nrdofHH + k2
!           ...Piola transformation
               u(1:3) = (dxdxi(1:3,1)*shapVV(1,k2)                 &
                      +  dxdxi(1:3,2)*shapVV(2,k2)                 &
                      +  dxdxi(1:3,3)*shapVV(3,k2))/rjac
               div_u  =  divVV(k2)/rjac
!
               k = nk(k1,n2)
               AP(k) = AP(k)  &
                     + ZIMG*omeg*(dq(1)*u(1)+dq(2)*u(2)+dq(3)*u(3) - q*div_u)*weight
            enddo
!     ...end of loop through H1 test functions
         enddo
!
!     ...loop through enriched H(div) test functions in the enriched space
         do k1 = 1,nrdofVV
            n1 = nrdofHH+k1
!        ...Piola transformation
            v(1:3) = (dxdxi(1:3,1)*shapVV(1,k1)                 &
                   +  dxdxi(1:3,2)*shapVV(2,k1)                 &
                   +  dxdxi(1:3,3)*shapVV(3,k1))/rjac
            div_v  =  divVV(k1)/rjac
!
!     --- Residual ---
!
            BloadHV(n1) = BloadHV(n1) +                                 &
                        + (ZIMG*OMEGA*(zu(1)*v(1)+zu(2)*v(2)+zu(3)*v(3))  &
                        -  zp*div_v                                     &
                        -  zf(2)*v(1)- zf(3)*v(2) - zf(4)*v(3))         &
                        * weight
!
!     --- Gram matrix ---
!
            do k2 = k1,nrdofVV
               n2 = nrdofHH+k2
!
!           ...Piola transformation
               u(1:3) = (dxdxi(1:3,1)*shapVV(1,k2)             &
                      +  dxdxi(1:3,2)*shapVV(2,k2)             &
                      +  dxdxi(1:3,3)*shapVV(3,k2))/rjac
               div_u  =  divVV(k2)/rjac
!
               k = nk(n1,n2)
               AP(k) = AP(k)        &
                     + ((omeg**2+ALPHA)*(v(1)*u(1)+v(2)*u(2)+v(3)*u(3))  &
                     + (div_v*div_u))*weight
            enddo
!     ...end of loop through H(div) test functions
         enddo
!  ...end of loop through integration points
      enddo
!
!-----------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L S                             |
!-----------------------------------------------------------------------
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
         INTEGRATION = NORD_ADD
         call set_2Dint_DPG(ftype,norderf, nint,tloc,wtloc)
         INTEGRATION = 0
!
!     ...loop through integration points
         do l=1,nint
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
!        ...compute approximate trace at the point
            zsolH(1) = ZERO
            do k=1,nrdofHi
               zsolH(1) = zsolH(1) + zdofH(1,k)*shapH(k)
            enddo
            zp = zsolH(1)
!
!        ...compute approximate flux at the point
            zvec(1:3) = ZERO
            do k=1,nrdofVi
               zvec(1:3) = zvec(1:3) + zdofV(1,k)*shapV(1:3,k)
            enddo
            zsolV(1,1:3) = (dxdxi(1:3,1)*zvec(1)                  &
                         +  dxdxi(1:3,2)*zvec(2)                  &
                         +  dxdxi(1:3,3)*zvec(3))/rjac
!
            zun = ZsolV(1,1)*rn(1) + ZsolV(1,2)*rn(2) + ZsolV(1,3)*rn(3)
!
!        ...impedance BC boundary
            if (ibc(ifc,2).eq.3) then
!
!           ...get the boundary source
               call getg(Mdle,x,rn,ibc(ifc,2), zg)
!
!           ...loop through H1 enriched test functions
               do k1=1,nrdofHH
!
!              ...value of the shape function at the point
                  q = shapHH(k1)
!
!              ...accumulate for the load vector
                  BloadHV(k1) = BloadHV(k1) + q*(zp-zg)*weight
               enddo
!
!           ...regular boundary
            else
!           ...loop through enriched H1 test functions
               do k1 = 1,nrdofHH
!
!              ...value of the shape function at the point
                  q  = shapHH(k1)
!
!              ...accumulate for the load vector
                  BloadHV(k1) = BloadHV(k1) + zun*q*weight
!
!           ...end of loop through H1 test functions
               enddo
            endif
!         
!           ...loop through H(div) enriched test functions
            do k1 = 1,nrdofVV
               n1 = nrdofHH+k1
!
!           ...normal component of the test function (Piola transformation at work!)
               vec(1:3) = dxdxi(1:3,1)*shapVV(1,k1)   &
                        + dxdxi(1:3,2)*shapVV(2,k1)   &
                        + dxdxi(1:3,3)*shapVV(3,k1)
               vec(1:3) = vec(1:3)/rjac
               vn = vec(1)*rn(1) + vec(2)*rn(2) + vec(3)*rn(3)
!
               BloadHV(n1) = BloadHV(n1) + zp*vn*weight
!        ...end of loop through H(div) test functions
            enddo
!     ...end of loop through integration points
         enddo
!  ...end of loop through faces
      enddo
!
!-----------------------------------------------------------------------
!
!  ...factorize the test stiffness matrix
      nrTEST = nrdofHH+nrdofVV
      call ZPPTRF('U', nrTEST, AP, info)
      if (info.ne.0) then
         write(*,*) 'elem_residual: AP info = ',info
         stop 1
      endif
!
      call ZTPTRS('U','C','N',nrTEST,1,AP,BloadHV,nrTest,info)
      if (info.ne.0) then
         write(*,*) 'elem: AP ZTPTRS: Mdle,info = ',Mdle,info
         stop 2
      endif
!
      Resid = ZERO
      do i=1,nrTest
         Resid = Resid + conjg(BloadHV(i))*BloadHV(i)
      enddo
!
      Nref_flag = 111
      deallocate(AP,BloadHV)
!
   end subroutine elem_residual

