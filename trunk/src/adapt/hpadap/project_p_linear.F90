!-----------------------------------------------------------------------
!> @brief Routine performs fine mesh to coarse mesh projections and p-refinement choice
!> @param[in]   Mdle           - index for element selected for refinement on the coarse mesh.
!> @param[in]   Flag_pref_loc  - flag stating that the element is p-refined or not during fine mesh formation
!> @param[out]  Error_org      - projection error for the coarse element
!> @param[out]  Rate_p         - maximum rate obtained via the anisotropic p-refinement.
!> @param[out]  Poly_flag      - polynomial refinement flag
!> @date May 2024
!-----------------------------------------------------------------------
subroutine project_p_linear(Mdle,Flag_pref_loc, Error_org,Rate_p,Poly_flag)
!
   use control
   use data_structure3D
   use element_data
   use parametersDPG
   use physics
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(in)  :: Flag_pref_loc
!
   real(8), intent(out) :: Error_org
   real(8), intent(out) :: Rate_p
   integer, intent(out) :: Poly_flag
!
   real(8), allocatable    :: ap(:)
   real(8), allocatable    :: zbload(:,:)
   integer, allocatable    :: mdle_sons(:)
!..stores the dofs corresponding to the interpolant on original element after projection from fine to coarse element.
   real(8) :: zdofQ(MAXEQNQ,MAXbrickQQ)
!..stores L2 dofs of the fine grid elements (sons of the marked element after the isotropic hp refinement).
   real(8) :: zdofQ_pp(MAXEQNQ,MAXbrickQQ)
!..stores the values of the solution for the original and sons.
   real(8) :: zvalQ(MAXEQNQ), zvalQpp(MAXEQNQ)
!
!..variables required for coarse element
!..number of shape functions in trial space for original element (in case of ultra weak formulation: only nrdofQ is useful)
   integer :: nrdofH,nrdofE,nrdofV,nrdofQ
!..number of shape functions in trial space for enriched sons (its same for all the sons)
   integer :: nrdofH_pp,nrdofE_pp,nrdofV_pp,nrdofQ_pp
!
!..L2 shape functions
   real(8) :: shapQ(MAXbrickQQ)
!..H1 shape functions
   real(8) :: shapH(MAXbrickHH), gradH(3,MAXbrickHH)
!..3D quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
!..geometry
   real(8) :: xi(3), x(3), dxdxi(3,3), dxidx(3,3), rn(3),xis(3)
   real(8) :: dxidt(3,2), dxdt(3,2), t(2)
!..geometry dof of sons
   real(8) :: xnod_pp(3,MAXbrickHH)
   integer :: etype
!..element order, orientation for edges and faces of the original element
   integer :: norder(19), norient_edge(12), norient_face(6)
   integer :: nint
!..element order, orientation for edges and faces for enriched sons
   integer :: norder_pp(19), norient_edge_pp(12), norient_face_pp(6)
   integer :: nint_pp, nrdof
   integer :: nord_glob
!..auxiliary variables
   real(8) :: wa,weight,rjac
   integer :: k1,k2,l,iflag,k,is,mdle_fine,j, first_son 
   integer :: iattr,icomp,ibeg
!..location of L2 variable in L2 solution array
   integer :: var_loc
!..number of sons of a refined element
   integer :: nr_mdle_sons 
   real(8) :: q1,q2,q
   real(8), allocatable :: weights_fine_store(:,:,:), quad_point_store(:,:,:)
   real(8), allocatable :: shap3DQ_fine_store(:,:,:),shap3DQ_coarse_store(:,:,:)
   integer, allocatable :: nint_pp_store(:)
!..options for polynomial selection process
   integer :: opt_poly_choice = 1
!
!..allocating memory to the p+1 order projection matrix
   etype = NODES(Mdle)%ntype
   call find_order(Mdle, norder)
   call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!
   allocate(ap(nrdofQ))
   allocate(zbload(nrdofQ,MAXEQNQ))
!
   if(etype .eq. MDLB) then
      nr_mdle_sons = 8
   endif
!
   allocate(mdle_sons(nr_mdle_sons))
   first_son = NODES(Mdle)%first_son
   do is = 1,nr_mdle_sons
      mdle_sons(is) = first_son + is - 1
   enddo
!
!.. looping over sons to compute the load vector of projection from fine to coarse
   ap = ZERO
   zbload = ZERO
!..allocating memory to arrays that will contain the quadrature and L2 shape function data
!  for the sons. We are storing these values so that we dont have to compute them again while
!  performing the fine to coarse projection.
   allocate(weights_fine_store(MAXNINT3ADD,2,nr_mdle_sons))
   allocate(quad_point_store(3,MAXNINT3ADD,nr_mdle_sons))
   allocate(nint_pp_store(nr_mdle_sons))
   allocate(shap3DQ_fine_store(MAXbrickQQ,MAXNINT3ADD,nr_mdle_sons))
   allocate(shap3DQ_coarse_store(MAXbrickQQ,MAXNINT3ADD,nr_mdle_sons))
!
   weights_fine_store = ZERO
   quad_point_store = ZERO
   nint_pp_store  =ZERO
   shap3DQ_coarse_store = ZERO
   shap3DQ_fine_store = ZERO
!
!..looping over the sons of the marked elements to perform the projection
   do is = 1,nr_mdle_sons
      mdle_fine = mdle_sons(is)
      etype = NODES(mdle_fine)%ntype
      call find_order(mdle_fine, norder_pp)
      call celndof(etype,norder_pp, nrdofH_pp,nrdofE_pp,nrdofV_pp,nrdofQ_pp)
      call find_orient(mdle_fine, norient_edge_pp,norient_face_pp)
      call nodcor(mdle_fine, xnod_pp)
      INTEGRATION = 1
      call set_3D_int_DPG(etype,norder_pp,norient_face_pp, nint_pp,xiloc,waloc)
      nint_pp_store(is) = nint_pp
!
!  ...extract the L2 dofs of the fine grid solution for the is^th son
      call solelm_L2(mdle_fine,zdofQ_pp)
      do l = 1,nint_pp
!     ...coordinates and weight of this integration point
         xi(1:3)=xiloc(1:3,l)
         quad_point_store(1:3,l,is) = xi(1:3)
         wa=waloc(l)
!     ...H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder_pp,norient_edge_pp,norient_face_pp, nrdof,shapH,gradH)
!     ...L2 shape function calls
         call shape3DQ(etype,xi,norder_pp, nrdof,shapQ)
         shap3DQ_fine_store(1:nrdofQ_pp,l,is) = shapQ(1:nrdofQ_pp)
!     ...geometry map
         call geom3D(mdle_fine,xi,xnod_pp,shapH,gradH,nrdofH_pp, x,dxdxi,dxidx,rjac,iflag)
         weight = rjac*wa
!     ...storing the weight and inverse of jacobian for sons
         weights_fine_store(l,1,is) = weight
         weights_fine_store(l,2,is) = rjac
!
         zvalQpp = ZERO
         do iattr = 1,NR_PHYSA
               if(D_TYPE(iattr) .eq. DISCON) then
                  ibeg = ADRES(iattr)
                  do icomp = 1, NR_COMP(iattr)
                     var_loc = ibeg + icomp
                     do k = 1,nrdofQ_pp
                           q = shapQ(k)/rjac
                           zvalQpp(var_loc) = zvalQpp(var_loc) + zdofQ_pp(var_loc,k) * q
                     enddo
                  enddo
               endif
         enddo
!     ...calling the map between son's master element and coarse element master element
         call fine_to_coarse_gp_map(is,xi,xis,etype)
         shapQ = ZERO
         call shape3DQ(etype,xis,norder_pp, nrdof,shapQ)
         shap3DQ_coarse_store(1:nrdofQ_pp,l,is) = shapQ(1:nrdofQ_pp)
!     ...scaling the jacobian to account for the isotropic refinement of coarse element
         rjac = rjac * real(nr_mdle_sons,8)
!
!     ...contribution to p+1 coarse projection matrix using subelement integration
!        for linear elements, the projection matrix is diagonal, hence,only storing
!     ...diagonal values
         do k1 = 1,nrdofQ_pp
               q1 = shapQ(k1)/rjac
               q2 = shapQ(k1)/rjac
               ap(k1) = ap(k1) + weight * q1 * q2
         enddo
!     ...computing  contributions to the projection load vector from each son (only for L2 variables)
         do iattr = 1,NR_PHYSA
               if(D_TYPE(iattr) .eq. DISCON) then
                  ibeg = ADRES(iattr)
                  do icomp = 1,NR_COMP(iattr)
                     var_loc = ibeg + icomp
                     do k = 1,nrdofQ_pp
                           q = shapQ(k)/rjac
                           zbload(k,var_loc) = zbload(k,var_loc) + weight *  q * zvalQpp(var_loc)
                     enddo
                  enddo
               endif
         enddo
      enddo
   enddo
!
   nord_glob = norder_pp(19)
!..call to the subroutine to perform the optimal p-only refinement on the coarse element
   call opt_polynomial_search_coarse_linear(Mdle,nr_mdle_sons,mdle_sons,nrdofQ,Flag_pref_loc,nord_glob, ap,zbload, &
                                                weights_fine_store,quad_point_store,shap3DQ_fine_store, &
                                                shap3DQ_coarse_store,nint_pp_store,&
                                                Error_org,Poly_flag,Rate_p)
!
!..deallocating the allocatable arrays
   deallocate(ap)
   deallocate(zbload)
   deallocate(mdle_sons)
   deallocate(weights_fine_store)
   deallocate(quad_point_store)
   deallocate(nint_pp_store)
   deallocate(shap3DQ_fine_store)
   deallocate(shap3DQ_coarse_store)
!
end subroutine project_p_linear
