!-----------------------------------------------------------------------
!> @brief Routine performs fine mesh to a h-refinement candidate element and computes
!!        optimal p-distribution for the h-ref candidate
!> @param[in]   Mdle            - middle node number of the coarse element
!> @param[in]   Flag_pref_loc   - flag stating that the element is p-refined or not during fine mesh formation
!> @param[in]   Kref_loc        - refinement flag for the h-ref candidate
!> @param[in]   Error_org       - projection error for the coarse element
!> @param[out]  Error_opt       - error for polynomial order that provides maximum error reduction rate
!> @param[out]  G_rate_max      - maximum error reduction rate
!> @param[out]  Rate_hcomp      - maximum error reduction rate for the competitive h-refinements
!> @param[out]  Nord_href       - polynomial order for guranteed rate
!> @param[out]  Nord_threshold  - polynomial distribution during the traversal along maximum error reduction path
!> @param[out]  Count_threshold - number of steps in maximum error reduction path, once the refinements becomes competitive
!> @param[out]  error_rate      - error rates while traversing the maximum error reduction path
!> @param[out]  Loc_max_rate    - index of the step which produced the best error reduction rate
!> @date May 2024
!-----------------------------------------------------------------------
subroutine elem_proj_h_linear(Mdle,Flag_pref_loc,Kref_loc,Error_org, Error_opt,&
                              G_rate_max,Rate_hcomp,Nord_href,&
                              Nord_threshold,Count_threshold,Err_reduction_rate,&
                              Loc_max_rate)
!
   use control
   use data_structure3D
   use element_data
   use parametersDPG
   use MPI 
!
   implicit none
!
   integer,    intent(in)  :: Mdle
   integer,    intent(in)  :: Flag_pref_loc
   integer,    intent(in)  :: Kref_loc
   real(8),    intent(in)  :: Error_org
!
   real(8),    intent(out) :: Error_opt
   real(8),    intent(out) :: G_rate_max
   real(8),    intent(out) :: Rate_hcomp
   integer,    intent(out) :: Nord_href(8)
   integer,    intent(out) :: Nord_threshold(100,8)
   integer,    intent(out) :: Count_threshold
   real(8),    intent(out) :: Err_reduction_rate(100)
   integer,    intent(out) :: Loc_max_rate
!
   integer  :: etype
!..element order, orientation for edges and faces of the original element
   integer  :: norder(19)
   integer  :: norder_pp(19), norient_edge_pp(12), norient_face_pp(6)
!.. number of shape functions in trial space for enriched sons (its same for all of them)
   integer  :: nrdofH_pp, nrdofE_pp, nrdofV_pp, nrdofQ_pp  
!.. L2 shape functions
   real(8)  :: shapQ(MAXbrickQQ)
!..H1 shape functions
   real(8)  :: shapH(MAXbrickHH), gradH(3,MAXbrickHH)
!..3D quadrature data
   real(8)  :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
!..geometry
   real(8)  :: xi(3), x(3), dxdxi(3,3), dxidx(3,3),xis(3)
!..geometry dof of sons
   real(8)  :: xnod_pp(3,MAXbrickHH)
!..to store extracted coefficeints from 
   real(8)  :: zdofQ_pp(MAXEQNQ,MAXbrickQQ), zdofH_pp(MAXEQNH,MAXbrickHH),zdofE_pp(MAXEQNE,MAXbrickEE), &
               zdofV_pp(MAXEQNV,MAXbrickVV)
!
!..to store the values of the solution for the original and sons respectively
   real(8)  :: zvalQpp(MAXEQNQ)
   integer  :: nr_subsons,hx,hy,hz
   integer  :: pxm,pym,pzm ! max poly order in each direction
!..projection matrix for the child of the coarse elements 
!  produced by anisotropic h ref
   real(8),    allocatable    :: subsons_ap(:,:,:)
!..load vector for the projection problem for the child of the coarse element
!  produced by anisotropic h ref
   real(8),    allocatable    :: subsons_zbload(:,:,:) 
!..stores the overlap of the fine grid with subsons
   integer,    allocatable    :: subsons_overlap(:,:)
!
   integer,    allocatable    :: mdle_sons(:)
   integer,    allocatable    :: nord_org_subsons(:)
   integer,    allocatable    :: nextract(:)
!
   real(8),    allocatable    :: awork(:,:),bwork(:,:),ap(:,:),zbload(:,:)
   real(8),    allocatable    :: subsons_awork(:,:,:),subsons_bwork(:,:,:)
!
   real(8)  :: proj_error_subson, proj_error, proj_error_net_subson
   integer  :: dof_diff(100)
!
   integer  :: poly_dist_array(100,8)
   real(8)  :: error_rate(100)
!
   integer :: nrdofgQ,nrdofmQ, nrdof_org
   integer :: nr_mdle_sons,mdle_fine,first_son
   integer :: nint_pp,nrdof,nord_glob, nord_org, polyflag
   integer, allocatable    :: nord_mep(:,:),nord_old(:)
   real(8), allocatable    :: error_subsons(:)
!     
   integer, allocatable    :: nord_max(:)
!..shape functions for the load vector computation in the projection problem
   real(8) :: q1,q2,q 
   real(8) :: max_error_subson, ratio_mep
!..store quadrature data for fine child elements for error computation
   real(8), allocatable :: weights_fine_store(:,:,:,:),quad_point_store(:,:,:,:)
   real(8), allocatable :: shap3DQ_fine_store(:,:,:,:),shap3DQ_coarse_store(:,:,:,:)
   integer, allocatable :: nint_pp_store(:,:) 
!..auxiliary variables
   integer :: j,iss,is,l,k,iflag,k1,k2,iso_p,nrdof_tmp
   integer :: Nref, local_order_check
   integer :: pxc,pyc,pzc, pxg,pyg,pzg
   integer :: iattr,icomp,ibeg
   integer :: var_loc
   real(8) :: wa, weight,rjac,g_rate_tmp
   integer,    allocatable :: subsons_nextract_prev(:,:)
!
!..options for polynomial selection process
   integer :: opt_poly_choice = 1
!
   etype = NODES(Mdle)%ntype
   call find_order(Mdle, norder)
!
   nord_glob = norder(19) !interior order for L2 variables
   if(etype .eq. MDLB) then
!
      if(Flag_pref_loc .eq. 1) then
         nord_org = nord_glob - 111
      else
         nord_org = nord_glob
      endif
!  ...dofs on the coarse element before hp refinement
      call ddecode(nord_org,pxm,pym,pzm)
      nrdof_org = pxm * pym * pzm
!  ...number of fine childs
      nr_mdle_sons = 8
!
      allocate(mdle_sons(nr_mdle_sons))
      call ddecode(nord_glob,pxg,pyg,pzg)
      nrdofgQ = pxg * pyg *pzg
      iso_p = MAX(pxg,pyg,pzg)
!
      call ddecode(Kref_loc,hx,hy,hz)
      nr_subsons = 2**(hx+hy+hz)
!
      allocate(subsons_ap(nrdofgQ,nrdofgQ,nr_subsons))
      subsons_ap = ZERO
      allocate(subsons_zbload(nrdofgQ,MAXEQNQ,nr_subsons))
      subsons_zbload = ZERO
      allocate(subsons_overlap(nr_subsons,8/nr_subsons))
      subsons_overlap = ZERO
!
      select case(Kref_loc)
         case(100)
            subsons_overlap(1,:) = (/1,4,5,8/)
            subsons_overlap(2,:) = (/2,3,6,7/)
         case(010)
            subsons_overlap(1,:) = (/1,2,5,6/)
            subsons_overlap(2,:) = (/3,4,7,8/)
         case(001)
            subsons_overlap(1,:) = (/1,2,3,4/)
            subsons_overlap(2,:) = (/5,6,7,8/)
         case(110)
            subsons_overlap(1,:) = (/1,5/)
            subsons_overlap(2,:) = (/2,6/)
            subsons_overlap(3,:) = (/3,7/)
            subsons_overlap(4,:) = (/4,8/)
         case(101)
            subsons_overlap(1,:) = (/1,4/)
            subsons_overlap(2,:) = (/2,3/)
            subsons_overlap(3,:) = (/6,7/)
            subsons_overlap(4,:) = (/5,8/)
         case(011)
            subsons_overlap(1,:) = (/1,2/)
            subsons_overlap(2,:) = (/3,4/)
            subsons_overlap(3,:) = (/7,8/)
            subsons_overlap(4,:) = (/5,6/)
         case(111)
            do j = 1,nr_subsons
               subsons_overlap(j,1) = j
            enddo
         case default
            write(*,*) " Not a valid kref"
      end select
!
      first_son = NODES(Mdle)%first_son
      do is = 1,nr_mdle_sons
         mdle_sons(is) = first_son + is - 1
      enddo
!     
      allocate(weights_fine_store(MAXNINT3ADD,2,8/nr_subsons,nr_subsons))
      allocate(quad_point_store(3,MAXNINT3ADD,8/nr_subsons,nr_subsons))
      allocate(nint_pp_store(8/nr_subsons,nr_subsons))
      allocate(shap3DQ_coarse_store(MAXbrickQQ,MAXNINT3ADD,8/nr_subsons,nr_subsons))
      allocate(shap3DQ_fine_store(MAXbrickQQ,MAXNINT3ADD,8/nr_subsons,nr_subsons))
!
      weights_fine_store = ZERO
      quad_point_store = ZERO
      nint_pp_store  =ZERO
      shap3DQ_coarse_store = ZERO
      shap3DQ_fine_store = ZERO
!
      do iss = 1,nr_subsons
         do is = 1,8/nr_subsons
!
            mdle_fine = mdle_sons(subsons_overlap(iss,is))
            etype = NODES(mdle_fine)%ntype
!
            call find_order(mdle_fine, norder_pp)
            call celndof(etype,norder_pp, nrdofH_pp,nrdofE_pp,nrdofV_pp,nrdofQ_pp)
!
            call find_orient(mdle_fine, norient_edge_pp,norient_face_pp)
            call nodcor(mdle_fine, xnod_pp)
            INTEGRATION = 1
            call set_3D_int_DPG(etype,norder_pp,norient_face_pp, nint_pp,xiloc,waloc)
            nint_pp_store(is,iss) = nint_pp
!        ...extract the coefficeints of the fine grid solution for the is^th son
            call solelm(mdle_fine,zdofH_pp,zdofE_pp,zdofV_pp,zdofQ_pp)
!
            do l = 1,nint_pp

               xi(1:3) = xiloc(1:3,l)
               quad_point_store(1:3,l,is,iss) = xi(1:3)
               wa = waloc(l)
                           !  ...H1 shape functions (for geometry)
               call shape3DH(etype,xi,norder_pp,norient_edge_pp,norient_face_pp, nrdof,shapH,gradH)
               !  ...L2 shape function calls
               call shape3DQ(etype,xi,norder_pp, nrdof,shapQ)
               shap3DQ_fine_store(1:nrdofQ_pp,l,is,iss) = shapQ(1:nrdofQ_pp)
               !  ...geometry map
               call geom3D(mdle_fine,xi,xnod_pp,shapH,gradH,nrdofH_pp, x,dxdxi,dxidx,rjac,iflag)

               weight = rjac*wa
               weights_fine_store(l,1,is,iss) = weight
               weights_fine_store(l,2,is,iss) = rjac
               zvalQpp = ZERO
!
               do iattr = 1,NR_PHYSA

                  if(D_TYPE(iattr).eq. DISCON) then
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
!           ...calling the map between iso-8 ref master element and candidate h-ref element's master element
               call fine_to_subson_gp_map(etype,Kref_loc,subsons_overlap(iss,is),xi,xis)
               shapQ = ZERO
               call shape3DQ(etype,xis,norder_pp,nrdof,shapQ)
               shap3DQ_coarse_store(1:nrdofQ_pp,l,is,iss) = shapQ(1:nrdofQ_pp)
!           ...scaling the jacobian for isotropic refinement of coarse element
               rjac = rjac * real(nr_mdle_sons/nr_subsons,8)
!
               do k1 = 1,nrdofQ_pp
                  q1 = shapQ(k1)/rjac
                  do k2 = 1,nrdofQ_pp
                        q2 = shapq(k2)/rjac
                        subsons_ap(k1,k2,iss) = subsons_ap(k1,k2,iss) + weight * q1 * q2
                  enddo
               enddo
!
               do iattr = 1,NR_PHYSA
                  if(D_TYPE(iattr).eq. DISCON) then
                        ibeg = ADRES(iattr)
                        do icomp = 1,NR_COMP(iattr)
                           var_loc = ibeg + icomp
                           do k = 1,nrdofQ_pp
                              q = shapQ(k)/rjac
                              subsons_zbload(k,var_loc,iss) = subsons_zbload(k,var_loc,iss) + weight * q * zvalQpp(var_loc)
                           enddo
                        enddo
                  endif
               enddo
            enddo
         enddo
      enddo
!
!  ...now solving the projection problem starting from base order on each subson by looping over them
      allocate(ap(nrdofgQ,nrdofgQ))
      allocate(zbload(nrdofgQ,MAXEQNQ))
!
      proj_error = 0.d0
      proj_error_subson = 0.d0
!
      allocate(nord_org_subsons(nr_subsons))
      allocate(nord_mep(nr_subsons,nr_subsons * (3 * iso_p)))
      allocate(nord_old(nr_subsons))
      allocate(nord_max(nr_subsons))
      allocate(error_subsons(nr_subsons))
      allocate(subsons_nextract_prev(nr_subsons,nrdofgQ))
      allocate(subsons_awork(nrdofgQ,nrdofgQ,nr_subsons))
      allocate(subsons_bwork(nrdofgQ,MAXEQNQ,nr_subsons))
!
      nord_org_subsons = ZERO
      nord_mep = ZERO
      nord_mep(1:nr_subsons,1) = 222
      nord_old = ZERO
      nord_max = ZERO
      error_subsons   = ZERO
      subsons_nextract_prev = ZERO
      subsons_awork = ZERO
      subsons_bwork = ZERO
!
      dof_diff = ZERO
      error_rate = ZERO
!
      poly_dist_array = ZERO
      Count_threshold = ZERO
      Err_reduction_rate = ZERO
      Nord_threshold = ZERO
      Loc_max_rate = ZERO
!           
      Nref = 1
!  ...ratio for selecting subsons to p-refine  
      ratio_mep = 0.7d0
!  ...maximum guranteed rate
      G_rate_max = 0.d0 
      do
!  ...solve over all subsons
         do iss = 1,nr_subsons
!        ...only solve if the subson has order changed
            if(nord_old(iss) .ne. nord_mep(iss,Nref)) then

               allocate(awork(nrdofgQ,nrdofgQ))
               allocate(bwork(nrdofgQ,MAXEQNQ))
               awork = ZERO
               bwork = ZERO
!           ...copying the projection matrix
               awork(1:nrdofgQ,1:nrdofgQ)  = subsons_awork(1:nrdofgQ,1:nrdofgQ,iss)
               bwork(1:nrdofgQ,1:MAXEQNQ)  = subsons_bwork(1:nrdofgQ,1:MAXEQNQ,iss)
!
               nord_org_subsons(iss) = nord_mep(iss,Nref)
               call ddecode(nord_org_subsons(iss),pxm,pym,pzm)
               nrdofmQ = pxm * pym * pzm
!           ...currently useless(nstep and Mblock variables are needed for pbi solver)
!
               allocate(nextract(nrdofmQ))
               nextract = ZERO        
               ap(1:nrdofgQ,1:nrdofgQ) = subsons_ap(1:nrdofgQ,1:nrdofgQ,iss)
               zbload(1:nrdofgQ,1:MAXEQNQ) = subsons_zbload(1:nrdofgQ,1:MAXEQNQ,iss)
!
               call extraction_vector_new(nord_old(iss),nord_org_subsons(iss),nord_glob,nrdofmQ,subsons_nextract_prev(iss,:),nextract)
               do iattr = 1,NRQVAR
                  do l = 1,nrdofmQ
                        bwork(l,iattr) = zbload(nextract(l),iattr)/ap(nextract(l),nextract(l))
                  enddo
               enddo
!
               subsons_nextract_prev(iss,1:nrdofmQ) = nextract(1:nrdofmQ)
               proj_error_subson = 0.d0
               proj_error_net_subson = 0.d0
!
               do iattr = 1,NRQVAR
                  call fine_to_subson_projection_error(Kref_loc,bwork(1:nrdofmQ,iattr),nextract,subsons_overlap(iss,:),&
                                                      nrdofmQ,nrdofgQ,iattr,Mdle,&
                                                      weights_fine_store(:,:,:,iss),quad_point_store(:,:,:,iss),nint_pp_store(:,iss),&
                                                      shap3DQ_fine_store(:,:,:,iss),shap3DQ_coarse_store(:,:,:,iss),proj_error_subson)
!
                  proj_error_net_subson = proj_error_net_subson + proj_error_subson
               enddo
!
               error_subsons(iss) = proj_error_net_subson
               subsons_awork(1:nrdofgQ,1:nrdofgQ,iss) = awork(1:nrdofgQ,1:nrdofgQ)
               subsons_bwork(1:nrdofgQ,1:MAXEQNQ,iss) = bwork(1:nrdofgQ,1:MAXEQNQ)
!
               deallocate(awork)
               deallocate(bwork)
               deallocate(nextract)
!
!
            endif
         enddo
         proj_error = sum(error_subsons)
!     ...conditional max
         max_error_subson = 0.d0
         nord_old(1:nr_subsons) = nord_mep(1:nr_subsons,Nref)
         nrdof_tmp = 0
!
         do iss = 1,nr_subsons
               if((error_subsons(iss) .ge. max_error_subson)   .and. (nord_mep(iss,Nref) .lt. nord_glob)) then
                  max_error_subson = error_subsons(iss)
               endif
               call ddecode(nord_old(iss),pxm,pym,pzm)
               nrdof_tmp = nrdof_tmp + pxm*pym*pzm
!
         enddo
!
         dof_diff(Nref) = nrdof_tmp - nrdof_org
!
!     ...increasing the order for max error in subsons with order  .le. nord_glob
         do iss = 1,nr_subsons
            nord_mep(iss,Nref + 1) = nord_mep(iss,Nref)
            if((error_subsons(iss) .ge. ratio_mep * max_error_subson) .and. (nord_mep(iss,Nref) .lt. nord_glob)) then
!
               call ddecode(nord_mep(iss,Nref),pxc,pyc,pzc)
               if(((pxg - pxc) .ge. 1) .and. ((pyg - pyc) .ge. 1) .and. ((pzg - pzc) .ge. 1)) then
!
                  allocate(awork(nrdofgQ,nrdofgQ))
                  allocate(bwork(nrdofgQ,MAXEQNQ))
                  nrdofmQ = pxc * pyc * pzc
!
                  awork = ZERO
                  bwork = ZERO
!
                  awork(1:nrdofgQ,1:nrdofgQ)           = subsons_awork(1:nrdofgQ,1:nrdofgQ,iss)
                  bwork(1:nrdofgQ,1:MAXEQNQ)           = subsons_bwork(1:nrdofgQ,1:MAXEQNQ,iss)
!
                  ap(1:nrdofgQ,1:nrdofgQ) = subsons_ap(1:nrdofgQ,1:nrdofgQ,iss)
                  zbload(1:nrdofgQ,1:MAXEQNQ) = subsons_zbload(1:nrdofgQ,1:MAXEQNQ,iss)
                  if(opt_poly_choice .eq. 1) then
                        call opt_polynomial_search_subson_linear(Kref_loc,Mdle,nrdofgQ,nord_mep(iss,Nref),nord_glob,subsons_overlap(iss,:),error_subsons(iss), &
                                                                ap,zbload,subsons_nextract_prev(iss,1:nrdofmQ),&
                                                                weights_fine_store(:,:,:,iss),quad_point_store(:,:,:,iss),nint_pp_store(:,iss),&
                                                                shap3DQ_fine_store(:,:,:,iss),shap3DQ_coarse_store(:,:,:,iss),polyflag)
                  endif             
                  nord_mep(iss, Nref + 1) = polyflag
                  deallocate(awork)
                  deallocate(bwork)
               else
                  nord_mep(iss, Nref + 1) = nord_glob
               endif
            endif
         enddo
!
         if(nrdof_tmp .gt. nrdof_org) then
               g_rate_tmp = (Error_org - proj_error)/abs(real(nrdof_tmp*NRQVAR,8) - real(nrdof_org*NRQVAR,8))
               error_rate(Nref) = g_rate_tmp  
               poly_dist_array(Nref,1:nr_subsons) = nord_old(1:nr_subsons)
               Count_threshold = Nref
         endif
!
         if(g_rate_tmp .ge. G_rate_max) then
               G_rate_max = g_rate_tmp
               Error_opt = proj_error
               nord_max(1:nr_subsons) = nord_old(1:nr_subsons)
         endif
         local_order_check = 0
         do iss = 1,nr_subsons
               call ddecode(nord_mep(iss,Nref+1),pxc,pyc,pzc)
               if((pxc .gt. pxg) .or. (pyc .gt. pyg) .or. (pzc .gt. pzg)) then
                  local_order_check = local_order_check + 1
               endif
         enddo
!
         if(local_order_check .gt. 0) then
               exit
         endif
!
         local_order_check = 0
         do iss = 1,nr_subsons
               call ddecode(nord_mep(iss,Nref),pxc,pyc,pzc)
               if((pxc .eq. pxg) .and. (pyc .eq. pyg) .and. (pzc .eq. pzg)) then
                  local_order_check = local_order_check + 1
               endif
         enddo
         if(local_order_check .eq. nr_subsons) then
               exit
         endif
! 
         Nref = Nref + 1
      enddo

      Nord_href(1:nr_subsons) = nord_max(1:nr_subsons)
      Err_reduction_rate(1:100) = error_rate(1:100)
      Nord_threshold(1:100,1:8) = poly_dist_array(1:100,1:8)
      Rate_hcomp = 0.d0
!
      do k = 1,Nref
         if((dof_diff(k)) .gt. 0) then
                  if(Rate_hcomp .le. error_rate(k)) then
                     Rate_hcomp = error_rate(k)
                     Loc_max_rate = k
                  endif
         endif
      enddo
!  ...deallocating allocatable arrays
      deallocate(subsons_ap)
      deallocate(subsons_zbload)
      deallocate(subsons_overlap)
      deallocate(mdle_sons)
      deallocate(nord_org_subsons)
      deallocate(ap)
      deallocate(zbload)
      deallocate(subsons_awork)
      deallocate(subsons_bwork)
      deallocate(nord_mep)
      deallocate(nord_old)
      deallocate(error_subsons)
      deallocate(nord_max)
      deallocate(weights_fine_store)
      deallocate(quad_point_store)
      deallocate(shap3DQ_fine_store)
      deallocate(shap3DQ_coarse_store)
      deallocate(nint_pp_store)
      deallocate(subsons_nextract_prev)
!
   endif
!
end subroutine elem_proj_h_linear