!-----------------------------------------------------------------------
!> @brief Routine performs fine mesh to coarse mesh projections and
!         p-refinement for an h-ref candidate
!> @param[in] Kref_loc              - h-ref flag of the h-ref candidate
!> @param[in] Mdle                  - Middle node number of the coarse element
!> @param[in] NrdofgQ               - number of L2 dofs for the isotropic p+1 order
!> @param[in] Nord_old              - previous order on the h-ref candidate    
!> @param[in] Nord_glob             - isotropic p+1 order
!> @param[in] Subson_overlap        - which isotropic fine child elements makes up the h-ref candidate
!> @param[in] Base_error            - projection error on the coarse element at the polynomial order
!                                     correspoding to the previous mesh
!> @param[in] Ap                    - Projection matrix
!> @param[in] Zbload                - load vector for the projection.
!> @param[in] Nextract_prev         - extraction vector correspoding to Nord_old
!> @param[in] Weights_fine_store    - weights and inverse of jacobian at each quadrature point
!!                                    for the sons
!> @param[in]  Quad_point_store     - coordinates of each quadrature point for the sons
!> @param[in]  Nint_pp_store        - number of quadrature points for the sons
!> @param[in]  Shap3DQ_fine_store   - values of L2 shape functions for the sons
!> @param[in]  Shap3DQ_coarse_store - values of L2 shape functions for the coarse element
!
!> @param[out] Polyflag             - optimal polynomial order for the given step
!> @date May 2024
!-----------------------------------------------------------------------
subroutine opt_polynomial_search_subson_linear(Kref_loc,Mdle,NrdofgQ,Nord_old,Nord_glob,&
                                              Subson_overlap,Base_error,Ap,Zbload,&
                                              Nextract_prev,Weights_fine_store,&
                                              Quad_point_store,Nint_pp_store,&
                                              Shap3DQ_fine_store,Shap3DQ_coarse_store,&
                                              Polyflag)
!
   use control
   use data_structure3D
   use element_data
   use parametersDPG
!
   implicit none
   integer,    intent(in)  :: Kref_loc
   integer,    intent(in)  :: Mdle
   integer,    intent(in)  :: NrdofgQ
   integer,    intent(in)  :: Nord_old
   integer,    intent(in)  :: Nord_glob
   real(8),    intent(in)  :: Base_error
   real(8),    intent(in)  :: Ap(NrdofgQ,NrdofgQ)
   real(8),    intent(in)  :: Zbload(NrdofgQ,MAXEQNQ)
   integer,    intent(in)  :: Nextract_prev(*)
   integer,    intent(in)  :: Subson_overlap(*)
   real(8),    intent(in)  :: Weights_fine_store(MAXNINT3ADD,2,*)
   real(8),    intent(in)  :: Quad_point_store(3,MAXNINT3ADD,*)
   real(8),    intent(in)  :: Shap3DQ_fine_store(MAXbrickQQ,MAXNINT3ADD,*)
   real(8),    intent(in)  :: Shap3DQ_coarse_store(MAXbrickQQ,MAXNINT3ADD,*)
   integer,    intent(in)  :: Nint_pp_store(*)
!
   integer,    intent(out) :: Polyflag
!
   real(8),    allocatable :: ap_subson(:,:)
   real(8),    allocatable :: zbload_subson(:,:)
   real(8),    allocatable :: bwork_subson(:,:), nextract_save_lvl(:,:)
   integer,    allocatable :: nextract_prev_subson(:)
   integer,    allocatable :: nextract_subson(:)
   integer  :: nrdofmQ,pxm,pym,pzm
   integer  :: max_loc(1)
   integer  :: etype
!..aux variables
   integer  ::  k,k1,k2,k3,j,l,iattr
   integer  ::  Nord_mod, order_add,Nord_prev, nrdof_old
   real(8)  ::  proj_error,rate_p,proj_error_net
   real(8)  ::  rate_max_lvl
   real(8)  ::  error_rates(6)
   integer  ::  poly_flags_lvl(6), poly_flag_chosen,counter
!
   allocate(ap_subson(NrdofgQ,NrdofgQ))
   allocate(zbload_subson(NrdofgQ,MAXEQNQ))
   allocate(bwork_subson(NrdofgQ,MAXEQNQ))
   ap_subson = ZERO
   zbload_subson =ZERO
   bwork_subson = ZERO
!
   ap_subson(1:NrdofgQ,1:NrdofgQ) = Ap(1:NrdofgQ,1:NrdofgQ)
   zbload_subson(1:NrdofgQ,1:MAXEQNQ) = Zbload(1:NrdofgQ,1:MAXEQNQ)
!
   etype = NODES(Mdle)%ntype
   if(etype .eq. MDLB) then
!
!  ...lvl 1: base order is already solved outside
      call ddecode(Nord_old,pxm,pym,pzm)
      nrdofmQ = pxm * pym * pzm
      nrdof_old = nrdofmQ
!
      allocate(nextract_save_lvl(NrdofgQ,3))
      allocate(nextract_prev_subson(NrdofgQ))
!
      nextract_save_lvl    = ZERO
      nextract_prev_subson = ZERO
      nextract_prev_subson(1:nrdofmQ) = Nextract_prev(1:nrdofmQ)
      Nord_mod = Nord_old
!
      counter = 0
      do k = 1,3
         order_add = int(10**(3-k))
         Nord_prev = Nord_mod
         Nord_mod  = Nord_mod + order_add
         call ddecode(Nord_mod,pxm,pym,pzm)
         nrdofmQ = pxm * pym * pzm
         allocate(nextract_subson(nrdofmQ))
         nextract_subson = ZERO
!
         call extraction_vector_new(Nord_prev,Nord_mod,Nord_glob,nrdofmQ,nextract_prev_subson,nextract_subson)
!
         do iattr = 1,NRQVAR
            do l = 1,nrdofmQ
               bwork_subson(l,iattr) = zbload_subson(nextract_subson(l),iattr)/ap_subson(nextract_subson(l),nextract_subson(l))
            enddo
         enddo
         proj_error_net = 0.d0
         do iattr = 1,NRQVAR
            call fine_to_subson_projection_error(Kref_loc,bwork_subson(1:nrdofmQ,iattr),nextract_subson,Subson_overlap,&
                                                nrdofmQ,NrdofgQ,iattr,Mdle,&
                                                Weights_fine_store,Quad_point_store,Nint_pp_store,&
                                                Shap3DQ_fine_store,Shap3DQ_coarse_store,proj_error)
            proj_error_net = proj_error_net + proj_error
         enddo
!
         nextract_save_lvl(1:nrdofmQ,k) = nextract_subson(1:nrdofmQ)
         error_rates(k)    = (Base_error - proj_error_net)/(abs(real(nrdofmQ*NRQVAR,8) - real(nrdof_old*NRQVAR,8)))
         poly_flags_lvl(k) =  Nord_mod
!
         Nord_mod = Nord_old
         counter = counter + 1
         deallocate(nextract_subson)
!
      enddo
!  ...selecting the first lvl candidate
      rate_max_lvl = maxval(error_rates(1:3))
      max_loc      = maxloc(error_rates(1:3))
      poly_flag_chosen = poly_flags_lvl(max_loc(1))
      Nord_mod = poly_flag_chosen
      call ddecode(Nord_mod,pxm,pym,pzm)
      nrdofmQ = pxm * pym * pzm
      nextract_prev_subson = ZERO
      nextract_prev_subson(1:nrdofmQ) = nextract_save_lvl(1:nrdofmQ,max_loc(1))
      nextract_save_lvl = ZERO
!  ...level 2
!
      do k = 1,3
         if(k .ne. max_loc(1)) then
            order_add = int(10**(3-k))
            Nord_prev = Nord_mod
            Nord_mod  = Nord_mod + order_add
            call ddecode(Nord_mod,pxm,pym,pzm)
            nrdofmQ = pxm * pym * pzm
            allocate(nextract_subson(nrdofmQ))
            nextract_subson = ZERO
!
            call extraction_vector_new(Nord_prev,Nord_mod,Nord_glob,nrdofmQ,nextract_prev_subson,nextract_subson)
!
            do iattr = 1,NRQVAR
               do l = 1,nrdofmQ
                  bwork_subson(l,iattr) = zbload_subson(nextract_subson(l),iattr)/ap_subson(nextract_subson(l),nextract_subson(l))
               enddo
            enddo
!
            proj_error_net = 0.d0
            do iattr = 1,NRQVAR
               call fine_to_subson_projection_error(Kref_loc,bwork_subson(1:nrdofmQ,iattr),nextract_subson,Subson_overlap,&
                                                    nrdofmQ,NrdofgQ,iattr,Mdle,&
                                                    Weights_fine_store,Quad_point_store,Nint_pp_store,&
                                                    Shap3DQ_fine_store,Shap3DQ_coarse_store,proj_error)
               proj_error_net = proj_error_net + proj_error
            enddo
!
            nextract_save_lvl(1:nrdofmQ,counter + 1 - 3) = nextract_subson(1:nrdofmQ)
            error_rates(counter + 1)    = (Base_error - proj_error_net)/(abs(real(nrdofmQ*NRQVAR,8) - real(nrdof_old*NRQVAR,8)))
            poly_flags_lvl(counter + 1) = Nord_mod
            Nord_mod = poly_flag_chosen
            counter = counter + 1
            deallocate(nextract_subson)
!
         endif
      enddo
!
!  ...selecting the second level candidate
      rate_max_lvl = maxval(error_rates(4:5))
      max_loc      = maxloc(error_rates(4:5))
      poly_flag_chosen = poly_flags_lvl(3 + max_loc(1))
!  ...final stage
      Nord_prev = poly_flag_chosen
      call ddecode(Nord_prev,pxm,pym,pzm)
      nrdofmQ = pxm * pym * pzm
      nextract_prev_subson = ZERO
      nextract_prev_subson(1:nrdofmQ) = nextract_save_lvl(1:nrdofmQ,max_loc(1))
      Nord_mod = Nord_old + 111
      call ddecode(Nord_mod,pxm,pym,pzm)
      nrdofmQ = pxm * pym * pzm
      allocate(nextract_subson(nrdofmQ))
      nextract_subson = ZERO
!
      call extraction_vector_new(Nord_prev,Nord_mod,Nord_glob,nrdofmQ,nextract_prev_subson,nextract_subson)

      do iattr = 1,NRQVAR
         do l = 1,nrdofmQ
            bwork_subson(l,iattr) = zbload_subson(nextract_subson(l),iattr)/ap_subson(nextract_subson(l),nextract_subson(l))
         enddo
      enddo
!
      proj_error_net = 0.d0
      do iattr = 1,NRQVAR
         call fine_to_subson_projection_error(Kref_loc,bwork_subson(1:nrdofmQ,iattr),nextract_subson,Subson_overlap,&
                                              nrdofmQ,NrdofgQ,iattr,Mdle,&
                                              Weights_fine_store,Quad_point_store,Nint_pp_store,&
                                              Shap3DQ_fine_store,Shap3DQ_coarse_store,proj_error)
         proj_error_net = proj_error_net + proj_error
      enddo
!
      error_rates(counter + 1)    = (Base_error - proj_error_net)/(abs(real(nrdofmQ*NRQVAR,8) - real(nrdof_old*NRQVAR,8)))
      poly_flags_lvl(counter + 1) = Nord_mod
      rate_p    = maxval(error_rates)
      max_loc   = maxloc(error_rates)
      Polyflag = poly_flags_lvl(max_loc (1))
      deallocate(nextract_subson)
!
      if(Nord_glob .lt. Polyflag) write(*,*) "Error = ",Polyflag, Nord_glob,max_loc(1)
    endif
!
end subroutine opt_polynomial_search_subson_linear
