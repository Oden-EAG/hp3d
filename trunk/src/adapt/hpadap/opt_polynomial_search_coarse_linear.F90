!-----------------------------------------------------------------------
!> @brief Routine performs fine mesh to coarse mesh projections and p-refinement choice
!> @param[in]  Mdle                 - middle node number of the coarse element
!> @param[in]  Nr_mdle_sons         - number of sons after isotropic refinement
!> @param[in]  Mdle_sons            - middle node numbers of the sons
!> @param[in]  NrdofQ               - number of L2 dofs for the coarse element
!> @param[in]  Nord_glob            - polynomial order of the coarse element
!> @param[in]  Ap                   - Diagonal of the projection matrix
!> @param[in]  Zbload               - load vector for the projection problem
!> @param[in]  Weights_fine_store   - weights and inverse of jacobian at each quadrature point
!!                                    for the sons
!> @param[in]  Quad_point_store     - coordinates of each quadrature point for the sons
!> @param[in]  Shap3DQ_fine_store   - values of L2 shape functions for the sons
!> @param[in]  Shap3DQ_coarse_store - values of L2 shape functions for the coarse element
!> @param[in]  Nint_pp_store        - number of quadrature points for the sons
!> @param[out] Error_org            - projection error for the coarse element
!> @param[out] Polyflag             - optimal polynomial order chosen for p-only refinement
!> @param[out] rate_p               - rate of error decrement for p-only refinement
!> @date May 2024
!-----------------------------------------------------------------------
subroutine opt_polynomial_search_coarse_linear(Mdle,Nr_mdle_sons,Mdle_sons,NrdofQ,Flag_pref_loc,&
                                              Nord_glob, Ap,Zbload,Weights_fine_store,&
                                              Quad_point_store,Shap3DQ_fine_store, &
                                              Shap3DQ_coarse_store,Nint_pp_store,&
                                              Error_org,Polyflag,rate_p)
!
   use data_structure3D
   use element_data
   use parametersDPG
   use physics
   use mpi_param, only: RANK
!
   implicit none
   integer,    intent(in)  :: Mdle
   integer,    intent(in)  :: Nr_mdle_sons
   integer,    intent(in)  :: Mdle_sons(Nr_mdle_sons)
   integer,    intent(in)  :: Nord_glob
   integer,    intent(in)  :: NrdofQ
   integer,    intent(in)  :: Flag_pref_loc
   real(8),    intent(in)  :: Ap(NrdofQ)
   real(8),    intent(in)  :: Zbload(NrdofQ,MAXEQNQ)
   real(8),    intent(in)  :: Weights_fine_store(MAXNINT3ADD,2,Nr_mdle_sons)
   real(8),    intent(in)  :: Quad_point_store(3,MAXNINT3ADD,Nr_mdle_sons)
   real(8),    intent(in)  :: Shap3DQ_fine_store(MAXbrickQQ,MAXNINT3ADD,Nr_mdle_sons)
   real(8),    intent(in)  :: Shap3DQ_coarse_store(MAXbrickQQ,MAXNINT3ADD,Nr_mdle_sons)
   integer,    intent(in)  :: Nint_pp_store(Nr_mdle_sons)
!
   real(8), intent(out)    :: Error_org
   integer, intent(out)    :: Polyflag
   real(8), intent(out)    :: rate_p
!..arrays containing the indices of the L2 shape functions generated using tensor product.
   integer, allocatable    :: nextract(:)
   integer, allocatable    :: nextract_prev(:)
!..variables for the function
   real(8), allocatable :: bwork(:,:), nextract_save_lvl(:,:)
   integer :: px,py,pz, nrdofgQ, nrdofmQ, nrdof_org
   integer :: nord_org, nord_mod, nord_prev
   integer :: order_add
   integer :: etype
   integer :: iattr
   integer :: max_loc(1)
   real(8) :: proj_error_net,proj_error, rate_max_lvl
   real(8) :: error_rates(6)
   integer :: poly_flags_lvl(6), poly_flag_chosen
!..auxiliary variables
   integer :: l,i,counter
!
   etype           = NODES(Mdle)%ntype
   error_rates     = ZERO
   poly_flags_lvl  = ZERO
!..cuurently only for Hexa elements
   if(etype .eq. MDLB) then
      call ddecode(Nord_glob,px,py,pz)
!  ...L2 dofs correspoding to isotropic p-refinement
      nrdofgQ = px * py * pz
!
      allocate(nextract_prev(nrdofgQ))
      allocate(nextract_save_lvl(nrdofgQ,3))
      allocate(bwork(nrdofgQ,MAXEQNQ))
!
      bwork = ZERO
!  ...if p-refinement done during fine mesh generation
      if(Flag_pref_loc .eq. 1) then 
         nord_org = Nord_glob - 111
      else
         nord_org = Nord_glob
      endif
!  ...traversing the set of candidate polynomial order starts at the order of the coarse element
!  ...level 0 (base polynomial order)
      call ddecode(nord_org,px,py,pz)
      nrdofmQ = px * py * pz
!  ...number of dofs for original order before iso hp-refinement
      nrdof_org = nrdofmQ
      allocate(nextract(nrdofmQ))
!  ...extraction vector corresponding to the current order
      nextract = ZERO
!  ...extraction vector corresponding to the previous order
      nextract_prev = ZERO
!  ...saves  the extraction vectors for increment (px,py,pz) -> (px+1,py,pz),(px,py+1,pz),(px,py,pz+1)
!  in different columns
      nextract_save_lvl = ZERO
!  ...calling the routine that computes the extraction vector: it utillizes the previous
!  extraction vector to build the current one.
      call extraction_vector_new(0,nord_org,Nord_glob,nrdofmQ,nextract_prev,nextract)
!  ...projection on coarse element (diagonal system)
      do iattr = 1,NRQVAR
         do l = 1,nrdofmQ
               bwork(l,iattr) = Zbload(nextract(l),iattr)/Ap(nextract(l))
         enddo
      enddo
!  ...storing the extraction vector of the base level (px,py,pz) into nextract_prev
      nextract_prev(1:nrdofmQ) = nextract(1:nrdofmQ)
      proj_error_net = 0.d0
!  ...computing the projection error by looping over all the L2 variables
      do iattr = 1,NRQVAR
         call fine_to_coarse_projection_error(Nr_mdle_sons,bwork(1:nrdofmQ,iattr),nextract,nrdofmQ, &
                                              nrdofgQ,Mdle,iattr,Mdle_sons,&
                                              Weights_fine_store,Quad_point_store,Nint_pp_store,&
                                              Shap3DQ_fine_store,Shap3DQ_coarse_store,proj_error)
         proj_error_net = proj_error_net + proj_error
      enddo
      deallocate(nextract)
!  ...projection error for the coarse element
      Error_org = proj_error_net
!  ...optimal p-only refinement search will be done only if the element was p-refined
!  during the fine mesh generation
      if(Flag_pref_loc .eq. 1) then
!     ...nord_mod : modified order,nord_prev : previous order
         nord_mod  = nord_org
         nord_prev = nord_org
!     ...level 1
         counter = 0
!     ...order of increase on x then y and then z on base order
         do i = 1,3
!
            order_add = int(10**(3-i))
            nord_mod  = nord_mod + order_add
!
            call ddecode(nord_mod,px,py,pz)
            nrdofmQ = px * py * pz
            allocate(nextract(nrdofmQ))
            nextract = ZERO
!
            call extraction_vector_new(nord_prev,nord_mod,Nord_glob,nrdofmQ,nextract_prev,nextract)
!
            do iattr = 1,NRQVAR
               do l = 1,nrdofmQ
                  bwork(l,iattr) = Zbload(nextract(l),iattr)/Ap(nextract(l))
               enddo
            enddo
!
            proj_error_net = 0.d0
            do iattr = 1,NRQVAR
               call fine_to_coarse_projection_error(Nr_mdle_sons,bwork(1:nrdofmQ,iattr),nextract,nrdofmQ, &
                                                    nrdofgQ,Mdle,iattr,Mdle_sons,&
                                                    Weights_fine_store,Quad_point_store,Nint_pp_store,&
                                                    Shap3DQ_fine_store,Shap3DQ_coarse_store,proj_error)

               proj_error_net = proj_error_net + proj_error
            enddo
!
            nextract_save_lvl(1:nrdofmQ,i) = nextract(1:nrdofmQ)
!        ...computing the error rates and saving the corresponding polynomial order
            error_rates(i)    = (error_org - proj_error_net)/abs(real(NRQVAR * nrdofmQ,8) - real(NRQVAR * nrdof_org,8))
            poly_flags_lvl(i) = nord_mod

            nord_mod = nord_org
            counter = counter + 1
            deallocate(nextract)
         enddo
!     ...selecting the first lvl candidate (one with the maximum error decrement rate)
         rate_max_lvl = maxval(error_rates(1:3))
         max_loc      = maxloc(error_rates(1:3))
         poly_flag_chosen = poly_flags_lvl(max_loc(1))

         nord_mod  = poly_flag_chosen
         nord_prev = poly_flag_chosen

         call ddecode(nord_mod,px,py,pz)
         nrdofmQ = px * py * pz

         nextract_prev            = ZERO
         nextract_prev(1:nrdofmQ) = nextract_save_lvl(1:nrdofmQ,max_loc(1))
         nextract_save_lvl        = ZERO
!
!..level 2: starting from order selected at the end of the level 1
         do i = 1,3
            if(i .ne. max_loc(1)) then

               order_add = int(10**(3-i))
               nord_mod  = nord_mod + order_add
               call ddecode(nord_mod,px,py,pz)
               nrdofmQ = px * py * pz
               allocate(nextract(nrdofmQ))
               nextract = ZERO
!
               call extraction_vector_new(nord_prev,nord_mod,Nord_glob,nrdofmQ,nextract_prev,nextract)
!
               do iattr = 1,NRQVAR
                  do l = 1,nrdofmQ
                        bwork(l,iattr) = Zbload(nextract(l),iattr)/Ap(nextract(l))
                  enddo
               enddo
!
               proj_error_net = 0.d0
               do iattr = 1,NRQVAR
                  call fine_to_coarse_projection_error(Nr_mdle_sons,bwork(1:nrdofmQ,iattr),nextract,nrdofmQ, &
                                                       nrdofgQ,Mdle,iattr,Mdle_sons,&
                                                       Weights_fine_store,Quad_point_store,Nint_pp_store,&
                                                       Shap3DQ_fine_store,Shap3DQ_coarse_store,proj_error)
                  proj_error_net = proj_error_net + proj_error
               enddo
!
               nextract_save_lvl(1:nrdofmQ,counter + 1 - 3) = nextract(1:nrdofmQ)
!
               error_rates(counter + 1)    = (error_org - proj_error_net)/abs(real(NRQVAR * nrdofmQ,8) - real(NRQVAR * nrdof_org,8))
               poly_flags_lvl(counter + 1) = nord_mod
!
               nord_mod = poly_flag_chosen
               counter = counter + 1
               deallocate(nextract)
            endif
         enddo
!
!     ...selecting the second level candidate
         rate_max_lvl = maxval(error_rates(4:5))
         max_loc      = maxloc(error_rates(4:5))
         poly_flag_chosen = poly_flags_lvl(3 + max_loc(1))
!     ...final level
         nord_prev = poly_flag_chosen
         call ddecode(nord_prev,px,py,pz)
         nrdofmQ = px * py * pz
!
         nextract_prev            = ZERO
         nextract_prev(1:nrdofmQ) = nextract_save_lvl(1:nrdofmQ,max_loc(1))
!
         nord_mod = nord_org + 111
         call ddecode(nord_mod,px,py,pz)
         nrdofmQ = px * py * pz
!
         allocate(nextract(nrdofmQ))
         nextract = ZERO
!
         call extraction_vector_new(nord_prev,nord_mod,Nord_glob,nrdofmQ,nextract_prev,nextract)
!
         do iattr = 1,NRQVAR
               do l = 1,nrdofmQ
                  bwork(l,iattr) = Zbload(nextract(l),iattr)/Ap(nextract(l))
               enddo
         enddo
!
         proj_error_net = 0.d0
         do iattr = 1,NRQVAR
               call fine_to_coarse_projection_error(Nr_mdle_sons,bwork(1:nrdofmQ,iattr),nextract,nrdofmQ, &
                                                    nrdofgQ,Mdle,iattr,Mdle_sons,&
                                                    Weights_fine_store,Quad_point_store,Nint_pp_store,&
                                                    Shap3DQ_fine_store,Shap3DQ_coarse_store,proj_error)
               proj_error_net = proj_error_net + proj_error
         enddo
         error_rates(counter + 1)    = (error_org - proj_error_net)/abs(real(NRQVAR * nrdofmQ,8) - real(NRQVAR * nrdof_org,8))
         poly_flags_lvl(counter + 1) = nord_mod

         rate_p    = maxval(error_rates)
         max_loc   = maxloc(error_rates)
         Polyflag = poly_flags_lvl(max_loc (1))
         deallocate(nextract)
      else
         rate_p = 0.0
         Polyflag = nord_org
      endif
    endif
end subroutine opt_polynomial_search_coarse_linear
