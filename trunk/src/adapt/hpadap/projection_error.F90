!-----------------------------------------------------------------------
!> @brief Routine computes projection error on the coarse element while projecting the reference solution
!> @param[in]  Nr_mdle_sons            - Number of isotropic hp refinement
!> @param[in]  Coeff                   - Coefficients of the shape functions for an L2 attr on the coarse element
!> @param[in]  Nextract                - Extraction vector corresponding to the coarse element
!> @param[in]  NrdofmQ                 - number of L2 dofs for the coarse element
!> @param[in]  NrdofgQ                 - number of L2 dofs correspoding isotropic(p+1) order
!> @param[in]  Mdle                    - Middle node number of the coarse element
!> @param[in]  Iattr                   - Index of the L2 attribute
!> @param[in]  Mdle_sons               - Array containing the middle node numbers of the isotropic childs
!> @param[in]  Weights_fine_store      - weights and inverse of jacobian at each quadrature point
!!                                       for the sons
!> @param[in]  Quad_point_store        - coordinates of each quadrature point for the sons
!> @param[in]  Nint_pp_store           - number of quadrature points for the sons
!> @param[in]  Shap3DQ_fine_store      - values of L2 shape functions for the sons
!> @param[in]  Shap3DQ_coarse_store    - values of L2 shape functions for the coarse element
!!
!> @param[out] Proj_error              - projection error
!> @date June 2024
!-----------------------------------------------------------------------
subroutine fine_to_coarse_projection_error(Nr_mdle_sons,Coeff,Nextract,NrdofmQ,NrdofgQ,Mdle,Iattr,&
                                           Mdle_sons,Weights_fine_store,Quad_point_store,&
                                           Nint_pp_store,Shap3DQ_fine_store,Shap3DQ_coarse_store,&
                                           Proj_error)
!
   use control
   use data_structure3D
   use element_data
   use parametersDPG
!
   implicit none
!
   integer,    intent(in)  :: NrdofmQ
   integer,    intent(in)  :: NrdofgQ
   integer,    intent(in)  :: nr_mdle_sons
   integer,    intent(in)  :: Mdle
   integer,    intent(in)  ::  Iattr
   integer,    intent(in)  :: Mdle_sons(nr_mdle_sons)
   real(8),    intent(in)  :: Coeff(NrdofmQ)
   integer,    intent(in)  :: Nextract(NrdofmQ)
   real(8),    intent(in)  :: Weights_fine_store(MAXNINT3ADD,2,nr_mdle_sons)
   real(8),    intent(in)  :: Quad_point_store(3,MAXNINT3ADD,nr_mdle_sons)
   real(8),    intent(in)  :: Shap3DQ_fine_store(MAXbrickQQ,MAXNINT3ADD,nr_mdle_sons)
   real(8),    intent(in)  :: Shap3DQ_coarse_store(MAXbrickQQ,MAXNINT3ADD,nr_mdle_sons)
   integer,    intent(in)  :: Nint_pp_store(nr_mdle_sons)
!
   real(8),    intent(out) :: Proj_error
!
!..to store the Coefficients of the basis functions of interpolant on original element after projection from fine to coarse element.
   real(8) :: zdofQ(MAXbrickQQ)
!..to store extracted Coefficeints
   real(8) :: zdofQ_pp(MAXEQNQ,MAXbrickQQ), zdofH_pp(MAXEQNH,MAXbrickHH),zdofE_pp(MAXEQNE,MAXbrickEE), &
              zdofV_pp(MAXEQNV,MAXbrickVV)
!..to store the values of the solution for the original and sons respectively
   real(8) :: zvalQ(MAXEQNQ), zvalQpp(MAXEQNQ)
!..variables required for coarse element
!.. number of shape functions in trial space for original element (in case of ultra weak formulation: only nrdofQ is useful)
   integer :: nrdofH ,nrdofE ,nrdofV ,nrdofQ
!.. number of shape functions in trial space for enriched sons (its same for all of them)
   integer :: nrdofH_pp, nrdofE_pp, nrdofV_pp, nrdofQ_pp
!.. L2 shape functions
   real(8) :: shapQ(MAXbrickQQ)
!..H1 shape functions
   real(8) :: shapH(MAXbrickHH), gradH(3,MAXbrickHH)
!..3D quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
!..geometry
   real(8) :: xi(3), x(3), dxdxi(3,3), dxidx(3,3), rn(3),xis(3)
   real(8) :: dxidt(3,2), dxdt(3,2), t(2)
!..geometry dof of original element
   real(8) :: xnod(3,MAXbrickHH)
!..geometry dof of sons
   real(8) :: xnod_pp(3,MAXbrickHH)
   integer :: etype
!..element order, orientation for edges and faces of the original element
   integer :: norder(19), norient_edge(12), norient_face(6)
   integer :: nint, nrdof
!..element order, orientation for edges and faces for enriched sons
   integer :: norder_pp(19), norient_edge_pp(12), norient_face_pp(6)
   integer :: nint_pp
!..auxiliary variables
   real(8) :: wa,weight,rjac,bjac,fval
   integer :: k1,k2,l,iflag,k,is,mdle_fine,j,order_add,Nord_mod    
   integer :: idec
   real(8) :: q
!
   etype = NODES(Mdle)%ntype
   call find_order(Mdle, norder)
   call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!
   Proj_error = 0.d0
   do is = 1,nr_mdle_sons
      mdle_fine = Mdle_sons(is)
      etype = NODES(mdle_fine)%ntype
      call find_order(mdle_fine, norder_pp)
      call celndof(etype,norder_pp, nrdofH_pp,nrdofE_pp,nrdofV_pp,nrdofQ_pp)
      call solelm_L2(mdle_fine,zdofQ_pp)
      nint_pp = Nint_pp_store(is)
      do l = 1,nint_pp
!     ...coordinates and weight of this integration point
         xi(1:3) = Quad_point_store(1:3,l,is)
         shapQ(1:nrdofQ_pp) = Shap3DQ_fine_store(1:nrdofQ_pp,l,is)
         weight = Weights_fine_store(l,1,is)
         rjac = Weights_fine_store(l,2,is)
         zvalQpp = ZERO
!
         do k = 1,nrdofQ_pp
               q = shapQ(k)/rjac
               zvalQpp(Iattr) = zvalQpp(Iattr) + zdofQ_pp(Iattr,k) * q
         enddo
!
         shapQ(1:nrdofQ_pp) = Shap3DQ_coarse_store(1:nrdofQ_pp,l,is)
!     ...scaling the jacobian for isotropic refinement of coarse element
         rjac = rjac * real(nr_mdle_sons,8)
!     ...recontructing the projection at the transformed gauss point for
!        the transformed coarse elements.
         zvalQ = ZERO
         do k = 1,NrdofmQ
               q = shapQ(Nextract(k))/rjac
               zvalQ(Iattr) = zvalQ(Iattr) + Coeff(k) * q
         enddo
         Proj_error = Proj_error + weight * (zvalQ(Iattr) - zvalQpp(Iattr))**2
      enddo
   enddo
!
end subroutine fine_to_coarse_projection_error
!
!-----------------------------------------------------------------------
!> @brief Routine computes projection error on a h-ref candidate while projecting the reference solution
!> @param[in]  Kref                    - kind of h-ref (e.g. for brick element 011,101,001,100,010,110,111)
!> @param[in]  Coeff                   - Coefficients of the shape functions for an L2 attr on the h-ref candidate
!> @param[in]  Nextract                - Extraction vector corresponding to the h-ref candidate
!> @param[in]  Overlap                 - contains which isotropic child elements forms the h-ref candidate
!> @param[in]  NrdofmQ                 - number of L2 dofs for the coarse element
!> @param[in]  NrdofgQ                 - number of L2 dofs correspoding isotropic(p+1) order
!> @param[in]  Iattr                   - Index of the L2 attribute
!> @param[in]  Mdle                    - Middle node number of the coarse element
!> @param[in]  Weights_fine_store      - weights and inverse of jacobian at each quadrature point
!!                                       for the sons
!> @param[in]  Quad_point_store        - coordinates of each quadrature point for the sons
!> @param[in]  Nint_pp_store           - number of quadrature points for the sons
!> @param[in]  Shap3DQ_fine_store      - values of L2 shape functions for the sons
!> @param[in]  Shap3DQ_coarse_store    - values of L2 shape functions for the coarse element
!!
!> @param[out] Proj_error              - projection error
!> @date June 2024
!-----------------------------------------------------------------------
subroutine fine_to_subson_projection_error(Kref,Coeff,Nextract,Overlap,NrdofmQ,NrdofgQ,Iattr,&
                                                Mdle,Weights_fine_store,Quad_point_store,&
                                                Nint_pp_store,Shap3DQ_fine_store,Shap3DQ_coarse_store,&
                                                Proj_error)
   use control
   use data_structure3D
   use element_data
   use parametersDPG
!
   implicit none
!
   integer,    intent(in)  ::  Kref
   integer,    intent(in)  ::  NrdofmQ
   integer,    intent(in)  ::  NrdofgQ
   integer,    intent(in)  ::  Iattr
   real(8),    intent(in)  :: Coeff(NrdofmQ)
   integer,    intent(in)  :: Nextract(NrdofmQ)
   real(8),    intent(in)  :: Weights_fine_store(MAXNINT3ADD,2,*)
   real(8),    intent(in)  :: Quad_point_store(3,MAXNINT3ADD,*)
   real(8),    intent(in)  :: Shap3DQ_fine_store(MAXbrickQQ,MAXNINT3ADD,*)
   real(8),    intent(in)  :: Shap3DQ_coarse_store(MAXbrickQQ,MAXNINT3ADD,*)
   integer,    intent(in)  :: Nint_pp_store(*)
   integer,    intent(in)  :: Mdle
   integer,    intent(in)  :: Overlap(*)
!
   real(8),    intent(out) ::  Proj_error
!
   integer :: norder_pp(19), norient_edge_pp(12), norient_face_pp(6)
   integer :: nrdofH_pp, nrdofE_pp, nrdofV_pp, nrdofQ_pp
!..to store extracted Coefficeints
   real(8) :: zdofQ_pp(MAXEQNQ,MAXbrickQQ), zdofH_pp(MAXEQNH,MAXbrickHH),zdofE_pp(MAXEQNE,MAXbrickEE), &
   zdofV_pp(MAXEQNV,MAXbrickVV)
!
!..to store the values of the solution for the original and sons respectively
   real(8) :: zvalQ(MAXEQNQ), zvalQpp(MAXEQNQ)
!.. L2 shape functions
   real(8) :: shapQ(MAXbrickQQ)
!..H1 shape functions
   real(8) :: shapH(MAXbrickHH), gradH(3,MAXbrickHH)
!..geometry dof of sons
   real(8) :: xnod_pp(3,MAXbrickHH)
!..3D quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
!..geometry
   real(8) :: xi(3), x(3), dxdxi(3,3), dxidx(3,3), rn(3),xis(3)
   real(8) :: dxidt(3,2), dxdt(3,2), t(2)
!..refinement indicators in x,y,z directions
   integer :: hx,hy,hz
   integer :: etype
   integer :: etype_coarse
!..auxiliary variables
   real(8) :: rjac,wa,weight,q
   integer :: first_son,nr_subsons,nr_mdle_sons,iss,iflag,nrdof
   integer :: is, mdle_fine,Overlap_count,nint_pp
   integer :: k,l
   integer,allocatable :: Mdle_olp_fine(:)
!
   etype_coarse = NODES(Mdle)%ntype
!
   if(etype_coarse .eq. MDLB) then
      nr_mdle_sons = 8
!
      call ddecode(Kref,hx,hy,hz)
      nr_subsons = 2**(hx+hy+hz)
      Overlap_count = nr_mdle_sons/nr_subsons
!
      allocate(Mdle_olp_fine(Overlap_count))
      first_son = NODES(Mdle)%first_son
!
      do is = 1,Overlap_count
         Mdle_olp_fine(is) = first_son + Overlap(is) - 1
      enddo
!
      Proj_error = 0.d0
      do is = 1,Overlap_count
!
         mdle_fine = Mdle_olp_fine(is)
         etype = NODES(mdle_fine)%ntype
         call find_order(mdle_fine, norder_pp)
         call celndof(etype,norder_pp, nrdofH_pp,nrdofE_pp,nrdofV_pp,nrdofQ_pp)
         call solelm_L2(mdle_fine,zdofQ_pp)
         nint_pp = Nint_pp_store(is)
         do l = 1,nint_pp
!
               xi(1:3) = Quad_point_store(1:3,l,is)
               shapQ(1:nrdofQ_pp) = Shap3DQ_fine_store(1:nrdofQ_pp,l,is)
               weight = Weights_fine_store(l,1,is)
               rjac = Weights_fine_store(l,2,is)
               zvalQpp = ZERO
!
               do k = 1,nrdofQ_pp
                  q = shapQ(k)/rjac
                  zvalQpp(Iattr) = zvalQpp(Iattr) + zdofQ_pp(Iattr,k) * q
               enddo
!
               shapQ = ZERO
               shapQ(1:nrdofQ_pp) = Shap3DQ_coarse_store(1:nrdofQ_pp,l,is)
!
!           ...scaling the jacobian for isotropic refinement of coarse element
               rjac = rjac * real(Overlap_count,8)
!
!           ...recontructing the projection at the transformed gauss point for
!              the transformed coarse elements.
               zvalQ = ZERO
               do k = 1,NrdofmQ
                  q = shapQ(Nextract(k))/rjac
                  zvalQ(Iattr) = zvalQ(Iattr) + Coeff(k) * q
               enddo
!
               Proj_error = Proj_error + weight * (zvalQ(Iattr) - zvalQpp(Iattr))**2
         enddo
      enddo
   endif
!
end subroutine fine_to_subson_projection_error