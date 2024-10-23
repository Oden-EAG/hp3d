!-----------------------------------------------------------------------
!> @brief Routine performs fine mesh element to candidate h-refined element projections
!> @param[in]   Mdle                -    element number
!> @param[in]   Flag_pref_loc       - flag indicating whether the element is p-refined or not
!> @param[in]   Error_org           - projection error at the original polynomial order on the coarse element
!> @param[in]   Rate_p              - highest error rate obtained from p-only refinement
!> @param[in]   Polyflag            - polynomial refinement flag associated with Rate_p
!> @param[out]  Elem_grate          - element guranteed rate 
!> @param[out]  Ref_indicator_flag  - refinement indicator flag array
!> @param[out]  Mep_Nord_href       - contains the polynomial flags for all the possible h-refinement candidates
!!                                    while traversing the maximum error reduction path
!> @param[out]  Elem_err_red_rate   - contains the error reduction rate while following maximum error reduction path
!!                                    for h-ref candidate who provides the maximum guranteed rate
!> @param[out]  Mep_count           - total number of steps taken while following the maximum error reduction path
!> @param[out]  Loc_max_rate        - maximum guranteed rate among all h-refinement candidates
!> @date May 2024
!-----------------------------------------------------------------------
subroutine project_h(Mdle,Flag_pref_loc,Error_org,Rate_p,Poly_flag,Istep, &
                     Elem_grate,Ref_indicator_flag, &
                     Mep_Nord_href,Elem_err_red_rate, &
                     Mep_count,Loc_max_rate)    
!
   use control
   use data_structure3D
   use element_data
   use parametersDPG
!
   implicit none
   integer,    intent(in)  :: Mdle
   integer,    intent(in)  :: Flag_pref_loc
   real(8),    intent(in)  :: Error_org
   real(8),    intent(in)  :: Rate_p
   integer,    intent(in)  :: Poly_flag
   integer,    intent(in)  :: Istep
!
   real(8),    intent(out) :: Elem_grate
   integer,    intent(out) :: Ref_indicator_flag(10)
   integer,    intent(out) :: Mep_Nord_href(100,8)
   real(8),    intent(out) :: Elem_err_red_rate(100)
   integer,    intent(out) :: Mep_count
   integer,    intent(out) :: Loc_max_rate
!
!..number of h-refinement options for an elemenet type
   integer :: ref_opts
! 
   integer :: k,j,m,l
   integer :: etype
!..projection error for each refinement options
   real(8), allocatable   ::  error_opt(:)
!..guranteed rate for each refinement option
   real(8), allocatable   ::  g_rate_ref(:)
!..refinement markers for each refinement option
   integer, allocatable   ::  kref_opts(:)
!..array containing the poly distributio for each ref opt for which we obtain guranteed rate
   integer, allocatable   ::  nord_max_href(:,:)
!..array to store error rate for competitive h-refinements
   real(8), allocatable   ::  rate_hcomp(:)
!..maximum rate over all refinement options
   real(8) :: max_rate_element
!..maximum rate for comp refs and error for correspoding comp-href
   real(8) :: max_rate_hcomp, error_max_rate
!
   integer, dimension(100,8,8)    ::  nord_threshold
   real(8), dimension(100,8)      ::  error_rate_ref
!..array keeping the count of configuration after N_coarse threshold is crossed.
   integer, allocatable           ::  count_ref(:)
!..array for location of max rate while maximum error reduction path
   integer, allocatable           ::  loc_max_rate_ref(:)
!
   integer :: max_rate_element_loc(1)
   integer :: nrdofgQ, nrdof_org
   integer :: nord_glob, nord_org ,px,py,pz
   integer :: norder(19)
!
   etype = NODES(Mdle)%ntype
!
   if(etype .eq. MDLB) then
      call find_order(Mdle, norder)
      nord_glob = norder(19) 
      call ddecode(nord_glob,px,py,pz)
      nrdofgQ = px * py *pz
!   
      if(Flag_pref_loc .eq. 1) then
         nord_org = norder(19) - 111
      else
         nord_org = norder(19)
      endif
!
      call ddecode(nord_org,px,py,pz)
      nrdof_org = px * py * pz
      ref_opts = 8
      allocate(error_opt(ref_opts))
      error_opt = ZERO
!
      allocate(kref_opts(ref_opts))
      kref_opts = ZERO
!
      allocate(g_rate_ref(ref_opts))
      g_rate_ref = ZERO
!
      allocate(rate_hcomp(ref_opts))
      rate_hcomp = ZERO
!
      allocate(nord_max_href(ref_opts,8))
      nord_max_href = ZERO
!
      allocate(count_ref(ref_opts))
      count_ref = ZERO
!
      allocate(loc_max_rate_ref(ref_opts))
      loc_max_rate_ref = ZERO

      m = 1
      do k = 0,1
         do j = 0,1
               do l = 0,1
                  kref_opts(m) = k * 10**(2) + j * 10 **1 + l * 10**0
                  m = m + 1
               enddo
         enddo
      enddo
!
      m = m - 1
      do j = 2,m
         call elem_proj_h_linear(Mdle,Flag_pref_loc,kref_opts(j),Error_org,error_opt(j),g_rate_ref(j),rate_hcomp(j),nord_max_href(j,1:8), &
                                 nord_threshold(1:100,1:8,j),count_ref(j),error_rate_ref(1:100,j),loc_max_rate_ref(j))
      enddo
!
      max_rate_element = maxval(g_rate_ref(2:m))
      error_max_rate  = maxval(error_opt(2:m))
      max_rate_element_loc =  maxloc(g_rate_ref(2:m))
      max_rate_hcomp = maxval(rate_hcomp(2:m))
!  ...selection between href or pref
      if(Flag_pref_loc .eq. 0) then
         Ref_indicator_flag(1) = 1
!  ...since 1 element of kref_opts means 000, hence we need to offset it by one. Thus,we have max_rate_element(1) + 1
         Ref_indicator_flag(2) = kref_opts(max_rate_element_loc(1)+1)
         Ref_indicator_flag(3:10) = nord_max_href(max_rate_element_loc(1) + 1,1:8)

         Mep_Nord_href(1:100,1:8) = nord_threshold(1:100,1:8,max_rate_element_loc(1) + 1)
         Loc_max_rate = loc_max_rate_ref(max_rate_element_loc(1) + 1)
         Mep_count = count_ref(max_rate_element_loc(1) + 1)
         Elem_err_red_rate(1:100) = error_rate_ref(1:100,max_rate_element_loc(1) + 1)
         Elem_grate = max_rate_element
      else
         if(max_rate_hcomp .ge. Rate_p) then
               Ref_indicator_flag(1) = 1
!           ...since 1 element of kref_opts means 000, hence we need to offset it by one. Thus,we have max_rate_element(1) + 1
               Ref_indicator_flag(2) = kref_opts(max_rate_element_loc(1)+1)
               Ref_indicator_flag(3:10) = nord_max_href(max_rate_element_loc(1) + 1,1:8)
               Mep_Nord_href(1:100,1:8) = nord_threshold(1:100,1:8,max_rate_element_loc(1) + 1)
               Loc_max_rate = loc_max_rate_ref(max_rate_element_loc(1) + 1)
               Mep_count = count_ref(max_rate_element_loc(1) + 1)
               Elem_err_red_rate(1:100) = error_rate_ref(1:100,max_rate_element_loc(1) + 1)
               Elem_grate = max_rate_element
         else
               Ref_indicator_flag(1) = 2
               Ref_indicator_flag(2) = 000
               Ref_indicator_flag(3) = Poly_flag !p-refinement by degree 1.
               Elem_grate = Rate_p
         endif
      endif
!  ...deallocating allocatable arrays
      deallocate(error_opt)
      deallocate(g_rate_ref)
      deallocate(kref_opts)
      deallocate(nord_max_href)
      deallocate(rate_hcomp)
      deallocate(count_ref)
      deallocate(loc_max_rate_ref)
   endif
end subroutine project_h
