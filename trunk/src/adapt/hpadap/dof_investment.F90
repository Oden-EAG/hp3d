!-----------------------------------------------------------------------
!> @brief Routine implements the maximum investment strategy
!> @param[in] Threshold_frac        -  fraction of guranteed rate fow which 
!                                      investment in terms of dof occurs. 
!> @param[in] Grate_mesh            -  guranteed rate for the element.
!> @param[in] Mep_Nord_href_elem    -  array containing polynomial orders for the
!                                      for the maximum error path of the  
!                                      h-ref candidate which provides the highest guranteed rate.
!> @param[in] Elem_err_red_rate     -  array containing the error rate for maximum error rate path 
!                                      for the h-ref candidate.
!> @param[in] Loc_max_rate          -  index of the polynomial order that provided the guranteed rate
!> @param[in] Mep_count             -  total number of steps in maximum error rate path.
!
!> @param[out]Ref_indicator_flags   -  Refinement flag array containing h-ref and associated 
!                                      polynomial order.
!---------------------------------------------------------------------------
subroutine dof_investment(Threshold_frac,Grate_mesh, Mep_Nord_href_elem, &
                          Elem_err_red_rate,Loc_max_rate,Mep_count,      &
                          Ref_indicator_flags)
!
   implicit none
!
   real(8), intent(in)                      :: Threshold_frac
   real(8), intent(in)                      :: Grate_mesh
   integer, dimension(100,8),  intent(in)   :: Mep_Nord_href_elem
   real(8), dimension(100),    intent(in)   :: Elem_err_red_rate
   integer,    intent(in)                   :: Loc_max_rate
   integer,    intent(in)                   :: Mep_count
!
   integer, dimension(10),     intent(out)  :: Ref_indicator_flags
!
!..aux variables
   real(8) :: threshold_rate
   integer :: i,idx
!
   threshold_rate = Grate_mesh * Threshold_frac
   idx = Loc_max_rate
!
   do i=Loc_max_rate,Mep_count
      if(Elem_err_red_rate(i) .ge. threshold_rate) then
         idx = i
      endif
   enddo
!
   Ref_indicator_flags(3:10) = Mep_Nord_href_elem(idx,1:8)
!   
end subroutine dof_investment