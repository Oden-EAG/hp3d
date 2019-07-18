!----------------------------------------------------------------------
!
!     routine:    zoltan_lb_param_block
!     purpose:    set up block partition load balancer
!
!     The following parameters [default values] are available:
!        none
!
!----------------------------------------------------------------------
subroutine zoltan_lb_param_block(Ierr_out)
   use zoltan_wrapper, only: zz,zoltan_w_handle_err
   use zoltan
   implicit none
   integer(Zoltan_int), intent(out) :: Ierr_out
   integer(Zoltan_int) :: ierr
   Ierr_out = ZOLTAN_OK
!
   ierr = Zoltan_Set_Param(zz,'LB_METHOD','BLOCK')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..block partitioner supports weights
   ierr = Zoltan_Set_Param(zz,'OBJ_WEIGHT_DIM','1')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
end subroutine zoltan_lb_param_block
!
!
!----------------------------------------------------------------------
!
!     routine:    zoltan_lb_param_random
!     purpose:    set up random partition load balancer
!
!     The following parameters [default values] are available:
!        RANDOM_MOVE_FRACTION [1.0]
!
!----------------------------------------------------------------------
subroutine zoltan_lb_param_random(Ierr_out)
   use zoltan_wrapper, only: zz,zoltan_w_handle_err
   use zoltan
   implicit none
   integer(Zoltan_int), intent(out) :: Ierr_out
   integer(Zoltan_int) :: ierr
   Ierr_out = ZOLTAN_OK
!
   ierr = Zoltan_Set_Param(zz,'LB_METHOD','RANDOM')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..random partitioner does not support weights
   ierr = Zoltan_Set_Param(zz,'OBJ_WEIGHT_DIM','0')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..RANDOM_MOVE_FRACTION
!     The fraction of objects to randomly move
!     1.0 = move all; 0.0 = move nothing
   ierr = Zoltan_Set_Param(zz,'RANDOM_MOVE_FRACTION','1.0')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
end subroutine zoltan_lb_param_random
!
!
!----------------------------------------------------------------------
!
!     routine:    zoltan_lb_param_rcb
!     purpose:    set up recursive coordinate bisection load balancer
!
!     The following parameters [default values] are available:
!        RCB_OVERALLOC [1.2]
!        RCB_REUSE [0]
!        RCB_OUTPUT_LEVEL [0]
!        CHECK_GEOM [1]
!        KEEP_CUTS [0]
!        AVERAGE_CUTS [0]
!        RCB_LOCK_DIRECTIONS [0]
!        REDUCE_DIMENSIONS [0]
!        DEGENERATE_RATIO [10]
!        RCB_SET_DIRECTIONS [0]
!        RCB_RECTILINEAR_BLOCKS [0]
!        RCB_RECOMPUTE_BOX [0]
!        OBJ_WEIGHTS_COMPARABLE [0]
!        RCB_MULTICRITERIA_NORM [1]
!        RCB_MAX_ASPECT_RATIO [10]
!
!----------------------------------------------------------------------
subroutine zoltan_lb_param_rcb(Ierr_out)
   use zoltan_wrapper, only: zz,zoltan_w_handle_err
   use zoltan
   implicit none
   integer(Zoltan_int), intent(out) :: Ierr_out
   integer(Zoltan_int) :: ierr
   Ierr_out = ZOLTAN_OK
!
   ierr = Zoltan_Set_Param(zz,'LB_METHOD','RCB')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..RCB partitioner supports weights
   ierr = Zoltan_Set_Param(zz,'OBJ_WEIGHT_DIM','1')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..RCB_OUTPUT_LEVEL
!     0: no output
!     1: print summary
!     2: print for each proc
   ierr = Zoltan_Set_Param(zz,'RCB_OUTPUT_LEVEL','1')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
end subroutine zoltan_lb_param_rcb
!
!
!----------------------------------------------------------------------
!
!     routine:    zoltan_lb_param_rib
!     purpose:    set up recursive intertial bisection load balancer
!
!     The following parameters [default values] are available:
!        RIB_OVERALLOC [1.2]
!        RIB_OUTPUT_LEVEL [0]
!        CHECK_GEOM [1]
!        KEEP_CUTS [0]
!        AVERAGE_CUTS [0]
!        REDUCE_DIMENSIONS [0]
!        DEGENERATE_RATIO [10]
!
!----------------------------------------------------------------------
subroutine zoltan_lb_param_rib(Ierr_out)
   use zoltan_wrapper, only: zz,zoltan_w_handle_err
   use zoltan
   implicit none
   integer(Zoltan_int), intent(out) :: Ierr_out
   integer(Zoltan_int) :: ierr
   Ierr_out = ZOLTAN_OK
!
   ierr = Zoltan_Set_Param(zz,'LB_METHOD','RIB')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..RIB partitioner supports weights
   ierr = Zoltan_Set_Param(zz,'OBJ_WEIGHT_DIM','1')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..RIB_OUTPUT_LEVEL
!     0: no output
!     1: print summary
!     2: print for each proc
   ierr = Zoltan_Set_Param(zz,'RIB_OUTPUT_LEVEL','1')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
end subroutine zoltan_lb_param_rib
!
!
!----------------------------------------------------------------------
!
!     routine:    zoltan_lb_param_hsfc
!     purpose:    set up Hilbert space filling curve load balancer
!
!     The following parameters [default values] are available:
!        KEEP_CUTS [0]
!        REDUCE_DIMENSIONS [0]
!        DEGENERATE_RATIO [10]
!
!----------------------------------------------------------------------
subroutine zoltan_lb_param_hsfc(Ierr_out)
   use zoltan_wrapper, only: zz,zoltan_w_handle_err
   use zoltan
   implicit none
   integer(Zoltan_int), intent(out) :: Ierr_out
   integer(Zoltan_int) :: ierr
   Ierr_out = ZOLTAN_OK
!
   ierr = Zoltan_Set_Param(zz,'LB_METHOD','HSFC')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..HSFC partitioner supports weights
   ierr = Zoltan_Set_Param(zz,'OBJ_WEIGHT_DIM','1')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
end subroutine zoltan_lb_param_hsfc
