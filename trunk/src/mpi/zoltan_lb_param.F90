!----------------------------------------------------------------------
!
!> @name    zoltan_lb_param_block
!> @brief   set up block partition load balancer
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
!> @name    zoltan_lb_param_random
!> @brief   set up random partition load balancer
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
!> @name    zoltan_lb_param_rcb
!> @brief   set up recursive coordinate bisection load balancer
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
   ierr = Zoltan_Set_Param(zz,'RCB_OUTPUT_LEVEL','0')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..RCB_LOCK_DIRECTIONS
!     0: don't lock directions
!     1: lock directions
!   ierr = Zoltan_Set_Param(zz,'RCB_LOCK_DIRECTIONS','1')
!   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!   if (ierr > Ierr_out) Ierr_out = ierr
!
!..RCB_SET_DIRECTIONS
!     0 = don't order cuts;
!     1 = xyz; 2 = xzy; 3 = yzx;
!     4 = yxz; 5 = zxy; 6 = zyx;
!   ierr = Zoltan_Set_Param(zz,'RCB_SET_DIRECTIONS','6')
!   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!   if (ierr > Ierr_out) Ierr_out = ierr
!
!..RCB_MAX_ASPECT_RATIO
!   ierr = Zoltan_Set_Param(zz,'RCB_MAX_ASPECT_RATIO','5')
!   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!   if (ierr > Ierr_out) Ierr_out = ierr
!
end subroutine zoltan_lb_param_rcb
!
!
!----------------------------------------------------------------------
!
!> @name    zoltan_lb_param_rib
!> @brief   set up recursive intertial bisection load balancer
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
   ierr = Zoltan_Set_Param(zz,'RIB_OUTPUT_LEVEL','0')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..REDUCE_DIMENSIONS
!   ierr = Zoltan_Set_Param(zz,'REDUCE_DIMENSIONS','1')
!   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!   if (ierr > Ierr_out) Ierr_out = ierr
!
end subroutine zoltan_lb_param_rib
!
!
!----------------------------------------------------------------------
!
!> @name    zoltan_lb_param_hsfc
!> @brief   set up Hilbert space filling curve load balancer
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
!
!
!----------------------------------------------------------------------
!
!> @name    zoltan_lb_param_hsfc
!> @brief   set up Hilbert space filling curve load balancer
!
!     The following parameters [default values] are available:
!        GRAPH_PACKAGE [PHG]
!     ParMETIS parameters:
!        LB_APPROACH [Repartition]
!        PARMETIS_METHOD [AdaptiveRepart]
!        PARMETIS_OUTPUT_LEVEL [0]
!        PARMETIS_COARSE_ALG [2]
!        PARMETIS_SEED [15]
!        PARMETIS_ITR [100]
!        USE_OBJ_SIZE [1]
!        CHECK_GRAPH [1]
!        SCATTER_GRAPH [1]
!        GRAPH_SYMMETRIZE [NONE]
!        GRAPH_SYM_WEIGHT [ADD]
!
!----------------------------------------------------------------------
subroutine zoltan_lb_param_graph(Ierr_out)
   use zoltan_wrapper, only: zz,zoltan_w_handle_err
   use zoltan
   implicit none
   integer(Zoltan_int), intent(out) :: Ierr_out
   integer(Zoltan_int) :: ierr
   Ierr_out = ZOLTAN_OK
!
   ierr = Zoltan_Set_Param(zz,'LB_METHOD','GRAPH')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..Graph partitioner supports weights
   ierr = Zoltan_Set_Param(zz,'EDGE_WEIGHT_DIM','1')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..GRAPH_PACKAGE
!     PHG     : use Zoltan's PHG hypergraph package
!     PARMETIS: use METIS/ParMETIS ordering
!     SCOTCH  : use Scotch/PT-Scotch ordering
   ierr = Zoltan_Set_Param(zz,'GRAPH_PACKAGE','PARMETIS')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..LB_APPROACH
!     PARTITION  : partition "from scratch", not taking into account
!                  current data distribution (static load balancing)
!     REPARTITION: partition but take into account current data
!                  distribution to keep data migration low
!                  (dynamic load balancing)
!     REFINE     : quickly improve the current data distribution
   ierr = Zoltan_Set_Param(zz,'LB_APPROACH','REPARTITION')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..PARMETIS_OUTPUT_LEVEL
!     0: no output
!     1: print timing info
   ierr = Zoltan_Set_Param(zz,'PARMETIS_OUTPUT_LEVEL','0')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..PARMETIS_ITR
!     Ratio of interprocessor communication time to redistribution time.
!     A high value will emphasize reducing the edge cut, while a small value
!     will try to keep the change in the new partition (distribution) small.
!     (100 to 1000 recommended)
   ierr = Zoltan_Set_Param(zz,'PARMETIS_ITR','100')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..CHECK_GRAPH
!     Level of error checking for graph input:
!        0 = no checking
!        1 = on-processor checking
!        2 = full checking (very slow)
   ierr = Zoltan_Set_Param(zz,'CHECK_GRAPH','1')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
!..SCATTER_GRAPH
!     Scatter graph data by distributing contiguous chunks of objects
!     to each processor before calling the partitioner.
!        0 = don't scatter
!        1 = scatter only if all objects are on a single processor
!        2 = scatter if at least one processor owns no objects
!        3 = always scatter
   ierr = Zoltan_Set_Param(zz,'SCATTER_GRAPH','2')
   call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   if (ierr > Ierr_out) Ierr_out = ierr
!
end subroutine zoltan_lb_param_graph
