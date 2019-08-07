!
#include "implicit_none.h"
!
!----------------------------------------------------------------------
!
!     module:              zoltan_wrapper
!     last modified:       July 2019
!
!----------------------------------------------------------------------
module zoltan_wrapper
!
   use data_structure3D
   use element_data
   use MPI
   use mpi_param
   use zoltan
!
   implicit none
!
   type(Zoltan_Struct), pointer, save :: zz
!
   logical, save :: ZOLTAN_IS_INIT = .false.
!
!..Load balancing strategy
!  0: nelcon (block partition via nelcon)
!  1: BLOCK  (block partition)
!  2: RANDOM (random partition)
!  3: RCB    (recursive coordinate bisection)
!  4: RIB    (recursive inertial bisection)
!  5: HSFC   (Hilbert space-filling curves)
   integer, save :: ZOLTAN_LB = 0
!
   contains
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_init
!     purpose:    initialize zoltan environment, and set parameters'
!
!     The following parameters [default values] are available:
!        (general)
!        NUM_GID_ENTRIES [1]
!        NUM_LID_ENTRIES [1]
!        DEBUG_LEVEL [1]
!        DEBUG_PROCESSOR [0]
!        DEBUG_MEMORY [1]
!        OBJ_WEIGHT_DIM [0]
!        EDGE_WEIGHT_DIM [0]
!        TIMER [wall]
!
!        (load balancing)
!        LB_METHOD [RCB]
!        LB_APPROACH [REPARTITION]
!        NUM_GLOBAL_PARTS [NUM_PROCS]
!        NUM_LOCAL_PARTS [1]
!        RETURN_LISTS [ALL]
!        REMAP [1]
!        IMBALANCE_TOL [1.1]
!        MIGRATE_ONLY_PROC_CHANGES [1]
!        AUTO_MIGRATE [FALSE]
!----------------------------------------------------------------------
   subroutine zoltan_w_init()
!
      integer(Zoltan_int)   :: ierr
      real   (Zoltan_float) :: ver
      character(len=4)      :: str
!
      if (ZOLTAN_IS_INIT) then
         write(*,*) 'zoltan_w_init: Zoltan has already been initialized.'
         return
      endif
!
!  ...initialize Zoltan environment
      ierr = Zoltan_Initialize(ver)
      call zoltan_w_handle_err(ierr,'Zoltan_Initialize')
      if (RANK .eq. ROOT) then
         write(*,100) 'Zoltan initialized sucessfully. Version = ', ver
      endif
  100 format(A,F5.2)
!
!  ...create Zoltan memory, and set default parameters
      zz => Zoltan_Create(MPI_COMM_WORLD)
      if (.not. associated(zz)) then
         write(*,*) 'zoltan_w_init: Fatal error in Zoltan_Create.'
         return
      endif
!
!  ...set global and local parts in partition
!     (should have been set by Zoltan_Create)
!      write(str,101) NUM_PROCS
!  101 format(I4)
!      ierr = Zoltan_Set_Param(zz,'NUM_GLOBAL_PARTS',str)
!      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!      ierr = Zoltan_Set_Param(zz,'NUM_LOCAL_PARTS','1')
!      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
!  ...set Zoltan query functions
!     Return dimensionality of the problem
      ierr = Zoltan_Set_Fn(zz,ZOLTAN_NUM_GEOM_FN_TYPE,zoltan_w_query_dim)
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Fn')
!     Return coordinates of an element
      ierr = Zoltan_Set_Fn(zz,ZOLTAN_GEOM_FN_TYPE,    zoltan_w_query_coords)
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Fn')
!     Return number of elements in subdomain
      ierr = Zoltan_Set_Fn(zz,ZOLTAN_NUM_OBJ_FN_TYPE, zoltan_w_query_subd)
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Fn')
!     Return global IDs of elements in subdomain
      ierr = Zoltan_Set_Fn(zz,ZOLTAN_OBJ_LIST_FN_TYPE,zoltan_w_query_elems)
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Fn')
!
!  ...set global IDs to one integer each (mdle)
      ierr = Zoltan_Set_Param(zz,'NUM_GID_ENTRIES','1')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
!  ...do not use local IDs
      ierr = Zoltan_Set_Param(zz,'NUM_LID_ENTRIES','0')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
!  ...set debug level for Zoltan library
      ierr = Zoltan_Set_Param(zz,'DEBUG_LEVEL','0')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
!  ...set weights to one floating point value per element
      ierr = Zoltan_Set_Param(zz,'OBJ_WEIGHT_DIM','0')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
!  ...set default timer to MPI_Wtime
      ierr = Zoltan_Set_Param(zz,'TIMER','wall')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
!  ...return only information about elements to be exported from a processor
!        ALL   : info about elements to be exported or imported
!        EXPORT: info about elements to be exported
!        IMPORT: info about elements to be imported
      ierr = Zoltan_Set_Param(zz,'RETURN_LISTS','ALL')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
      ZOLTAN_IS_INIT = .true.
!
   end subroutine zoltan_w_init
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_finalize
!     purpose:    close zoltan environment
!----------------------------------------------------------------------
   subroutine zoltan_w_finalize()
!
      if (.not. ZOLTAN_IS_INIT) then
         write(*,*) 'zoltan_w_finalize: Zoltan has not been initialized.'
         return
      endif
!
      call Zoltan_Destroy(zz)
!
      ZOLTAN_IS_INIT = .false.
!
   end subroutine zoltan_w_finalize
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_set_lb
!     purpose:    set load balancing strategy
!----------------------------------------------------------------------
   subroutine zoltan_w_set_lb(LB)
      integer, intent(in) :: LB
      integer(Zoltan_int) :: ierr
      ierr = ZOLTAN_OK
      select case(LB)
         case(0); ierr = Zoltan_Set_Param(zz,'LB_METHOD','NONE')
         case(1); call zoltan_lb_param_block(ierr)
         case(2); call zoltan_lb_param_random(ierr)
         case(3); call zoltan_lb_param_rcb(ierr)
         case(4); call zoltan_lb_param_rib(ierr)
         case(5); call zoltan_lb_param_hsfc(ierr)
         case default
            write(*,*) 'zoltan_w_set_lb: invalid param LB =', LB
            return
      end select
      call zoltan_w_handle_err(ierr,'zoltan_w_set_lb')
      ZOLTAN_LB = LB
   end subroutine zoltan_w_set_lb
!
!----------------------------------------------------------------------
!     function:   zoltan_w_query_dim
!                 Zoltan query function (ZOLTAN_NUM_GEOM_FN_TYPE)
!     purpose :   returns the number of dimensions of the problem
!----------------------------------------------------------------------
   function zoltan_w_query_dim(Dat, Ierr)
      integer(Zoltan_int) :: zoltan_w_query_dim
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(out) :: Ierr
      zoltan_w_query_dim = NDIMEN
      Ierr = ZOLTAN_OK
   end function zoltan_w_query_dim
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_query_coords
!                 Zoltan query function (ZOLTAN_GEOM_FN_TYPE)
!     purpose:    returns the coordinates of an element
!----------------------------------------------------------------------
   subroutine zoltan_w_query_coords(Dat,NrGIDs,NrLIDs,GID,LID, Coords,Ierr)
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int)   , intent(in)  :: NrGIDs,NrLIDs
      integer(Zoltan_int)   , intent(in)  :: GID(*)
      integer(Zoltan_int) :: LID(1) ! dummy declaration, do not use
      integer(Zoltan_double), intent(out) :: Coords(*)
      integer(Zoltan_int)   , intent(out) :: Ierr
      integer :: mdle,i,nrv
      real*8  :: x(NDIMEN), xnod(NDIMEN,8)
      mdle = GID(1)
      call nodcor_vert(mdle, xnod)
      nrv = nvert(NODES(mdle)%type)
      x(1:3) = 0.d0
      do i = 1,nrv
         x(1:3) = x(1:3) + xnod(1:3,i)
      enddo
      Coords(1:3) = x(1:3) / nrv
      !write(*,200) 'Mdle = ', mdle,', Coords = ', x(1:3)/nrv
  200 format(A,I5,A,3F6.2)
      Ierr = ZOLTAN_OK
   end subroutine zoltan_w_query_coords
!
!----------------------------------------------------------------------
!     function:   zoltan_w_query_subd
!                 Zoltan query function (ZOLTAN_NUM_OBJ_FN_TYPE)
!     purpose :   returns number of elements in subdomain
!----------------------------------------------------------------------
   function zoltan_w_query_subd(Dat, Ierr)
      integer(Zoltan_int) :: zoltan_w_query_subd
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(out) :: Ierr
      integer(Zoltan_int) :: nreles_subd
      integer :: iel,mdle,subd
      Nreles_subd = 0_Zoltan_int
      do iel=1,NRELES
         mdle = ELEM_ORDER(iel)
         call get_subd(mdle, subd)
         if (subd .eq. RANK) then
            Nreles_subd = Nreles_subd + 1_Zoltan_int
         endif
      enddo
      zoltan_w_query_subd = Nreles_subd
      Ierr = ZOLTAN_OK
   end function zoltan_w_query_subd
!
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_query_elems
!                 Zoltan query function (ZOLTAN_OBJ_LIST_FN_TYPE)
!     purpose:    returns elements within subdomain, and their
!                 respective weights for load balancing
!----------------------------------------------------------------------
   subroutine zoltan_w_query_elems(Dat,NrGIDs,NrLIDs, GIDs,LIDs, &
                                             Wgt_dim, Wgts,Ierr)
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(in)  :: NrGIDs,NrLIDs
      integer(Zoltan_int), intent(out) :: GIDs(*)
      integer(Zoltan_int) :: LIDs(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(in)  :: Wgt_dim
      real(Zoltan_float) , intent(out) :: Wgts(*)
      integer(Zoltan_int), intent(out) :: Ierr
      integer :: iel,j,mdle,subd
      j = 0
      do iel=1,NRELES
         mdle = ELEM_ORDER(iel)
         call get_subd(mdle, subd)
         if (subd .eq. RANK) then
            j = j + 1
            GIDs(j) = mdle
            if (Wgt_dim > 0) Wgts(j) = 1.0_Zoltan_float
         endif
      enddo
      Ierr = ZOLTAN_OK
   end subroutine zoltan_w_query_elems
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_partition
!     purpose:    perform partitioning of current domain
!----------------------------------------------------------------------
   subroutine zoltan_w_partition(Mdle_subd)
      integer, intent(out) :: Mdle_subd(NRELES)
      logical :: changes,print_stats
      integer(Zoltan_int) :: ierr,nrGIDs,nrLIDs,nrImp,nrExp
      integer(Zoltan_int),pointer,dimension(:) :: impGIDs,impLIDs,impProcs,impParts
      integer(Zoltan_int),pointer,dimension(:) :: expGIDs,expLIDs,expProcs,expParts
      integer :: iel,i,mdle,subd,src,count,ierr_mpi
!
!  ...find new partition
      ierr = Zoltan_LB_Partition(zz, changes,nrGIDs,nrLIDs,                   &
                                     nrImp,impGIDs,impLIDs,impProcs,impParts, &
                                     nrExp,expGIDs,expLIDs,expProcs,expParts )
      call zoltan_w_handle_err(ierr,'Zoltan_LB_Partition')
!
      !write(*,*) 'zoltan_w_partition:'
      !write(*,300) '   changes  = ', changes
      !write(*,301) '   nrImp    = ', nrImp
      !write(*,301) '   nrExp    = ', nrExp
      !if (nrImp > 0) write(*,310) '   impProcs = ', impProcs
      !if (nrExp > 0) write(*,320) '   expProcs = ', expProcs
  300 format(A,L5)
  301 format(A,I5)
  310 format(A,<nrImp>I5)
  320 format(A,<nrExp>I5)
!
!  ...collect new subdomains for mesh
      Mdle_subd(1:NRELES) = 0
!
      count = 1
      do iel=1,NRELES
         mdle = ELEM_ORDER(iel)
         call get_subd(mdle, subd)
         src = subd
         if (RANK .eq. src) then
!        ...determine if this mdle node is exported to new processor
            do i=1,nrExp
               if (mdle .eq. expGIDs(i)) subd = expProcs(i)
            enddo
         endif
         call MPI_BCAST (subd,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr_mpi)
         Mdle_subd(iel) = subd
      enddo
      !write(*,*) 'Suggested new partition is as follows'
      !write(*,330) Mdle_subd
  330 format(<NRELES>I4)
!
!  ...deallocate internal data structures
      ierr = Zoltan_LB_Free_Part(impGIDs,impLIDs,impProcs,impParts)
      call zoltan_w_handle_err(ierr,'Zoltan_LB_Free_Part')
      ierr = Zoltan_LB_Free_Part(expGIDs,expLIDs,expProcs,expParts)
      call zoltan_w_handle_err(ierr,'Zoltan_LB_Free_Part')
!
   end subroutine zoltan_w_partition
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_eval
!     purpose:    evaluates the quality of the current partition
!----------------------------------------------------------------------
   subroutine zoltan_w_eval()
      logical :: print_stats
      integer(Zoltan_int) :: ierr
      !if (RANK .eq. ROOT) write(*,*) 'Evaluating quality of current partition...'
      print_stats=.true.
      ierr = Zoltan_LB_Eval(zz,print_stats)
      call zoltan_w_handle_err(ierr,'Zoltan_LB_Eval')
!
   end subroutine zoltan_w_eval
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_handle_err
!     purpose:    handle error code returned by a Zoltan function
!----------------------------------------------------------------------
   subroutine zoltan_w_handle_err(Ierr,Str)
      integer(Zoltan_int), intent(in) :: Ierr
      character(len=*)   , intent(in) :: Str
      if (Ierr .ne. ZOLTAN_OK) then
         select case(Ierr)
            case(ZOLTAN_WARN)  ; write(*,400) Str,Ierr,' (ZOLTAN_WARN)'
            case(ZOLTAN_FATAL) ; write(*,400) Str,Ierr,' (ZOLTAN_FATAL)'
            case(ZOLTAN_MEMERR); write(*,400) Str,Ierr,' (ZOLTAN_MEMERR)'
         end select
     400 format(' zoltan_w_handle_err: ',A,': Ierr = ',I2,A)
      endif
   end subroutine zoltan_w_handle_err
!
!
end module zoltan_wrapper
