!
#include "typedefs.h"
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
!  6: GRAPH  (Graph partitioners: ParMETIS,PT-Scotch)
!  7: FIBER  (Custom partitioner for waveguide geometry)
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
!      character(len=4)      :: str
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
!      ierr = Zoltan_Set_Fn(zz,ZOLTAN_GEOM_FN_TYPE,zoltan_w_query_coords)
      ierr = Zoltan_Set_Fn(zz,ZOLTAN_GEOM_MULTI_FN_TYPE, zoltan_w_query_coords_omp)
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Fn')
!     Return number of elements in subdomain
      ierr = Zoltan_Set_Fn(zz,ZOLTAN_NUM_OBJ_FN_TYPE, zoltan_w_query_subd)
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Fn')
!     Return global IDs of elements in subdomain
      ierr = Zoltan_Set_Fn(zz,ZOLTAN_OBJ_LIST_FN_TYPE,zoltan_w_query_elems)
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Fn')
!     Return number of element neighbors for a given element
!      ierr = Zoltan_Set_Fn(zz,ZOLTAN_NUM_EDGES_FN_TYPE,zoltan_w_query_num_neig)
      ierr = Zoltan_Set_Fn(zz,ZOLTAN_NUM_EDGES_MULTI_FN_TYPE,zoltan_w_query_num_neig_omp)
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Fn')
!     Return lists of neighbors, their owners, and optionally face weights
!     for elements sharing a face with a given element
!      ierr = Zoltan_Set_Fn(zz,ZOLTAN_EDGE_LIST_FN_TYPE,zoltan_w_query_neig_list)
      ierr = Zoltan_Set_Fn(zz,ZOLTAN_EDGE_LIST_MULTI_FN_TYPE,zoltan_w_query_neig_list_omp)
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
!  ...set object weights to one floating point value per element
      ierr = Zoltan_Set_Param(zz,'OBJ_WEIGHT_DIM','1')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
!  ...set edge weights to one floating point value per neighbor
      ierr = Zoltan_Set_Param(zz,'EDGE_WEIGHT_DIM','1')
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
         case(6); call zoltan_lb_param_graph(ierr)
         case(7); ierr = Zoltan_Set_Param(zz,'LB_METHOD','NONE')
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
      real(8) :: x(3), xnod(3,8)
      mdle = GID(1)
      call nodcor_vert(mdle, xnod)
      nrv = nvert(NODES(mdle)%ntype)
      x(1:3) = 0.d0
      do i = 1,nrv
         x(1:3) = x(1:3) + xnod(1:3,i)
      enddo
      Coords(1:3) = x(1:3) / nrv
#if DEBUG_MODE
      write(*,200) 'Mdle = ', mdle,', Coords = ', x(1:3)/nrv
  200 format(A,I5,A,3F6.2)
#endif
      Ierr = ZOLTAN_OK
   end subroutine zoltan_w_query_coords
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_query_coords_omp
!                 Zoltan query function (ZOLTAN_GEOM_MULTI_FN_TYPE)
!     purpose:    returns the coordinates of each element in list
!----------------------------------------------------------------------
   subroutine zoltan_w_query_coords_omp(Dat,NrGIDs,NrLIDs,NumObj, &
                                            GIDs,LIDs,NumDim, Coords,Ierr)
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int)   , intent(in)  :: NrGIDs,NrLIDs,NumObj
      integer(Zoltan_int)   , intent(in)  :: GIDs(*)
      integer(Zoltan_int) :: LIDs(1) ! dummy declaration, do not use
      integer(Zoltan_int)   , intent(in)  :: NumDim
      integer(Zoltan_double), intent(out) :: Coords(*)
      integer(Zoltan_int)   , intent(out) :: Ierr
      integer :: mdle,i,k,nrv
      real(8) :: x(3), xnod(3,8)
!
!$OMP PARALLEL DO PRIVATE(mdle,i,nrv,x,xnod)
      do k = 1,NumObj
         mdle = GIDs(k)
         call nodcor_vert(mdle, xnod)
         nrv = nvert(NODES(mdle)%ntype)
         x(1:3) = 0.d0
         do i = 1,nrv
            x(1:3) = x(1:3) + xnod(1:3,i)
         enddo
         i = (k-1)*NumDim
         Coords(i+1:i+3) = x(1:3) / nrv
#if DEBUG_MODE
         !$OMP CRITICAL
         write(*,200) 'Mdle = ', mdle,', Coords = ', x(1:3)/nrv
     200 format(A,I5,A,3F6.2)
         !$OMP END CRITICAL
#endif
      enddo
!$OMP END PARALLEL DO
!
      Ierr = ZOLTAN_OK
   end subroutine zoltan_w_query_coords_omp
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
      zoltan_w_query_subd = NRELES_SUBD
      Ierr = ZOLTAN_OK
   end function zoltan_w_query_subd
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
      integer :: nrdofb(NR_PHYSA)
      integer :: iel,mdle
!$OMP PARALLEL DO PRIVATE(mdle)
      do iel=1,NRELES_SUBD
         mdle = ELEM_SUBD(iel)
         GIDs(iel) = mdle
         if (Wgt_dim > 0) then
            call zoltan_w_get_nrdofb(Mdle, nrdofb)
            Wgts(iel) = sum(nrdofb)
         endif
      enddo
!$OMP END PARALLEL DO
      Ierr = ZOLTAN_OK
   end subroutine zoltan_w_query_elems
!
!----------------------------------------------------------------------
!     function:   zoltan_w_query_num_neig
!                 Zoltan query function (ZOLTAN_NUM_EDGES_FN_TYPE)
!     purpose :   returns number of neighbors of an element
!----------------------------------------------------------------------
   function zoltan_w_query_num_neig(Dat,NrGIDs,NrLIDs,GID,LID, Ierr)
      integer(Zoltan_int) :: zoltan_w_query_num_neig
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(in)  :: NrGIDs,NrLIDs
      integer(Zoltan_int), intent(in)  :: GID(1)
      integer(Zoltan_int) :: LID(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(out) :: Ierr
      integer :: mdle,i,j,num_neig
      integer :: neig_list(4,6)
      mdle = GID(1)
      call zoltan_w_find_neig(mdle, neig_list)
      num_neig = 0
      do i=1,6
         do j=1,4
            if (neig_list(j,i) .ne. 0) num_neig = num_neig + 1
         enddo
      enddo
      zoltan_w_query_num_neig = num_neig
      Ierr = ZOLTAN_OK
   end function zoltan_w_query_num_neig
!
!----------------------------------------------------------------------
!   subroutine:   zoltan_w_query_num_neig_omp
!                 Zoltan query function (ZOLTAN_NUM_EDGES_MULTI_FN_TYPE)
!     purpose :   returns number of neighbors of each element in list
!----------------------------------------------------------------------
   subroutine zoltan_w_query_num_neig_omp(Dat,NrGIDs,NrLIDs,NumObj,GIDs,LIDs, NumNeig,Ierr)
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(in)  :: NrGIDs,NrLIDs,NumObj
      integer(Zoltan_int), intent(in)  :: GIDs(*)
      integer(Zoltan_int) :: LIDs(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(out) :: NumNeig(*)
      integer(Zoltan_int), intent(out) :: Ierr
      integer :: mdle,i,j,k,num_neig
      integer :: neig_list(4,6)
!$OMP PARALLEL DO PRIVATE(i,j,k,mdle,neig_list,num_neig)
      do k = 1,NumObj
         mdle = GIDs(k)
         call zoltan_w_find_neig(mdle, neig_list)
         num_neig = 0
         do i=1,6
            do j=1,4
               if (neig_list(j,i) .ne. 0) num_neig = num_neig + 1
            enddo
         enddo
         NumNeig(k) = num_neig
      enddo
!$OMP END PARALLEL DO
      Ierr = ZOLTAN_OK
   end subroutine zoltan_w_query_num_neig_omp
!
!----------------------------------------------------------------------
!   subroutine:   zoltan_w_query_neig_list
!                 Zoltan query function (ZOLTAN_EDGE_LIST_MULTI_FN_TYPE)
!     purpose :   returns list of neighbors of an element, and their
!                 respective face (connectivity) weights
!----------------------------------------------------------------------
   subroutine zoltan_w_query_neig_list(Dat,NrGIDs,NrLIDs,GID,LID, Neig_GIDs,Neig_subd,  &
                                                         Wgt_dim, Wgts,Ierr)
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(in)  :: NrGIDs,NrLIDs
      integer(Zoltan_int), intent(in)  :: GID(1)
      integer(Zoltan_int) :: LID(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(out) :: Neig_GIDs(*)
      integer(Zoltan_int), intent(out) :: Neig_subd(*)
      integer(Zoltan_int), intent(in)  :: Wgt_dim
      real(Zoltan_float) , intent(out) :: Wgts(*)
      integer(Zoltan_int), intent(out) :: Ierr
      integer :: mdle,i,j,k,num_neig,neig,subd
      integer :: ndofH,ndofE,ndofV,ndofQ
      integer :: neig_list(4,6)
      integer :: iprint
      iprint=0
      mdle = GID(1)
      call zoltan_w_find_neig(mdle, neig_list)
      num_neig = 0
      do i=1,6
         do j=1,4
            neig = neig_list(j,i)
            if (neig .ne. 0) then
               num_neig = num_neig + 1
               call get_subd(neig, subd)
               Neig_GIDs(num_neig) = neig
               Neig_subd(num_neig) = subd
               if (Wgt_dim > 0) then
                  call ndof_nod(NODES(neig)%ntype,NODES(neig)%order, ndofH,ndofE,ndofV,ndofQ)
                  Wgts(num_neig) = 0.0_Zoltan_float
                  do k=1,NR_PHYSA
                     if (.not. PHYSAm(k)) cycle
                     select case(DTYPE(k))
                        case('contin') ; Wgts(num_neig) = Wgts(num_neig) + ndofH*NR_COMP(k)
                        case('tangen') ; Wgts(num_neig) = Wgts(num_neig) + ndofE*NR_COMP(k)
                        case('normal') ; Wgts(num_neig) = Wgts(num_neig) + ndofV*NR_COMP(k)
                     end select
                  enddo
#if DEBUG_MODE
                  if (Wgts(num_neig) .lt. 1.0_Zoltan_float) then
                     write(*,*) 'zoltan_w_query_neig_list: Wgts(num_neig)',Wgts(num_neig)
                     Wgts(num_neig) = 1.0_Zoltan_float
                  endif
#endif
               endif
#if DEBUG_MODE
               if (Is_inactive(neig)) then
                  write(*,*) 'zoltan_w_query_neig_list: neig not active'; stop
               endif
               if (subd .lt. 0 .or. subd .ge. NUM_PROCS) then
                  write(*,*) 'zoltan_w_query_neig_list: neig invalid subd'; stop
               endif
               if (neig .eq. mdle) then
                  write(*,*) 'zoltan_w_query_neig_list: neig .eq. mdle'; stop
               endif
#endif
            endif
         enddo
      enddo
#if DEBUG_MODE
      if (iprint .eq. 1) then
         write(*,500) 'Neighbors of mdle = ', mdle
         do i=1,num_neig
            write(*,510) '  neig, subd = ', Neig_GIDs(i), Neig_subd(i)
         enddo
500      format(A,I10)
510      format(A,I10,I4)
      endif
#endif
      Ierr = ZOLTAN_OK
   end subroutine zoltan_w_query_neig_list
!
!----------------------------------------------------------------------
!   subroutine:   zoltan_w_query_neig_list_omp
!                 Zoltan query function (ZOLTAN_EDGE_LIST_FN_TYPE)
!     purpose :   returns list of neighbors of each element in list,
!                 and their respective face (connectivity) weights
!----------------------------------------------------------------------
   subroutine zoltan_w_query_neig_list_omp(Dat,NrGIDs,NrLIDs,NumObj,GIDs,LIDs, &
                                           NumEdges, Neig_GIDs,Neig_subd,      &
                                           Wgt_dim,  Wgts,Ierr)
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(in)  :: NrGIDs,NrLIDs,NumObj
      integer(Zoltan_int), intent(in)  :: GIDs(*)
      integer(Zoltan_int) :: LIDs(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(in)  :: NumEdges(*)
      integer(Zoltan_int), intent(out) :: Neig_GIDs(*)
      integer(Zoltan_int), intent(out) :: Neig_subd(*)
      integer(Zoltan_int), intent(in)  :: Wgt_dim
      real(Zoltan_float) , intent(out) :: Wgts(*)
      integer(Zoltan_int), intent(out) :: Ierr
      integer :: mdle,i,j,k,l,n,num_neig,neig,subd,nsum,ndofH,ndofE,ndofV,ndofQ
      integer :: neig_list(4,6)
      integer :: noffset(NumObj)
      integer :: iprint
      iprint=0
      noffset = 0; nsum = 0
      do k=1,NumObj
         noffset(k) = nsum
         nsum = nsum + NumEdges(k)
      enddo
!$OMP PARALLEL DO PRIVATE(i,j,k,l,n,mdle,neig_list,num_neig,neig,subd,ndofH,ndofE,ndofV,ndofQ)
      do k=1,NumObj
         mdle = GIDs(k)
         call zoltan_w_find_neig(mdle, neig_list)
         num_neig = 0
         do i=1,6
            do j=1,4
               neig = neig_list(j,i)
               if (neig .ne. 0) then
                  num_neig = num_neig + 1
                  call get_subd(neig, subd)
                  n = noffset(k) + num_neig
                  Neig_GIDs(n) = neig
                  Neig_subd(n) = subd
                  if (Wgt_dim > 0) then
                     call ndof_nod(NODES(neig)%ntype,NODES(neig)%order, ndofH,ndofE,ndofV,ndofQ)
                     Wgts(n) = 0.0_Zoltan_float
                     do l=1,NR_PHYSA
                        if (.not. PHYSAm(l)) cycle
                        select case(DTYPE(l))
                           case('contin') ; Wgts(n) = Wgts(n) + ndofH*NR_COMP(l)
                           case('tangen') ; Wgts(n) = Wgts(n) + ndofE*NR_COMP(l)
                           case('normal') ; Wgts(n) = Wgts(n) + ndofV*NR_COMP(l)
                        end select
                     enddo
#if DEBUG_MODE
                     if (Wgts(n) .lt. 1.0_Zoltan_float) then
                        !$OMP CRITICAL
                        write(*,*) 'zoltan_w_query_neig_list_omp: Wgts(n)',Wgts(n)
                        !$OMP END CRITICAL
                        Wgts(n) = 1.0_Zoltan_float
                     endif
#endif
                  endif
#if DEBUG_MODE
                  if (Is_inactive(neig)) then
                     !$OMP CRITICAL
                     write(*,*) 'zoltan_w_query_neig_list_omp: neig not active'; stop
                     !$OMP END CRITICAL
                  endif
                  if (subd .lt. 0 .or. subd .ge. NUM_PROCS) then
                     !$OMP CRITICAL
                     write(*,*) 'zoltan_w_query_neig_list_omp: neig invalid subd'; stop
                     !$OMP END CRITICAL
                  endif
                  if (neig .eq. mdle) then
                     !$OMP CRITICAL
                     write(*,*) 'zoltan_w_query_neig_list_omp: neig .eq. mdle'; stop
                     !$OMP END CRITICAL
                  endif
#endif
               endif
            enddo
         enddo
#if DEBUG_MODE
         if (iprint .eq. 1) then
            !$OMP CRITICAL
            write(*,600) 'Neighbors of mdle = ', mdle
            do i=1,num_neig
               write(*,610) '  neig, subd = ', Neig_GIDs(noffset(k)+i), Neig_subd(noffset(k)+i)
            enddo
            !$OMP END CRITICAL
        600 format(A,I10)
        610 format(A,I10,I4)
         endif
#endif
      enddo
!$OMP END PARALLEL DO
      Ierr = ZOLTAN_OK
   end subroutine zoltan_w_query_neig_list_omp
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
      character(16) :: fmt
!  ...find new partition
      ierr = Zoltan_LB_Partition(zz, changes,nrGIDs,nrLIDs,                   &
                                     nrImp,impGIDs,impLIDs,impProcs,impParts, &
                                     nrExp,expGIDs,expLIDs,expProcs,expParts )
      call zoltan_w_handle_err(ierr,'Zoltan_LB_Partition')
!
      print_stats = .false.
      if (print_stats) then
         write(*,*) 'zoltan_w_partition:'
         write(*,300) '   changes  = ', changes
         write(*,301) '   nrImp    = ', nrImp
         write(*,301) '   nrExp    = ', nrExp
!
         if (nrImp > 0) then
            write(fmt,'("(A,",I0,"I5)")') nrImp
            write(6,fmt) '   impProcs = ', impProcs
         endif
!
         if (nrExp > 0) then
            write(fmt,'("(A,",I0,"I5)")') nrExp
            write(6,fmt) '   expProcs = ', expProcs
         endif
!
  300    format(A,L5)
  301    format(A,I5)
      endif
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
!--------------------------------------------------------------------------
!> @brief Finds neighbors (up to 4, for an h4 refined face) across faces
!! @param[in]  Mdle      - middle node
!! @param[out] Neig_list - neighbors
!!
!> @date Mar 2023
!--------------------------------------------------------------------------
   subroutine zoltan_w_find_neig(Mdle, Neig_list)
!
      integer,                 intent(in)  :: Mdle
      integer, dimension(4,6), intent(out) :: Neig_list
!
      integer, dimension(27) :: nodesl,norientl
      integer, dimension(2)  :: neig,nsid_list,norient_list
      integer :: i,j,k,l,nod,nodson,nrneig,mdle_neig,ntype
!
!  ...initialize
      Neig_list(1:4,1:6)=0
!
!---------------------------------------------------
! Step 0: initial mesh elements
!---------------------------------------------------
      ntype=NODES(Mdle)%ntype
      if (is_root(Mdle)) then
         do i=1,nface(ntype)
            mdle_neig = ELEMS(Mdle)%neig(i)
!
            if (mdle_neig .eq. 0) then
               cycle
            elseif (NODES(mdle_neig)%ref_kind .eq. 0) then
               Neig_list(1,i) = ELEMS(Mdle)%neig(i)
            else
               nod = ELEMS(Mdle)%nodes(nvert(ntype)+nedge(ntype)+i)
!
               l = 0
               do j=1,NODES(nod)%nr_sons
                  nodson = Son(nod,j)
!
                  select case(NODES(nodson)%ntype)
                  case(MDLQ,MDLT)
                     l = l+1
                  case default
                     cycle
                  end select
!
                  call neig_face(nodson, nrneig,neig,nsid_list,norient_list)
                  do k=1,nrneig
                     if (NODES(neig(k))%father .eq. mdle_neig) Neig_list(l,i) = neig(k)
                  enddo
               enddo
            endif
         enddo
         return
      endif
!
!---------------------------------------------------
! Step 1: use neig_face
!---------------------------------------------------
      call elem_nodes(Mdle, nodesl,norientl)
      do i=1,nface(ntype)
         nod = nodesl(nvert(ntype)+nedge(ntype)+i)
         call neig_face(nod, nrneig,neig,nsid_list,norient_list)
         select case (nrneig)
            case(2)
               ! pick the other one
               if (Mdle.eq.neig(1)) then
                  mdle_neig = neig(2)
               else
                  mdle_neig = neig(1)
               endif
            case default
               cycle
         end select
!
         if (NODES(nod)%ref_kind .eq. 0) then
            Neig_list(1,i) = mdle_neig
         else
            l = 0
            do j=1,NODES(nod)%nr_sons
!
               nodson = Son(nod,j)
!
               select case(NODES(nodson)%ntype)
               case(MDLQ,MDLT)
                  l = l+1
               case default
                  cycle
               end select
!
               call neig_face(nodson, nrneig,neig,nsid_list,norient_list)
               do k=1,nrneig
                  if (NODES(neig(k))%father .eq. mdle_neig) Neig_list(l,i) = neig(k)
               enddo
            enddo
         endif
      enddo
!
   end subroutine zoltan_w_find_neig
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_get_nrdofb
!     purpose:    computes number of bubble dofs to be used for weights
!                 in graph partitioning
!----------------------------------------------------------------------
subroutine zoltan_w_get_nrdofb(Mdle, Nrdofb)
   integer, intent(in)  :: Mdle
   integer, intent(out) :: Nrdofb(NR_PHYSA)
   integer :: i,ntype,nord,nrdofH,nrdofE,nrdofV,nrdofQ
   integer :: ncase(NR_PHYSA)
!
   Nrdofb(1:NR_PHYSA)=0
!
!..set node to be middle node number (global id in NODES)
   call decod(NODES(Mdle)%case,2,NR_PHYSA, ncase)
!
!..number of dofs for a SINGLE H1,H(curl),H(div),L2 component
   ntype = NODES(Mdle)%ntype
   call ndof_nod(ntype,nord, nrdofH,nrdofE,nrdofV,nrdofQ)
!
!..number of dofs for ALL H1,H(curl),H(div),L2 components
   do i=1,NR_PHYSA
!  ...skip if variable is not supported on middle node
      if (ncase(i) .eq. 0) cycle
!  ...skip if variable is deactivated
      if (.not. PHYSAm(i)) cycle
!  ...skip interface variables
      if (PHYSAi(i))       cycle
!
      select case(DTYPE(i))
         case('contin') ; Nrdofb(i)=nrdofH*NR_COMP(i)
         case('tangen') ; Nrdofb(i)=nrdofE*NR_COMP(i)
         case('normal') ; Nrdofb(i)=nrdofV*NR_COMP(i)
         case('discon') ; Nrdofb(i)=nrdofQ*NR_COMP(i)
      end select
   enddo
!
#if DEBUG_MODE
   if (sum(Nrdofb) .lt. 1) then
      write(*,*) 'zoltan_w_get_nrdofb: sum(Nrdofb)=',sum(Nrdofb)
      Nrdofb(1) = 1
   endif
#endif
!
end subroutine zoltan_w_get_nrdofb
!
!
end module zoltan_wrapper
