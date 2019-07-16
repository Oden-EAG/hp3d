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
   logical, save :: ZOLTAN_IS_INIT = .false.
!
   type(Zoltan_Struct), pointer, save :: zz
!
   contains
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_init
!     purpose:    initialize zoltan environment, and set parameters
!----------------------------------------------------------------------
   subroutine zoltan_w_init()
!
      integer(Zoltan_int)   :: ierr
      real   (Zoltan_float) :: ver
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
  100 format(/,A,F5.2)
!
!  ...create Zoltan memory, and set default parameters
      zz => Zoltan_Create(MPI_COMM_WORLD)
      if (.not. associated(zz)) then
         write(*,*) 'zoltan_w_init: Fatal error in Zoltan_Create.'
         return
      endif
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
!  ...set debug level for Zoltan library
      ierr = Zoltan_Set_Param(zz,'debug_level','4')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
!  ...set global IDs to one integer each (mdle)
      ierr = Zoltan_Set_Param(zz,'NUM_GID_ENTRIES','1')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
!  ...set weights to one floating point value per element
      ierr = Zoltan_Set_Param(zz,'OBJ_WEIGHT_DIM','1')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
!  ...omit local IDs
      ierr = Zoltan_Set_Param(zz,'NUM_LID_ENTRIES','0')
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
!     purpose:
!----------------------------------------------------------------------
   subroutine zoltan_w_query_coords(Dat,NrGIDs,NrLIDs,GID,LID, Coords,Ierr)
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int)   , intent(in)  :: NrGIDs,NrLIDs
      integer(Zoltan_int)   , intent(in)  :: GID(*)
      integer(Zoltan_int) :: LID(1) ! dummy declaration, do not use
      integer(Zoltan_double), intent(out) :: Coords(*)
      integer(Zoltan_int)   , intent(out) :: Ierr
      integer :: mdle,i,nrv
      real*8  :: x(NDIMEN), xnod(NDIMEN,MAXbrickH)
      mdle = GID(1)
      call nodcor_vert(mdle, xnod)
      nrv = nvert(NODES(mdle)%type)
      x(1:3) = 0.d0
      do i = 1,nrv
         x(1:3) = x(1:3) + xnod(1:3,i)
      enddo
      Coords(1:3) = x(1:3) / nrv
      Ierr = ZOLTAN_OK
   end subroutine zoltan_w_query_coords
!
!----------------------------------------------------------------------
!     function:   zoltan_w_query_subd
!                 Zoltan query function (ZOLTAN_NUM_OBJ_FN_TYPE)
!     purpose :   returns numbers of element in subdomain
!----------------------------------------------------------------------
   function zoltan_w_query_subd(Dat, Ierr)
      integer(Zoltan_int) :: zoltan_w_query_subd
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int), intent(out) :: Ierr
      integer(Zoltan_int) :: nreles_subd
      integer :: iel,mdle,subd
      mdle = 0; Nreles_subd = 0_Zoltan_int
      do iel=1,NRELES
         call nelcon(mdle, mdle)
         call get_subd(mdle, subd)
         if (subd .eq. RANK) then
            Nreles_subd = Nreles_subd + 1_Zoltan_int
         endif
      enddo
      zoltan_w_query_subd = Nreles_subd
      Ierr = ZOLTAN_OK
   end function zoltan_w_query_subd
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_query_elems
!                 Zoltan query function (ZOLTAN_OBJ_LIST_FN_TYPE)
!     purpose:    returns elements within subdomain, and their
!                 respective weights
!----------------------------------------------------------------------
   subroutine zoltan_w_query_elems(Dat,NrGIDs,NrLIDs, GIDs,LIDs, &
                                             Wgt_dim, Wgts,Ierr)
      integer(Zoltan_int) :: Dat(1) ! dummy declaration, do not use
      integer(Zoltan_int)  , intent(in)  :: NrGIDs,NrLIDs
      integer(Zoltan_int)  , intent(out) :: GIDs(*)
      integer(Zoltan_int) :: LIDs(1) ! dummy declaration, do not use
      integer(Zoltan_int)  , intent(in)  :: Wgt_dim
      integer(Zoltan_float), intent(out) :: Wgts(*)
      integer(Zoltan_int)  , intent(out) :: Ierr
      integer :: iel,j,mdle,subd
      mdle = 0; j = 0
      do iel=1,NRELES
         call nelcon(mdle, mdle)
         call get_subd(mdle, subd)
         if (subd .eq. RANK) then
            j = j + 1
            GIDs(j) = mdle
            Wgts(j) = 1.0_Zoltan_float
         endif
      enddo
      Ierr = ZOLTAN_OK
   end subroutine zoltan_w_query_elems
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_handle_err
!     purpose:    handle error code returned by a Zoltan function
!----------------------------------------------------------------------
   subroutine zoltan_w_handle_err(Ierr,Str)
      integer(Zoltan_int), intent(in) :: Ierr
      character(len=*)   , intent(in) :: Str
      if (Ierr .ne. ZOLTAN_OK) then
         write(*,*) str,': Ierr = ', Ierr
      endif
      ! ZOLTAN_WARN
      ! ZOLTAN_FATAL
      ! ZOLTAN_MEMERR
   end subroutine zoltan_w_handle_err
!
!
end module zoltan_wrapper
