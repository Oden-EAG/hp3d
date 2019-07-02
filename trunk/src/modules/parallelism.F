c----------------------------------------------------------------------
c
c   module name        - module parallelism
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jun 07
c
c   purpose            - module stores basic parallel execution params
c                        like processor ranks and number of processors
c
c
c----------------------------------------------------------------------
c
      module parallelism
c
      save

c-------------------------------------------------------------------
c     GENERAL INFORMATION
c-------------------------------------------------------------------
c  ...rank of this processor
      integer :: MYRANK
c  ...rank of this processor converted to a string
      character(len=4) :: PRINTRANK
c  ...total number of processors
      integer :: NRPROC
c  ...maximum number of processors handled
      integer, parameter :: MAXPROC = 256
c  ...0/1 graphics will show full domain, or only one subdomain
      integer :: RESTRICT_GRAPHICS

c  ...buffer size use for messaging, for loop involving all processors
      integer, parameter :: BUFFER_SIZE = 1000000
c  ...max allowed message size
      integer, parameter :: MAX_MSG_SIZE = 100000
c  ...tag of the first part of the message, containing size of the buffer
      integer, parameter :: MSG_SIZE_TAG = 1
c  ...integer denoting end of the buffer      
      integer,parameter   :: END_OF_BUFFER = -1

c  ...turn on/off checkpoints for parallel adaptivity
      integer, parameter :: PARALLEL_CHECKPOINTS = 0

c  ...total parallel execution time
      double precision :: DTIME_PARALLEL
      integer :: ITIME_PARALLEL
c
c-------------------------------------------------------------------
c     DATA FOR DOMAIN REDISTRIBUTION
c-------------------------------------------------------------------

      type elem_parallel_data_type
        integer :: processor_owner
      endtype elem_parallel_data_type

      type(elem_parallel_data_type), allocatable :: 
     .  ELEMS_PARALLEL_DATA(:)

      end module parallelism
