      module communicators
      use parallelism , only : MAXPROC
      save
c  ...group involving only current process
      integer :: MYRANK_GROUP(1:MAXPROC)
c  ...communicator involving only current process
      integer :: COMM_MYRANK(1:MAXPROC)
      integer :: COMM_SHIFTED
      end
