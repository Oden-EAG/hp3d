module mesh_partition
  implicit none
      integer, pointer     :: vwgt=>null(), vsize=>null()!, options=>null()
      integer, allocatable :: eptr(:),eind(:),epart(:),npart(:),options(:)
      real(8), pointer     :: tpwgts=>null()

  contains

    subroutine mesh_partition_0
      use geometry_polydpg

      implicit none

      ! integer, pointer     :: vwgt=>null(), vsize=>null()!, options=>null()
      ! integer, allocatable :: eptr(:),eind(:),epart(:),npart(:),options(:)
      ! real(8), pointer     :: tpwgts=>null()
      integer :: nfile_elem,nfile_node,ne,nn,ncommon,npel,objval,nparti

    ! Open output files
      nfile_elem=1
      open(unit=nfile_elem,file='output/epart', &
              form='formatted',access='sequential',status='unknown')

      nfile_node=2
      open(unit=nfile_node,file='output/npart', &
              form='formatted',access='sequential',status='unknown')
    ! For tetrahedral meshes, each pair of elements sharing a face have 3 common nodes
      ncommon = 3

      ne=2 ! NUMC
      nn=6 ! NUMV
      npel=4

      allocate(eptr(ne+1),eind(ne*npel),epart(ne),npart(nn),options(0:40))

      eptr=(/0,4,8/)
      eind=(/0,1,2,3,1,4,5,2/)

      nparti=2 ! NUMC/8

      call METIS_SetDefaultOptions(options)
      options(17) = 0
    ! ! METIS mesh partitioning option array
    !   options(METIS_OPTION_OBJTYPE)   = METIS_OBJTYPE_CUT
    !   options(METIS_OPTION_NUMBERING) = 1
      call METIS_PartMeshDual(ne,nn,eptr,eind,vwgt,vsize,ncommon,nparti,tpwgts,options,objval,epart,npart)
      ! call METIS_PartMeshDual(ne,nn,eptr,eind,(/1,1/),(/1,1/),ncommon,nparti,(/0.5d0,0.5d0/),options,objval,epart,npart)

      write(nfile_elem,*) epart
      write(nfile_node,*) npart

      write(*,*) epart
      write(*,*) npart
    ! close output files
      close(nfile_elem)
      close(nfile_node)

    end subroutine mesh_partition_0

end module mesh_partition