!> Purpose : recover the mesh from history file
subroutine reproduce_refinements
  use data_structure3D
  implicit none
  integer, allocatable ::  nhistory(:,:)
  integer :: mdle, kref, ic, nr_ref

  read(NHIST,*) mdle,kref
  ic=0
  do while (mdle.ne.0)
     ic=ic+1
     read(NHIST,*) mdle,kref
  enddo
  nr_ref = ic

  allocate(nhistory(2,nr_ref))

  rewind(NHIST)
  do ic=1,nr_ref
     read(NHIST,*) nhistory(1:2,ic)
  enddo

  rewind(NHIST)
  do ic=1,nr_ref
     call break(nhistory(1,ic),nhistory(2,ic))
  enddo

  call close
  deallocate(nhistory)
end subroutine reproduce_refinements

