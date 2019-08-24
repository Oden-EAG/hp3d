!--------------------------------------------------------------------
!> Purpose : routine performs a global isotropic h-refinement
!
!> @date Aug 2019
!--------------------------------------------------------------------
subroutine global_href_omp
!
   use error
   use data_structure3D
   use environment , only : QUIET_MODE
   use mpi_param   , only : RANK,ROOT
   use MPI         , only : MPI_COMM_WORLD
!
   implicit none
!
   integer :: nr_elem,mdle,i,kref,ierr
   real(8) :: MPI_Wtime,start_time,end_time
!
!..OpenMP auxiliary variables
   integer :: nrleafs,nrsons
   integer, allocatable :: elem_leaf,elem_kref
!
!--------------------------------------------------------------------
!
   allocate(elem_leaf(NRELES)); elem_leaf = 0
   allocate(elem_kref(NRELES)); elem_kref = 0
   nrleafs = 0
!
!$OMP PARALLEL PRIVATE(mdle,nrsons)
!$OMP DO REDUCTION(sum:nrleafs)
   do i=1,NRELES
      mdle = ELEM_ORDER(i)
      if (is_leaf(mdle)) then
         nrleafs = nrleafs + 1
         elem_leaf(i) = 1
         call get_isoref(mdle, elem_kref(mdle))
         ! get nrsons from breaking this mdle node
         ! --> this part must be sequential (depends on breaking previous elements)
      endif
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
!  maybe mark node ownership?
!
!  Q: does it matter in which order the elements are broken??
!     if not, compute index offsets in each MPI subdomain separately,
!     and do a reduction/prefix sum computation over indices based on
!     the MPI RANK
!
!..collect elements
   nr_elem = NRELES
!
!..break the elements
!  TODO OpenMP parallelize
   do i=1,nr_elem
      mdle = ELEM_ORDER(i)
      if (is_leaf(mdle)) then
         call get_isoref(mdle, kref)
         call break(mdle,kref)
      endif
   enddo
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   call refresh
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
!
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,100) ' # of elements broken         = ', nr_elem
      write(*,200) ' # of current elements, nodes = ', NRELES, NRNODS
      write(*,300) end_time-start_time
   endif
 100 format(A,I8)
 200 format(A,I8,', ',I9)
 300 format(' refresh    : ',f12.5,'  seconds')
!
end subroutine global_href_omp
!
!
!--------------------------------------------------------------------
!> Purpose : routine performs a global anisotropic h-refinement
!            note: works only for meshes with hexas/prisms,
!                  and breaking along the same axis
!                  (no mesh irregularity checks are done here)
!
!> Arguments:
!     in:    Krefxy -   1 (break in xy)  - hexa/prism mesh
!            Krefz  -   1 (break in  z)  - hexa/prism mesh
!
!> @date Aug 2019
!--------------------------------------------------------------------
subroutine global_href_aniso_omp(Krefxy,Krefz)
!
   use error
   use data_structure3D
   use environment , only : QUIET_MODE
   use mpi_param   , only : RANK,ROOT
   use MPI         , only : MPI_COMM_WORLD
!
   implicit none
!
   integer, intent(in) :: Krefxy,Krefz
!
   integer :: nr_elem,mdle,i,ierr
   integer :: kref_mdlb,kref_mdlp
   real(8) :: MPI_Wtime,start_time,end_time
!
!..check if valid refinement
   kref_mdlb = Krefxy*110 + Krefz
   kref_mdlp = Krefxy*10  + Krefz
!
!..collect elements
   nr_elem = NRELES
!
!..break the elements
!  TODO OpenMP parallelize
   do i=1,nr_elem
      mdle = ELEM_ORDER(i)
      if (is_leaf(mdle)) then
         select case(NODES(mdle)%type)
            case('mdlb'); call break(mdle,kref_mdlb)
            case('mdlp'); call break(mdle,kref_mdlp)
            case default
               write(*,*) 'global_href_aniso: unexpected node type. stop.'
               stop 1
         end select
      endif
   enddo
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   call refresh
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
!
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,100) ' # of elements broken         = ', nr_elem
      write(*,200) ' # of current elements, nodes = ', NRELES, NRNODS
      write(*,300) end_time-start_time
   endif
 100 format(A,I8)
 200 format(A,I8,', ',I9)
 300 format(' refresh    : ',f12.5,'  seconds')
!
end subroutine global_href_aniso_omp
