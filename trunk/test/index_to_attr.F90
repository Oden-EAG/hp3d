!> @brief Tests index_to_attr subroutine
!> @date Sep 2023
program test_index_to_attr
!
   use mpi_param
   use mpi_wrapper
   use physics
!
   implicit none
!
   integer :: NPASS
!
   integer :: attr,comp,index
!
   QUIET_MODE = .true.
!
!..initialize MPI environment
   call mpi_w_init
!
!..initialize physics, geometry, etc.
   call initialize
!
!..this test is set up for sequential computation
   if (NUM_PROCS > 1) then
      write(*,*) 'test_index_to_attr: NUM_PROCS>1'
      NPASS = 1; goto 99
   endif
!
   NPASS = 1
!
!..physics attributes components
!  attr | comp | index | description
!     1 |  1-2 |   1-2 | H1
!     2 |  1-2 |   3-4 | H(curl)
!     3 |  1-2 |   5-6 | H(div)
!     4 |  1-2 |   7-8 | L2
!
   index = 1
   call index_to_attr(index, attr,comp)
   if (attr .ne. 1) NPASS = 0
   if (comp .ne. 1) NPASS = 0
!
   index = 2
   call index_to_attr(index, attr,comp)
   if (attr .ne. 1) NPASS = 0
   if (comp .ne. 2) NPASS = 0
!
   index = 3
   call index_to_attr(index, attr,comp)
   if (attr .ne. 2) NPASS = 0
   if (comp .ne. 1) NPASS = 0
!
   index = 4
   call index_to_attr(index, attr,comp)
   if (attr .ne. 2) NPASS = 0
   if (comp .ne. 2) NPASS = 0
!
   index = 5
   call index_to_attr(index, attr,comp)
   if (attr .ne. 3) NPASS = 0
   if (comp .ne. 1) NPASS = 0
!
   index = 6
   call index_to_attr(index, attr,comp)
   if (attr .ne. 3) NPASS = 0
   if (comp .ne. 2) NPASS = 0
!
   index = 7
   call index_to_attr(index, attr,comp)
   if (attr .ne. 4) NPASS = 0
   if (comp .ne. 1) NPASS = 0
!
   index = 8
   call index_to_attr(index, attr,comp)
   if (attr .ne. 4) NPASS = 0
   if (comp .ne. 2) NPASS = 0
!
 99 continue
!
!..finalize MPI environment
   call mpi_w_finalize
!
   if (NPASS.ne.1) stop 1
!
end program test_index_to_attr
!
!----------------------------------------------------------------------
!> @brief initialization
subroutine initialize
!
   use environment
   use control
   use GMP
   use data_structure3D
   use refinements
!
   implicit none
!
   logical :: q
!
   q = QUIET_MODE
   QUIET_MODE = .true.
!
#if HP3D_USE_OPENMP
   call omp_set_num_threads(1)
#endif
!
!..initialize refinements arrays
   call init_refinements('../files/ref')
!
!..generate constraint arrays
   call init_cnstr
!
!..read control file
   call read_control('../files/control/control_0')
!
!..set GMP parameters
!                                 NDIM //    MANDIM //
   call set_gmp_parameters(         3 ,        3 , &
!                                MAXSU //     MAXNP //     MAXNC //
                                    0 ,        8 ,          12 , &
!                                MAXTR //     MAXRE //     MAXBT //
                                    0 ,        6 ,           0 , &
!                                MAXHE //     MAXTE //     MAXPY //
                                    1 ,        0 ,           0)
!
!..set hp3D parameters
!                      NRCOMS, NRRHS
   call set_parameters(     1,     1)
!
!..read geometry file
   call read_geometry('../files/mesh/hexa_orient_0')
!
!..read physics file and generate mesh data structures
   call hp3gen('../files/physics/physics_2')
!
   QUIET_MODE = q
!
end subroutine initialize
!
!----------------------------------------------------------------------
!> @brief initialization
subroutine set_initial_mesh(Nelem_order)
!
   use data_structure3D
!
   implicit none
!
   integer :: i,iel
   integer :: ibc(6,NRINDEX)
!
   integer, intent(out) :: Nelem_order(NRELIS)
!
   do iel=1,NRELIS
!
!  ...set physics
      ELEMS(iel)%nrphysics = 4
      allocate(ELEMS(iel)%physics(4))
      ELEMS(iel)%physics(1) ='ucont'
      ELEMS(iel)%physics(2) ='etang'
      ELEMS(iel)%physics(3) ='vnorm'
      ELEMS(iel)%physics(4) ='qdisc'
!
!  ...set initial mesh element order
      Nelem_order(iel) = 111
!
!  ...set boundary condition flags: 0 - no BC ; 1 - Dirichlet
      ibc(1:6,1) = 1
      ibc(1:6,2:NRINDEX) = 0
      
!  ...allocate BC flags (one per attribute component)
      allocate(ELEMS(iel)%bcond(NRINDEX))
!
!  ...for each component, encode face BC into a single BC flag
      do i=1,NRINDEX
         call encodg(ibc(1:6,i),10,6, ELEMS(iel)%bcond(i))
      enddo
!
   enddo
!
end subroutine set_initial_mesh
!
!----------------------------------------------------------------------
!> @brief dummy routine
subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
   ValH = 0; DvalH = 0
   ValE = 0; DvalE = 0
   ValV = 0; DvalV = 0
!
end subroutine dirichlet
!----------------------------------------------------------------------
!> @brief dummy routine
subroutine celem(Mdle,Idec,Nrdofs,Nrdofm,Nrdofc,Nodm, &
                 NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,Bload,Astif)
!
end subroutine celem
