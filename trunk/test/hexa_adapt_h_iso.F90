!> @brief Tests h-adaptive refinements for hexa mesh
!> @details Randomly isotropically h-refines a distributed
!!          hexa mesh repeatedly and verifies that mesh
!!          consistency is maintained after each iteration.
!> @date Mar 2023
program test_hexa_adapt_h_iso
!
   use data_structure3D
   use environment
   use mpi
   use mpi_param
   use mpi_wrapper
   use par_mesh
!
   implicit none
!
   integer :: NPASS
!
!..percentage of elements to refine per iteration
   real(8), parameter :: PERC = 0.2d0
!..number of iterations
   integer, parameter :: NITR = 5
!..max number of randomly refined elements
   integer, parameter :: NREF_MAX = 1000
!
   QUIET_MODE = .true.
!
!..initialize MPI environment
   call mpi_w_init
!
!..initialize physics, geometry, etc.
   call initialize
!
!..perform one uniform h-refinement
   call global_href
   call update_gdof
   call update_Ddof
!
   if (DISTRIBUTED) then
      EXCHANGE_DOF = .false.
      call par_verify
      call zoltan_w_set_lb(0)
      call distr_mesh
   endif
!
   call hexa_adapt_h_iso
!
!..finalize MPI environment
   call mpi_w_finalize
!
   if (NPASS.ne.1) stop 1
!
   contains
!
!----------------------------------------------------------------------
! !> @brief Refines mesh randomly
! !> @date Mar 2023
!----------------------------------------------------------------------
   subroutine hexa_adapt_h_iso
!
      integer :: idx,icnt,itr,ierr,kref,mdle,nel,nref
      real(8) :: x
!
      integer :: mdle_ref(NREF_MAX)
!
      integer :: iprint
      iprint=0
!
      NPASS = 1
!
      mdle_ref = 0; icnt = 0
!
      if (iprint.ge.1 .and. RANK.eq.ROOT) then
         write(*,110) ' NITR,PERC = ',NITR,PERC
     110 format(' hexa_adapt_h_iso: ',A,I3,',',F5.2)
      endif
!
      do itr=1,NITR
!
         nref = max(1,int(NRELES*PERC))
         if (iprint.ge.1 .and. RANK.eq.ROOT) then
            write(*,120) '  itr,nref = ',itr,nref
        120 format(' hexa_adapt_h_iso: ',A,I3,',',I5)
         endif
!
         ! iterate nref times
         nel = NRELES
         do while (nref.gt.0)
!
!        ...random number generate from 0.0 to 1.0
            if (RANK.eq.ROOT) then
               mdle = 1
               do while (.not.Is_leaf(mdle))
                  call random_seed
                  call random_number(x)
!              ...pick one mdle from list
                  idx  = min(int(nel*x+1),nel)
                  mdle = ELEM_ORDER(idx)
               enddo
            endif
!
            if (NUM_PROCS.gt.1) then
               call MPI_BCAST(mdle,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD, ierr)
            endif
!
!        ...keep track of refinements
            icnt = icnt+1
            if (icnt.le.NREF_MAX) then
               mdle_ref(icnt) = mdle
            else
               if (RANK.eq.ROOT) then
                  write(*,130) '  itr,NITR,NREF_MAX = ',itr,NITR,NREF_MAX
              130 format(' hexa_adapt_h_iso: ',A,I5,',',I5,',',I5)
               endif
               exit
            endif
!
            if (iprint.gt.1 .and. RANK.eq.ROOT) then
               write(*,140) '       mdle = ',mdle
           140 format(' hexa_adapt_h_iso: ',A,I6)
            endif
!
            call get_isoref(mdle, kref)
            call refine(mdle, kref)
!
            nref = nref-1
!
!     ...end loop over nref
         enddo
!
         call close_mesh
!
         if (DISTRIBUTED) then
!        ...verify mesh consistency
!           TODO: change par_verify to return flag whether it passed
!                 (instead of ending computation with error msg)
            call par_verify
            if (NPASS.ne.1) then
               ! TEST FAILS
               if (RANK.eq.ROOT) then
                  write(*,150) '      icnt = ',icnt
                  do idx=1,icnt
                     write(*,150) '  mdle_ref = ',mdle_ref(idx)
                  enddo
              150 format(' hexa_adapt_h_iso: ',A,I6)
               endif
               return
            endif
!        ...repartition
            if (itr .le. 2) call distr_mesh
         endif
!
         if (icnt.gt.NREF_MAX) exit
!
!  ...end loop over itr
      enddo
!
   end subroutine hexa_adapt_h_iso
!
!
end program test_hexa_adapt_h_iso
!
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
   use zoltan_wrapper
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
   call hp3gen('../files/physics/physics_0')
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
