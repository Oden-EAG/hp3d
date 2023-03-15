!> @date Mar 2023
program test_hexa_adapt_h
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
   real(8) :: per
   integer :: nitr,npass,nref
   
   integer :: idx,icnt,itr,kref,mdle
   real(8) :: x
!
   real(8), parameter :: perc  = 0.1d0
   integer, parameter :: nritr = 1
!
   integer :: mdle_ref(1000)
!
!..initialize MPI environment
   call mpi_w_init
!
!..initialize physics, geometry, etc.
   call initialize
!
!..
   QUIET_MODE = .false.
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
   mdle_ref = 0; icnt = 0
!
   do itr=1,nritr
!
      nref = int(NRELES*Per)
!
      ! iterate nref times
!
!  ...random number generate from 0.0 to 1.0
      if (RANK.eq.ROOT) then
         call random_seed
         call random_number(x)
      endif
      ! if (num_procs>1) broadcast from root

!  ...pick one mdle from list
      idx  = min(int(NRELES*x+1),NRELES)
      mdle = ELEM_ORDER(idx)
!
      icnt = icnt+1
      mdle_ref(icnt) = mdle
!
      write(*,*) 'call refine'
      write(*,*) 'mdle=',mdle
      call get_isoref(mdle, kref)
      call refine(mdle, kref)
!
      write(*,*) 'call close_mesh'
      call close_mesh
      write(*,*) 'returned close_mesh'
!
      if (DISTRIBUTED) then
         call par_verify
      endif
!
      call update_gdof
      call update_Ddof
!
   enddo
!
!..
!
!..percentage of elements to refine per iteration
!   per  = 0.1
!..number of refinements (iterations)
!   nitr = 5
!
!   write(*,*) 'test not working at the moment (needs revision)'
!   goto 99
!
!   call random_refine_test(per,nitr,npass)
!   if (npass .eq. 1) then
!      write(*,*) 'test_random_refine PASSED.'
!   else
!      write(*,*) 'test_random_refine FAILED.'
!   endif
!
!   99 continue
!
!..finalize MPI environment
   call mpi_w_finalize
!
end program test_hexa_adapt_h
!
!----------------------------------------------------------------------
! !> @brief Refines mesh randomly
! !> @date Mar 2023
!----------------------------------------------------------------------
   subroutine test_random_refine_aux(Per,Nitr, Npass)
!
      use data_structure3D
!
      implicit none
!
      real(8), intent(in)  :: Per
      integer, intent(in)  :: Nitr
      integer, intent(out) :: Npass
!
      integer :: kref_prism(3) = (/11,10,1/)
      integer :: mdle_list(NRELES)
!
      integer :: i,iel,icnt,idx,kref
      integer :: mdle,nelts,nref
      real(8) :: x
!
      integer :: iprint
      iprint=0
!
      Npass = 1
!
      do i=1,Nitr
!
!  .....collect all mdle
        mdle = 0
        do iel=1, NRELES
          call nelcon(mdle, mdle)
          mdle_list(iel) = mdle
        enddo
!
!  .....set number of refinements to be done
        icnt = 0
        kref = 0
        x = 0.0
        nref = int(NRELES*Per)
        nelts = NRELES
!
        if (iprint.eq.1) then
           write(*,7010) i, nelts, nref 
 7010      format('test_random_refine_aux: i=',i3,' nelts=',i6,' nref=',i6)
        endif

        do while (nref.gt.0)

!  .......random number generage from 0.0 to 1.0
          call random_seed
          call random_number(x)

!  .......pick one mdle from list
          idx  = int(nelts*x+1)
          mdle = mdle_list(idx)
          if (NODES(mdle)%ref_kind.eq.0) then
            select case (NODES(mdle)%ntype)
            case (MDLN)
              call get_isoref(mdle, kref)
              call refine(mdle, kref)
            case (MDLP)
              kref = kref_prism(mod(int(x*100),3)+1)
              call refine(mdle, kref)
            case default
              write(*,*) 'test_random_refine_aux: MDLE type'
              stop
            end select
            nref = nref - 1
            if (iprint.eq.1) then
              write(*,7000) mdle, S_Type(NODES(mdle)%ntype), kref
 7000         format('mdle =',i6,' ', a5, 'kref=', i3)
            endif
          endif

!  .......loop exit condition to prevent infinite loop
          icnt = icnt + 1
          if (icnt.gt.(nelts*1000)) then
            exit
          endif
        enddo
!
        if (.false.) Npass = 0
!
      enddo
!
   end subroutine test_random_refine_aux

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
!                        NRCOMS // MAXNRHS //
   call set_parameters(      1 ,        1 ,  &
!                       MAXEQNH // MAXEQNE // MAXEQNV // MAXEQNQ //
                             1 ,        1,         1,         1)
!
!..read geometry file
   call read_geometry('../files/mesh/hexa_orient_0')
!
!..read physics file and generate mesh data structures
   call hp3gen('../files/physics/physics_0')
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
      
   enddo
!
end subroutine set_initial_mesh
!
!----------------------------------------------------------------------
!> @brief dummy routine
subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
   zvalH = 0; zdvalH = 0
   zvalE = 0; zdvalE = 0
   zvalV = 0; zdvalV = 0
!
end subroutine dirichlet
!----------------------------------------------------------------------
!> @brief dummy routine
subroutine celem(Mdle,Idec,Nrdofs,Nrdofm,Nrdofc,Nodm, &
                 NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,Bload,Astif)
!
end subroutine celem
