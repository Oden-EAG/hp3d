!> @brief Tests datastructure dumpout and dumpin debug routines
!> @details Tests whether exporting the datastructure to file
!!          and then importing retains all DOFs correctly
!> @date Sep 2023
program test_dump_pois
!
   use control
   use data_structure3D
   use environment
   use mpi_wrapper
   use par_mesh
   use stc
!
   implicit none
!
   logical, external :: dnear
!
   integer :: NPASS
!
   QUIET_MODE = .true.
!
!..initialize MPI environment
   call mpi_w_init
!
#if HP3D_COMPLEX
   write(*,*) 'test_dump_pois: HP3D_COMPLEX'
   NPASS = 1; goto 99
#endif
!
!..initialize physics, geometry, etc.
   call initialize
!
!..this test is set up for sequential computation
   if (NUM_PROCS > 1) then
      write(*,*) 'test_dump_pois: NUM_PROCS>1'
      NPASS = 1; goto 99
   endif
!
   call dump_pois
!
 99 continue
!
!..finalize MPI environment
   call mpi_w_finalize
!
   if (NPASS.ne.1) stop 1
!
   contains
!
!----------------------------------------------------------------------
!> @brief Computes H1 error for Poisson Galerkin problem with
!!        manufactured sinusoidal solution before dumping out
!!        datastructure and after dumping in again
!> @date Sep 2023
!----------------------------------------------------------------------
   subroutine dump_pois
!
      real(8) :: errorH,errorE,errorV,errorQ
      real(8) :: rnormH,rnormE,rnormV,rnormQ
!
      real(8) :: err_curr,err_dump
!
      character(len=15) :: dump_file
!
      integer :: iflag(1),iel
!
      integer :: iprint
      iprint=0
!
      NPASS = 1
!
!  ...FLAGS
      ISTC_FLAG = .true.
      STORE_STC = .true.
      HERM_STC  = .false.
!
      iflag(1) = 1
!
!  ...refine
      call global_href
!
!  ...update geometry and Dirichlet DOFs
      call update_gdof
      call update_Ddof
!
!  ...solve
      call mumps_sc('G')
!
!  ...compute H1 error
      err_curr = 0.d0
      do iel=1,NRELES
         call element_error(ELEM_ORDER(iel),iflag,          &
                            errorH,errorE,errorV,errorQ,    &
                            rnormH,rnormE,rnormV,rnormQ)
         err_curr = err_curr + errorH
      enddo
      err_curr = sqrt(err_curr)
!
!  ...dumpout datastructure
      dump_file = 'dump_pois_file'
      call dumpout_hp3d(dump_file)
!
!  ...refine the mesh to change current DOFs
      call global_href
!
!  ...dumpin datastructure
      call dumpin_hp3d(dump_file,Delete_file=.true.)
!
!  ...compute the error again
      err_dump = 0.d0
      do iel=1,NRELES
         call element_error(ELEM_ORDER(iel),iflag,          &
                            errorH,errorE,errorV,errorQ,    &
                            rnormH,rnormE,rnormV,rnormQ)
         err_dump = err_dump + errorH
      enddo
      err_dump = sqrt(err_dump)
!
      if (.not. dnear(err_curr,err_dump)) then
         NPASS = 0
         iprint = 1
         goto 80
      endif
!
!  ...update geometry and Dirichlet DOFs
      call update_gdof
      call update_Ddof
!
!  ...solve
      call mumps_sc('G')
!
!  ...compute the error again after resolving the problem
      err_dump = 0.d0
      do iel=1,NRELES
         call element_error(ELEM_ORDER(iel),iflag,          &
                            errorH,errorE,errorV,errorQ,    &
                            rnormH,rnormE,rnormV,rnormQ)
         err_dump = err_dump + errorH
      enddo
      err_dump = sqrt(err_dump)
!
      if (.not. dnear(err_curr,err_dump)) then
         NPASS = 0
         iprint = 1
      endif
!
   80 continue
!
      if (iprint.eq.1) then
         write(*,210) err_curr,err_dump
     210 format(  'dump_pois: err_curr, err_dump = ',I10,' ,',ES12.4)
      endif
!
   end subroutine dump_pois
!
!
end program test_dump_pois
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
   call set_parameters(     2,     1)
!
!..read geometry file
   call read_geometry('../files/mesh/hexa_orient_0')
!
!..read physics file and generate mesh data structures
   call hp3gen('../files/physics/physics_H')
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
   if (NRINDEX .ne. 1) then
      write(*,*) 'dump_pois: NRINDEX = ',NRINDEX
      stop
   endif
!
   do iel=1,NRELIS
!
!  ...set physics
      ELEMS(iel)%nrphysics = 1
      allocate(ELEMS(iel)%physics(1))
      ELEMS(iel)%physics(1) ='field'
!
!  ...set initial mesh element order
      Nelem_order(iel) = 222
!
!  ...set boundary condition flags: 0 - no BC ; 1 - Dirichlet
      ibc(1:6,1) = 1
      
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
!> @brief Assembles element matrices
subroutine elem(Mdle, Itest,Itrial)
!
   use data_structure3D
   use assembly, only: ALOC,BLOC
!
   implicit none
!
   integer ,intent(in)  :: Mdle
   integer ,intent(out) :: Itest(1),Itrial(1)
!
   integer :: norder(19)
   integer :: nrdofH,nrdofE,nrdofV,nrdofQ
!
!..activate physics variable (H1) for assembly
   Itest(1) = 1; Itrial(1) = 1
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..find number of dof for each energy space supported by the element
   call celndof(NODES(Mdle)%ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!
!..call element integration routine
   call elem_poisson(Mdle,nrdofH, ALOC(1,1)%array,BLOC(1)%array)
!
end subroutine elem
!
!----------------------------------------------------------------------
!> @brief Auxiliary routine for assembling element matrices
subroutine elem_poisson(Mdle,Nrdof, Zaloc,Zbloc)
!
   use data_structure3D
   use element_data
   use parameters
   use mpi_param, only: RANK
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(in)  :: Nrdof
   real(8), intent(out) :: Zaloc(Nrdof,Nrdof), Zbloc(Nrdof)
!
   real(8) :: rjac, fval, wa, weight, q, p
   integer :: etype, iflag, nrv, nre, nrf
   integer :: nrdofH, nint, k1, k2, l
   integer :: norder(19), norient_edge(12), norient_face(6)
   real(8) :: xnod(3,MAXbrickH)
   real(8) :: xi(3), x(3), dxdxi(3,3), dxidx(3,3)
   real(8) :: shapH(MAXbrickH), gradH(3,MAXbrickH)
   real(8) :: xiloc(3,MAX_NINT3), waloc(MAX_NINT3)
   real(8) :: dq(3), dp(1:3)
!----------------------------------------------------------------------
!
   Zaloc = ZERO; Zbloc = ZERO
!
   etype = NODES(Mdle)%ntype
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
   call find_order(Mdle, norder)
   call find_orient(Mdle, norient_edge,norient_face)
   call nodcor(Mdle, xnod)
   call set_3D_int(etype,norder,norient_face, nint,xiloc,waloc)
!
!..loop over integration points
   do l=1,nint
      xi(1:3)=xiloc(1:3,l)
      wa=waloc(l)
!  ...H1 shape functions
      call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, x,dxdxi,dxidx,rjac,iflag)
!  ...integration weight
      weight = rjac*wa
!  ...get the RHS
      call getf(Mdle,x, fval)
!  ...loop through H1 test functions
      do k1=1,nrdofH
!     ...Piola transformation
         q = shapH(k1)
         dq(1:3) = gradH(1,k1)*dxidx(1,1:3) &
                 + gradH(2,k1)*dxidx(2,1:3) &
                 + gradH(3,k1)*dxidx(3,1:3)
!     ...accumulate for the load vector
         Zbloc(k1) = Zbloc(k1) + q*fval*weight
!     ...loop through H1 trial functions
         do k2=1,nrdofH
!        ...Piola transformation
            p = shapH(k2)
            dp(1:3) = gradH(1,k2)*dxidx(1,1:3) &
                    + gradH(2,k2)*dxidx(2,1:3) &
                    + gradH(3,k2)*dxidx(3,1:3)
!        ...accumulate for the stiffness matrix
            Zaloc(k1,k2) = Zaloc(k1,k2) + weight * ( &
                           dq(1)*dp(1) + dq(2)*dp(2) + dq(3)*dp(3))
!     ...end of loop through H1 trial functions
         enddo
!  ...end of loop through H1 test functions
      enddo
!..end of loop through integration points
   enddo
!
end subroutine elem_poisson

!
!----------------------------------------------------------------------
!> @brief Returns Dirichlet data
subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!
   use parameters, only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ
!
   implicit none
!
   real(8), intent(in)  :: X(3)
   integer, intent(in)  :: Icase,Mdle
!
   real(8),dimension(  MAXEQNH    ) ::   ValH
   real(8),dimension(  MAXEQNH,3  ) ::  DvalH
   real(8),dimension(  MAXEQNH,3,3) :: d2valH
   real(8),dimension(3,MAXEQNE    ) ::   ValE
   real(8),dimension(3,MAXEQNE,3  ) ::  DvalE
   real(8),dimension(3,MAXEQNE,3,3) :: d2valE
   real(8),dimension(3,MAXEQNV    ) ::   ValV
   real(8),dimension(3,MAXEQNV,3  ) ::  DvalV
   real(8),dimension(3,MAXEQNV,3,3) :: d2valV
   real(8),dimension(  MAXEQNQ    ) ::   valQ
   real(8),dimension(  MAXEQNQ,3  ) ::  dvalQ
   real(8),dimension(  MAXEQNQ,3,3) :: d2valQ
!
   call exact(X,Icase, ValH,DvalH,d2valH,ValE,DvalE,d2valE, &
                       ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
!
end subroutine dirichlet
!
!----------------------------------------------------------------------
!> @brief Returns exact solution data
subroutine exact(X,Icase, ValH,DvalH,D2valH, &
                          ValE,DvalE,D2valE, &
                          ValV,DvalV,D2valV, &
                          ValQ,DvalQ,D2valQ)
!
   use control   , only : NEXACT
   use parameters, only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ
!
   implicit none
!
   real(8), intent(in)  :: X(3)
   integer, intent(in)  :: Icase
!
   real(8),dimension(  MAXEQNH    ) ::   ValH
   real(8),dimension(  MAXEQNH,3  ) ::  DvalH
   real(8),dimension(  MAXEQNH,3,3) :: D2valH
   real(8),dimension(3,MAXEQNE    ) ::   ValE
   real(8),dimension(3,MAXEQNE,3  ) ::  DvalE
   real(8),dimension(3,MAXEQNE,3,3) :: D2valE
   real(8),dimension(3,MAXEQNV    ) ::   ValV
   real(8),dimension(3,MAXEQNV,3  ) ::  DvalV
   real(8),dimension(3,MAXEQNV,3,3) :: D2valV
   real(8),dimension(  MAXEQNQ    ) ::   ValQ
   real(8),dimension(  MAXEQNQ,3  ) ::  dvalQ
   real(8),dimension(  MAXEQNQ,3,3) :: D2valQ
!
   real(8), parameter   :: PI = 4.d0*datan(1.d0)
!
   real(8) :: x1,x2,x3
   x1 = X(1); x2 = X(2); x3 = X(3)
!
!..make sure exact solution is available
   if (NEXACT.ne.1) then
      write(*,*) 'dump_pois: NEXACT =',NEXACT
      stop
   endif
!
!..smooth sin sinh z solution on unit cube
!  non-vanishing boundary only on x=1,y=1,z=1
   ValH(1)    = dsin(x1*PI)*dsinh(x2*PI)*(x3**2)
!
   DValH(1,1) = PI*dcos(x1*PI)*dsinh(x2*PI)*(x3**2)
   DValH(1,2) = PI*dsin(x1*PI)*dcosh(x2*PI)*(x3**2)
   DValH(1,3) = 2.*dsin(x1*PI)*dsinh(x2*PI)*x3
!
   D2valH(1,1,1) = -PI*PI*dsin(x1*PI)*dsinh(x2*PI)*(x3**2)
   D2valH(1,1,2) =  PI*PI*dcos(x1*PI)*dcosh(x2*PI)*(x3**2)
   D2valH(1,1,3) =  2.*PI*dcos(x1*PI)*dsinh(x2*PI)*x3
   D2valH(1,2,1) =  PI*PI*dcos(x1*PI)*dcosh(x2*PI)*(x3**2)
   D2valH(1,2,2) =  PI*PI*dsin(x1*PI)*dsinh(x2*PI)*(x3**2)
   D2valH(1,2,3) =  2.*PI*dsin(x1*PI)*dcosh(x2*PI)*x3
   D2valH(1,3,1) =  2.*PI*dcos(x1*PI)*dsinh(x2*PI)*x3
   D2valH(1,3,2) =  2.*PI*dsin(x1*PI)*dcosh(x2*PI)*x3
   D2valH(1,3,3) =  2.*   dsin(x1*PI)*dsinh(x2*PI)
!
end subroutine exact
!
!----------------------------------------------------------------------
!> @brief Returns load data for known solution
subroutine getf(Mdle,X, Fval)
!
   use parameters, only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ
!
   implicit none
!
   integer, intent(in)  :: Mdle
   real(8), intent(in)  :: X(3)
   real(8), intent(out) :: Fval
!
   real(8),dimension(  MAXEQNH    ) ::   valH
   real(8),dimension(  MAXEQNH,3  ) ::  dvalH
   real(8),dimension(  MAXEQNH,3,3) :: d2valH
   real(8),dimension(3,MAXEQNE    ) ::   valE
   real(8),dimension(3,MAXEQNE,3  ) ::  dvalE
   real(8),dimension(3,MAXEQNE,3,3) :: d2valE
   real(8),dimension(3,MAXEQNV    ) ::   valV
   real(8),dimension(3,MAXEQNV,3  ) ::  dvalV
   real(8),dimension(3,MAXEQNV,3,3) :: d2valV
   real(8),dimension(  MAXEQNQ    ) ::   valQ
   real(8),dimension(  MAXEQNQ,3  ) ::  dvalQ
   real(8),dimension(  MAXEQNQ,3,3) :: d2valQ
!
   integer :: icase
   icase = 1
!
   call exact(X,icase, ValH,DvalH,d2valH,ValE,DvalE,d2valE, &
                       ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
!
   Fval = -(d2valH(1,1,1)+d2valH(1,2,2)+d2valH(1,3,3))
!
end subroutine getf
!
!----------------------------------------------------------------------
!> @brief dummy routine
subroutine celem(Mdle,Idec,Nrdofs,Nrdofm,Nrdofc,Nodm, &
                 NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,Bload,Astif)
!
   write(*,*) 'dump_pois: celem'
   stop
!
end subroutine celem
