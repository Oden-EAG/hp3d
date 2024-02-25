!> @brief Tests polynomial solution for Maxwell Galerkin problem
!> @details Tests exact representation of a polynomial solution
!!          for the variational curl-curl Maxwell problem
!> @date Mar 2023
program test_poly_maxw
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
   integer :: NPASS
!
!..number of uniform h-refinements
   integer, parameter :: nref = 2
!
   QUIET_MODE = .true.
!
!..initialize MPI environment
   call mpi_w_init
!
#if C_MODE==0
   write(*,*) 'test_poly_maxw: C_MODE=0'
   NPASS = 1; goto 99
#endif
!
!..initialize physics, geometry, etc.
   call initialize
!
!..this test is set up for sequential computation
   if (NUM_PROCS > 1) then
      write(*,*) 'test_poly_maxw: NUM_PROCS>1'
      NPASS = 1; goto 99
   endif
!
   call poly_maxw(nref)
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
! !> @brief Computes H(curl) error for Maxwell Galerkin problem with
!!          manufactured polynomial solution
! !> @date Mar 2023
!----------------------------------------------------------------------
   subroutine poly_maxw(Nref)
!
      use assembly_sc, only: NRDOF_TOT
!
      integer, intent(in) :: Nref
!
      real(8) :: errorH,errorE,errorV,errorQ
      real(8) :: rnormH,rnormE,rnormV,rnormQ
!
      real(8) :: err_curr
!
      integer :: iflag(1),iel,i
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
!  ...do refinements and solve
      do i=1,Nref
!     ...refine
         call global_href
!     ...update geometry and Dirichlet DOFs
         call update_gdof
         call update_Ddof
!     ...solve
         call mumps_sc('G')
!
!     ...compute H(curl) error
         err_curr = 0.d0
         do iel=1,NRELES
            call element_error(ELEM_ORDER(iel),iflag,          &
                               errorH,errorE,errorV,errorQ,    &
                               rnormH,rnormE,rnormV,rnormQ)
            err_curr = err_curr + errorE
         enddo
         err_curr = sqrt(err_curr)
!
         if (err_curr > 1.0d-14) then
            NPASS = 0
            iprint = 1
         endif
!
         if (iprint.eq.1) then
            write(*,210) NRDOF_TOT,err_curr
        210 format(  'poly_maxw: dof_curr, err_curr = ',I10,' ,',ES12.4)
         endif
!
      enddo
!
   end subroutine poly_maxw
!
!
end program test_poly_maxw
!
#include "typedefs.h"
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
   call hp3gen('../files/physics/physics_E')
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
      write(*,*) 'poly_maxw: NRINDEX = ',NRINDEX
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
   call elem_maxwell(Mdle,nrdofE, ALOC(1,1)%array,BLOC(1)%array)
!
end subroutine elem
!
!----------------------------------------------------------------------
!> @brief Auxiliary routine for assembling element matrices
subroutine elem_maxwell(Mdle,Nrdof, Zaloc,Zbloc)
!
   use data_structure3D
   use element_data
   use parameters
   use mpi_param, only: RANK
!
   implicit none
!
   integer   , intent(in)  :: Mdle
   integer   , intent(in)  :: Nrdof
   complex(8), intent(out) :: Zaloc(Nrdof,Nrdof), Zbloc(Nrdof)
!
   real(8), parameter :: PI = 4.d0*datan(1.d0)
   real(8), parameter :: OMEGA   = 1.d0 * PI
   real(8), parameter :: EPSILON = 1.d0
   real(8), parameter :: MU      = 1.d0
   real(8), parameter :: SIGMA   = 0.d0
!
!-------------------------------------------------------------------------
!
   complex(8) :: zJ(3), za, zb
!
   real(8) :: rjac, wa, weight
   integer :: iflag, etype, nrdofH, nrdofE, nint, k1, k2, l
!
   integer :: norder(19), norient_edge(12), norient_face(6)
   real(8) :: xnod(3,MAXbrickH)
   real(8) :: xi(3), x(3), dxdxi(3,3), dxidx(3,3)
   real(8) :: shapH(  MAXbrickH), gradH(3,MAXbrickH)
   real(8) :: shapE(3,MAXbrickE), curlE(3,MAXbrickE)
   real(8) :: xiloc(3,MAX_NINT3), waloc(MAX_NINT3)
   real(8) :: E(3), CE(3), F(3), CF(3)
!
!----------------------------------------------------------------------
!
   Zaloc = ZERO; Zbloc = ZERO
!
   etype = NODES(Mdle)%ntype
   call find_order(Mdle, norder)
   call find_orient(Mdle, norient_edge,norient_face)
   call nodcor(Mdle, xnod)
!
   call set_3D_int(etype,norder,norient_face, nint,xiloc,waloc)
!
!..loop over integration points
   do l=1,nint
!  ...coordinates and weight of this integration point
      xi(1:3)=xiloc(1:3,l)
      wa=waloc(l)
!  ...H1 shape functions (for geometry)
      call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!  ...H(curl) shape functions
      call shape3DE(etype,xi,norder,norient_edge,norient_face, nrdofE,shapE,curlE)
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, x,dxdxi,dxidx,rjac,iflag)
!  ...integration weight
      weight = rjac*wa
!  ...get the RHS
      call getf(Mdle,x, zJ)
!
!  ...loop through H(curl) test functions
      do k1=1,nrdofE
!
!     ...Piola transformation
         F(1:3) = shapE(1,k1)*dxidx(1,1:3) &
                + shapE(2,k1)*dxidx(2,1:3) &
                + shapE(3,k1)*dxidx(3,1:3)
         CF(1:3) = dxdxi(1:3,1)*curlE(1,k1) &
                 + dxdxi(1:3,2)*curlE(2,k1) &
                 + dxdxi(1:3,3)*curlE(3,k1)
         CF(1:3) = CF(1:3)/rjac
!
!     ...accumulate for the load vector: (-iωJ,F)
         za = F(1)*zJ(1) + F(2)*zJ(2) + F(3)*zJ(3)
         Zbloc(k1) = Zbloc(k1) - ZIMG*OMEGA*za*weight
!
!     ...loop through H(curl) trial functions
         do k2=1,nrdofE
!
!        ...Piola transformation
            E(1:3) = shapE(1,k2)*dxidx(1,1:3) &
                   + shapE(2,k2)*dxidx(2,1:3) &
                   + shapE(3,k2)*dxidx(3,1:3)
            CE(1:3) = dxdxi(1:3,1)*curlE(1,k2) &
                    + dxdxi(1:3,2)*curlE(2,k2) &
                    + dxdxi(1:3,3)*curlE(3,k2)
            CE(1:3) = CE(1:3)/rjac
!
!        ...accumulate for the stiffness matrix: ((1/μ)curl E, curl F)-((ω^2ε-iωσ)E, F)
            za = (CE(1)*CF(1) + CE(2)*CF(2) + CE(3)*CF(3)) / MU
            zb = (OMEGA*OMEGA*EPSILON - ZIMG*OMEGA*SIGMA)*(E(1)*F(1) + E(2)*F(2) + E(3)*F(3))
            Zaloc(k1,k2) = Zaloc(k1,k2) + (za-zb)*weight
!
!     ...end of loop through H(curl) trial functions
         enddo
!  ...end of loop through H(curl) test functions
      enddo
!..end of loop through integration points
   enddo
!
end subroutine elem_maxwell
!
!----------------------------------------------------------------------
!> @brief Returns Dirichlet data
subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!
   use parameters, only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
!
   implicit none
!
   real(8), intent(in)  :: X(3)
   integer, intent(in)  :: Icase,Mdle
!
   VTYPE,dimension(  MAXEQNH    ) ::   ValH
   VTYPE,dimension(  MAXEQNH,3  ) ::  DvalH
   VTYPE,dimension(  MAXEQNH,3,3) :: d2valH
   VTYPE,dimension(3,MAXEQNE    ) ::   ValE
   VTYPE,dimension(3,MAXEQNE,3  ) ::  DvalE
   VTYPE,dimension(3,MAXEQNE,3,3) :: d2valE
   VTYPE,dimension(3,MAXEQNV    ) ::   ValV
   VTYPE,dimension(3,MAXEQNV,3  ) ::  DvalV
   VTYPE,dimension(3,MAXEQNV,3,3) :: d2valV
   VTYPE,dimension(  MAXEQNQ    ) ::   valQ
   VTYPE,dimension(  MAXEQNQ,3  ) ::  dvalQ
   VTYPE,dimension(  MAXEQNQ,3,3) :: d2valQ
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
   use parameters, only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
!
   implicit none
!
   real(8), intent(in)  :: X(3)
   integer, intent(in)  :: Icase
!
   VTYPE,dimension(  MAXEQNH    ) ::   ValH
   VTYPE,dimension(  MAXEQNH,3  ) ::  DvalH
   VTYPE,dimension(  MAXEQNH,3,3) :: D2valH
   VTYPE,dimension(3,MAXEQNE    ) ::   ValE
   VTYPE,dimension(3,MAXEQNE,3  ) ::  DvalE
   VTYPE,dimension(3,MAXEQNE,3,3) :: D2valE
   VTYPE,dimension(3,MAXEQNV    ) ::   ValV
   VTYPE,dimension(3,MAXEQNV,3  ) ::  DvalV
   VTYPE,dimension(3,MAXEQNV,3,3) :: D2valV
   VTYPE,dimension(  MAXEQNQ    ) ::   ValQ
   VTYPE,dimension(  MAXEQNQ,3  ) ::  dvalQ
   VTYPE,dimension(  MAXEQNQ,3,3) :: D2valQ
!
   real(8) :: x1,x2,x3
   x1 = X(1); x2 = X(2); x3 = X(3)
!
!..make sure exact solution is available
   if (NEXACT.ne.1) then
      write(*,*) 'poly_maxw: NEXACT =',NEXACT
      stop
   endif
!
!..Field values (polynomial solution with homogeneous tangential BC field)
   ValE(1,1) =  X(2)*(1.d0-X(2))*X(3)*(1.d0-X(3))
   ValE(2,1) =  X(2)*X(1)*(1.d0-X(1))*X(3)*(1.d0-X(3))
   ValE(3,1) =  X(1)*(1.d0-X(1))*X(2)*(1.d0-X(2))
!
!..1st order derivatives
!  grad(E_x)
   DvalE(1,1,1) = ZERO !...derivative wrt x
   DvalE(1,1,2) = (1.d0-2.d0*X(2))*X(3)*(1.d0-X(3)) !...derivative wrt y
   DvalE(1,1,3) = X(2)*(1.d0-X(2))*(1.d0-2.d0*X(3)) !...derivative wrt z
!  grad(E_y)
   DvalE(2,1,1) = X(2)*(1.d0-2.d0*X(1))*X(3)*(1.d0-X(3))
   DvalE(2,1,2) = X(1)*(1.d0-X(1))*X(3)*(1.d0-X(3))
   DvalE(2,1,3) = X(2)*X(1)*(1.d0-X(1))*(1.d0-2.d0*X(3))
!  grad(E_z)
   DvalE(3,1,1) = (1.d0-2.d0*X(1))*X(2)*(1.d0-X(2))
   DvalE(3,1,2) = X(1)*(1.d0-X(1))*(1.d0-2.d0*X(2))
   DvalE(3,1,3) = ZERO
!
!..2nd order derivatives
!  grad grad E_x
   D2valE(1,1,1,1) = ZERO !...derivative wrt xx
   D2valE(1,1,1,2) = ZERO !...derivative wrt xy
   D2valE(1,1,1,3) = ZERO !...derivative wrt xz
   D2valE(1,1,2,1) = D2valE(1,1,1,2) !...derivative wrt yx
   D2valE(1,1,2,2) = -2.d0*X(3)*(1.d0-X(3))  !...derivative wrt yy
   D2valE(1,1,2,3) = (1.d0-2.d0*X(2))*(1.d0-2.d0*X(3)) !...derivative wrt yz
   D2valE(1,1,3,1) = D2valE(1,1,1,3) !...derivative wrt zx
   D2valE(1,1,3,2) = D2valE(1,1,2,3) !...derivative wrt zy
   D2valE(1,1,3,3) = X(2)*(1.d0-X(2))*(-2.d0) !...derivative wrt zz
!  grad grad E_y
   D2valE(2,1,1,1) = X(2)*(-2.d0)*X(3)*(1.d0-X(3))
   D2valE(2,1,1,2) = (1.d0-2.d0*X(1))*X(3)*(1.d0-X(3))
   D2valE(2,1,1,3) = X(2)*(1.d0-2.d0*X(1))*(1.d0-2.d0*X(3))
   D2valE(2,1,2,1) = D2valE(2,1,1,2)
   D2valE(2,1,2,2) = ZERO
   D2valE(2,1,2,3) = X(1)*(1.d0-X(1))*(1.d0-2.d0*X(3))
   D2valE(2,1,3,1) = D2valE(2,1,1,3)
   D2valE(2,1,3,2) = D2valE(2,1,2,3)
   D2valE(2,1,3,3) = X(2)*X(1)*(1.d0-X(1))*(-2.d0)
!  grad grad E_z
   D2valE(3,1,1,1) = (-2.d0)*X(2)*(1.d0-X(2))
   D2valE(3,1,1,2) = (1.d0-2.d0*X(1))*(1.d0-2.d0*X(2))
   D2valE(3,1,1,3) = ZERO
   D2valE(3,1,2,1) = D2valE(3,1,1,2)
   D2valE(3,1,2,2) = X(1)*(1.d0-X(1))*(-2.d0)
   D2valE(3,1,2,3) = ZERO
   D2valE(3,1,3,1) = D2valE(3,1,1,3)
   D2valE(3,1,3,2) = D2valE(3,1,2,3)
   D2valE(3,1,3,3) = ZERO
!
end subroutine exact
!
!----------------------------------------------------------------------
!> @brief Returns load data for known solution
subroutine getf(Mdle,X, zJ)
!
   use parameters, only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZIMG
!
   implicit none
!
   integer, intent(in)  :: Mdle
   real(8), intent(in)  :: X(3)
   VTYPE  , intent(out) :: zJ(3)
!
   VTYPE,dimension(  MAXEQNH    ) ::   valH
   VTYPE,dimension(  MAXEQNH,3  ) ::  dvalH
   VTYPE,dimension(  MAXEQNH,3,3) :: d2valH
   VTYPE,dimension(3,MAXEQNE    ) ::   valE
   VTYPE,dimension(3,MAXEQNE,3  ) ::  dvalE
   VTYPE,dimension(3,MAXEQNE,3,3) :: d2valE
   VTYPE,dimension(3,MAXEQNV    ) ::   valV
   VTYPE,dimension(3,MAXEQNV,3  ) ::  dvalV
   VTYPE,dimension(3,MAXEQNV,3,3) :: d2valV
   VTYPE,dimension(  MAXEQNQ    ) ::   valQ
   VTYPE,dimension(  MAXEQNQ,3  ) ::  dvalQ
   VTYPE,dimension(  MAXEQNQ,3,3) :: d2valQ
!
   complex(8) :: za(3),zb(3)
   complex(8) :: CE_x(3), CE_y(3), CE_z(3), CC_E(3)
!
   real(8), parameter :: PI = 4.d0*datan(1.d0)
   real(8), parameter :: OMEGA   = 1.d0 * PI
   real(8), parameter :: EPSILON = 1.d0
   real(8), parameter :: MU      = 1.d0
   real(8), parameter :: SIGMA   = 0.d0
!
   integer :: icase
   icase = 1
!
   call exact(X,icase, ValH,DvalH,d2valH,ValE,DvalE,d2valE, &
                       ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
!
!..1st order derivatives of curl(E)
!  grad(curl(E)_x)
   CE_x(1) = D2valE(3,1,2,1) - D2valE(2,1,3,1) ! (curl E)_1,x
   CE_x(2) = D2valE(3,1,2,2) - D2valE(2,1,3,2) ! (curl E)_1,y
   CE_x(3) = D2valE(3,1,2,3) - D2valE(2,1,3,3) ! (curl E)_1,z
!  grad(curl(E)_y)
   CE_y(1) = D2valE(1,1,3,1) - D2valE(3,1,1,1) ! (curl E)_2,x
   CE_y(2) = D2valE(1,1,3,2) - D2valE(3,1,1,2) ! (curl E)_2,y
   CE_y(3) = D2valE(1,1,3,3) - D2valE(3,1,1,3) ! (curl E)_2,z
!  grad(curl(E)_z)
   CE_z(1) = D2valE(2,1,1,1) - D2valE(1,1,2,1) ! (curl E)_3,x
   CE_z(2) = D2valE(2,1,1,2) - D2valE(1,1,2,2) ! (curl E)_3,y
   CE_z(3) = D2valE(2,1,1,3) - D2valE(1,1,2,3) ! (curl E)_3,z
!
!..curl curl E
   CC_E(1) = CE_z(2) - CE_y(3)
   CC_E(2) = CE_x(3) - CE_z(1)
   CC_E(3) = CE_y(1) - CE_x(2)
!
!..compute impressed current from exact solution
!  -iωJ = curl((1/μ)curl E) - (ω^2ε-iωσ)E
   za = CC_E / MU
   zb = (OMEGA*OMEGA*EPSILON - ZIMG*OMEGA*SIGMA) * ValE(1:3,1)
   zJ = (za-zb) / (-ZIMG * OMEGA)
!
end subroutine getf
!
!----------------------------------------------------------------------
!> @brief dummy routine
subroutine celem(Mdle,Idec,Nrdofs,Nrdofm,Nrdofc,Nodm, &
                 NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,Bload,Astif)
!
   write(*,*) 'poly_maxw: celem'
   stop
!
end subroutine celem
