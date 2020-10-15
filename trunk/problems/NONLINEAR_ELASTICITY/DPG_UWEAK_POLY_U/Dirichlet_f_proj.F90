



subroutine Dirichlet_f_proj(mdlf,nrv_f,xverts_f,nord,nrdofU,nrdofF,zdofU,zdofF)

  use control,       only: NEXACT, INTEGRATION
  use parameters,    only: MAXTRIAQ,MAXTRIAH
  use parametersdpg, only: MAXNINT2ADD, NORD_ADD

  implicit none

  integer                   ,      intent(in ) :: nrv_f,nord,mdlf
  real*8, dimension(3,nrv_f),      intent(in ) :: xverts_f
  integer,                         intent(out) :: nrdofU,nrdofF
  real*8, dimension(3,MAXTRIAH),   intent(out) :: zdofU
  real*8, dimension(3,MAXTRIAQ),   intent(out) :: zdofF

  real*8                       :: rjac_f,b_area,r_f,wa,weight
  real*8, dimension(3)         :: b0,b1,b2,x_aux,fn,x
  real*8, dimension(3,3)       :: b_aff
  real*8, dimension(2,2)       :: dxdxi_bf
  real*8, dimension(2)         :: t
  real*8, dimension(2,(nrv_f-2)*MAXNINT2ADD) :: tloc
  real*8, dimension((nrv_f-2)*MAXNINT2ADD)   :: wt
  real*8, dimension(3,(nrv_f-2)*MAXNINT2ADD) :: xloc

  ! solution variables
  real*8, dimension(3)         :: u,sigma_n
  real*8, dimension(3,3)       :: gradu
  real*8, dimension(3,3)       :: epsilon
  real*8, dimension(3,3)       :: sigma
  real*8, dimension(3)         :: divsigma
  real*8, dimension(MAXTRIAH)  :: shapH_f
  real*8, dimension(3*MAXTRIAH):: rhsU
  real*8, dimension(2,MAXTRIAH):: gradH_f
  real*8, dimension(3*MAXTRIAH*(3*MAXTRIAH+1)/2) :: ProjU
  real*8, dimension(MAXTRIAQ)  :: shapQ_f
  real*8, dimension(3*MAXTRIAQ):: rhsF
  real*8, dimension(3*MAXTRIAQ*(3*MAXTRIAQ+1)/2) :: ProjF

  integer, dimension(5) :: nordf
  integer, dimension(4) :: norie

  integer :: k1,k2,ipt,nint,i,jv,info,info1,icomp,jcomp,m1,m2,k,nrdofQ_f,nrdofH_f , imat
! statement function for packed format index for symmetric matrices
  integer :: nk
  nk(k1,k2) = (k2-1)*k2/2+k1

! get normal vector
  call face_normal(xverts_f(:,1:3),fn)
! get face affine coordinates points b0,b1,b2
  call poly_affine_2d(xverts_f,nrv_f,b0,b1,b2)
  b_aff(:,1) = b0
  b_aff(:,2) = b1
  b_aff(:,3) = b2
! get area of triangle b0b1b2
  call trian_centr_area(b0,b1,b2,x_aux,b_area)
  rjac_f = 2.d0*b_area
! 
! get quadrature points for polygon using the affine coordinates
  INTEGRATION = max(0, NORD_ADD - 1)
  call set_quadrature_affine_face_DPG(mdlf,nord,nrv_f,xverts_f,b_aff,nint,tloc,wt)
  INTEGRATION = 0
! get physical coordinates for integration points so that exact solution can be evaluated 
  xloc = 0.d0
  call poly_affine_2d_to_x(tloc(1:2,1:nint),nint,b_aff,xloc)
! 
! set up face order in usual structure
  nordf = 0
  nordf(1:4) = nord
! set orientation vector to zero
  norie = 0
! 
! ***************************************************************
!                  PROJECTION OF TRACE \hat u
! ***************************************************************
!  initialize projection matrix and rhs
  ProjU = 0.d0
  rhsU = 0.d0
!    .....loop through face integration points
  do ipt=1,nint
!
!    .......face coordinates
    t(1:2) = tloc(1:2,ipt); wa = wt(ipt); x(1:3) = xloc(1:3,ipt)
!
!   evaluate shape functions of trace space (polynomials of order p)
    call shape2DH('tria',t,nordf, norie, nrdofH_f,shapH_f,gradH_f)
! 
!   compute exact soution
    call find_material(x,fn,imat)
    call elast_solution(imat,X, u,gradu,epsilon,sigma,divsigma)

    weight = wa
!
!...OUTER loop through Face H1 shape functions
    do k1=1,nrdofH_f
      do jcomp=1,3
        m1 = (k1-1)*3+jcomp
!
!    ...INNER loop through trial dofs
        do k2=k1,nrdofH_f
          do icomp=1,3
            m2 = (k2-1)*3+icomp
            if (icomp.eq.jcomp) then
              k = nk(m1,m2)
              ProjU(k) = ProjU(k) + shapH_f(k2)*shapH_f(k1)*weight
              ! write(*,*) 'Dirichlet_f_proj :   m1,m2=',m1,m2
              ! write(*,*) 'Dirichlet_f_proj :   k1,k2=',k1,k2
              ! write(*,*) 'Dirichlet_f_proj :   jcomp,icomp=',jcomp,icomp
              ! write(*,*) 'Dirichlet_f_proj :   shapQ_f(k2),shapQ_f(k1)=',shapQ_f(k2),shapQ_f(k1)
            endif
          enddo
        enddo

        rhsU(m1) = rhsU(m1) + u(jcomp)*shapH_f(k1)*weight

      enddo
    enddo
! integration point loop ends
  enddo

! apply Cholesky Factorization to Proj
  call DPPTRF('U',3*nrdofH_f,ProjU,info)
  if (info.ne.0) then
    write(*,*) 'Dirichlet_f_proj - ProjU Gram Cholesky Factorization: info = ',info
    stop
  endif
! 
! ...Proj^-1 * rhs
  call DPPTRS('U',3*nrdofH_f,1,ProjU,rhsU,3*MAXTRIAH,info1)
  if (info1.ne.0) then
    write(*,*) 'Dirichlet_f_proj - ProjU System Solution after Cholesky: info1 = ',info1 ; stop
  endif


! ***************************************************************
!             PROJECTION OF NORMAL FLUX \hat Sigma_n
! ***************************************************************
!  initialize projection matrix and rhs
  ProjF = 0.d0
  rhsF = 0.d0
!    .....loop through face integration points
  do ipt=1,nint
!
!    .......face coordinates
    t(1:2) = tloc(1:2,ipt); wa = wt(ipt); x(1:3) = xloc(1:3,ipt)
!
!   evaluate shape functions of normal flux space (polynomials of order p-1)
    call shape2DQ('tria',t,nordf, nrdofQ_f,ShapQ_f)
! 
!   compute exact soution
    call find_material(x,fn,imat)
    call elast_solution(imat,X, u,gradu,epsilon,sigma,divsigma)
! 
    sigma_n = matmul(sigma,fn)
! 
! ..Change coordinates so the shape functions are on the physical element
!   L2 (trial)
    shapQ_f(1:nrdofQ_f) = shapQ_f(1:nrdofQ_f)/rjac_f
!
    weight = wa
!
!...OUTER loop through Face L2 shape functions
    do k1=1,nrdofQ_f
      do jcomp=1,3
        m1 = (k1-1)*3+jcomp
!
!    ...INNER loop through trial dofs
        do k2=k1,nrdofQ_f
          do icomp=1,3
            m2 = (k2-1)*3+icomp
            if (icomp.eq.jcomp) then
              k = nk(m1,m2)
              ProjF(k) = ProjF(k) + shapQ_f(k2)*shapQ_f(k1)*weight
          ! write(*,*) 'Dirichlet_f_proj :   m1,m2=',m1,m2
          ! write(*,*) 'Dirichlet_f_proj :   k1,k2=',k1,k2
          ! write(*,*) 'Dirichlet_f_proj :   jcomp,icomp=',jcomp,icomp
          ! write(*,*) 'Dirichlet_f_proj :   shapQ_f(k2),shapQ_f(k1)=',shapQ_f(k2),shapQ_f(k1)        
            endif
          enddo        
        enddo
! 
        rhsF(m1) = rhsF (m1) + sigma_n(jcomp)*shapQ_f(k1)*weight
      enddo
! 
    enddo
! 
! integration point loop ends
  enddo
! 
!! apply Cholesky Factorization to Proj
  call DPPTRF('U',3*nrdofQ_f,ProjF,info)
  if (info.ne.0) then
    write(*,*) 'Dirichlet_f_proj - ProjF Gram Cholesky Factorization: info = ',info
    stop
  endif
! 
! ...Proj^-1 * rhs
  call DPPTRS('U',3*nrdofQ_f,1,ProjF,rhsF,3*MAXTRIAQ,info1)
  if (info1.ne.0) then
    write(*,*) 'Dirichlet_f_proj - ProjF System Solution after Cholesky: info1 = ',info1 ; stop
  endif
! 
! Write the output
! number of dofs
  nrdofU = nrdofH_f; nrdofF = nrdofQ_f
! values of dofs
  zdofU = 0.d0; zdofF = 0.d0

  do k1 = 1,nrdofH_f
    k=3*(k1-1)
    zdofU(1:3,k1) = rhsU(k+1:k+3)
  enddo
  do k1 = 1,nrdofQ_f
    k=3*(k1-1)
    zdofF(1:3,k1) = rhsF(k+1:k+3)
  enddo

end subroutine

subroutine update_Ddof_poly      
      use data_structure3D_poly
      use connectivity_poly
      use geometry_polydpg
      use physics
      use parameters
#include "syscom.blk"
      integer, allocatable :: nverts_f(:),ibface(:)
      real*8, allocatable :: xverts_f(:,:)
      real*8, dimension(3) :: fn
      real*8, dimension(3,MAXTRIAH) :: zdofU
      real*8, dimension(3,MAXTRIAQ) :: zdofF
! c-----------------------------------------------------------------------
      iprint=0
! c-----------------------------------------------------------------------
! c
! c  Update   F A C E   dof for Dirichlet nodes
! c
      if (iprint.ge.1) write(*,*) 'update_Ddof_poly: MID-FACE NODES...'

! c
      allocate(ibface(NUMF))
      ibface = 0
      nr_bf = 0

      do nf=1,NUMF
! c  .......get global node number
        mdlf = NRELES + NUMV + NUME + nf        
!       check if face actually corresponds to an active node
        if (NODES(mdlf)%act.ne.0) then
          ibc = NODES(mdlf)%bcond
          if (ibc.gt.0) then
            nr_bf = nr_bf + 1
            ibface(nr_bf) = mdlf
          endif
        endif
      enddo

!$OMP PARALLEL
!$OMP DO PRIVATE(mdlf,ibc,nord,nrv_f,nverts_f,jv,xverts_f,fn,nrdofU,nrdofF,zdofU,zdofF) &
!$OMP SCHEDULE(DYNAMIC)
      do nf = 1, nr_bf
! c  .......get global node number
        mdlf = ibface(nf)
! !       check if face actually corresponds to an active node
!         if (NODES(mdlf)%act.eq.0) cycle
        ! 
        ibc = NODES(mdlf)%bcond
        nord= NODES(mdlf)%order
        ! if (iprint.ge.1) write(*,*) 'nf,mdlf,ibc=',nf,mdlf,ibc
        ! call find_ndof_poly(mdlf, ndofH,ndofE,ndofV,ndofQ)
! c     -- Face L2 --
        if (ibc.gt.0) then
          call face_vert_list(mdlf,nrv_f,nverts_f)
          allocate(xverts_f(3,nrv_f))
          do jv=1,nrv_f
            xverts_f(1:3,jv) = NODES(nverts_f(jv))%coord(1:3,1)
          enddo
          ! deallocate(nverts_f)
  !       get face normal in element coordinates
          call face_normal(xverts_f(:,1:3),fn)

          call Dirichlet_f_proj(mdlf,nrv_f,xverts_f,nord,nrdofU,nrdofF,zdofU,zdofF)
          deallocate(xverts_f,nverts_f)
        endif
        select case (ibc)

        case (1)
          ! first index in zdofF is for variable (every attribute component)
          ! second index in zdofF is for first degree of freedom: constants
          ! PRESCRIBED DISPLACEMENTS: 0.d0 ON INDICATED FACES
          
          NODES(mdlf)%zdofU(1:3,1:nrdofU) = zdofU(1:3,1:nrdofU)
          NODES(mdlf)%index(1:3)          = 11
          
        case (2)
          ! PRESCRIBED TRACTION: -pressure, normal to surface          
          ! NODES(mdlf)%zdofF(4:6,1) = -1.d0*fn(1:3)
          NODES(mdlf)%zdofF(1:3,1:nrdofF) = zdofF(1:3,1:nrdofF)
          NODES(mdlf)%index(4:6)          = 7

        case (3)
          ! CAUCHY: prescribed DISPLACEMENTS and TRACTION
          ! NODES(mdlf)%zdofF(1:3,1) = 0.d0
          NODES(mdlf)%zdofU(1:3,1:nrdofU) = zdofU(1:3,1:nrdofU)
          ! NODES(mdlf)%zdofF(4:6,1) = 0.d0!-1.d0*fn(1:3)
          NODES(mdlf)%zdofF(1:3,1:nrdofF) = zdofF(1:3,1:nrdofF)
          NODES(mdlf)%index(1:3)          = 11
          NODES(mdlf)%index(4:6)          = 7
        case (5)
          ! SYMMETRY B.C. on plane xy
          NODES(mdlf)%zdofU(  3,1:nrdofU) = 0.d0
          NODES(mdlf)%zdofF(1:2,1:nrdofF) = 0.d0
          NODES(mdlf)%index(  3)          = 11
          NODES(mdlf)%index(4:5)          = 7
        case (6)
          ! SYMMETRY B.C. on plane yz
          NODES(mdlf)%zdofU(  1,1:nrdofU) = 0.d0
          NODES(mdlf)%zdofF(2:3,1:nrdofF) = 0.d0
          NODES(mdlf)%index(  1)          = 11
          NODES(mdlf)%index(5:6)          = 7
        case (7)
          ! SYMMETRY B.C. on plane zx
          NODES(mdlf)%zdofU(  2,1:nrdofU) = 0.d0
          NODES(mdlf)%zdofF(1  ,1:nrdofF) = 0.d0
          NODES(mdlf)%zdofF(3  ,1:nrdofF) = 0.d0
          NODES(mdlf)%index(2)            = 11
          NODES(mdlf)%index(4)            = 7
          NODES(mdlf)%index(6)            = 7
        case default
          NODES(mdlf)%zdofU = ZERO
          NODES(mdlf)%zdofF = ZERO
        end select
!        PROJECTION BASED INTERPOLATION - DISABLED
! if (associated(NODES(mdlf)%zdofF)) then
        !   call dhpface_poly(mdlf,NODES(mdlf)%zdofF)               
        ! endif
! c
!         if (iprint.eq.4) then
!           call find_ndof_poly(mdlf, ndofH,ndofE,ndofV,ndofQ)
!           do k=1,ndofV
!             write(*,7012) k,NODES(mdlf)%zdofV(1:NRFVAR,k)
! #if C_MODE
!  7012       format(i2,2x,30(/,5(2e12.5,2x)))
! #else
!  7012       format(i2,2x,30(/,10(e12.5,2x)))
! #endif
!           enddo
!         endif
      enddo
!$OMP END DO
!$OMP END PARALLEL

      deallocate(ibface)
! c
      if (iprint.ge.1) write(*,*) 'update_Ddof_poly: DONE...'

end subroutine update_Ddof_poly