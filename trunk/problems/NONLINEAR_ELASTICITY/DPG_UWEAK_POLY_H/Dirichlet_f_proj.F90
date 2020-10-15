
subroutine Dirichlet_f_proj(mdlf,nrv_f,xverts_f,nord,nrdofU,nrdofF,zdofU,zdofF)

  use control,       only: NEXACT, INTEGRATION
  use parameters,    only: MAXTRIAQ,MAXTRIAH
  use parametersdpg, only: MAXNINT2ADD, NORD_ADD
  use geometry_polydpg
  use data_structure3D_poly
  use gauss_quadrature

  implicit none

  integer                   ,      intent(in ) :: nrv_f,nord,mdlf
  real*8, dimension(3,3),          intent(in ) :: xverts_f
  integer,                         intent(out) :: nrdofU,nrdofF
  real*8, dimension(3,MAXTRIAH),   intent(out) :: zdofU
  real*8, dimension(3,MAXTRIAQ),   intent(out) :: zdofF

  real*8, dimension(MAXTRIAH,3) :: zdofU_tr
  
  real*8                       :: rjac_f,b_area,r_f,wa,weight , scale , xi , ejac , rjac
  real*8, dimension(3)         :: b0,b1,b2,x_aux,fn,x, gradu_n , rt
  real*8, dimension(3,3)       :: b_aff , dxdxi, dxidx
  real*8, dimension(3,2)       :: covariant
  real*8, dimension(2,2)       :: dxdxi_bf
  real*8, dimension(2)         :: t
  real*8, dimension(2,MAXNINT2ADD) :: tloc
  real*8, dimension(MAXNINT2ADD)   :: wt
  real*8, dimension(3,MAXNINT2ADD) :: xloc

  ! solution variables
  real*8, dimension(3)         :: u,sigma_n,dudt,dudt_diff,dudt_proj
  real*8, dimension(3,3)       :: gradu, gradu_proj, gradu_t, gradu_t_diff
  real*8, dimension(3,3)       :: epsilon
  real*8, dimension(3,3)       :: sigma
  real*8, dimension(3)         :: divsigma
  real*8, dimension(MAXTRIAH)  :: shapH_f, dshapHdt
  real*8, dimension(MAXTRIAH,3):: rhsU
  real*8, dimension(2,MAXTRIAH):: gradH_f
  real*8, dimension(3,MAXTRIAH):: gradH_t
  real*8, dimension(MAXTRIAH*(MAXTRIAH+1)/2) :: ProjU
  real*8, dimension(MAXTRIAQ)  :: shapQ_f
  real*8, dimension(MAXTRIAQ,3):: rhsF
  real*8, dimension(MAXTRIAQ*(MAXTRIAQ+1)/2) :: ProjF

  integer, dimension(5) :: nordf
  integer, dimension(4) :: norie

  integer :: k1,k2,ipt,nint,i,jv,info,info1,icomp,jcomp,m1,m2,k,nrdofQ_f,nrdofH_f , imat,         &
             nrebubU,nrfbubU, ie, nint1
! statement function for packed format index for symmetric matrices
  integer :: nk
  nk(k1,k2) = (k2-1)*k2/2+k1

  if (nrv_f.gt.3) then
    write(*,*) 'Dirichlet_f_proj: This problem only supports triangular faces; mdlf, nrv_f=', mdlf,nrv_f
    stop
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up face order in usual structure
  nordf = 0
  nordf(1:4) = nord
! set orientation vector to zero
  norie = 0

  norie(1:3) = NODES(mdlf)%face_orient(1:3)
! in this case we drop the edge orientation convention for general 
! polygons (where we had CCW as the positive orientation). Here,
! we want v0->v2 to be the positive orientation for edge no. 3 so that
! we can naturally use the orientation embedded shape functions. Then,
  norie(3) = 1 - norie(3)

  ! nedges(1:3)= NODES(mdlf)%faces(1:3)

  nrdofU = (nord+1)*(nord+2)/2
  nrebubU = max(0,nord-1)
  nrfbubU = max(0,(nord-1)*(nord-2)/2)

  ! call face_vert_list(mdlf,nrv,nverts)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! initialize temporary(transposed) dofs
  zdofU_tr = 0.d0
  ! xnod(1:3,1:3) = xverts_f


! get normal vector
  call face_normal(xverts_f(:,1:3),fn)
! get face affine coordinates points b0,b1,b2
  call poly_affine_2d(xverts_f,nrv_f,b0,b1,b2)  ! when face is a triangle, b0,b1,b2 are explicitly vertices 1,2,3, respectively.
  b_aff(:,1) = b0
  b_aff(:,2) = b1
  b_aff(:,3) = b2
! get area of triangle b0b1b2
  call trian_centr_area(b0,b1,b2,x_aux,b_area)
  rjac_f = 2.d0*b_area


! determine approximate scale by computig dist(b_2 , b_1 )
  call norm( b_aff(:,3)-b_aff(:,2) , scale )

! construct linear transformation associated to affine map. Interpreted as classic jacobian
  dxdxi(:,1) = b1 - b0
  dxdxi(:,2) = b2 - b0
  dxdxi(:,3) = fn*scale
! get inverse and determinant of linear transformation
  call geom(dxdxi,dxidx,rjac,info)
  covariant = transpose(dxidx(1:2,1:3))

! 
! get quadrature points for polygon using the affine coordinates
  INTEGRATION = max(0, NORD_ADD - 1)
  call set_quadrature_affine_face_DPG(mdlf,nord,nrv_f,xverts_f,b_aff,nint,tloc,wt)
  INTEGRATION = 0
! get physical coordinates for integration points so that exact solution can be evaluated 
  xloc = 0.d0
  call poly_affine_2d_to_x(tloc(1:2,1:nint),nint,b_aff,xloc)
! 
! 
! ***************************************************************
!                  PROJECTION OF TRACE \hat u
! ***************************************************************
! 
! H1 PB-INTERPOLATION ON VERTICES
  do jv = 1,3
    x(1:3) = xverts_f(1:3,jv)
    !   compute exact soution
    call find_material(x,fn,imat)
    call elast_solution(imat,X, u,gradu,epsilon,sigma,divsigma)
    ! vertices degrees of freedom equal u(x_vert)
    do icomp = 1,3
      zdofU_tr(jv,icomp) = u(icomp)
    enddo
  enddo

! H1 PB-INTERPOLATION ON EDGES
if (nrebubU.gt.0) then

  ProjU = 0.d0
  rhsU = 0.d0

!  ...loop through element edges
  do ie=1,3
!   tangential vector
    select case(ie)
    case(1)
      rt = b1 - b0
    case(2)
      rt = b2 - b1
    case(3)
      rt = b2 - b0
    end select
    call norm(rt,ejac)
    rt = rt / ejac
!
!  # integration points
    nint1 = nord + 1
!
!  loop through integration points
    do ipt=1,nint1
!
!  ...edge coordinate
      xi = XIGAUS1(ipt,nint1); wa = WAGAUS1(ipt,nint1)
      select case(ie)
      case(1)
        t = (/        xi , 0.d0 /)
      case(2)
        t = (/ 1.d0 - xi ,   xi /)
      case(3)
        t = (/      0.d0 ,   xi /)
      end select
!
!  ...determine face H1 shape functions
      shapH_f = 0.d0; gradH_f = 0.d0
      call shape2DH('tria',t,nordf,norie, nrdofH_f,shapH_f,gradH_f)
! ... use covariant map to get gradient of physical shape functions
      gradH_t = matmul(covariant  ,gradH_f)
!  ...dot product with unit tangent vector to get physical tangential 
!     derivative of shape functions
      dshapHdt(:) = rt(1)*gradH_t(1,:)                               &
                  + rt(2)*gradH_t(2,:)                               &
                  + rt(3)*gradH_t(3,:)
!  ...multiply dshapHdt with PB-INTERPOLATION coefficients obtained so far
      dudt_proj = matmul(dshapHdt,zdofU_tr)
!  ...we need to do the same for the exact solution u
!  ...get physical coordinates
      x = b0 + matmul(dxdxi(:,1:2),t)
!  ...compute exact soution
      call find_material(x,fn,imat)
      call elast_solution(imat,X, u,gradu,epsilon,sigma,divsigma)
!  ...physical tangential derivative of exact solution u
      dudt = matmul(gradu,rt)
!  ...compute difference
      dudt_diff = dudt-dudt_proj
!  ...line integral weight
      weight = ejac*wa
!  ...loop over U components
      do icomp = 1,3
!    ...loop over edge bubble shape functions
        do k1=1,nrebubU
!      ...retrieve shape function index          
          m1=(ie-1)*nrebubU+k1
          rhsU(m1,icomp) = rhsU(m1,icomp) + dudt_diff(icomp)*dshapHdt(3+m1) * weight
        enddo
      enddo
      do k1=1,nrebubU
        m1=(ie-1)*nrebubU+k1
        do k2=k1,nrebubU
          m2=(ie-1)*nrebubU+k2
          k = nk(m1,m2)
          ProjU(k) = ProjU(k) + dshapHdt(3+m1)*dshapHdt(3+m2)*weight
        enddo
      enddo
!  .....end of loop through integration points
    enddo
!  ...end of loop through edges
  enddo
! apply Cholesky Factorization to Proj
  call DPPTRF('U',3*nrebubU,ProjU,info)
  if (info.ne.0) then
    write(*,*) 'Dirichlet_f_proj - ProjU Gram Cholesky Factorization (edges): info = ',info
    stop
  endif
! 
! ...Proj^-1 * rhs
  call DPPTRS('U',3*nrebubU,3,ProjU,rhsU,MAXTRIAH,info1)
  if (info1.ne.0) then
    write(*,*) 'Dirichlet_f_proj - ProjU System Solution after Cholesky (edges): info1 = ',info1 ; stop
  endif
!
  do icomp=1,3
    zdofU_tr(4:3+3*nrebubU,icomp) = rhsU(1:3*nrebubU,icomp)
  enddo
!
endif
!
!
! H1 PB-INTERPOLATION ON FACES
if (nrfbubU.gt.0) then

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
    shapH_f = 0.d0; gradH_f = 0.d0
    call shape2DH('tria',t,nordf, norie, nrdofH_f,shapH_f,gradH_f)

    gradH_t = matmul(covariant  ,gradH_f)
!
    gradu_proj = matmul(gradH_t,zdofU_tr) ! this is already the tangential gradient
! 
!   compute exact soution
    call find_material(x,fn,imat)
    call elast_solution(imat,X, u,gradu,epsilon,sigma,divsigma)

    gradu_n = matmul(gradu,fn)
    do icomp=1,3
      ! call scalar_product(fn,gradu_proj(:,icomp),gradu_p_n)
      ! gradu_t_p(:,icomp) = gradu_proj(:,icomp) - fn*gradu_p_n      
      gradu_t(:,icomp) = gradu(icomp,:) - fn*gradu_n(icomp)
    enddo
    gradu_t_diff = gradu_t - gradu_proj! - gradu_t_p
!
    ! call poly_x_to_affine_2d(gradu_t_diff,3,b_aff,gradu_t_aff)
!
!
    weight = wa
!
!...OUTER loop through Face H1 shape functions
    do k1=1,nrfbubU
!     
      m1 = nrdofU - nrfbubU + k1
      do jcomp=1,3
        ! rhsU(k1,jcomp) = rhsU(k1,jcomp) + ( gradu_t_aff(1,jcomp)*gradH_f(1,k1) +            &
        !                                     gradu_t_aff(2,jcomp)*gradH_f(2,k1)   ) * weight
        rhsU(k1,jcomp) = rhsU(k1,jcomp) + ( gradu_t_diff(1,jcomp)*gradH_t(1,m1) +            &
                                            gradu_t_diff(2,jcomp)*gradH_t(2,m1) +            &
                                            gradu_t_diff(3,jcomp)*gradH_t(3,m1)   ) * weight
      enddo
!    ...INNER loop through trial dofs
      do k2=k1,nrfbubU
        k = nk(k1,k2)
        m2 = nrdofU - nrfbubU + k2
        ! ProjU(k) = ProjU(k) + ( gradH_f(1,k2)*gradH_f(1,k1) +             &
        !                         gradH_f(2,k2)*gradH_f(2,k1)   ) * weight
        ProjU(k) = ProjU(k) + ( gradH_t(1,m2)*gradH_t(1,m1) +             &
                                gradH_t(2,m2)*gradH_t(2,m1) +             &
                                gradH_t(3,m2)*gradH_t(3,m1)   ) * weight
            ! write(*,*) 'Dirichlet_f_proj :   m1,m2=',m1,m2
            ! write(*,*) 'Dirichlet_f_proj :   k1,k2=',k1,k2
            ! write(*,*) 'Dirichlet_f_proj :   jcomp,icomp=',jcomp,icomp
            ! write(*,*) 'Dirichlet_f_proj :   shapQ_f(k2),shapQ_f(k1)=',shapQ_f(k2),shapQ_f(k1)
      enddo
    enddo
! integration point loop ends
  enddo

! apply Cholesky Factorization to Proj
  call DPPTRF('U',nrfbubU,ProjU,info)
  if (info.ne.0) then
    write(*,*) 'Dirichlet_f_proj - ProjU Gram Cholesky Factorization (faces): info = ',info
    stop
  endif
! 
! ...Proj^-1 * rhs
  call DPPTRS('U',nrfbubU,3,ProjU,rhsU,MAXTRIAH,info1)
  if (info1.ne.0) then
    write(*,*) 'Dirichlet_f_proj - ProjU System Solution after Cholesky (faces): info1 = ',info1 ; stop
  endif

  do icomp=1,3
    zdofU_tr(4+3*nrebubU:nrdofH_f,icomp) = rhsU(1:nrfbubU,icomp)
  enddo

endif
! Write the output for U variables
! number of dofs
  ! nrdofU = nrdofH_f;
! values: return result in expected array shape
  zdofU = 0.d0
  zdofU(1:3,1:MAXTRIAH) = transpose(zdofU_tr)

! ***************************************************************
!             PROJECTION OF NORMAL FLUX \hat Sigma_n
! ***************************************************************
!  initialize projection matrix and rhs

  ! nrdofF = (nord  )*(nord+1)/2
  ProjF = 0.d0
  rhsF = 0.d0
!    .....loop through face integration points
  do ipt=1,nint
!
!    .......face coordinates
    t(1:2) = tloc(1:2,ipt); wa = wt(ipt); x(1:3) = xloc(1:3,ipt)
!
!   evaluate shape functions of normal flux space (polynomials of order p-1)
    shapQ_f = 0.d0
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
        rhsF(k1,jcomp) = rhsF (k1,jcomp) + sigma_n(jcomp)*shapQ_f(k1)*weight
        ! m1 = (k1-1)*3+jcomp
      enddo
!    ...INNER loop through trial dofs
      do k2=k1,nrdofQ_f
              k = nk(k1,k2)
              ProjF(k) = ProjF(k) + shapQ_f(k2)*shapQ_f(k1)*weight
          ! write(*,*) 'Dirichlet_f_proj :   m1,m2=',m1,m2
          ! write(*,*) 'Dirichlet_f_proj :   k1,k2=',k1,k2
          ! write(*,*) 'Dirichlet_f_proj :   jcomp,icomp=',jcomp,icomp
          ! write(*,*) 'Dirichlet_f_proj :   shapQ_f(k2),shapQ_f(k1)=',shapQ_f(k2),shapQ_f(k1)        
      enddo              
    enddo
! 
! integration point loop ends
  enddo
! 
!! apply Cholesky Factorization to Proj
  call DPPTRF('U',nrdofQ_f,ProjF,info)
  if (info.ne.0) then
    write(*,*) 'Dirichlet_f_proj - ProjF Gram Cholesky Factorization: info = ',info
    stop
  endif
! 
! ...Proj^-1 * rhs
  call DPPTRS('U',nrdofQ_f,3,ProjF,rhsF,MAXTRIAQ,info1)
  if (info1.ne.0) then
    write(*,*) 'Dirichlet_f_proj - ProjF System Solution after Cholesky: info1 = ',info1 ; stop
  endif
! 
! Write the output for F variables
! number of dofs
  nrdofF = nrdofQ_f
! values: return result in expected array shape
  zdofF = transpose(rhsF)
! 
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

