!--------------------------------------------------------------------------
!> Purpose : return stiffness matrix and load vector for element
!!
!! @param[in]  Mdle   - an element (middle node) number
!! @param[out] Nrdof  - number of dof for a single component
!! @param[out] Itest  - index for assembly
!! @param[out] Itrial - index for assembly
!--------------------------------------------------------------------------
!
subroutine elem(Mdle, Itest,Itrial)

  use physics   , only : NR_PHYSA
  use data_structure3D_poly
!--------------------------------------------------------------------------
  implicit none
  integer,                    intent(in)  :: Mdle
  integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!--------------------------------------------------------------------------

  Itest (1:NR_PHYSA) = 1
  Itrial(1:NR_PHYSA) = 1

  call elem_DPG_UWEAK_poly(Mdle)

!
end subroutine elem


!------------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for primal DPG elasticity problem
!! @param[in]  Mdle      - middle node number
!------------------------------------------------------------------------------------------
!
!                |   \tau \in H(div)^3   |       v \in (H1)^3     |
!
!                   - <\hat u,(\tau n)>  +            0
!  +                        0            + - <\hat (\sigma n),v>
!  + \int_\Omega [   u \cdot div(\tau)   +            0           ]
!  + \int_\Omega [    A \sigma : \tau    +    \sigma : grad(v)    ]
!  + \int_\Omega [     \omega : \tau     +            0           ]
!  = \int_\Omega [          0            +        f \cdot v       ]
!
!------------------------------------------------------------------------------------------
!
subroutine elem_DPG_UWEAK_poly(Mdle)

      use uweak_module_poly
      use connectivity_poly
      use control, only: INTEGRATION
      use parameters
      use parametersDPG
      use data_structure3D_poly
      use element_data_poly
      use isotropic_elast_material
      use physics   , only : NR_PHYSA
      use assembly_poly  , only : ALOC,BLOC,NR_RHS
      use common_prob_data, only: SYMMETRY_TOL, TEST_NORM
!------------------------------------------------------------------------------------------
      implicit none
      integer,                                  intent(in)  :: Mdle
!  ...element and face type
      character(len=4) :: etype,ftype
!
!  ...number of topological entities (vertices,edges,faces)
      integer :: nrv,nre,nrf
!
!  ...element and face order, enriched order
      integer, dimension(19) :: norder
      integer, dimension(5)  :: nordf
      integer                :: nordP
!
!  ...node edge and face orientations
      integer, dimension(12) :: nedge_orient
      integer, dimension(6)  :: nface_orient
!
!  ...SHAPE FUNCTIONS
!     H1  (geometry and trial)
      real*8, dimension(  MAXtetraH)  :: shapH
      real*8, dimension(3,MAXtetraH)  :: gradH
      integer                         :: nrdofH
!     H(div)  (trial)
      real*8, dimension(3,MAXtetraV)  :: shapV
      real*8, dimension(  MAXtetraV)  :: divV
      real*8, dimension(  MAXtetraV)  :: shapV_n
      integer                         :: nrdofV
!     Face L2 (trial)      
      real*8, dimension(MAXtriaQ   )  :: shapQ_f
!     L2  (trial)
      real*8, dimension(  MAXtetraQ)  :: shapQ
      integer                         :: nrdofQ
!     H1   (test)
      real*8, dimension(  MAXtetraHH) :: shapHH
      real*8, dimension(3,MAXtetraHH) :: gradHH
      integer                         :: nrdofHH
!     H(div)  (test)
      real*8, dimension(3,MAXtetraVV) :: shapVV
      real*8, dimension(  MAXtetraVV) :: divVV
      real*8, dimension(  MAXtetraVV) :: shapVV_n
      integer                         :: nrdofVV
!
!
!  ...geometry
      real*8, dimension(3,MAXtetraH) :: xnod
      real*8, dimension(3)           :: xi,x,rn,vn1,vn2,vn3,vn4,xg_e,fn,xg_f,vt1,vt2,vt3
      real*8, dimension(3,3)         :: dxdxi_e,dxidx_e,Qrot_e,Qrot_f,maptet,maptetinv
      real*8, dimension(2,2)         :: dxdxi_f,dxidx_f,dxdxi_bf,maptri
      real*8, dimension(2)           :: t
      real*8, dimension(3,2)         :: dxidt,dxdt,rt
      integer                        :: nsign
!
!  ...tensors in physical coordinates
      real*8, dimension(3,3,3,3) :: A
!
!  ...source term (don't need Neumann term)
      real*8, dimension(3,MAXNRHS) :: fval
!
!  ...3D quadrature data
      real*8, dimension(3,MAXNINT3ADD) :: xiloc
      real*8, dimension(MAXNINT3ADD)   :: wxi
!
!  ...2D quadrature data for boundary terms
      real*8, dimension(2,MAXNINT2ADD) :: tloc
      real*8, dimension(MAXNINT2ADD)   :: wt
!
!  ...miscellaneous
      integer :: i,j,k,l,m,n,k1,k2,k3,k4,k5,m1,m2,m3,m4,m5,n1,n2,n3,ipt,ifc,  &
                 icomp,jcomp,nint,iprint,iflag,info,info1,info2,info3,  &
                 kH,kV,kQ,lH,lV,lQ,kmin,kmax,enrdof,       &
                 Nrv_f,mdlf,jv,jf,loc,nrdofQ_f, &
                 gdump,bdump1,bdump2,bdump3,bdump4,bdump5,bdump6,bdump7

      integer, dimension(NR_PHYSA) :: ndofphysics
      real*8  :: weight,wa,rjac,brjac,tmp,diffmax,dmax,rjac_e,rjac_f,rjac_bf,r_e, &
                 maptetdet,maptridet,Area_f,r_f
      
      integer, allocatable :: nfaces(:),Norientf(:),Nedges(:),Nverts(:),nverts_f(:)
      real*8, allocatable :: Xvertl(:,:),xverts_f(:,:),Rotverts_f(:,:),Gramtmp(:,:)
!
!  ...LAPACK stuff
      character uplo,transa,transb
! NOTE: nk is a "statement function"
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
      iprint=0

!     Compute jacobian for the bounding element w.r.t master tetrahedron
      etype = 'tetr'
      ftype = 'tria'
      vn1(1:3) = BOUND_TETRA(1:3,1)
      vn2(1:3) = BOUND_TETRA(1:3,2)
      vn3(1:3) = BOUND_TETRA(1:3,3)
      vn4(1:3) = BOUND_TETRA(1:3,4)
      call tetra_affine_map(vn1,vn2,vn3,vn4,dxdxi_e)
      call geom(dxdxi_e,dxidx_e,rjac_e,iflag)
!     Compute jacobian for the bounding face w.r.t master triangle
      vt1=0.d0; vt1(1:2) = BOUND_TRIAN(1:2,1)
      vt2=0.d0; vt2(1:2) = BOUND_TRIAN(1:2,2)
      vt3=0.d0; vt3(1:2) = BOUND_TRIAN(1:2,3)
      call trian_affine_map(vt1(1:2),vt2(1:2),vt3(1:2),dxdxi_bf)
      rjac_bf = dxdxi_bf(1,1)*dxdxi_bf(2,2)-dxdxi_bf(1,2)*dxdxi_bf(2,1)
      ! dxidx_f = reshape(                                           &
      !   (/dxdxi_f(2,2),-dxdxi_f(2,1),-dxdxi_f(1,2),dxdxi_f(1,1)/),   &
      !   (/2,2/))
      ! dxidx_f = dxidx_f / rjac_f

!
!  ...element type
      ! etype = NODES(Mdle)%type
      ! nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)

!     retrieve element node data
      call elem_nodes_poly(Mdle,Nfaces,Norientf,Nrf,Nedges,Nre,Nverts,Nrv)
! !     allocate array for vertices and copy coordinates
!       allocate(xnod(3,nrv))
!       do jv = 1,nrv
!         xnod(1:3,jv) = NODES(Nverts(jv))%coord(1:3,1)
!       enddo


      allocate(Xvertl(3,Nrv))
!     compute centroid, radius, rotation (unitary) matrix and vertices in local coord
      call elem_vertex_coordl(Mdle,Nfaces(1:Nrf),Norientf(1:Nrf),Nrf,Nverts(1:Nrv),Nrv, &
                              Xg_e,R_e,Qrot_e,Xvertl)
!     save the first vertex of the polyhedron as vertex 1 of all subtetrahedra
      vn1(1:3) = Xvertl(1:3,1)
!     compute actual jacobian, its inverse and det: Jac <--- R_e * Qrot_e * Jac
      dxdxi_e = R_e * matmul( Qrot_e , dxdxi_e)
      dxidx_e = 1 / R_e * matmul( dxidx_e, transpose(Qrot_e))
      rjac_e = R_e**3.d0 * rjac_e

      rjac_bf = rjac_bf * R_e**2.d0

!     get order, in the norder structure of a uniform tetrahedral element
      norder = 0
      norder(1:15) = NODES(Mdle)%order      
!  ...set the enriched order of appoximation
      nordP = NODES(Mdle)%order + NORD_ADD

!    ...set up the element quadrature
      INTEGRATION = NORD_ADD
      call set_3Dint(etype,norder, nint,xiloc,wxi)
      INTEGRATION = 0

!  ...initialize the enriched local element stiffness matrices and load vectors
      EnrTraceDispl=ZERO; EnrTraceStress=ZERO; EnrFieldDispl=ZERO; EnrFieldStress=ZERO; EnrFieldOmega=ZERO
      EnrLoad=ZERO
! !  ...initialize the Gram matrix
      Gram=ZERO
!
!     loop over faces
      do jf = 1,nrf
        mdlf = Nfaces(jf)
        call face_vert_list(mdlf,nrv_f,nverts_f)
        allocate(xverts_f(3,nrv_f))
        do jv=1,Nrv_f
          xverts_f(1:3,jv) = NODES(nverts_f(jv))%coord(1:3,1)-Xg_e(1:3)
        enddo
!       transform vertices to element coordinate system
        xverts_f = matmul(transpose(Qrot_e),xverts_f)
        xverts_f = xverts_f / R_e
!       get face normal in element coordinates
        call face_normal(xverts_f(:,1:3),fn)
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!      
!       Check if vertex 1 is the face vertex list.
!       If not, the face is a valid pyramid/tetrahedron base.
        call locate(Nverts(1),nverts_f,nrv_f,loc)
        if (loc.eq.0) then
!         complete subtetrahedron vertex list vn2,3,4
          vn2(1:3) = xverts_f(1:3,1)
!         loop on subtetrahedra
          do kv = 1,Nrv_f-2
!         check if face normal locally goes outward (Norientf = 0)
            if(Norientf(jf).eq.0) then
              vn3(1:3) = xverts_f(1:3,kv+1)
              vn4(1:3) = xverts_f(1:3,kv+2)
            else
!           if normal goes inward swap order of vertices 3 and 4            
              vn3(1:3) = xverts_f(1:3,kv+2)
              vn4(1:3) = xverts_f(1:3,kv+1)
            endif

            call tetra_affine_map(vn1,vn2,vn3,vn4,maptet)

            call geom(maptet,maptetinv,maptetdet,iflag)
!  
!        ...loop through integration points
            do ipt=1,nint
              xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt)*maptetdet
!  
!             to compute physical coordinate and evaluate functions:
!
!             first, find xi in the local coordinates
              xi = matmul(maptet,xi) + vn1
!       
!        .....Compute shape functions needed for test/trial field variables and geometry
!             L2 (field trial)
              call shape3Q(etype,xi,norder, nrdofQ,shapQ)              
!             H1 (test)
              call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)              
!             H(div) (test)
              call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!
!       .....Change coordinates so the shape functions are on the physical element
!             L2 (trial)
              shapQ(1:nrdofQ) = shapQ(1:nrdofQ)/rjac_e
!             H1 (test)
              do k=1,nrdofHH
                gradHH(1:3,k) = gradHH(1,k)*dxidx_e(1,1:3)  &
                              + gradHH(2,k)*dxidx_e(2,1:3)  &
                              + gradHH(3,k)*dxidx_e(3,1:3)
              enddo
!             H(div) (test)
              do k=1,nrdofVV
                shapVV(1:3,k) = dxdxi_e(1:3,1)*shapVV(1,k)  &
                              + dxdxi_e(1:3,2)*shapVV(2,k)  &
                              + dxdxi_e(1:3,3)*shapVV(3,k)
              enddo
              shapVV(1:3,1:nrdofVV) = shapVV(1:3,1:nrdofVV)/rjac_e
              divVV(1:nrdofVV) = divVV(1:nrdofVV)/rjac_e
!
!             second: descale, rotate and translate
              x = xi * r_e
              x = matmul(Qrot_e,x)
              x = x + xg_e

!        .....integration weight
              weight = wa*rjac_e
!
!        .....compute the compliance tensor
              call getA(x, A)
!
!        .....get the source term
              call getf(Mdle,x, fval)
!
!
!    P A R T  1 : go through \tau\in H(div)^3 test space (this fills the first set of rows)
!
!
!        .....FIRST OUTER loop through enriched H(div) dofs
              do k1=1,nrdofVV
!        .......OUTER loop through components
                do jcomp=1,3
                  m1 = (k1-1)*3+jcomp
!     
!                 E N R I C H E D   L O A D   V E C T O R
!
!   0
!
!           G R A M   M A T R I X
!
                  select case(TEST_NORM)
!
!   (\tau_2,\tau)+(div(\tau_2),div(tau))
!
                  case(2)
                    do m2=m1,3*nrdofVV
                      icomp = mod(m2-1,3)+1
                      if (icomp.eq.jcomp) then
                        k = nk(m1,m2)
                        k2 = int((m2-1)/3)+1
                        Gram(k) = Gram(k)  &
                                + ( divVV(k1)*divVV(k2)  &
                                  + shapVV(1,k1)*shapVV(1,k2)  &
                                  + shapVV(2,k1)*shapVV(2,k2)  &
                                  + shapVV(3,k1)*shapVV(3,k2) )*weight
                      endif
                    enddo
                  end select
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   C A U C H Y    S T R E S S   )
!
!   A \sigma : \tau
!
!          ( sigma1   sigma4  sigma5 )
!  sigma = ( sigma4   sigma2  sigma6 )
!          ( sigma5   sigma6  sigma3 )
!
!  .........INNER loop through trial dofs for Cauchy stress
                  do k3=1,nrdofQ
                    do icomp=1,6
                      m3 = (k3-1)*6+icomp
                      if (icomp.eq.1) then
                        tmp = 0.d0
                        do n=1,3
                          tmp = tmp  &
                              + A(1,1,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                        enddo
                      elseif (icomp.eq.2) then
                        tmp = 0.d0
                        do n=1,3
                          tmp = tmp  &
                              + A(2,2,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                        enddo
                      elseif (icomp.eq.3) then
                        tmp = 0.d0
                        do n=1,3
                          tmp = tmp  &
                              + A(3,3,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                        enddo
                      elseif (icomp.eq.4) then
                        tmp = 0.d0
                        do n=1,3
                          tmp = tmp  &
                              + A(1,2,jcomp,n)*shapQ(k3)*shapVV(n,k1) &
                              + A(2,1,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                        enddo
                      elseif (icomp.eq.5) then
                        tmp = 0.d0
                        do n=1,3
                          tmp = tmp  &
                              + A(1,3,jcomp,n)*shapQ(k3)*shapVV(n,k1) &
                              + A(3,1,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                        enddo
                      elseif (icomp.eq.6) then
                        tmp = 0.d0
                        do n=1,3
                          tmp = tmp  &
                              + A(2,3,jcomp,n)*shapQ(k3)*shapVV(n,k1) &
                              + A(3,2,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                        enddo
                      endif
                      EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                            + tmp*weight
                    enddo
                  enddo
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   u \cdot div(\tau)
!
!  .........INNER loop through trial dofs for displacement
                  do k4=1,nrdofQ
                    do icomp=1,3
                      m4 = (k4-1)*3+icomp
                      if (icomp.eq.jcomp) then
                        EnrFieldDispl(m1,m4) = EnrFieldDispl(m1,m4)  &
                                             + shapQ(k4)*divVV(k1)*weight
                      endif
                    enddo
                  enddo
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                   (   L A G R A N G E    M U L T I P L I E R   )
!
!   + omega : \tau
!
!          (   0     omega3  -omega2 )    1
!  omega = (-omega3    0      omega1 ) = --- ( grad(u)-grad(u)^T)
!          ( omega2  -omega1    0    )    2
!
!  .........INNER loop through trial dofs for omega (antisymmetric part of displacement gradient)
                  do k5=1,nrdofQ
                    do icomp=1,3
                      m5 = (k5-1)*3+icomp
                      if (icomp.eq.1) then
                        if (jcomp.eq.2) then
                          EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                               + shapQ(k5)*shapVV(3,k1)*weight
                        elseif (jcomp.eq.3) then
                          EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                               - shapQ(k5)*shapVV(2,k1)*weight
                        endif
                      elseif (icomp.eq.2) then
                        if (jcomp.eq.1) then
                          EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                               - shapQ(k5)*shapVV(3,k1)*weight
                        elseif (jcomp.eq.3) then
                          EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                               + shapQ(k5)*shapVV(1,k1)*weight
                        endif
                      elseif (icomp.eq.3) then
                        if (jcomp.eq.1) then
                          EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                               + shapQ(k5)*shapVV(2,k1)*weight
                        elseif (jcomp.eq.2) then
                          EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                               - shapQ(k5)*shapVV(1,k1)*weight
                        endif
                      endif
                    enddo
                  enddo
!  .....END OUTER LOOP through test stresses
                enddo
              enddo
!
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows)
!
!        .....SECOND OUTER loop through enriched H1 dofs
              do k1=1,nrdofHH
!        .......OUTER loop through components
                do jcomp=1,3
            ! counter of row
                  m1 = 3*nrdofVV+(k1-1)*3+jcomp
!
!
!            E N R I C H E D   L O A D   V E C T O R
!
!   f \cdot v
!
                  EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS) &
                                 + fval(jcomp,1:NR_RHS)*shapHH(k1)*weight
!
!            G R A M   M A T R I X
!
                  select case(TEST_NORM)
!
!   (v_2,v) + (grad(v_2),grad(v))
!
                  case(2)
                    do m2=m1,3*nrdofVV+3*nrdofHH
                      icomp = mod(m2-1,3)+1
                      if (icomp.eq.jcomp) then
                        k = nk(m1,m2)
                        k2 = int((m2-3*nrdofVV-1)/3)+1
                        Gram(k) = Gram(k)  &
                                + ( shapHH(k1)*shapHH(k2)  &
                                  + gradHH(1,k1)*gradHH(1,k2)  &
                                  + gradHH(2,k1)*gradHH(2,k2)  &
                                  + gradHH(3,k1)*gradHH(3,k2) )*weight
                      endif
                    enddo
                  end select
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   C A U C H Y    S T R E S S   )
!
!   + \sigma : grad(v)
!
!          ( sigma1   sigma4  sigma5 )
!  sigma = ( sigma4   sigma2  sigma6 )
!          ( sigma5   sigma6  sigma3 )
!
!  .........INNER loop through trial dofs for Cauchy stress
                  do k3=1,nrdofQ
                    do icomp=1,6
                      m3 = (k3-1)*6+icomp
                      if (icomp.eq.1) then
                        if (jcomp.eq.1) then
                          EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                                + shapQ(k3)*gradHH(1,k1)*weight
                        endif
                      elseif (icomp.eq.2) then
                        if (jcomp.eq.2) then
                          EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                                + shapQ(k3)*gradHH(2,k1)*weight
                        endif
                      elseif (icomp.eq.3) then
                        if (jcomp.eq.3) then
                          EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                                + shapQ(k3)*gradHH(3,k1)*weight
                        endif
                      elseif (icomp.eq.4) then
                        if (jcomp.eq.1) then
                          EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                                + shapQ(k3)*gradHH(2,k1)*weight
                        elseif (jcomp.eq.2) then
                          EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                                + shapQ(k3)*gradHH(1,k1)*weight
                        endif
                      elseif (icomp.eq.5) then
                        if (jcomp.eq.1) then
                          EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                                + shapQ(k3)*gradHH(3,k1)*weight
                        elseif (jcomp.eq.3) then
                          EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                                + shapQ(k3)*gradHH(1,k1)*weight
                        endif
                      elseif (icomp.eq.6) then
                        if (jcomp.eq.2) then
                          EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                                + shapQ(k3)*gradHH(3,k1)*weight
                        elseif (jcomp.eq.3) then
                          EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                                + shapQ(k3)*gradHH(2,k1)*weight
                        endif
                      endif
                    enddo
                  enddo
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   0
!
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                   (   L A G R A N G E    M U L T I P L I E R   )
!
!   0
!
!
!        .......END OUTER LOOP
                enddo
              enddo
!           end of integration point loop
            enddo


!         end of sub tetrahedra loop
          enddo
!
        endif
!
!-----------------------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L                                           |
!-----------------------------------------------------------------------------------
!
!       compute face centroid in element coord, compute radius and 
!       transform vertices to face coordinates
        call face_centroid(Xverts_f,Nrv_f,Xg_f,Area_f)
        call polyhedron_radius(xverts_f,nrv_f,xg_f,r_f)
!       update face jacobian
        rjac_f = rjac_bf * r_f**2.d0
        do jv = 1, nrv_f
          xverts_f(1:3,jv) = xverts_f(1:3,jv) - xg_f(1:3)
        enddo
        allocate(Rotverts_f(3,Nrv_f))
        call face_rotation_to_2D(Xverts_f,Nrv_f,Fn,Rotverts_f,Qrot_f)
        ! Rotverts_f = matmul(transpose(Qrot_f),Rotverts_f) / r_f
        Rotverts_f = Rotverts_f / r_f
!       2d shape functions jacobian is not updtaed since is not needed here.
!       However, if we require Piola maps on faces, we will need to fix this !!!
!
!       set up face order in usual structure
        nordf = 0
        nordf(1:4) = NODES(Mdle)%order
!  .....set up the face quadrature
        INTEGRATION = NORD_ADD
        call set_2Dint(ftype,nordf, nint,tloc,wt)
        INTEGRATION = 0
!       define first vertex of subtriangle always as vertex 1 of current face
        vt1(1:3) = Rotverts_f(1:3,1)
!
!   ....transform normal vector to physical coordinates, to multiply by shape funs
        fn = matmul(Qrot_e,fn)
!       check if face normal locally goes outward (Norientf = 0)
        if(Norientf(jf).ne.0) then
!       if normal goes inward, correct it
          fn = -1.d0*fn
        endif
        write(*,*)'elem: mdle, jf,fn=',mdle,jf,fn
!       loop on subtriangles
        do kv = 1, nrv_f-2
          vt2(1:3) = Rotverts_f(1:3,kv+1)
          vt3(1:3) = Rotverts_f(1:3,kv+2)

          call trian_affine_map(vt1(1:2),vt2(1:2),vt3(1:2),maptri)
          maptridet = maptri(1,1)*maptri(2,2)-maptri(1,2)*maptri(2,1)
!    .....loop through face integration points
          do ipt=1,nint
!
!    .......face coordinates
            t(1:2) = tloc(1:2,ipt); wa = wt(ipt)*maptridet
!
!           in order to compute physical coordinate and evaluate functions:
!
!           first, find t in the current face coordinates
            t = matmul(maptri,t) + vt1(1:2)
!           evaluate trial (2d face) and test (3d element) functions
            call shape2DQ(ftype,t,nordf, nrdofQ_f,ShapQ_f)
!           get element coordinates xi for face coordinates t
            xi(1:3) = (Qrot_f(1:3,1)*t(1)+Qrot_f(1:3,2)*t(2))*r_f+xg_f(1:3)
!           H1 (test)
            call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!           H(div) (test)
            call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!       .....Change coordinates so the shape functions are on the physical element
!           L2 (trial)
            shapQ_f(1:nrdofQ_f) = shapQ_f(1:nrdofQ_f)/rjac_f
!
!           H(div) (test)
            do k=1,nrdofVV
              shapVV(1:3,k) = dxdxi_e(1:3,1)*shapVV(1,k)  &
                            + dxdxi_e(1:3,2)*shapVV(2,k)  &
                            + dxdxi_e(1:3,3)*shapVV(3,k)
              shapVV_n(k) = shapVV(1,k)*fn(1)  &
                          + shapVV(2,k)*fn(2)  &
                          + shapVV(3,k)*fn(3)
            enddo

            shapVV_n(1:nrdofVV) = shapVV_n(1:nrdofVV)/rjac_e
!
!           get physical coordinates x: descale, rotate and translate
            x = xi * r_e
            x = matmul(Qrot_e,x)
            x = x + xg_e

            weight = wa*rjac_f
!
!
!    P A R T  1 : go through \tau\in H(div)^3 test space (this fills the first set of rows)
!
!
!  .......OUTER loop through enriched H(div) test functions
            do k1=1,nrdofVV
              do jcomp=1,3
                m1 = (k1-1)*3+jcomp
!
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   - <\hat u,(\tau n)>
!
!  ...........INNER loop through \hat{u} bdry trial dofs
                do k2=1,nrdofQ_f
                  do icomp=1,3
                    m2 = (k2-1)*3+icomp + 3*nrdofQ_f*(jf-1)
                    if (icomp.eq.jcomp) then
                      EnrTraceDispl(m1,m2) = EnrTraceDispl(m1,m2)  &
                                           - shapQ_f(k2)*shapVV_n(k1)*weight
                      ! write(*,*) 'elem    EnrTraceDispl integr:   m1,m2=',m1,m2
                      ! write(*,*) 'elem    EnrTraceDispl integr:   k1,k2=',k1,k2
                      ! write(*,*) 'elem    EnrTraceDispl integr:   jcomp,icomp=',jcomp,icomp
                      ! write(*,*) 'elem    EnrTraceDispl integr:   shapQ_f(k2),shapVV_n(k1)=',shapQ_f(k2),shapVV_n(k1)

                    endif
                  enddo
                enddo
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   C A U C H Y    S T R E S S   )
!
!   0
!
!
!  .........END OUTER LOOP
              enddo
            enddo
!
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows)
!
!  .......SECOND OUTER loop through enriched H1 test function
            do k1=1,nrdofHH
!  .........OUTER loop through components
              do jcomp=1,3
              ! counter of row
                m1 = 3*nrdofVV+(k1-1)*3+jcomp
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   0
!
!
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   - <\hat \sigma,v>
!
!  ...........INNER loop through H(div) bdry trial dofs
                do k3=1,nrdofQ_f
                  do icomp=1,3
                    m3 = (k3-1)*3+icomp + 3*nrdofQ_f*(jf-1)
                    if (icomp.eq.jcomp) then
                      EnrTraceStress(m1,m3) = EnrTraceStress(m1,m3)  &
                                            - shapQ_f(k3)*shapHH(k1)*weight
                    endif
                  enddo
                enddo
!
!    .........END OUTER LOOP
              enddo
            enddo
!    .....end of loop over integration points
          enddo
!   ....end sub triangle loop
        enddo
!
        deallocate(xverts_f,nverts_f,Rotverts_f)
!     end of face loop
      enddo
!
!
        enrdof = 3*nrdofVV+3*nrdofHH
        write(*,*) 'Gram       mdle = ',mdle
        do i=1,18
          write(*,6000) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 6000     format('i = ',i3,'  ',25e12.5)
        enddo
        if (mdle.eq.1) then
        gdump=75
        open(unit=gdump,file='output/gram', &
          form='formatted',access='sequential',status='unknown')
        do i=1,enrdof
          write(gdump,5999)Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 5999     format(75(e22.15,","))
        enddo
        close(gdump)
        bdump1=76
        open(unit=bdump1,file='output/b_u', &
          form='formatted',access='sequential',status='unknown')
        do i=1,enrdof
          write(bdump1,5998) EnrFieldDispl(i,1:3*nrdofQ)
 5998     format(3(e22.15,","))        
        enddo
        close(bdump1)
        bdump2=77
        open(unit=bdump2,file='output/b_sigma', &
          form='formatted',access='sequential',status='unknown')
        do i=1,enrdof
          write(bdump2,5997) EnrFieldStress(i,1:6*nrdofQ)
 5997     format(6(e22.15,","))
        enddo
        close(bdump2)
        bdump3=78
        open(unit=bdump3,file='output/b_omega', &
          form='formatted',access='sequential',status='unknown')
        do i=1,enrdof
          write(bdump3,5998) EnrFieldOmega(i,1:3*nrdofQ) 
        enddo
        close(bdump3)
        bdump4=79
        open(unit=bdump4,file='output/b_uhat', &
          form='formatted',access='sequential',status='unknown')
        do i=1,enrdof
          write(bdump4,5996) EnrTraceDispl(i,1:3*nrdofQ_f*Nrf)
 5996     format(27(e22.15,","))        
        enddo
        close(bdump4)
        bdump5=80
        open(unit=bdump5,file='output/b_sigmahat', &
          form='formatted',access='sequential',status='unknown')
        do i=1,enrdof
          write(bdump5,5996) EnrTraceStress(i,1:3*nrdofQ_f*Nrf)
        enddo
        close(bdump5)

      endif
        ! call pause











      if (iprint.eq.1) then
        write(*,*) 'Gram = '
        do i=1,25
          write(*,6000) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 ! 6000     format('i = ',i3,'  ',25e12.5)
        enddo

        call pause

        write(*,*) 'EnrFieldDispl = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrFieldDispl(i,1:3*nrdofQ)
 6001     format('i = ',i4,'  ',15(/,10e12.5))
        enddo

        call pause

        write(*,*) 'EnrFieldStress = '
        do i=1,3*nrdofVV
          write(*,6001) i,EnrFieldStress(i,1:6*nrdofQ)
        enddo

        call pause

        do i=1+3*nrdofVV,3*nrdofVV+3*nrdofHH
          write(*,6001) i,EnrFieldStress(i,1:6*nrdofQ)
        enddo

        call pause

        write(*,*) 'EnrFieldOmega = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrFieldOmega(i,1:3*nrdofQ)
        enddo

        call pause

        write(*,*) 'EnrTraceDispl = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrTraceDispl(i,1:3*nrdofQ_f*Nrf)
        enddo

        call pause

        write(*,*) 'EnrTraceStress = '
        do i=3*nrdofVV,3*nrdofVV+3*nrdofHH
          write(*,6001) i,EnrTraceStress(i,1:3*nrdofQ_f*Nrf)
        enddo

        call pause

      endif

!
!-----------------------------------------------------------------------------------
!     D P G    L O C A L    A S S E M B L Y                                        |
!-----------------------------------------------------------------------------------
!
!  ...Compact enriched number of rows (total enriched test dof)
      enrdof = 3*nrdofVV+3*nrdofHH
!
!  ...factor the Gram matrix
      uplo = 'U'


      ! allocate(Gramtmp(enrdof,enrdof))
 !      Gramtmp= ZERO
 !      call DTPTTR(uplo,enrdof,Gram,Gramtmp,enrdof,info)
 !      call DPOTRF(uplo,enrdof,Gramtmp,enrdof,info)

 !      gdump=75
 !        open(unit=gdump,file='output/gram_fact', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(gdump,5999)Gramtmp(i,1:75)
 ! ! 5990     format(75(e12.5,","))
 !        enddo
 !        close(gdump)

 !      call DTRTTP(uplo,enrdof,Gramtmp,enrdof,Gram,info)


      call DPPTRF(uplo, enrdof, Gram, info)

 !      gdump=75
 !        open(unit=gdump,file='output/gram', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(gdump,5999)Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 ! ! 5999     format(75(e12.5,","))
 !        enddo
 !        close(gdump)
      ! if (info.ne.0) then
        write(*,*) 'elem_POLYDPG, mdle,enrdof=',mdle,enrdof
        write(*,*) 'elem_POLYDPG Gram Cholesky Factorization: info = ',info
        ! stop
      ! endif
!
!  ...Create vector of indices with dof of each physics variable (in the same order of physics file)
      ndofphysics = (/3*nrdofQ_f*Nrf,3*nrdofQ_f*Nrf,3*nrdofQ,6*nrdofQ,3*nrdofQ/)
!
!  ...Construct EnrEverything by packing all Enriched Matrices in the order H1,Hcurl,Hdiv,L2
!     (in the same order of physics file) and load
      EnrEverything=ZERO
!  ...EnrTraceDispl (Face L2)
      kmin=0
      kmax=kmin+ndofphysics(1)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrTraceDispl(1:enrdof,1:(kmax-kmin))
!  ...EnrTraceStress (Face L2)
      kmin=kmax
      kmax=kmin+ndofphysics(2)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrTraceStress(1:enrdof,1:(kmax-kmin))
!  ...EnrFieldDispl (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(3)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrFieldDispl(1:enrdof,1:(kmax-kmin))
!  ...EnrFieldStress (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(4)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrFieldStress(1:enrdof,1:(kmax-kmin))
!  ...EnrFieldOmega (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(5)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrFieldOmega(1:enrdof,1:(kmax-kmin))
!  ...EnrLoad
      kmin=kmax
      kmax=kmin+NR_RHS
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrLoad(1:enrdof,1:(kmax-kmin))
!
!  ...Save copy of EnrStiffness which implies deleting the load part from EnrEverything
      EnrStiffness=ZERO
      EnrStiffness(1:enrdof,1:kmin)=EnrEverything(1:enrdof,1:kmin)
!
      bdump7=82
      open(unit=bdump7,file='output/EnrEverything', &
          form='formatted',access='sequential',status='unknown')
        do i=1,75
          write(bdump7,5995) EnrEverything(i,1:67)
 ! 5995     format(66(e12.5,","))        
        enddo
      close(bdump7)
! !  ...save copies of enriched stiffness matrices
!       EnrTraceDisplc=EnrTraceDispl; EnrTraceStressc=EnrTraceStress
!       EnrFieldDisplc=EnrFieldDispl; EnrFieldStressc=EnrFieldStress
!       EnrFieldOmegac=EnrFieldOmega
!
!  ...G^-1 * EnrEverything
      call DPPTRS(uplo,enrdof,kmax,Gram,EnrEverything(:,1:kmax),3*MAXtetraVV+3*MAXtetraHH,info1)
      if (info1.ne.0) then
        write(*,*) 'elem_DPG_UWEAK: info1 = ',info1 ; stop
      endif
!
!  ...Build full DPG matrix (stiffness + load) in one go
      FullDPG = ZERO
      transa = 'T'
      transb = 'N'
      m = kmin
      n = kmax
      k = enrdof
      l = 3*MAXtetraVV+3*MAXtetraHH
      call DGEMM(transa,transb,m,n,k,1.d0,EnrStiffness,l,EnrEverything,l,0.d0,FullDPG,6*MAXtriaQ*MAX_NRF+12*MAXtetraQ)
!
!     ULTIMATE DPG LOAD VECTORS AND STIFFNESS MATRICES THROUGH STATIC CONDENSATION
!
!  ...Populate the assembly matrices accordingly
      n1=0
      do i=1,NR_PHYSA
        n=ndofphysics(i)
        n2=n1+n
        m1=0
        do j=1,NR_PHYSA
          m=ndofphysics(j)
          m2=m1+m
!  .......First initialize
          ALOC(i,j)%array = ZERO
!  .......Then populate
          ALOC(i,j)%array(1:n,1:m) = FullDPG(n1+1:n2,m1+1:m2)
          m1=m2
        enddo
        m2=m1+NR_RHS
!  .....First initialize
        BLOC(i)%array = ZERO
!  .....Then populate
        BLOC(i)%array(1:n,1:NR_RHS) = FullDPG(n1+1:n2,m1+1:m2)
        n1=n2
      enddo
!
      bdump6=81
      open(unit=bdump6,file='output/FullDPG', &
          form='formatted',access='sequential',status='unknown')
        do i=1,66
          write(bdump6,5995) FullDPG(i,1:66)
 5995     format(66(e22.15,","))        
        enddo
      close(bdump6)
      
      deallocate(Nfaces,Norientf,Nedges,Nverts,Xvertl)
!
end subroutine elem_DPG_UWEAK_poly