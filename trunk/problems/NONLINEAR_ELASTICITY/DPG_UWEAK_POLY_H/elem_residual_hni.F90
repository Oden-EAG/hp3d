!--------------------------------------------------------------------------
!> Purpose : return stiffness matrix and Residual vector for element
!!
!! @param[in]  Mdle      - an element (middle node) number
!! @param[out] Resid     - element residual (squared)
!! @param[out] Nref_flag - suggested h-refinement flag
!--------------------------------------------------------------------------
!
      subroutine elem_residual(Mdle, Resid)
!
      use uweak_module_poly
      use connectivity_poly
      use control, only : INTEGRATION
      use parametersDPG
      use element_data_poly
      use data_structure3D_poly
      use isotropic_elast_material_compo
      use common_prob_data, only: SYMMETRY_TOL, TEST_NORM, IP
      use geometry_polydpg,only: ALLC_DOM
!------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in)  :: Mdle
      real*8,  intent(out) :: Resid
      ! integer, intent(out) :: Nref_flag
!------------------------------------------------------------------------------------------
!  ...MATRICES
!     stiffnes and load matrices for the enriched test space
  ! real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),(3*MAXtriaH*MAX_NRFC)) :: EnrTraceDispl!,EnrTraceDisplc
  ! real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),3*MAXtriaQ*MAX_NRFC) :: EnrTraceStress!,EnrTraceStressc
  ! real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),3*MAXtetraQ) :: EnrFieldDispl!,EnrFieldDisplc
  ! real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),6*MAXtetraQ) :: EnrFieldStress!,EnrFieldStressc
  ! real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),3*MAXtetraQ) :: EnrFieldOmega!,EnrFieldOmegac
  ! real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),MAXNRHS_MOD) :: EnrLoad
  ! real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),                                               &
  !                   3*(MAXtriaQ+MAXtriaH)*MAX_NRFC+12*MAXtetraQ) :: EnrStiffness
  ! real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),                                               &
  !                   3*(MAXtriaQ+MAXtriaH)*MAX_NRFC+12*MAXtetraQ+MAXNRHS_MOD) :: EnrEverything
  ! real*8, dimension((3*(MAXtriaQ+MAXtriaH)*MAX_NRFC+12*MAXtetraQ),                             &
  !                   3*(MAXtriaQ+MAXtriaH)*MAX_NRFC+12*MAXtetraQ+MAXNRHS_MOD) :: FullDPG
!     Gram matrix for the local Riesz matrix in LAPACK format
  real*8, dimension((3*MAXtetraVV+3*MAXtetraHH)*(3*MAXtetraVV+3*MAXtetraHH+1)/2) :: Gram

  real*8, dimension(3*MAXtetraVV+3*MAXtetraHH) :: diag_pc

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
      integer, dimension(4)  :: norie
!
!  ...SHAPE FUNCTIONS
!     Discont Face H1 (trial)      
      real*8, dimension(MAXtriaH   )  :: shapH_f
      real*8, dimension(2,MAXtriaH )  :: gradH_f
!     Face L2 (trial)      
      real*8, dimension(MAXtriaQ   )  :: shapQ_f
!     L2  (trial)
      real*8, dimension(  MAXtetraQ)  :: shapQ
      integer                         :: nrdofQ
!     H1  (geometry and trial)
      real*8, dimension(  MAXtetraH)  :: shapH
      real*8, dimension(3,MAXtetraH)  :: gradH
      integer                         :: nrdofH
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
      integer, dimension(MAXtetraHH) ::degH
      integer, dimension(MAXtetraVV) ::degV
      integer, dimension(MAXtetraQ ) ::degQ
!
!  ...geometry
      real*8, dimension(3,MAXtetraH) :: xnod
      real*8, dimension(3)           :: xi,x,rn,fn,x_aux,b0,b1,b2,a0,a1,a2,a3,fn_a
      real*8, dimension(3,3)         :: dxdxi_e,dxidx_e,Qrot_e,Qrot_f,maptet,maptetinv,b_aff
      real*8, dimension(2,2)         :: dxdxi_f,dxidx_f,dxdxi_bf,maptri
      real*8, dimension(2)           :: t
      real*8, dimension(3,2)         :: dxidt,dxdt,rt,dxi_eta
      real*8, dimension(3,4)         :: a_aff
      integer                        :: nsign
!
!  ...Resid vector for the enriched space
      real*8, dimension(3*MAXtetraVV+3*MAXtetraHH) :: EnrResid,EnrResidc
!
!  ...tensors in physical coordinates
      real*8, dimension(3,3,3,3) :: A,AA
!
!  ...temporary variables
      real*8                 :: tmp
!
!  ...source term (don't need Neumann term)
      real*8, dimension(3,MAXNRHS) :: fval
!
!  ...3D quadrature data
      real*8, allocatable :: xiloc(:,:),xloc(:,:)
      real*8, allocatable :: wxi(:)
!
!  ...2D quadrature data for boundary terms
      real*8, allocatable :: tloc(:,:)
      real*8, allocatable :: wt(:)
!
!  ...approximate solution
      real*8, dimension(MAXEQNH,MAXtetraH) :: dofH
      real*8, dimension(MAXEQNE,MAXtetraE) :: dofE
      real*8, dimension(MAXEQNV,MAXtetraV) :: dofV
      real*8, dimension(MAXEQNQ,MAXtetraQ) :: dofQ
      real*8, dimension(  MAXEQNH  )       :: solH
      real*8, dimension(  MAXEQNH,3)       :: dsolH
      real*8, dimension(3,MAXEQNE  )       :: solE
      real*8, dimension(3,MAXEQNE  )       :: curlE
      real*8, dimension(3,MAXEQNV  )       :: solV
      real*8, dimension(  MAXEQNV  )       :: divV
      real*8, dimension(  MAXEQNQ  )       :: solQ
!     displacement
      real*8, dimension(3)                 :: u,uhat
!     stress
      real*8, dimension(3,3)               :: sigma, omega
      real*8, dimension(3)                 :: sigma_n
!
!  ...miscellaneous
      integer :: i,j,k,m,n,k1,k2,m1,m2,ipt,icomp,jcomp,ijcomp,kcomp,lcomp,klcomp,    &
                 nint,ifc,iprint,iload,iflag,info,info1,info2,info3,nordtmp, enrdof,       &
                 Nrv_f,mdlf,jv,jf,loc,nrdofQ_f,nrdofH_f, idom
      real*8  :: weight,wa,rjac,brjac,diff,alpha

      real*8  :: diffmax,dmax,rjac_e,rjac_f,rjac_bf,r_e, &
                 maptetdet,maptridet,Area_f,r_f, dir,b_area,a_vol
      real*8, dimension(9) :: faceint
      integer, allocatable :: nfaces(:),norientf(:),nedges(:),nverts(:),nverts_f(:)
      real*8, allocatable :: xverts(:,:),xverts_f(:,:),Gramtmp(:,:)

!     For solution evaluation
      real*8, dimension(MAXEQNQ,MAXtetraQ):: zdofQ
      real*8, dimension(MAXEQNU,MAXtriaH) :: zdofU
      real*8, dimension(MAXEQNF,MAXtriaQ) :: zdofF
      real*8, dimension(  MAXEQNQ  )      :: zsolQ
      real*8, dimension(  MAXEQNU  )      :: zsolU
      real*8, dimension(  MAXEQNF  )      :: zsolF
      real*8                              :: DDOT
      integer :: ndofH,ndofE,ndofV,ndofU,ndofF,ndofQ, ivar, nord_mdlf, &
                 nord_add_local,nord_add_ini , gdump, rdump , imat
!  ...LAPACK stuff
      character uplo,transa,transb
!
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1

!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
!  ...Load number for which to calculate residual from
      iload = 1
      alpha = 1d-4
!
      iprint = 0

      etype = 'tetr'
      ftype = 'tria'
!     Retrieve solved MDLE node degrees of freedom
      call find_ndof_poly(Mdle, ndofH,ndofE,ndofV,ndofQ)
! 
      zdofQ = 0.d0
      zdofQ(1:NRQVAR,1:ndofQ) = NODES(mdle)%zdofQ(1:NRQVAR,1:ndofQ)   
!     retrieve element node data
      call elem_nodes_poly(Mdle,Nfaces,Norientf,Nrf,Nedges,Nre,Nverts,Nrv)
! !     allocate array for vertices 
      allocate(xverts(3,Nrv))
!     get vertices coordinates
      do jv = 1, nrv
        xverts(1:3,jv) = NODES(nverts(jv))%coord(1:3,1)
      enddo
!
!     get element affine coordinates points a0,a1,a2,a3
      call poly_affine_3d(xverts,nrv,a0,a1,a2,a3)
!     construct array a_aff with a0,a1,a2,a3 as columns
      a_aff(:,1) = a0; a_aff(:,2) = a1; a_aff(:,3) = a2; a_aff(:,4) = a3
! ! computational of jacobian determinant
! ! get normal of face a1a2a3
!       call face_normal(a_aff(:,2:4),fn_a)
! ! get tetrahedron volume
!       call tetra_centr_vol(a0,a1,a2,a3,fn_a,x_aux,a_vol)
! ! correct
!       rjac_e = 6.d0*a_vol
! construct linear transformation associated to affine map. Interpreted as classic jacobian
      dxdxi_e(:,1) = a1 - a0
      dxdxi_e(:,2) = a2 - a0
      dxdxi_e(:,3) = a3 - a0
! get inverse and determinant of linear transformation
      call geom(dxdxi_e,dxidx_e,rjac_e,iflag)
!
      if (iflag.ne.0) then
        write(*,8999) Mdle,rjac
 8999   format('Negative Jacobian!Mdle,rjac=',i8,2x,e12.5)
        stop
      endif
!     get order, in the norder structure of a uniform tetrahedral element
      norder = 0
      norder(1:15) = NODES(Mdle)%order      
!
!
        idom = ALLC_DOM(mdle)
        select case(idom)
        case(1)
          imat = 1
        case default
          imat = 2
        end select


!     set nord_add_local
      nord_add_ini = NORD_ADD
      if (VARIABLE_DP.eq.1) then
        call correct_dp(Mdle,nord_add_ini,NODES(Mdle)%order,nrf,nord_add_local)
      ! write(*,*) 'mdle, nord_add_local = ',Mdle,nord_add_local
      else
        nord_add_local = NORD_ADD
      endif
! !     dp higher in residual computation
!       nord_add_local = nord_add_local  + 1
!
!  ...set the enriched order of appoximation
      nordP = min(NODES(Mdle)%order + nord_add_local , MAXPP)
! 
!  ...set up the element quadrature
      INTEGRATION = max(0, nord_add_local)
      allocate(xiloc(3,(Nrf-3)*(MAX_NRVF-2)*MAXNINT2ADD),wxi((Nrf-3)*(MAX_NRVF-2)*MAXNINT2ADD))
      call set_quadrature_elem_HNI_DPG(mdle,Nrf,NODES(Mdle)%order,a_aff,nint,xiloc,wxi)
      INTEGRATION = 0
      allocate(xloc(3,nint))
      xloc = 0.d0
      call poly_affine_3d_to_x(xiloc(:,1:nint),nint,a_aff,xloc)
      if (iprint.eq.1) then
        write(*,7020) xnod(1,1:8),xnod(2,1:8),xnod(3,1:8)
 7020   format('elem_residual: xnod  = ',8(f8.3,2x),  &
         2(  /,'                       ',8(f8.3,2x)))
        write(*,7030) dofH(1,1:8),dofV(1,1:6),dofQ(1:2,1)
 7030   format('elem_residual: dofH(1,1:8) = ',8(e12.5,2x),  &
             /,'               dofV(1,1:6) = ',6(e12.5,2x),  &
             /,'               dofQ(1:2,1) = ',2(e12.5,2x))
      endif
!
!  ...initialize the enriched local element residual vector
      EnrResid = ZERO
! !  ...initialize the Gram matrix
      Gram = ZERO
!
!-----------------------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L  ( O N L Y    G R A M    B Y    H N I )     |
!-----------------------------------------------------------------------------------
!
      do ipt=1,nint
!              
        xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt); x(1:3) = xloc(1:3,ipt)
!       
!        .....Compute shape functions needed for test/trial field variables and geometry
!             L2 (field trial)
        call shap3Q_mono(xi,NODES(Mdle)%order, nrdofQ,shapQ,degQ)              
!             H1 (test)
        call shap3HH_mono(xi,nordP, nrdofHH,shapHH,gradHH,degH)              
!             H(div) (test)
        call shap3VV_mono(xi,nordP, nrdofVV,shapVV,divVV,degV)
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
!        .....integration weight
        weight = wa *rjac_e
!
!       Evaluate the computed approximate solution
        zsolQ(1:MAXEQNQ)=ZERO
        do k=1,nrdofQ
          do ivar=1,MAXEQNQ
            zsolQ(ivar) = zsolQ(ivar) + zdofQ(ivar,k)*shapQ(k)
          enddo
        enddo
!       store with physically meaning variable names
        u(1:3)        = zsolQ(1:3)
        ! write(*,*) 'u=',u(1),u(2),u(3)
        ! 
        sigma(1,1:3)  = (/  zsolQ(4)  ,  zsolQ(7)  ,  zsolQ(8)  /)
        sigma(2,1:3)  = (/  zsolQ(7)  ,  zsolQ(5)  ,  zsolQ(9)  /)
        sigma(3,1:3)  = (/  zsolQ(8)  ,  zsolQ(9)  ,  zsolQ(6)  /)
        ! 
        omega(1,1:3)  = (/      0.d0  ,  zsolQ(12) , -zsolQ(11) /)
        omega(2,1:3)  = (/ -zsolQ(12) ,      0.d0  ,  zsolQ(10) /)
        omega(3,1:3)  = (/  zsolQ(11) , -zsolQ(10) ,      0.d0  /)
        ! write(*,*) 'sigma='
        ! write(*,*) sigma(1,:)
        ! write(*,*) sigma(2,:)
        ! write(*,*) sigma(3,:)
        ! write(*,*) 'omega='
        ! write(*,*) omega(1,:)
        ! write(*,*) omega(2,:)
        ! write(*,*) omega(3,:)
        ! call find_material(x,(/0.d0,0.d0,0.d0/),imat)
        idom = ALLC_DOM(mdle)
        select case(idom)
        case(1)
          imat = 1
        case default
          imat = 2
        end select
! 
!  .....compute the compliance tensor
        call getA(imat,x, A)
!
!  .....need this for the adjoint graph norm
        call getAA(imat,X, AA)
!
!  .....get the source term
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
!
!           G R A M   M A T R I X
!
            select case(TEST_NORM)
!
!   ADJOINT GRAPH
!   (A:\tau_2,A:\tau)+(div(\tau_2),div(tau))+(\tau_2,\tau)+alpha*(\tau_2,\tau)
!
            case(1)
              do m2=m1,3*nrdofVV
                icomp = mod(m2-1,3)+1
                k = nk(m1,m2)
                k2 = int((m2-1)/3)+1
                do kcomp=1,3
                  do lcomp=1,3
                    if ((icomp.eq.jcomp.and.kcomp.eq.lcomp).or. &
                        (icomp.eq.kcomp.and.jcomp.eq.lcomp).or. &
                        (icomp.eq.lcomp.and.jcomp.eq.kcomp))    &
                      Gram(k) = Gram(k)  &
                              + ( AA(icomp,kcomp,jcomp,lcomp)                    &
                                 *shapVV(lcomp,k1)                               &
                                 *shapVV(kcomp,k2)/real(3+degV(k1)+degV(k2),8) )   &
                                *weight
                  enddo
                enddo
                if (icomp.eq.jcomp) then
                  Gram(k) = Gram(k)  + divVV(k1)*divVV(k2)*weight
                  Gram(k) = Gram(k)  &
                          + (1.d0+alpha)*(                          &
                            + shapVV(1,k1)*shapVV(1,k2)             &
                            + shapVV(2,k1)*shapVV(2,k2)             &
                            + shapVV(3,k1)*shapVV(3,k2) )           &
                            / real(3+degV(k1)+degV(k2),8)  *  weight
                endif
              enddo
!
!   (A:\tau_2, \grad v)
!
              do m2=3*nrdofVV+1,3*nrdofVV+3*nrdofHH
                icomp = mod(m2-1,3)+1
                k = nk(m1,m2)
                k2 = int((m2-3*nrdofVV-1)/3)+1
                do kcomp=1,3
                  do lcomp=1,3
                      if ((icomp.eq.jcomp.and.kcomp.eq.lcomp).or. &
                        (icomp.eq.kcomp.and.jcomp.eq.lcomp).or. &
                        (icomp.eq.lcomp.and.jcomp.eq.kcomp))    &
                        Gram(k) = Gram(k)  &
                                + ( A(icomp,kcomp,jcomp,lcomp)                 &
                                   *shapVV(lcomp,k1)*gradHH(kcomp,k2) )        &
                                    / real(3+degV(k1)+degH(k2)-1,8)         *weight
                  enddo
                enddo
              enddo
!
!   MATHEMATICIAN'S
!   (\tau_2,\tau)+(div(\tau_2),div(tau))
!
            case(2)
              do m2=m1,3*nrdofVV
                icomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( divVV(k1)*divVV(k2)                     &
                              /real(3+degV(k1)-1+degV(k2)-1,8)        &
                            +(shapVV(1,k1)*shapVV(1,k2)               &
                            + shapVV(2,k1)*shapVV(2,k2)               &
                            + shapVV(3,k1)*shapVV(3,k2))              &
                            /real(3+degV(k1)+degV(k2),8)      )*weight
                endif
              enddo
            end select
!  .....END OUTER LOOP through test stresses
          enddo
        enddo
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows)
!
!  .....SECOND OUTER loop through enriched H1 dofs
        do k1=1,nrdofHH
!  .......OUTER loop through components
          do icomp=1,3
            ! counter for part 2
            m1 = 3*nrdofVV+(k1-1)*3+icomp
! !
!            G R A M   M A T R I X
!
            select case(TEST_NORM)
!
!   ADJOINT GRAPH
!   (\grad(v_1),A:\tau_2+\grad(v_2))+alpha*(v_1,v_2)
!
            case(1)

              do m2=m1,3*nrdofVV+3*nrdofHH
                jcomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-3*nrdofVV-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( alpha*shapHH(k1)*shapHH(k2)           &
                             /real(3+degH(k1)+degH(k2),8)           &
                            +(gradHH(1,k1)*gradHH(1,k2)             &
                            + gradHH(2,k1)*gradHH(2,k2)             &
                            + gradHH(3,k1)*gradHH(3,k2))            &
                            /real(3+degH(k1)-1+degH(k2)-1,8))*weight
                endif
              enddo
!
!   MATHEMATICIAN'S
!   (v_2,v) + (grad(v_2),grad(v))
!
            case(2)
              do m2=m1,3*nrdofVV+3*nrdofHH
                jcomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-3*nrdofVV-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( shapHH(k1)*shapHH(k2)                 &
                             /real(3+degH(k1)+degH(k2),8)           &
                            +(gradHH(1,k1)*gradHH(1,k2)             &
                            + gradHH(2,k1)*gradHH(2,k2)             &
                            + gradHH(3,k1)*gradHH(3,k2))            &
                            /real(3+degH(k1)-1+degH(k2)-1,8))*weight
                endif
              enddo
            end select
!
!  ......END OUTER LOOP
          enddo
        enddo
!
!  ...end of loop through integration points
      enddo
!
!-----------------------------------------------------------------------------------
!    R E S I D U A L    V E C T O R    ( B Y   S U B T E S S E L L A T I O N )     |
!-----------------------------------------------------------------------------------
!
!
      deallocate(xiloc,wxi,xloc,xverts)
!  ...set up the element quadrature
      INTEGRATION = max(0,(nord_add_local+1)/2)
      allocate(xiloc(3,(Nrf-3)*(MAX_NRVF-2)*MAXNINT3ADD),wxi((Nrf-3)*(MAX_NRVF-2)*MAXNINT3ADD))
      call set_quadrature_affine_elem_mod(mdle,Nrf,NODES(Mdle)%order,a_aff,nint,xiloc,wxi)
      INTEGRATION = 0
      allocate(xloc(3,nint))
      xloc = 0.d0
      call poly_affine_3d_to_x(xiloc(:,1:nint),nint,a_aff,xloc)
!
      do ipt=1,nint
!              
        xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt); x(1:3) = xloc(1:3,ipt)
!       
!        .....Compute shape functions needed for test/trial field variables and geometry
!             L2 (field trial)
        call shap3Q_mono(xi,NODES(Mdle)%order, nrdofQ,shapQ,degQ)              
!             H1 (test)
        call shap3HH_mono(xi,nordP, nrdofHH,shapHH,gradHH,degH)              
!             H(div) (test)
        call shap3VV_mono(xi,nordP, nrdofVV,shapVV,divVV,degV)
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
!        .....integration weight
        weight = wa !*rjac_e
!
!       Evaluate the computed approximate solution
        zsolQ(1:MAXEQNQ)=ZERO
        do k=1,nrdofQ
          do ivar=1,MAXEQNQ
            zsolQ(ivar) = zsolQ(ivar) + zdofQ(ivar,k)*shapQ(k)
          enddo
        enddo
!       store with physically meaning variable names
        u(1:3)        = zsolQ(1:3)
        ! write(*,*) 'u=',u(1),u(2),u(3)
        ! 
        sigma(1,1:3)  = (/  zsolQ(4)  ,  zsolQ(7)  ,  zsolQ(8)  /)
        sigma(2,1:3)  = (/  zsolQ(7)  ,  zsolQ(5)  ,  zsolQ(9)  /)
        sigma(3,1:3)  = (/  zsolQ(8)  ,  zsolQ(9)  ,  zsolQ(6)  /)
        ! 
        omega(1,1:3)  = (/      0.d0  ,  zsolQ(12) , -zsolQ(11) /)
        omega(2,1:3)  = (/ -zsolQ(12) ,      0.d0  ,  zsolQ(10) /)
        omega(3,1:3)  = (/  zsolQ(11) , -zsolQ(10) ,      0.d0  /)
        ! write(*,*) 'sigma='
        ! write(*,*) sigma(1,:)
        ! write(*,*) sigma(2,:)
        ! write(*,*) sigma(3,:)
        ! write(*,*) 'omega='
        ! write(*,*) omega(1,:)
        ! write(*,*) omega(2,:)
        ! write(*,*) omega(3,:)
        ! call find_material(x,(/0.d0,0.d0,0.d0/),imat)
        idom = ALLC_DOM(mdle)
        select case(idom)
        case(1)
          imat = 1
        case default
          imat = 2
        end select
! 
!  .....compute the compliance tensor
        call getA(imat,x, A)
!
!  .....need this for the adjoint graph norm
        call getAA(imat,X, AA)
!
!  .....get the source term
        call getf(Mdle,x, fval)
!
        ! u= 0.d0
        ! sigma = 0.d0
        ! omega = 0.d0
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
!           E N R I C H E D   R E S I D U A L   V E C T O R
!
!   ( A:\sigma + \omega , \tau ) + ( u , div(\tau) ) +( \omega, \tau )
!
            tmp=0.d0
            do lcomp=1,3
              do kcomp=1,3
                do icomp=1,3
                tmp = tmp  &
                    + A(icomp,kcomp,lcomp,jcomp)*sigma(icomp,kcomp)*shapVV(lcomp,k1)
                enddo
              enddo
            enddo
! 
            tmp = tmp + u(jcomp)*divVV(k1)            
            do lcomp=1,3
              ! if (jcomp.eq.lcomp) then
              !   tmp = tmp + u(lcomp)*divVV(k1)
              ! endif
              tmp = tmp + omega(jcomp,lcomp)*shapVV(lcomp,k1)
            enddo
            EnrResid(m1) = EnrResid(m1) + tmp*weight
!  .....END OUTER LOOP through test stresses
          enddo
        enddo
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows)
!
!  .....SECOND OUTER loop through enriched H1 dofs
        do k1=1,nrdofHH
!  .......OUTER loop through components
          do icomp=1,3
            ! counter for part 2
            m1 = 3*nrdofVV+(k1-1)*3+icomp
!
!           E N R I C H E D   R E S I D U A L   V E C T O R
!
!   + (\sigma, grad(v) ) - (f,v)
!
            tmp=0.d0
            do kcomp=1,3
              tmp = tmp + sigma(icomp,kcomp)*gradHH(kcomp,k1)
            enddo
            tmp = tmp - fval(icomp,1)*shapHH(k1)
            EnrResid(m1) = EnrResid(m1) + tmp*weight
!
!  ......END OUTER LOOP
          enddo
        enddo
!
!  ...end of loop through integration points
      enddo
!-----------------------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L                                           |
!-----------------------------------------------------------------------------------
!
!  ...loop through element faces
      do jf=1,nrf
        mdlf = Nfaces(jf)

!       Retrieve computed MDLF node degrees of freedom
        nord_mdlf= NODES(mdlf)%order
        ndofU = (nord_mdlf+1)*(nord_mdlf+2)/2  
        ndofF = (nord_mdlf  )*(nord_mdlf+1)/2  
        zdofF = 0.d0
        zdofF(1:NRFVAR,1:ndofF) = NODES(mdlf)%zdofF(1:NRFVAR,1:ndofF)
        zdofU = 0.d0
        zdofU(1:NRUVAR,1:ndofU) = NODES(mdlf)%zdofU(1:NRUVAR,1:ndofU)

        call face_vert_list(mdlf,nrv_f,nverts_f)
        allocate(xverts_f(3,nrv_f))
        do jv=1,Nrv_f
          xverts_f(1:3,jv) = NODES(nverts_f(jv))%coord(1:3,1)
        enddo
!
! get face affine coordinates points b0,b1,b2
        call poly_affine_2d(xverts_f,nrv_f,b0,b1,b2)
        b_aff(:,1) = b0
        b_aff(:,2) = b1
        b_aff(:,3) = b2
! 
!       set up face order in usual structure
        nordf = 0
        nordf(1:4) = NODES(Mdlf)%order
        norie = 0
        if (nrv_f.eq.3) then
          norie(1:3) = NODES(Mdlf)%face_orient(1:3)
          norie(3) = 1 - norie(3)
        ! else
        !   write(*,*) 'elem_poly_hni: face is not a triangle! mdle,jf=',mdle,jf
        !   stop
        endif
!  .....set up the face quadrature
        allocate(tloc(2,(nrv_f-2)*MAXNINT2ADD),wt((nrv_f-2)*MAXNINT2ADD))
! get quadrature points for polygon using the affine coordinates
        INTEGRATION = max(0, (nord_add_local+1)/2)
        call set_quadrature_affine_face(mdlf,NODES(Mdlf)%order,nrv_f,xverts_f,b_aff,nint,tloc,wt)
        INTEGRATION = 0
! transform integration points' face affine coordinates to element affine coordinates
        xiloc = 0.d0
        dxi_eta = 0.d0
        call poly_affine_2d_to_3d(tloc,nint,b_aff,a_aff,xiloc,dxi_eta)
! get area of triangle b0b1b2
        call trian_centr_area(b0,b1,b2,x_aux,b_area)
        rjac_f = 2.d0*b_area
! 
        ! write(*,*) 'elem_poly:   mdlf,rjac_f=',mdlf,rjac_f
!       get face normal in physical coordinates
        call face_normal(xverts_f(:,1:3),fn)
! 
!       check if face normal locally goes outward (Norientf = 0). Also correct sign of integral (dir)
        dir = 1.d0
        if(Norientf(jf).ne.0) then
!       if normal fn goes inward, correct it
          ! fn = -1.d0*fn
          dir = -1.d0
        endif
!
!  .....loop through face integration points
        do ipt=1,nint
!    .......quadrature point in face coordinates
          t(1:2) = tloc(1:2,ipt); wa = wt(ipt); xi(1:3) = xiloc(1:3,ipt)
!           evaluate trial (2d face) and test (3d element) functions
!         Discont Face H1
          call shape2DH(ftype,t,nordf, norie,nrdofH_f,shapH_f,gradH_f)
!         Face L2
          call shape2DQ(ftype,t,nordf, nrdofQ_f,ShapQ_f)
!           H1 (test)
          call shap3HH_mono(xi,nordP, nrdofHH,shapHH,gradHH,degH)
!           H(div) (test)
          call shap3VV_mono(xi,nordP, nrdofVV,shapVV,divVV,degV)
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
          weight = wa !*rjac_f
!
!     evaluate the approximate solution for \hat u
          zsolU(1:MAXEQNU)=ZERO
          do k=1,nrdofH_f
            do ivar=1,MAXEQNU
              zsolU(ivar) = zsolU(ivar) + zdofU(ivar,k)*shapH_f(k)
            enddo
          enddo
!     evaluate the approximate solution for \hat sigma_n
          zsolF(1:MAXEQNF)=ZERO
          do k=1,nrdofQ_f
            do ivar=1,MAXEQNF
              zsolF(ivar) = zsolF(ivar) + zdofF(ivar,k)*shapQ_f(k)
            enddo
          enddo
!      relabel solution
          uhat(1:3)    = zsolU(1:3)
          sigma_n(1:3) = zsolF(1:3)
          ! write(*,*) 'u=',uhat
          ! write(*,*) 's=',sigma_n

          ! uhat = 0.d0
          ! sigma_n = 0.d0
!
!
!    P A R T  1 : go through \tau\in H(div)^3 test space (this fills the first set of rows)
!
!
 ! .......OUTER loop through enriched test functions
          do k1=1,nrdofVV
            do kcomp=1,3
              m1 = (k1-1)*3+kcomp
!
!
!            E N R I C H E D   L O A D   V E C T O R
!
!   - <\hat u,(\tau n)>
!
              EnrResid(m1) = EnrResid(m1)  &
                           - uhat(kcomp)*shapVV_n(k1)*weight*dir
!
!  .......OUTER loop
            enddo
          enddo
! !
! !    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows)
! !
! !  .......SECOND OUTER loop through enriched H1 test function
          do k1=1,nrdofHH
            do kcomp=1,3
              ! counter of row
              m1 = 3*nrdofVV+(k1-1)*3+kcomp
!
!
!            E N R I C H E D   T R A C E   V E C T O R
!
!   - <sigma_n,v>
!
              EnrResid(m1) = EnrResid(m1)  &
                           - sigma_n(kcomp)*shapHH(k1)*weight*dir
!
!  .......OUTER loop
            enddo
          enddo
!
!  .....end of loop over integration points
        enddo
!
        deallocate(xverts_f,nverts_f,tloc,wt)
!  ...end of loop over faces
      enddo
!
!-----------------------------------------------------------------------------------
!     D P G    L O C A L    A S S E M B L Y                                        |
!-----------------------------------------------------------------------------------
!
 !      enrdof = 3*nrdofVV+3*nrdofHH
 ! !      if (mdle.eq.6) then
 !        gdump=75
 !        open(unit=gdump,file='output/gram1_hni', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(gdump,5999) Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 ! 5999     format(528(e22.15,","))
 !        enddo
 !        close(gdump)
 !        rdump=74
 !        open(unit=rdump,file='output/res1_hni', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(rdump,5998) Enrresid(i)
 ! 5998     format(e22.15,",")
 !        enddo
 !        close(rdump)
 !      endif

      enrdof = 3*nrdofVV+3*nrdofHH
!  ...factor the Gram matrix
!  ... with diagonal preconditioning
      diag_pc = ZERO
      do i=1,enrdof
        j = nk(i,i)
        diag_pc(i) = 1.d0/sqrt(Gram(j))
      enddo

      do i=1,enrdof
        do j=i,enrdof
          k = nk(i,j)
          Gram(k) = diag_pc(i)*Gram(k)*diag_pc(j)
        enddo
      enddo

      uplo = 'U'
      call DPPTRF(uplo, enrdof, Gram, info)
      if (info.ne.0) then
        write(*,*) 'elem_residual: info = ',info ; stop
      endif
!
!  ...save copies of the RHS to compute later the residual
      EnrResidc = EnrResid

      do i=1,enrdof
        EnrResid(i) = diag_pc(i)*EnrResid(i)
      enddo      
!
!  ...compute the product of inverted test Gram matrix with RHS,
!     EnrResid is overwritten with the solution
      call DPPTRS(uplo, enrdof, 1, Gram, EnrResid, 3*MAXtetraVV+3*MAXtetraHH, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem_residual: info1 = ',info1
        stop
      endif

      do i=1,enrdof
        EnrResid(i) = diag_pc(i)*EnrResid(i)
      enddo   

!
!  ...compute the residual
      Resid = DDOT(enrdof,EnrResidc,1,EnrResid,1)
!
!-----------------------------------------------------------------------------------
!             R E F I N E M E N T  F L A G S                                       |
!-----------------------------------------------------------------------------------
!  ...determines what refinement flag to use
!  ...if isotropic h refinements
      ! call get_isoref(Mdle, Nref_flag)
      ! Nref_flag = 110
!  ...if anisotropic h refinements -> missing
!
      deallocate(Nfaces,Norientf,Nedges,Nverts)
      deallocate(xiloc,wxi)
!
      end subroutine elem_residual

!--------------------------------------------------------------------------
!> Purpose : returns global residual
!!
!--------------------------------------------------------------------------
!
      subroutine compute_residual
!
      use data_structure3D_poly
      use environment, only : QUIET_MODE
!------------------------------------------------------------------------------------------
      implicit none
!
!  ...middle node
      integer :: mdle
!
!  ...residual
      real*8 :: residual

      real*8,allocatable  :: resid_sq(:)
!
!  ...rate
      real*8 :: rate
!
!  ...number of dof for higher order node
      integer :: ndofH,ndofE,ndofV,ndofQ
      integer :: nrdof_total
!
!  ...refinement flag
      ! integer :: nref_flag
!
!  ...field variables flag
      integer, dimension(NR_PHYSA) :: nflag
!
!  ...visitation flag
      integer, save :: ivis = 0
!
!  ...total number of dof for the old mesh
      integer, save :: nrdof_total_old
!
!  ...residual for the old mesh
      real*8,  save :: residual_old
!
!  ...residuals and rates to display
      real*8 , dimension(2,10), save :: rwork
!
!  ...number of dof to display
      integer , dimension(10), save :: iwork
!
!  ...miscellaneous
      integer :: i,iel,iphys,iprint, ref_dump
      integer,allocatable :: ref_flag(:)
      real*8 :: resid_sq_max, tol_res_sq
!------------------------------------------------------------------------------------------
!
      iprint=0

      allocate(resid_sq(NRELES))
      resid_sq = 0.d0
!
! !  ...field variables flag
!       nflag(1:NR_PHYSA)=(/0,0,1,1,1/)
!  ...compute total residual and number of dof
      nrdof_total = NRDOFSU + NRDOFSF + NRDOFSQ
      ! residual = 0.d0
      ! mdle = 0
!$OMP PARALLEL
!$OMP DO &
!$OMP SCHEDULE(DYNAMIC)
      do iel=1,NRELES
        ! call nelcon(mdle, mdle)
        call elem_residual(iel, resid_sq(iel))
 !        if (iprint.eq.1) then
 !          write(*,7010) iel, mdle, resid
 ! 7010     format('compute_residual: iel, mdle = ',2i5,  &
 !                 ' element residual = ',e12.5)
        ! endif
        ! call find_ndof(mdle, ndofH,ndofE,ndofV,ndofQ)
! !  .....subtract the bubbles from trace variables
!         do iphys=1,NR_PHYSA
! !      ...skip if field variable, otherwise remove middle dof and
! !         leave trace dof only
!           if (nflag(iphys).eq.1) cycle
!           select case(DTYPE(iphys))
!           case('contin'); nrdof_total=nrdof_total-ndofH*NR_COMP(iphys)
!           case('tangen'); nrdof_total=nrdof_total-ndofE*NR_COMP(iphys)
!           case('normal'); nrdof_total=nrdof_total-ndofV*NR_COMP(iphys)
!           case('discon'); nrdof_total=nrdof_total-ndofQ*NR_COMP(iphys)
!           end select
        ! enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
!
      residual = dsqrt(sum(resid_sq))

!
!  ...compute rate
      rate = 0.d0
      if (ivis.ne.0) then
        if (nrdof_total.gt.nrdof_total_old) then
          rate = log(residual_old/residual)  &
                /log(float(nrdof_total_old)/float(nrdof_total))
        endif
      endif
!
!  ...save current data
      ivis = ivis+1
      nrdof_total_old = nrdof_total
      residual_old = residual
!
!  ...store data to display
      iwork(ivis) = nrdof_total
      rwork(1,ivis) = residual
      rwork(2,ivis) = rate
!
!  ...display the convergence history
      if (.NOT. QUIET_MODE) then
        write(*,*)''
        write(*,*)'         -- Residual Report --'
        write(*,7100)
 7100   format(' Mesh  ','  Nrdof  ', ' Residual   ','     Rate ')
        do i=1,ivis
          write(*,7110)i,iwork(i),rwork(1:2,i)
 7110     format(i3,4x,i7,2x,e12.5,2x,f8.3)
        enddo
        write(*,*)''
      endif

      tol_res_sq = 0.0625d0
      resid_sq_max = maxval(resid_sq)

      ref_dump = 18
      open(unit=ref_dump,file='output/elem_ref_flag', &
          form='formatted',access='sequential',status='unknown')
      
      allocate(ref_flag(NRELES))
      ref_flag = 0      
      do iel=1,NRELES
        if (resid_sq(iel).ge.tol_res_sq*resid_sq_max) ref_flag(iel) = 1
        write(ref_dump,*) iel,ref_flag(iel), resid_sq(iel)/resid_sq_max
 ! 7531 format(i6,i2)
      enddo

      deallocate(resid_sq)
      deallocate(ref_flag)
      close(ref_dump)
!
!
      end subroutine compute_residual
