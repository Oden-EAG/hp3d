!--------------------------------------------------------------------------
!> Purpose : return stiffness matrix and Residual vector for element
!!
!! @param[in]  Mdle      - an element (middle node) number
!! @param[out] Resid     - element residual (squared)
!! @param[out] Nref_flag - suggested h-refinement flag
!--------------------------------------------------------------------------
!
      subroutine elem_residual(Mdle, Resid,Nref_flag)
!
      use control, only : INTEGRATION
      ! use uweak_symm_module, only : Gram
      use parametersDPG
      use element_data
      use data_structure3D
      use isotropic_elast_material
      use common_prob_data, only: SYMMETRY_TOL, TEST_NORM, IP
!------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in)  :: Mdle
      real*8,  intent(out) :: Resid
      integer, intent(out) :: Nref_flag
!------------------------------------------------------------------------------------------

!     Gram matrix for the local Riesz matrix in LAPACK format
      real*8, dimension((3*MAXbrickVV+3*MAXbrickHH)*(3*MAXbrickVV+3*MAXbrickHH+1)/2) :: Gram

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
      real*8, dimension(  MAXbrickH)  :: shapH
      real*8, dimension(3,MAXbrickH)  :: gradH
      integer                         :: nrdofH
!     H1   (test)
      real*8, dimension(  MAXbrickHH) :: shapHH
      real*8, dimension(3,MAXbrickHH) :: gradHH
      integer                         :: nrdofHH
!     H(div)  (test)
      real*8, dimension(3,MAXbrickVV) :: shapVV
      real*8, dimension(  MAXbrickVV) :: divVV
      real*8, dimension(  MAXbrickVV) :: shapVV_n
      integer                         :: nrdofVV
!
!  ...geometry
      real*8, dimension(3,MAXbrickH) :: xnod
      real*8, dimension(3)           :: xi,x,rn
      real*8, dimension(3,3)         :: dxdxi,dxidx
      real*8, dimension(2)           :: t
      real*8, dimension(3,2)         :: dxidt,dxdt,rt
      integer                        :: nsign
!
!  ...Resid vector for the enriched space
      real*8, dimension(3*MAXbrickVV+3*MAXbrickHH) :: EnrResid,EnrResidc
!
!  ...tensors in physical coordinates
      real*8, dimension(3,3,3,3) :: A,AA,Symm
!
!  ...temporary variables
      real*8                 :: tmp
      ! real*8, dimension(3)   :: tmpDivTau1,tmpDivTau2,tmpTau_n
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
!  ...approximate solution
      real*8, dimension(MAXEQNH,MAXbrickH) :: dofH
      real*8, dimension(MAXEQNE,MAXbrickE) :: dofE
      real*8, dimension(MAXEQNV,MAXbrickV) :: dofV
      real*8, dimension(MAXEQNQ,MAXbrickQ) :: dofQ
      real*8, dimension(  MAXEQNH  )       :: solH
      real*8, dimension(  MAXEQNH,3)       :: dsolH
      real*8, dimension(3,MAXEQNE  )       :: solE
      real*8, dimension(3,MAXEQNE  )       :: curlE
      real*8, dimension(3,MAXEQNV  )       :: solV
      real*8, dimension(  MAXEQNV  )       :: divV
      real*8, dimension(  MAXEQNQ  )       :: solQ
!     displacement
      real*8, dimension(3)                 :: u, uhat
!     stress
      real*8, dimension(3,3)               :: sigma, omega
      real*8, dimension(3)                 :: sigma_n
!
!  ...miscellaneous
      integer :: i,j,k,m,n,k1,k2,m1,m2,ipt,icomp,jcomp,ijcomp,kcomp,lcomp,klcomp,  &
                 nint,ifc,iprint,iload,iflag,info,info1,info2,info3,nordtmp
      real*8  :: weight,wa,rjac,brjac,diff,DDOT,l2StressWeight, alpha
!
!  ...LAPACK stuff
      character uplo
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1

!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
!  ...Load number for which to calculate residual from
      iload = 1
      l2StressWeight = 1d-4

      alpha = 1.d0
!
      select case(Mdle)
      case(1)
        iprint=1
      case default
        iprint=0
      end select
!
!  ...element type
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
!  ...order of approximation, orientations, geometry dof's (don't need bc flags)
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
!
!  ...set the enriched order of approximation
      nordtmp = NORD_ADD
      ! nordtmp = 4 - IP !max(NORD_ADD,2)
      select case(etype)
      case('mdlb')        ; nordP = NODES(Mdle)%order + nordtmp*111
      case('mdln','mdld') ; nordP = NODES(Mdle)%order + nordtmp*1
      case('mdlp')        ; nordP = NODES(Mdle)%order + nordtmp*11
      end select
!
!  ...determine solution dof
      call solelm(Mdle, dofH,dofE,dofV,dofQ)
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
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!
!  ...set up the element quadrature
      INTEGRATION = nordtmp
      call set_3dint_DPG(etype,norder, nint,xiloc,wxi)
      INTEGRATION = 0
!
!  ...Auxiliary tensors
      call getSymm(Symm)
!
!  ...loop through integration points
      do ipt=1,nint
        xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt)
!
!  .....Compute shape functions needed for test/trial field variables and geometry
!       H1 (geometry)
        call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                     nrdofH,shapH,gradH)
!       H1 (test)
        call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!       H(div) (test)
        call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!  .....geometry map
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
                    x,dxdxi,dxidx,rjac,iflag)
        if (iflag.ne.0) then
          write(*,1000) Mdle,rjac
 1000     format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
          stop
        endif
!
!  .....Change coordinates so the test functions are on the physical element
!       H1 (test)
        do k=1,nrdofHH
          gradHH(1:3,k) = gradHH(1,k)*dxidx(1,1:3)  &
                        + gradHH(2,k)*dxidx(2,1:3)  &
                        + gradHH(3,k)*dxidx(3,1:3)
        enddo
!       H(div) (test)
        do k=1,nrdofVV
          shapVV(1:3,k) = dxdxi(1:3,1)*shapVV(1,k)  &
                        + dxdxi(1:3,2)*shapVV(2,k)  &
                        + dxdxi(1:3,3)*shapVV(3,k)
        enddo
        shapVV(1:3,1:nrdofVV) = shapVV(1:3,1:nrdofVV)/rjac
        divVV(1:nrdofVV) = divVV(1:nrdofVV)/rjac
!
!  .....integration weight
        weight = wa*rjac
!
!  .....compute the approximate solution on the PHYSICAL element
        call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                     x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)
!       store with physically meaning variable names
        u(1:3)        = solQ(1:3)
        ! write(*,*) 'u=',u(1),u(2),u(3)
        ! 
        sigma(1,1:3)  = (/  solQ(4)  ,  solQ(7)  ,  solQ(8)  /)
        sigma(2,1:3)  = (/  solQ(7)  ,  solQ(5)  ,  solQ(9)  /)
        sigma(3,1:3)  = (/  solQ(8)  ,  solQ(9)  ,  solQ(6)  /)
        ! 
        omega(1,1:3)  = (/     0.d0  ,  solQ(12) , -solQ(11) /)
        omega(2,1:3)  = (/ -solQ(12) ,     0.d0  ,  solQ(10) /)
        omega(3,1:3)  = (/  solQ(11) , -solQ(10) ,     0.d0  /)
!
!  .....compute the compliance tensor
        call getA(x, A)
!
!  .....need this for the adjoint graph norm
        call getAA(X, AA)
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
!
!           G R A M   M A T R I X
!
            select case(TEST_NORM)
!
!   ADJOINT GRAPH
!   (A:\tau_2,A:\tau)+(div(\tau_2),div(tau))+(\tau_2,\tau)
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
                                + ( AA(icomp,kcomp,jcomp,lcomp)                 &
                                   *shapVV(lcomp,k1)*shapVV(kcomp,k2) )*weight
                  enddo
                enddo
                if (icomp.eq.jcomp) then
                  Gram(k) = Gram(k)  + divVV(k1)*divVV(k2)*weight
                  Gram(k) = Gram(k)  &
                          + (1.d0+alpha)*(               &
                            + shapVV(1,k1)*shapVV(1,k2)  &
                            + shapVV(2,k1)*shapVV(2,k2)  &
                            + shapVV(3,k1)*shapVV(3,k2) )*weight
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
                                   *shapVV(lcomp,k1)*gradHH(kcomp,k2) )*weight
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
                          + ( divVV(k1)*divVV(k2)  &
                            + shapVV(1,k1)*shapVV(1,k2)  &
                            + shapVV(2,k1)*shapVV(2,k2)  &
                            + shapVV(3,k1)*shapVV(3,k2) )*weight
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
!            G R A M   M A T R I X
!
            select case(TEST_NORM)
!
!   (\grad(v_1),A:\tau_2+\grad(v_2))+alpha*(v_1,v_2)
!
            case(1)
              do m2=m1,3*nrdofVV+3*nrdofHH
                jcomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-3*nrdofVV-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( alpha*shapHH(k1)*shapHH(k2)  &
                            + gradHH(1,k1)*gradHH(1,k2)  &
                            + gradHH(2,k1)*gradHH(2,k2)  &
                            + gradHH(3,k1)*gradHH(3,k2) )*weight
                endif
              enddo
!
!   (v_2,v) + (grad(v_2),grad(v))
!
            case(2)
              do m2=m1,3*nrdofVV+3*nrdofHH
                jcomp = mod(m2-1,3)+1
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
!  ......END OUTER LOOP
          enddo
        enddo
!
!  ...end of loop through integration points
      enddo
!
!-----------------------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L                                           |
!-----------------------------------------------------------------------------------
!
!  ...loop through element faces
      do ifc=1,nrf
!
!  .....sign factor to determine the OUTWARD normal unit vector
        nsign = nsign_param(etype,ifc)
!
!  .....face type
        ftype = face_type(etype,ifc)
!
!  .....face order of approximation
        call face_order(etype,ifc,norder, nordf)
!
!  .....set up the face quadrature
        INTEGRATION = nordtmp
        call set_2dint_DPG(ftype,nordf, nint,tloc,wt)
        INTEGRATION = 0
!
!  .....loop through face integration points
        do ipt=1,nint
!
!  .......face coordinates
          t(1:2) = tloc(1:2,ipt); wa = wt(ipt)
!
!  .......master element coordinates using face parameterization
          call face_param(etype,ifc,t, xi,dxidt)
!
!  .......Compute shape functions needed for test variables and geometry
!         H1 (geometry)
          call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                     nrdofH,shapH,gradH)
!
!         H1 (test)
          call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!         H(div) (test)
          call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!
!  .......geometry map
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,  &
                       x,dxdxi,dxidx,rjac,dxdt,rn,brjac)


!         H(div) (test)
          do k=1,nrdofVV
            shapVV(1:3,k) = dxdxi(1:3,1)*shapVV(1,k)  &
                          + dxdxi(1:3,2)*shapVV(2,k)  &
                          + dxdxi(1:3,3)*shapVV(3,k)
            shapVV_n(k) = shapVV(1,k)*rn(1)  &
                        + shapVV(2,k)*rn(2)  &
                        + shapVV(3,k)*rn(3)
          enddo

          shapVV_n(1:nrdofVV) = shapVV_n(1:nrdofVV)/rjac

          weight = wa*brjac
!
!  .....compute the approximate solution on the PHYSICAL element
          call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                       x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)
          uhat(1:3)    = solH(1:3)
          sigma_n(1:3) = solV(1,1:3)*rn(1)+solV(2,1:3)*rn(2)+solV(3,1:3)*rn(3)
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
                           - uhat(kcomp)*shapVV_n(k1)*weight
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
                           - sigma_n(kcomp)*shapHH(k1)*weight
!
!  .......OUTER loop
            enddo
          enddo
!
!  .....end of loop over integration points
        enddo
!
!  ...end of loop over faces
      enddo
!
!-----------------------------------------------------------------------------------
!     D P G    L O C A L    A S S E M B L Y                                        |
!-----------------------------------------------------------------------------------
!
!  ...factor the Gram matrix
      uplo = 'U'
      call DPPTRF(uplo, 3*nrdofVV + 3*nrdofHH, Gram, info)
      if (info.ne.0) then
        write(*,*) 'elem_residual: info = ',info ; stop
      endif
!
!  ...save copies of the RHS to compute later the residual
      EnrResidc = EnrResid
!
!  ...compute the product of inverted test Gram matrix with RHS,
!     EnrResid is overwritten with the solution
      call DPPTRS(uplo, 3*nrdofVV + 3*nrdofHH, 1, Gram, EnrResid, 3*MAXbrickVV+3*MAXbrickHH, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem_residual: info1 = ',info1
        stop
      endif
!
!  ...compute the residual
      Resid = DDOT(3*nrdofVV + 3*nrdofHH,EnrResidc,1,EnrResid,1)
!
!-----------------------------------------------------------------------------------
!             R E F I N E M E N T  F L A G S                                       |
!-----------------------------------------------------------------------------------
!  ...determines what refinement flag to use
!  ...if isotropic h refinements
      call get_isoref(Mdle, Nref_flag)
      ! Nref_flag = 110
!  ...if anisotropic h refinements -> missing
!
!
      end subroutine elem_residual

!--------------------------------------------------------------------------
!> Purpose : returns global residual
!!
!--------------------------------------------------------------------------
!
      subroutine compute_residual
!
      use data_structure3D
      use environment, only : QUIET_MODE
!------------------------------------------------------------------------------------------
      implicit none
!
!  ...middle node
      integer :: mdle
!
!  ...residual
      real*8 :: resid,residual
!
!  ...rate
      real*8 :: rate
!
!  ...number of dof for higher order node
      integer :: ndofH,ndofE,ndofV,ndofQ
      integer :: nrdof_total
!
!  ...refinement flag
      integer :: nref_flag
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
      integer :: i,iel,iphys,iprint
!------------------------------------------------------------------------------------------
!
      iprint=0
!
!  ...field variables flag
      nflag(1:NR_PHYSA)=(/0,0,1,1,1/)
!  ...compute total residual and number of dof
      nrdof_total = NRDOFSH + NRDOFSE + NRDOFSV + NRDOFSQ
      residual = 0.d0
      mdle = 0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_residual(mdle, resid,nref_flag)
        if (iprint.eq.1) then
          write(*,7010) iel, mdle, resid
 7010     format('compute_residual: iel, mdle = ',2i5,  &
                 ' element residual = ',e12.5)
        endif
        residual = residual + resid
        call find_ndof(mdle, ndofH,ndofE,ndofV,ndofQ)
!  .....subtract the bubbles from trace variables
        do iphys=1,NR_PHYSA
!      ...skip if field variable, otherwise remove middle dof and
!         leave trace dof only
          if (nflag(iphys).eq.1) cycle
          select case(DTYPE(iphys))
          case('contin'); nrdof_total=nrdof_total-ndofH*NR_COMP(iphys)
          case('tangen'); nrdof_total=nrdof_total-ndofE*NR_COMP(iphys)
          case('normal'); nrdof_total=nrdof_total-ndofV*NR_COMP(iphys)
          case('discon'); nrdof_total=nrdof_total-ndofQ*NR_COMP(iphys)
          end select
        enddo
      enddo
      residual = sqrt(residual)
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
        write(*,*)'         -- Error Report --'
        write(*,7100)
 7100   format(' Mesh  ','  Nrdof  ', ' Residual   ','     Rate ')
        do i=1,ivis
          write(*,7110)i,iwork(i),rwork(1:2,i)
 7110     format(i3,4x,i7,2x,e12.5,2x,f8.3)
        enddo
        write(*,*)''
      endif
!
!
      end subroutine compute_residual
