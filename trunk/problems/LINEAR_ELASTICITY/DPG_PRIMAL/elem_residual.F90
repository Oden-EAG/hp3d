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
      use control, only: INTEGRATION
      use primal_module
      use parametersDPG
      use element_data
      use data_structure3D
      use isotropic_elast_material
      use common_prob_data, only: TEST_NORM, IP
!------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in)  :: Mdle
      real*8,  intent(out) :: Resid
      integer, intent(out) :: Nref_flag
!------------------------------------------------------------------------------------------
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
!  ...geometry
      real*8, dimension(3,MAXbrickH) :: xnod
      real*8, dimension(3)           :: xi,x,rn
      real*8, dimension(3,3)         :: dxdxi,dxidx
      real*8, dimension(2)           :: t
      real*8, dimension(3,2)         :: dxidt,dxdt,rt
!
!  ...SHAPE FUNCTIONS
!     H1  (geometry)
      real*8, dimension(  MAXbrickH)  :: shapH
      real*8, dimension(3,MAXbrickH)  :: gradH
      integer                         :: nrdofH
!     H1  (test)
      real*8, dimension(  MAXbrickHH) :: shapHH
      real*8, dimension(3,MAXbrickHH) :: gradHH
      integer                         :: nrdofHH
!
!  ...Resid vector for the enriched space
      real*8, dimension(3*MAXbrickHH) :: EnrResid,EnrResidc
!     Gram matrix for the local Riesz matrix in LAPACK format
      real*8, dimension(3*MAXbrickHH*(3*MAXbrickHH+1)/2) :: Gram
!
!  ...elasticities in master and physical coordinates,
      real*8, dimension(3,3,3,3) :: C,SIP,CC,Symm
!
!  ...source term, Neumann term
      real*8, dimension(3,MAXNRHS) :: fval
!
!  ...3D quadrature data
      real*8, dimension(3,MAX_NINT3) :: xiloc
      real*8, dimension(MAX_NINT3)   :: wxiloc
!
!  ...2D quadrature data for boundary terms
      real*8, dimension(2,MAXNINT2ADD) :: tloc
      real*8, dimension(MAXNINT2ADD)   :: wtloc
!
!  ...directional contributions to element residual
      real*8, dimension(0:4) :: residd
      real*8, dimension(3)   :: nref
!
!  ...error representation function
      real*8, dimension(3)   :: psi
      real*8, dimension(3,3) :: gradpsi
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
!     flux
      real*8, dimension(3)                  :: sigma_n
! !
! !  ...exact solution (for debugging)
!       real*8,    dimension(  MAXEQNH    ) ::   valH
!       real*8,    dimension(  MAXEQNH,3  ) ::  dvalH
!       real*8,    dimension(  MAXEQNH,3,3) :: d2valH
!       real*8,    dimension(3,MAXEQNE    ) ::   valE
!       real*8,    dimension(3,MAXEQNE,3  ) ::  dvalE
!       real*8,    dimension(3,MAXEQNE,3,3) :: d2valE
!       real*8,    dimension(3,MAXEQNV    ) ::   valV
!       real*8,    dimension(3,MAXEQNV,3  ) ::  dvalV
!       real*8,    dimension(3,MAXEQNV,3,3) :: d2valV
!       real*8,    dimension(  MAXEQNQ    ) ::   valQ
!       real*8,    dimension(  MAXEQNQ,3  ) ::  dvalQ
!       real*8,    dimension(  MAXEQNQ,3,3) :: d2valQ
!
!  ...miscellaneous
      integer :: i,j,k,l,m,n,mm,nn,k1,k2,k3,m1,m2,m3,ipt,icomp,jcomp,      &
                 nint,ifc,nsign,iprint,iload,iflag,info,info1,info2,info3, &
                 nordtmp,nrHbub,nrEbub,nrVbub,nrQbub
      real*8  :: weight,wa,rjac,bjac,tmp,diff
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
!
      ! select case(Mdle)
      ! case(1)
      !   iprint=0
      ! case default
      !   iprint=0
      ! end select
      iprint=0
!
!  ...element type
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
!  ...order of approximation, orientations, geometry dof's
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
!
!  ...set the enriched order of appoximation
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
        write(*,7030) dofH(1,1:8),dofV(1,1:6)
 7030   format('elem_residual: dofH(1,1:8) = ',8(e12.5,2x),  &
             /,'               dofV(1,1:6) = ',6(e12.5,2x))
      endif
!
!  ...initialize the enriched local element residual vector
      EnrResid = ZERO
!  ...initialize the Gram matrix
      Gram = ZERO

      call getSymm(Symm)
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!
!  ...use the enriched order to set the quadrature
      INTEGRATION = nordtmp
      call set_3dint_DPG(etype,norder, nint,xiloc,wxiloc)
      INTEGRATION = 0

      do ipt=1,nint
        xi(1:3) = xiloc(1:3,ipt); wa = wxiloc(ipt)
!
!  .....Compute shape functions needed for test field variables and geometry
!       H1 (geometry)
        call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                     nrdofH,shapH,gradH)
!       H1 (test)
        call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!  .....geometry map
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
                    x,dxdxi,dxidx,rjac,iflag)
!
!  .....integration weight
        weight = rjac*wa
!
!  .....Change coordinates so the shape functions are on the physical element
        do k=1,nrdofHH
          gradHH(1:3,k) = gradHH(1,k)*dxidx(1,1:3)  &
                        + gradHH(2,k)*dxidx(2,1:3)  &
                        + gradHH(3,k)*dxidx(3,1:3)
        enddo
!
!  .....compute the approximate solution on the PHYSICAL element
        call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                     x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)
!
!  .....compute the stiffness tensor
        call getC(x, C)
!
!  .....get the source term
        call getf(Mdle,x, fval)
!
!  .....OUTER loop through enriched H1 test functions
        do k1=1,nrdofHH
          do icomp=1,3
            m1 = (k1-1)*3+icomp
!
!           E N R I C H E D   R E S I D U A L   V E C T O R
!
            tmp = ZERO
            ! Bu part
            do m=1,3; do mm=1,3; do nn=1,3
              tmp = tmp  &
                  + C(m,mm,icomp,nn)*dsolH(m,mm)*gradHH(nn,k1)
            enddo; enddo; enddo
            ! Bu-l part
            tmp = tmp  &
                - fval(icomp,iload)*shapHH(k1)
            ! accumulate
            EnrResid(m1) = EnrResid(m1) + tmp*weight
!
!           G R A M   M A T R I X
!
            select case(TEST_NORM)
! GALERKIN
!   (C:grad(v_2),grad(v)) + (v_2,v)
            case(1)
              do m2=m1,3*nrdofHH
                jcomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-1)/3)+1
                  do m=1,3; do n=1,3
                  tmp = tmp  &
                      + C(jcomp,n,icomp,m)*gradHH(n,k2)*gradHH(m,k1)
                  enddo; enddo
                  Gram(k) = Gram(k)  &
                          + ( tmp + shapHH(k1)*shapHH(k2) )*weight
                endif
              enddo
! 
! MATHEMATICIANS
!   (grad(v_2),grad(v)) + (v_2,v)
            case(2)
              do m2=m1,3*nrdofHH
                jcomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( shapHH(k1)*shapHH(k2)  &
                            + gradHH(1,k1)*gradHH(1,k2)  &
                            + gradHH(2,k1)*gradHH(2,k2)  &
                            + gradHH(3,k1)*gradHH(3,k2) )*weight
                endif
              enddo
! PSEUDO-STRAIN
!   (C:eps(v_2),C:eps(v)) + (v_2,v)
              case(3)
                do m2=m1,3*nrdofHH
                  jcomp = mod(m2-1,3)+1
                  if (icomp.eq.jcomp) then
                    k = nk(m1,m2)
                    k2 = int((m2-1)/3)+1
                    do m=1,3; do n=1,3
                    tmp = tmp  &
                        + CC(jcomp,n,icomp,m)*gradHH(n,k2)*gradHH(m,k1)
                    enddo; enddo
                    Gram(k) = Gram(k)  &
                            + ( tmp + shapHH(k1)*shapHH(k2) )*weight
                  endif
                enddo
! PSEUDO-STRESS
!   (eps(v_2),eps(v)) + (v_2,v)
            case(4)
                do m2=m1,3*nrdofHH
                  jcomp = mod(m2-1,3)+1
                  if (icomp.eq.jcomp) then
                    k = nk(m1,m2)
                    k2 = int((m2-1)/3)+1
                    do m=1,3; do n=1,3
                    tmp = tmp  &
                        + Symm(jcomp,n,icomp,m)*gradHH(n,k2)*gradHH(m,k1)
                    enddo; enddo
                    Gram(k) = Gram(k)  &
                            + ( tmp + shapHH(k1)*shapHH(k2) )*weight
                  endif
                enddo
            end select
!  .....OUTER loop
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
!  .....set 2D quadrature
        INTEGRATION = nordtmp
        call set_2dint_DPG(ftype,nordf, nint,tloc,wtloc)
        INTEGRATION = 0
!
!  .....loop through integration points
        do ipt=1,nint
!
!  .......face coordinates
          t(1:2) = tloc(1:2,ipt); wa = wtloc(ipt)
!
!  .......master element coordinates using face parameterization
          call face_param(etype,ifc,t, xi,dxidt)
!
!  .......Compute shape functions needed for test field variables and geometry
!         H1 (geometry)
          call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                       nrdofH,shapH,gradH)
!         H1 (test)
          call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)

!  .......geometry
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,  &
                       x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
          weight = wa*bjac
!
!  .......compute the approximate solution on the PHYSICAL element
          call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                       x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)
!
!  .......compute approximate flux vector
          sigma_n(1:3) = solV(1,1:3)*rn(1)  &
                       + solV(2,1:3)*rn(2)  &
                       + solV(3,1:3)*rn(3)
 !          call exact(x,Mdle,  &
 !                    valH,dvalH,d2valH, valE,dvalE,d2valE,  &
 !                    valV,dvalV,d2valV, valQ,dvalQ,d2valQ)
 !          valVn = valV(1,1)*rn(1)+valV(2,1)*rn(2)+valV(3,1)*rn(3)
 !          if (iprint.eq.1) then
 !            write(*,8040) if,x(1:3)
 ! 8040       format('elem_residual: if = ',i6,' x = ',3f8.3)
 !            write(*,8050) solVn
 ! 8050       format('  approximate flux = ',e20.12)
 !            write(*,8060) valVn
 ! 8060       format('  exact flux       = ',e20.12)
 !          endif
!
!  .......OUTER loop through enriched H1 test functions
          do k1=1,nrdofHH
            do icomp=1,3
              m1 = (k1-1)*3+icomp
!
!  ...........accumulate for the residual vector
              EnrResid(m1) = EnrResid(m1)  &
                           - sigma_n(icomp)*shapHH(k1)*weight
!  .......OUTER loop
            enddo
          enddo
!
!  .....end of loop over integration points
        enddo
        ! if (iprint.eq.1) call pause
      enddo
      if (iprint.ge.1) then
        write(*,7015) EnrResid(1:nrdofHH)
 7015   format('elem_residual: FINAL EnrResid = ',10(/,10(e12.5,2x)))
        call pause
      endif
!
!-----------------------------------------------------------------------------------
!     D P G    L O C A L    A S S E M B L Y                                        |
!-----------------------------------------------------------------------------------
!
!  ...factorize the test stiffness matrix
      uplo = 'U'
      call DPPTRF(uplo, 3*nrdofHH, Gram, info)
      if (info.ne.0) then
        write(*,*) 'elem_residual: info = ',info
        stop 1
      endif
!
!  ...save copies of the RHS to compute later the residual
      EnrResidc = EnrResid
!
!  ...compute the product of inverted test Gram matrix with RHS,
!     EnrResid is overwritten with the solution
      call DPPTRS(uplo, 3*nrdofHH, 1, Gram, EnrResid, 3*MAXbrickHH, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem_residual: info1 = ',info1
        stop 1
      endif
!
!  ...compute the residual
      Resid = 0.d0
      do k=1,3*nrdofHH
        Resid = Resid  &
              + EnrResidc(k)*EnrResid(k)
      enddo
! !
! !-----------------------------------------------------------------------------------
! !     E R R O R    R E P R E S E N T A T I O N    F U N C T I O N                  |
! !-----------------------------------------------------------------------------------
! !
! !
! !  ...recompute the element residual through direct integration to
! !     establish anistropy flags
!       residd(0:3) = 0.d0
!       do ipt=1,nint
!         xi(1:3) = xiloc(1:3,ipt); wa = wxiloc(ipt)
! !
! !  .....H1 shape functions
!         call shape3H(etype,xi,norder,nedge_orient,nface_orient,  &
!                      nrdofH,shapH,gradH)
! !
! !  .....discontinuous H1 shape functions
!         call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
! !
! !  .....geometry map
!         call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
!                     x,dxdxi,dxidx,rjac,iflag)
!         if (iflag.ne.0) then
!           write(*,1000) Mdle,rjac
!  1000     format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
!           stop
!         endif
! !
! !  .....integration weight
!         weight = rjac*wa
! !
! !  .....Change coordinates so the shape functions are on the physical element
!         do k=1,nrdofHH
!           gradHH(1:3,k) = gradHH(1,k)*dxidx(1,1:3)  &
!                         + gradHH(2,k)*dxidx(2,1:3)  &
!                         + gradHH(3,k)*dxidx(3,1:3)
!         enddo
! !
! !  .....compute the error representation function (for only one load vector)
!         psi(1:3) = 0.d0; gradpsi(1:3,1:3) = 0.d0
!         do k1=1,nrdofHH
!           do icomp=1,3
!             m1 = (k1-1)*3+icomp
! !
!             psi(icomp)         = psi(icomp)         + EnrResid(m1)*shapHH(k1)
!             gradpsi(icomp,1:3) = gradpsi(icomp,1:3) + EnrResid(m1)*gradHH(1:3,k1)
!           enddo
!         enddo
! !
! !  .....compute the stiffness tensor
!         call getC(x, C)
! !
! !  .....compute the semi inner product tensor (depends upon chosen norm)
!         call getSemiIPTensor(x, SIP)

!         !  a weird way of partitioning the H1 norm of psi...
!         do k=1,3; do l=1,3
!           if (k.eq.l) then
!             do m=1,3; do n=1,3
!               residd(k) = residd(k)  &
!                         + SIP(m,k,n,k)*gradpsi(m,k)*gradpsi(n,k)*weight
!             enddo; enddo
!           else
!             do m=1,3; do n=1,3
!               residd(0) = residd(0)  &
!                         + SIP(m,k,n,l)*gradpsi(m,k)*gradpsi(n,l)*weight
!             enddo; enddo
!           endif
!         enddo; enddo
!         do icomp=1,3
!           residd(0) = residd(0) + psi(icomp)**2*weight
!         enddo
! !
! !  ...end of loop through integration points
!       enddo
! !
! !  ...a test
!       diff = residd(0)+residd(1)+residd(2)+residd(3) - Resid
!      if ((abs(diff).gt.1.d-8*abs(Resid)).and.(abs(diff).gt.1.d-14)) then
!        write(*,*) 'Resid = ', Resid
!        write(*,*) 'sum(residd) = ', residd(0)+residd(1)+residd(2)+residd(3)
!        write(*,*) 'residd = ', residd(0:3)
!        call pause
!      endif
! !
! !  ...determine the refinement flag (I see no reason to change this)
!       select case(etype)
!       case('mdlb')
!         if (residd(0).lt..1d0*Resid) then
!           nref(1:3) = 1
!           do i=1,3
!             if (residd(i).lt..1d0*Resid) nref(i)=0
!           enddo
!           Nref_flag = nref(1)*100+nref(2)*10+nref(3)
!         else
!           Nref_flag = 111
!         endif
!       case('mdln','mdld')
!         Nref_flag = 1
!       case('mdlp')
!         if (residd(0).lt..1d0*Resid) then
!           nref(1:2) = 1
!           if (residd(1)+residd(2).lt..2d0*Resid) nref(1)=0
!           if (residd(3).lt..1d0*Resid) nref(2)=0
!           Nref_flag = nref(1)*10+nref(2)
!         else
!           Nref_flag = 111
!         endif
!       end select
!
!-----------------------------------------------------------------------
!
      if (iprint.ge.1) then
        write(*,7010) Mdle, Resid
 7010   format('elem_residual: Mdle, Resid = ',i5,3x,e12.5)
        call pause
      endif
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
      nflag(1:NR_PHYSA)=(/1,0/)
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
