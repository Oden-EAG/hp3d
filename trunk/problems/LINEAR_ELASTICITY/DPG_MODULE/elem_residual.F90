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
      use parametersDPG
      use element_data
      use data_structure3D
      use isotropic_elast_material
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
!  ...derivatives wrt physical coordinates, flux
      real*8, dimension(3)           :: dv1,dv2,vec
!
!  ...SHAPE FUNCTIONS
!     H1
      real*8, dimension(  MAXbrickH) :: shapH
      real*8, dimension(3,MAXbrickH) :: gradH
      integer                        :: nrdofH
!     discontinuous H1
      real*8, dimension(  MAXbrickHH) :: shapHH
      real*8, dimension(3,MAXbrickHH) :: gradHH
      integer                         :: nrdofHH
!     H(div)
      real*8, dimension(3,MAXbrickV) :: shapV
      real*8, dimension(  MAXbrickV) :: divV ! never used
      integer                        :: nrdofV
!  ...flux
      real*8, dimension(3)            :: sigma,sigma_n
!
!  ...MATRICES
!     stiffnes matrix for the local Riesz matrix in LAPACK format
      real*8, dimension(9*MAXbrickHH*(3*MAXbrickHH+1)/2) :: Gram
!
!     Resid vector for the enriched space
      real*8, dimension(3*MAXbrickHH) :: EnrResid,EnrResidc
!
!  ...elasticities in master and physical coordinates,
      real*8, dimension(3,3,3,3) :: C,Cnew,SIP,SIPNew
      real*8, dimension(3,3)     :: kwave
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
      integer, dimension(MAXEQNH,MAXbrickH) :: dofH
      integer, dimension(MAXEQNE,MAXbrickE) :: dofE
      integer, dimension(MAXEQNV,MAXbrickV) :: dofV
      integer, dimension(MAXEQNQ,MAXbrickQ) :: dofQ
      real*8,  dimension(MAXEQNH)           :: solH
      real*8,  dimension(MAXEQNH,3)         :: dsolH
      real*8,  dimension(3)                 :: solV
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
      integer :: i,j,k,l,m,n,mm,nn,k1,k2,k3,m1,m2,m3,ipt,icomp,jcomp,  &
                 nint,ifc,nsign,iprint,iload,iflag,info,info1,info2,info3
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
!  ...Load number which to calculate residual from
      iload = 1
!
      select case(Mdle)
      case(1)
        iprint=0
      case default
        iprint=0
      end select
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
      select case(etype)
      case('mdlb'); nordP = NODES(Mdle)%order+NORD_ADD*111
      case('mdln','mdld'); nordP = NODES(Mdle)%order+NORD_ADD
      case('mdlp'); nordP = NODES(Mdle)%order+NORD_ADD*11
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
             /,'               dofV(1,1:8) = ',6(e12.5,2x))
      endif
!
!  ...initialize the enriched local element residual vector
      EnrResid = ZERO
!  ...initialize the Gram matrix
      Gram = ZERO
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!
!  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD
      call set_3Dint(etype,norder, nint,xiloc,wxiloc)
      INTEGRATION = 0

      do ipt=1,nint
        xi(1:3) = xiloc(1:3,ipt); wa = wxiloc(ipt)
!
!  .....H1 shape functions
        call shape3H(etype,xi,norder,nedge_orient,nface_orient,  &
                     nrdofH,shapH,gradH)
!
!  .....discontinuous H1 shape functions
        call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
                    x,dxdxi,dxidx,rjac,iflag)
!
!  .....integration weight
        weight = rjac*wa
!
!  .....compute the approximate solution on the master element
        solH(1:MAXEQNH) = ZERO; dsolH(1:MAXEQNH,1:3) = ZERO
        do k=1,nrdofH
          solH(1:MAXEQNH) = solH(1:MAXEQNH)  &
                          + dofH(1:MAXEQNH,k)*shapH(k)
          do m=1,3
            dsolH(1:MAXEQNH,m) = dsolH(1:MAXEQNH,m)  &
                               + dofH(1:MAXEQNH,k)*gradH(m,k)
          enddo
        enddo
!
 !        call exact(x,Mdle, valH,dvalH,d2valH, valE,dvalE,d2valE,  &
 !                           valV,dvalV,d2valV, valQ,dvalQ,d2valQ)
 !        if (iprint.eq.1) then
 !          write(*,8010) Mdle,x(1:3)
 ! 8010     format('elem_residual: Mdle = ',i6,',  x = ',3f8.3)
 !          write(*,8020) solH, gradH(1:3)
 ! 8020     format('  approximate solution and gradient = ',  &
 !                 e20.12,3x,3(e20.12,2x))
 !          write(*,8030) valH(1),dvalH(1,1:3)
 ! 8030     format('  exact solution and gradient       = ',  &
 !                 e20.12,3x,3(e20.12,2x))
 !        endif
!
!  .....compute the elastictiy tensor
        call getC(X, C)
!
! !  .....compute the semi inner product tensor (depends upon chosen norm)
!         call getSemiIPTensor(x, SIP)
!
!  .....precompute tensors in master coordinates
        CNew = ZERO; kwave = ZERO; SIPNew = ZERO
!
! IDEA: Could maybe optimize this further using symmetries of C
!
! !  .....loop through coordinates
!         do i=1,3; do j=1,3
!           do k=1,3; do l=1,3
! !  .........loop to change variables
!             do m=1,3; do n=1,3
!               CNew(i,j,k,l) = CNew(i,j,k,l)  &
!                             + C(i,m,k,n)*dxidx(j,m)*dxidx(l,n)
!               SIPNew(i,j,k,l) = SIPNew(i,j,k,l)  &
!                               + SIP(i,m,k,n)*dxidx(j,m)*dxidx(l,n)
!             enddo; enddo
!           enddo; enddo
!           if (i.eq.j) kwave(i,j) = - RHO*OMEGA**2
!         enddo; enddo
!
!  .....get the source term
        call getf(Mdle,x, fval)
!
!  .....OUTER loop through enriched H1 test functions
        do k1=1,nrdofHH
          do icomp=1,MAXEQNH
            m1 = (k1-1)*MAXEQNH+icomp
!
!           E N R I C H E D   R E S I D U A L   V E C T O R
!
            tmp = ZERO
            ! Bu part
            do m=1,3; do mm=1,3; do nn=1,3
              tmp = tmp  &
                  + Cnew(m,mm,icomp,nn)*dsolH(m,mm)*gradHH(nn,k1)
            enddo; enddo; enddo
            ! Bu-l part
            tmp = tmp  &
                - fval(icomp,iload)*shapHH(k1)
            ! accumulate
            EnrResid(m1) = EnrResid(m1) + tmp*weight
!
!  .........INNER loop through enriched H1 trial functions
            do m2=m1,3*nrdofHH
              k2 = int((m2-1)/3)+1
              jcomp = mod(m2-1,3)+1
              !  index must be given in triangular format for LAPACK
              k = nk(m1,m2)
!
!             G R A M   M A T R I X
!
              tmp = ZERO
!             SEMI INNER PRODUCT PART
              do m=1,3; do n=1,3
                tmp = tmp  &
                    + SIPNew(icomp,m,jcomp,n)*gradHH(m,k1)*gradHH(n,k2)
              enddo; enddo
!             L2 PART
              if (icomp.eq.jcomp) tmp = tmp + shapHH(k1)*shapHH(k2)
!             accumulate
              Gram(k) = Gram(k) + tmp*weight
!  .........INNER loop
            enddo
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
        INTEGRATION = NORD_ADD
        call set_2Dint(ftype,nordf, nint,tloc,wtloc)
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
!  .......H1 shape functions (for geometry)
          call shape3H(etype,xi,norder,nedge_orient,nface_orient,  &
                       nrdofH,shapH,gradH)
!
!  .......discontinuous H1 shape functions
          call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!  .......H(div) shape functions (for fluxes)
          call shape3V(etype,xi,norder,nface_orient,  &
                       nrdofV,shapV,divV)
!  .......geometry
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,  &
                       x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
          weight = wa*bjac
!
!  .......compute approximate flux vector
          do icomp=1,MAXEQNH
            !  compute i-th H(div) extension of solution (a vector)
            solV(1:3) = 0.d0
            do k1=1,nrdofV
              solV(1:3) = solV(1:3)  &
                        + dofV(icomp,m1)*shapV(1:3,k1)
            enddo
            !  change to physical coordinates
            sigma(1:3) = dxdxi(1:3,1)*solV(1)  &
                       + dxdxi(1:3,2)*solV(2)  &
                       + dxdxi(1:3,3)*solV(3)
            sigma(1:3) = sigma(1:3)/rjac
            !  compute flux
            sigma_n(icomp) = sigma(1)*rn(1)+sigma(2)*rn(2)+sigma(3)*rn(3)
          enddo
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
        if (iprint.eq.1) call pause
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
        write(*,*) 'elem_dpgH1: info = ',info
        stop1
      endif
!
!  ...save copies of the RHS to compute later the residual
      EnrResidc = EnrResid
!
!  ...compute the product of inverted test Gram matrix with RHS,
!     EnrResid is overwritten with the solution
      call DPPTRS(uplo, 3*nrdofHH, 1, Gram, EnrResid, 3*MAXbrickHH, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem_dpgH1: info1 = ',info1
        stop1
      endif
!
!  ...compute the residual
      Resid = 0.d0
      do k=1,3*nrdofHH
        Resid = Resid  &
              + EnrResidc(k)*EnrResid(k)
      enddo
!
!-----------------------------------------------------------------------------------
!     E R R O R    R E P R E S E N T A T I O N    F U N C T I O N                  |
!-----------------------------------------------------------------------------------
!
!
!  ...recompute the element residual through direct integration to
!     establish anistropy flags
      residd(0:3) = 0.d0
      do ipt=1,nint
        xi(1:3) = xiloc(1:3,ipt); wa = wxiloc(ipt)
!
!  .....H1 shape functions
        call shape3H(etype,xi,norder,nedge_orient,nface_orient,  &
                     nrdofH,shapH,gradH)
!
!  .....discontinuous H1 shape functions
        call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!  .....geometry map
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
                    x,dxdxi,dxidx,rjac,iflag)
        if (iflag.ne.0) then
          write(*,1000) Mdle,rjac
 1000     format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
          stop
        endif
!
!  .....integration weight
        weight = rjac*wa
!
!  .....compute the error representation function (for only one load vector)
        psi(1:3) = 0.d0; gradpsi(1:3,1:3) = 0.d0
        do k1=1,nrdofHH
          do icomp=1,3
            m1 = (k1-1)*3+icomp
!
            psi(icomp)         = psi(icomp)         + EnrResid(m1)*shapHH(k1)
            gradpsi(icomp,1:3) = gradpsi(icomp,1:3) + EnrResid(m1)*gradHH(1:3,k1)
          enddo
        enddo
!
!  .....compute the elastictiy tensor
        call getC(x, C)
!
! !  .....compute the semi inner product tensor (depends upon chosen norm)
!         call getSemiIPTensor(x, SIP)
!
!  .....precompute in master coordinates
        SIPNew = ZERO
! !  .....loop through coordinates
!         do i=1,3; do j=1,3; do k=1,3; do l=1,3
! !  .......loop to change variables
!           do m=1,3; do n=1,3
!             SIPNew(i,j,k,l) = SIPNew(i,j,k,l)  &
!                             + SIP(i,m,k,n)*dxidx(j,m)*dxidx(l,n)
!           enddo; enddo
!         enddo; enddo; enddo; enddo
!
        !  a weird way of partitioning the H1 norm of psi...
        do k=1,3; do l=1,3
          if (k.eq.l) then
            do m=1,3; do n=1,3
              residd(k) = residd(k)  &
                        + SIPNew(m,k,n,k)*gradpsi(m,k)*gradpsi(n,k)*weight
            enddo; enddo
          else
            do m=1,3; do n=1,3
              residd(0) = residd(0)  &
                        + SIPNew(m,k,n,l)*gradpsi(m,k)*gradpsi(n,l)*weight
            enddo; enddo
          endif
        enddo; enddo
        do icomp=1,3
          residd(0) = residd(0) + psi(icomp)**2*weight
        enddo
!
!  ...end of loop through integration points
      enddo
!
!  ...a test
      diff = residd(0)+residd(1)+residd(2)+residd(3) - Resid
      if (abs(diff).gt.1.d-8*abs(Resid)) then
        write(*,*) 'Resid = ',Resid,  &
                    residd(0)+residd(1)+residd(2)+residd(3)
        write(*,*) 'residd = ',residd(0:3)
        call pause
      endif
!
!  ...determine the refinement flag (I see no reason to change this)
      select case(etype)
      case('mdlb')
        if (residd(0).lt..1d0*Resid) then
          nref(1:3) = 1
          do i=1,3
            if (residd(i).lt..1d0*Resid) nref(i)=0
          enddo
          Nref_flag = nref(1)*100+nref(2)*10+nref(3)
        else
          Nref_flag = 111
        endif
      case('mdln','mdld')
        Nref_flag = 1
      case('mdlp')
        if (residd(0).lt..1d0*Resid) then
          nref(1:2) = 1
          if (residd(1)+residd(2).lt..2d0*Resid) nref(1)=0
          if (residd(3).lt..1d0*Resid) nref(2)=0
          Nref_flag = nref(1)*10+nref(2)
        else
          Nref_flag = 111
        endif
      end select
      write(*,*) 'residd = ',residd(0:3)
      write(*,*) 'Mdle,Nref_flag = ',Mdle,Nref_flag
!!!      call pause
!
!-----------------------------------------------------------------------
!
      if (iprint.ge.1) then
        write(*,7010) Mdle, Resid
 7010   format('elem_residual: Mdle, Resid = ',i5,3x,e12.5)
        call pause
      endif
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
      integer :: i,iel,iprint
!------------------------------------------------------------------------------------------
!
      iprint=0
!
!  ...compute total residual and number of dof
      nrdof_total = NRDOFSH+NRDOFSV; residual = 0.d0
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
!
!  .....subtract the middle node H(div) dof from the global count
        nrdof_total = nrdof_total - ndofV
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
      rwork(1,ivis) = residual; rwork(2,ivis) = rate
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






