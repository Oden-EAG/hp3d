!--------------------------------------------------------------------------
!> Purpose : return stiffness matrix and Residual vector for element
!!
!! @param[in]  Mdle      - an element (middle node) number
!! @param[out] Resid     - element residual (squared)
!! @param[out] Nref_flag - suggested h-refinement flag
!--------------------------------------------------------------------------!------------------------------------------------------------------------------------------
!
!     r = || \sigma - C(\grad(u)) ||^2 + || f + div(\sigma) ||^2
!
!------------------------------------------------------------------------------------------
!
      subroutine elem_residual(Mdle, Resid,Nref_flag)
!
      use control, only: INTEGRATION
      use parametersDPG
      use element_data
      use data_structure3D
      use isotropic_elast_material
      use common_prob_data, only: IP
!------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in)  :: Mdle
      real*8,  intent(out) :: Resid
      integer, intent(out) :: Nref_flag
!------------------------------------------------------------------------------------------
!  ...element and face type
      character(len=4) :: etype
!
!  ...element and face order, enriched order
      integer, dimension(19) :: norder
!
!  ...node edge and face orientations
      integer, dimension(12) :: nedge_orient
      integer, dimension(6)  :: nface_orient
!
!  ...geometry
      real*8, dimension(3,MAXbrickH) :: xnod
      real*8, dimension(3)           :: xi,x
      real*8, dimension(3,3)         :: dxdxi,dxidx
!
!  ...SHAPE FUNCTIONS
!     H1  (geometry)
      real*8, dimension(  MAXbrickH) :: shapH
      real*8, dimension(3,MAXbrickH) :: gradH
      integer                        :: nrdofH
!
!  ...elasticities in master and physical coordinates,
      real*8, dimension(3,3,3,3) :: C
!
!  ...source term, Neumann term
      real*8, dimension(3,MAXNRHS) :: fval
!
!  ...3D quadrature data
      real*8, dimension(3,MAX_NINT3) :: xiloc
      real*8, dimension(MAX_NINT3)   :: wxiloc
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
!
!  ...miscellaneous
      integer :: m,n,ipt,icomp,jcomp,nint,iprint,iload,iflag
      real*8  :: weight,wa,rjac,tmp
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
!  ...Load number for which to calculate residual from
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
!
!  ...order of approximation, orientations, geometry dof's
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
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
!  ...set up the element quadrature
      INTEGRATION = NORD_ADD
      ! INTEGRATION = 4 - IP
      call set_3Dint(etype,norder, nint,xiloc,wxiloc)
      INTEGRATION = 0
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!
      Resid = 0.d0
!
!  ...compute the element residual using the formula
      do ipt=1,nint
        xi(1:3) = xiloc(1:3,ipt); wa = wxiloc(ipt)
!
!  .....H1 shape functions
        call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                     nrdofH,shapH,gradH)
!
!  .....compute the approximate solution on the PHYSICAL element
        call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                     x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)
!
!  .....geometry map
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
                    x,dxdxi,dxidx,rjac,iflag)
!
!  .....integration weight
        weight = rjac*wa
!
!  .....compute the stiffness tensor
        call getC(x, C)
!
!  .....get the source term
        call getf(Mdle,x, fval)
!
!   || \sigma - C(\grad(u)) ||^2
!
        do icomp=1,3
          do jcomp=1,3

            tmp = solV(jcomp,icomp)
            do m=1,3; do n=1,3
              tmp = tmp - C(icomp,jcomp,m,n)*dsolH(m,n)
            enddo; enddo

            !  accumulate
            Resid = Resid + tmp**2*weight

          enddo
        enddo
!
!   || f + div(\sigma) ||^2
!
        tmp=0.d0
        do icomp=1,3

          tmp = fval(icomp,iload) + divV(icomp)

          !  accumulate
          Resid = Resid + tmp**2*weight

        enddo
!
!  ...end of loop through integration points
      enddo

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
      nflag(1:NR_PHYSA)=(/1,1/)
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
