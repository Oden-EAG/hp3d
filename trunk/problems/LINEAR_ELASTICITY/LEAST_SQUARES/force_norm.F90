!--------------------------------------------------------------------------
!> Purpose : returns global force norm in V'
!!
!--------------------------------------------------------------------------
!
      subroutine force_norm
!
      use data_structure3D
      use environment, only : QUIET_MODE
!------------------------------------------------------------------------------------------
      implicit none
!
!  ...middle node
      integer :: mdle
!
!  ...force
      real*8 :: frce,force
!
!  ...visitation flag
      integer, save :: ivis = 0
!
!  ...force to display
      real*8 , dimension(20), save :: rwork
!
!  ...miscellaneous
      integer :: i,iel,iprint
!------------------------------------------------------------------------------
!
!  ...printing flag for debugging
      iprint=0
!
      force = 0.d0
      mdle = 0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_force
        if (iprint.eq.1) then
          write(*,7010) iel, mdle, frce
 7010     format('force_norm: iel, mdle = ',2i5,  &
                 ' element force = ',e12.5)
        endif
        force = force + frce
      enddo
      force = sqrt(force)
!
!  ...save current data
      ivis = ivis+1
!
!  ...store data to display
      rwork(ivis) = force
!
!  ...display the force_norm history
      if (.NOT. QUIET_MODE) then
        write(*,*)''
        write(*,*)'         -- Report --'
        write(*,7100)
 7100   format(' Mesh  ','Force Norm ')
        do i=1,ivis
          write(*,7110)i,rwork(i)
 7110     format(i3,4x,e12.5)
        enddo
        write(*,*)''
      endif
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      contains
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  !--------------------------------------------------------------------------
  !> Purpose : return the force for element
  !--------------------------------------------------------------------------
  !
        subroutine elem_force
  !
        use control, only: INTEGRATION
        use parametersDPG
        use element_data
        use data_structure3D
        use common_prob_data, only: TEST_NORM, IP
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
        real*8, dimension(  MAXbrickH)  :: shapH
        real*8, dimension(3,MAXbrickH)  :: gradH
        integer                         :: nrdofH
  !
  !  ...source term, Neumann term
        real*8, dimension(3,MAXNRHS) :: fval
  !
  !  ...3D quadrature data
        real*8, dimension(3,MAX_NINT3) :: xiloc
        real*8, dimension(MAX_NINT3)   :: wxiloc
  !
  !  ...miscellaneous
        integer :: ipt,icomp,nint,iprint,iload,iflag
        real*8  :: weight,wa,rjac
  !
  !-----------------------------------------------------------------------------------
  !      I N I T I A L I Z A T I O N                                                 |
  !-----------------------------------------------------------------------------------
  !
  !  ...printing flag for debugging
        iprint=0
  !
  !  ...Load number for which to calculate force from
        iload = 1
  !
  !  ...element type
        etype = NODES(mdle)%type
  !
  !  ...order of approximation, orientations, geometry dof's
        call find_order (mdle, norder)
        call find_orient(mdle, nedge_orient,nface_orient)
        call nodcor     (mdle, xnod)
  !
        if (iprint.eq.1) then
          write(*,7020) xnod(1,1:8),xnod(2,1:8),xnod(3,1:8)
   7020   format('elem_force: xnod  = ',8(f8.3,2x),  &
           2(  /,'                       ',8(f8.3,2x)))
        endif
  !
  !-----------------------------------------------------------------------------------
  !      E L E M E N T    I N T E G R A L                                            |
  !-----------------------------------------------------------------------------------
  !
  !  ...set up the element quadrature
        INTEGRATION = NORD_ADD
        ! INTEGRATION = 4 - IP
        call set_3Dint(etype,norder, nint,xiloc,wxiloc)
        INTEGRATION = 0
  !
  !  ...initialize force calculation
        frce = 0.d0
  !
  !  ...integrate
        do ipt=1,nint
          xi(1:3) = xiloc(1:3,ipt); wa = wxiloc(ipt)
  !
  !  .....Compute shape functions needed for test field variables and geometry
  !       H1 (geometry)
          call shape3H(etype,xi,norder,nedge_orient,nface_orient,  &
                       nrdofH,shapH,gradH)
  !
  !  .....geometry map
          call geom3D(mdle,xi,xnod,shapH,gradH,nrdofH,  &
                      x,dxdxi,dxidx,rjac,iflag)
  !
  !  .....integration weight
          weight = rjac*wa
  !
  !  .....get the source term
          call getf(mdle,x, fval)
  !
          if (iprint.ge.2) then
            write(*,*) 'elem_force: x = ', x
            write(*,*) 'elem_force: fval = ', fval(1:3,iload) 
          endif
  !
  !  .....|| f ||^2
          do icomp=1,3
            frce = frce + fval(icomp,iload)**2*weight
          enddo
  !
  !  ...end of loop through integration points
        enddo
  !
        if (iprint.ge.1) then
          write(*,7010) mdle, frce
   7010   format('elem_force: Mdle, frce = ',i5,3x,e12.5)
          call pause
        endif
  !
  !
        end subroutine elem_force
!
      end subroutine force_norm
