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
        use primal_module, only : Gram
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
        integer                :: nordP
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
  !     H1  (test)
        real*8, dimension(  MAXbrickHH) :: shapHH
        real*8, dimension(3,MAXbrickHH) :: gradHH
        integer                         :: nrdofHH
  !
  !  ...Resid vector for the enriched space
        real*8, dimension(3*MAXbrickHH) :: EnrForce,EnrForcec
  !
  !  ...source term, Neumann term
        real*8, dimension(3,MAXNRHS) :: fval
  !
  !  ...3D quadrature data
        real*8, dimension(3,MAX_NINT3) :: xiloc
        real*8, dimension(MAX_NINT3)   :: wxiloc
  !
  !  ...miscellaneous
        integer :: i,j,k,l,m,n,k1,k2,m1,m2,ipt,icomp,jcomp,  &
                   nint,iprint,iload,iflag,info,info1,nordtmp
        real*8  :: weight,wa,rjac
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
  !  ...set the enriched order of appoximation
        nordtmp = NORD_ADD
        ! nordtmp = 4 - IP
        select case(etype)
        case('mdlb')        ; nordP = NODES(mdle)%order + nordtmp*111
        case('mdln','mdld') ; nordP = NODES(mdle)%order + nordtmp*1
        case('mdlp')        ; nordP = NODES(mdle)%order + nordtmp*11
        end select
  !
        if (iprint.eq.1) then
          write(*,7020) xnod(1,1:8),xnod(2,1:8),xnod(3,1:8)
   7020   format('elem_force: xnod  = ',8(f8.3,2x),  &
           2(  /,'                       ',8(f8.3,2x)))
        endif
  !
  !  ...initialize the enriched local element residual vector
        EnrForce = ZERO
  !  ...initialize the Gram matrix
        Gram = ZERO
  !
  !-----------------------------------------------------------------------------------
  !      E L E M E N T    I N T E G R A L                                            |
  !-----------------------------------------------------------------------------------
  !
  !  ...use the enriched order to set the quadrature
        INTEGRATION = nordtmp
        call set_3Dint(etype,norder, nint,xiloc,wxiloc)
        INTEGRATION = 0
  !
        do ipt=1,nint
          xi(1:3) = xiloc(1:3,ipt); wa = wxiloc(ipt)
  !
  !  .....Compute shape functions needed for test field variables and geometry
  !       H1 (geometry)
          call shape3H(etype,xi,norder,nedge_orient,nface_orient,  &
                       nrdofH,shapH,gradH)
  !       H1 (test)
          call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
  !
  !  .....geometry map
          call geom3D(mdle,xi,xnod,shapH,gradH,nrdofH,  &
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
  !  .....get the source term
          call getf(mdle,x, fval)
  !
          if (iprint.ge.2) then
            write(*,*) 'elem_force: x = ', x
            write(*,*) 'elem_force: fval = ', fval(1:3,iload)
          endif
  !
  !  .....OUTER loop through enriched H1 test functions
          do k1=1,nrdofHH
            do icomp=1,3
              m1 = (k1-1)*3+icomp
  !
  !           E N R I C H E D   R E S I D U A L   V E C T O R
  !
              ! accumulate
              EnrForce(m1) = EnrForce(m1) + fval(icomp,iload)*shapHH(k1)*weight
  !
  !           G R A M   M A T R I X
  !
              select case(TEST_NORM)
  !
  !   (v_2,v) + (grad(v_2),grad(v))
  !
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
              end select
  !  .....OUTER loop
            enddo
          enddo
  !
  !  ...end of loop through integration points
        enddo
  !
        if (iprint.ge.2) then
          write(*,7015) EnrForce(1:nrdofHH)
   7015   format('elem_force: FINAL EnrForce = ',10(/,10(e12.5,2x)))
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
          write(*,*) 'elem_force: info = ',info
          stop 1
        endif
  !
  !  ...save copies of the RHS to compute later the force
        EnrForcec = EnrForce
  !
  !  ...compute the product of inverted test Gram matrix with RHS,
  !     EnrForce is overwritten with the solution
        call DPPTRS(uplo, 3*nrdofHH, 1, Gram, EnrForce, 3*MAXbrickHH, info1 )
        if (info1.ne.0) then
          write(*,*) 'elem_force: info1 = ',info1
          stop 1
        endif
  !
  !  ...compute the force
        frce = 0.d0
        do k=1,3*nrdofHH
          frce = frce  &
                + EnrForcec(k)*EnrForce(k)
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