!---------------------------------------------------------------------------------------
!>  Purpose : compute element contributions to the total error
!!
!>  @date : Apr 2016
!---------------------------------------------------------------------------------------
!
subroutine custom_error(Mdle, Err,SolNorm)
!
      use element_data
      use data_structure3D
#if TRANSVERSE_MODE
      use transverse_isotropic_elast_material
#else
      use isotropic_elast_material
#endif
      use control, only: INTEGRATION
!------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in)  :: Mdle
      real*8 , intent(out) :: Err
      real*8 , intent(out) :: SolNorm
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
!
!  ...node edge and face orientations
      integer, dimension(12) :: nedge_orient
      integer, dimension(6)  :: nface_orient
!
!  ...SHAPE FUNCTIONS
!     H1 (geometry)
      real*8, dimension(  MAXbrickH)    :: shapH
      real*8, dimension(3,MAXbrickH)    :: gradH
      integer                           :: nrdofH
!
!  ...geometry
      real*8, dimension(3,MAXbrickH) :: xnod
      real*8, dimension(3)           :: xi,x
      real*8, dimension(3,3)         :: dxdxi,dxidx
!
!  ...stiffness tensors in master coordinates and physical coordinates
      real*8, dimension(3,3,3,3) :: C
!
!  ...temporary variables
      real*8 :: tmp1,tmp2
!
!  ...3D quadrature data
      real*8, dimension(3,MAX_NINT3+2) :: xiloc
      real*8, dimension(MAX_NINT3+2)   :: wxi
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
!  ...solution variables
      real*8, dimension(  MAXEQNH    ) ::   valH
      real*8, dimension(  MAXEQNH,3  ) ::  dvalH
      real*8, dimension(  MAXEQNH,3,3) :: d2valH
      real*8, dimension(3,MAXEQNE    ) ::   valE
      real*8, dimension(3,MAXEQNE,3  ) ::  dvalE
      real*8, dimension(3,MAXEQNE,3,3) :: d2valE
      real*8, dimension(3,MAXEQNV    ) ::   valV
      real*8, dimension(3,MAXEQNV,3  ) ::  dvalV
      real*8, dimension(3,MAXEQNV,3,3) :: d2valV
      real*8, dimension(  MAXEQNQ    ) ::   valQ
      real*8, dimension(  MAXEQNQ,3  ) ::  dvalQ
      real*8, dimension(  MAXEQNQ,3,3) :: d2valQ
!
      real*8, dimension(  MAXEQNH,3  ) ::  gradDiff
!
!  ...miscellaneous
      integer :: i,m,n,icomp,jcomp,nint,ipt,ifc,iprint,iload,iflag
      real*8  :: weight,wa,rjac
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
!  ...Load number
      iload = 1
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
!  ...determine solution dof
      call solelm(Mdle, dofH,dofE,dofV,dofQ)
!
!  ...initialize error and norm
      Err     = 0.d0
      SolNorm = 0.d0
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!
!  ...set up the element quadrature
      INTEGRATION = 3
      call set_3Dint(etype,norder, nint,xiloc,wxi)
      INTEGRATION = 0
!
!  ...loop through integration points
      do ipt=1,nint
        xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt)
!
!  .....Compute shape functions needed for trial field variables and geometry
!       H1 (geometry)
        call shape3H(etype,xi,norder,nedge_orient,nface_orient,  &
                     nrdofH,shapH,gradH)
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
!  .....get the solution
        call exact(x,1, valH,dvalH,d2valH,  &
                        valE,dvalE,d2valE,  &
                        valV,dvalV,d2valV,  &
                        valQ,dvalQ,d2valQ)
!
!  .....integration weight
        weight = wa*rjac
!
!  .....compute the approximate solution on the PHYSICAL element
        call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                     x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)
!
        call getC(X, C)
!
!  .....gradDiff
        gradDiff = dsolH-dvalH
!
!  .....accumulate errors
!
!  ||C:eps(u)||**2
!
        tmp1=0.d0
        do icomp=1,3; do jcomp=1,3
          tmp2=0.d0
          do m=1,3; do n=1,3
            tmp2 = tmp2  &
                 + C(icomp,jcomp,m,n)*gradDiff(m,n)
          enddo; enddo
          tmp1 = tmp1 + tmp2**2
        enddo; enddo
!
        Err = Err + tmp1*weight
!
!  .....accumulate for norm of solution
!
!  ||C:eps(u)||**2
!
        tmp1=0.d0
        do icomp=1,3; do jcomp=1,3
          tmp2=0.d0
          do m=1,3; do n=1,3
            tmp2 = tmp2  &
                 + C(icomp,jcomp,m,n)*dvalH(m,n)
          enddo; enddo
          tmp1 = tmp1 + tmp2**2
        enddo; enddo
  !
        SolNorm = SolNorm + tmp1*weight
!
!  ...end of loop through integration points
      enddo

end subroutine custom_error