!---------------------------------------------------------------------------------------
!>  Purpose : compute element contributions to the total error
!!
!>  @date : Apr 2016
!---------------------------------------------------------------------------------------
!
subroutine stress_error(Mdle, Err,SolNorm)
!
      use element_data
      use data_structure3D
      use sheathed_isotropic_materials
      use control, only: INTEGRATION
!------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in)  :: Mdle
      real*8 , intent(out) :: Err
      real*8 , intent(out) :: SolNorm
!------------------------------------------------------------------------------------------
!  ...element and face type
      integer :: etype
!
!  ...number of topological entities (vertices,edges,faces)
      integer :: nrv,nre,nrf
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
!  ...temporary variables
      real*8 :: tmp
!
!  ...3D quadrature data
      real*8, dimension(3,MAX_NINT3) :: xiloc
      real*8, dimension(MAX_NINT3)   :: wxi
!
!  ...stiffness tensor
      real*8, dimension(3,3,3,3) :: C
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
!  ...miscellaneous
      integer :: k,l,m,n,nint,ipt,iprint,iload,iflag,ndom
      real*8  :: weight,wa,rjac
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
!  ...Load number
      iload = 1
!
      select case(Mdle)
      case(1)
        iprint=1
      case default
        iprint=0
      end select
!
!  ...element type
      etype = NODES(Mdle)%ntype
      nrv = NVERT(etype); nre = NEDGE(etype); nrf = NFACE(etype)
!
!  ...order of approximation, orientations, geometry dof's (don't need bc flags)
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
      call find_domain(Mdle, ndom)
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
!  .....compute the APPROXIMATE solution on the PHYSICAL element
        call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                     x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)
!
!  .....get the EXACT solution
        call exact(x,1, valH,dvalH,d2valH,  &
                        valE,dvalE,d2valE,  &
                        valV,dvalV,d2valV,  &
                        valQ,dvalQ,d2valQ)
!
!  .....Jacobian
        call geom(dxdxi, dxidx,rjac,iflag)
        if (iflag /= 0) then
          write(*,9997) mdle,ndom,rjac
9997      format(' element_error: mdle,ndom,rjac = ',i8,2x,i2,2x,e12.5)
        endif
!
!  .....compute the stiffness tensor
        call getC(x,ndom, C)
!
!  .....integration weight
        weight = wa*rjac
!
!  .....accumulate errors
        dsolH = dsolH-dvalH
        tmp=0.d0;
        do k=1,3; do l=1,3; do m=1,3; do n=1,3
        tmp = tmp  &
            + C(k,l,m,n)*dsolH(k,l)*dsolH(m,n)
        enddo; enddo; enddo; enddo
        Err = Err + tmp*weight

        tmp=0.d0;
        do k=1,3; do l=1,3; do m=1,3; do n=1,3
        tmp = tmp  &
            + C(k,l,m,n)*dvalH(k,l)*dvalH(m,n)
        enddo; enddo; enddo; enddo
        SolNorm = SolNorm + tmp*weight
!
!  ...end of loop through integration points
      enddo

end subroutine stress_error
