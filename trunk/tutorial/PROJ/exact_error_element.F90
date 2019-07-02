subroutine exact_error_element(Mdle, Derr,Dnorm)
  use parameters
  use data_structure3D
  use proj
  
  implicit none
  !
  ! ** Arguments
  integer,                        intent(in)    :: Mdle
  real*8, dimension(MAXEQN_PROB), intent(inout) :: Dnorm, Derr

  ! ** Locals
  ! * order and solution
  integer, dimension(19)          :: norder
  real*8,  dimension(3,MAXbrickH) :: xnod

  real*8 :: &
       zdofH(MAXEQNH,MAXbrickH),zdofE(MAXEQNE,MAXbrickE), &
       zdofV(MAXEQNV,MAXbrickV),zdofQ(MAXEQNQ,MAXbrickQ)
  ! * node connectivity
  integer :: nedge_orient(12), nface_orient(6)

  ! * shape function
  real*8 :: &
       shapH(MAXbrickH),dshapH(3,MAXbrickH), &
       dshapHx(3,MAXbrickH)

  ! * geometry
  real*8 :: xi(3),x(3),dxdxi(3,3),dxidx(3,3)

  ! * quadrature
  real*8 :: xiloc(3,maxint),wxi(maxint)

  real*8 :: &
       zsolH(  MAXEQNH),zdsolH(  MAXEQNH,3)

  real*8 :: &
       zvalH(  MAXEQNH),zdvalH(   MAXEQNH,3  ), &
       zd2valH(  MAXEQNH,3,3), &
       zvalE(3,MAXEQNE),zdvalE( 3,MAXEQNE,3  ), &
       zcvalE( 3,MAXEQNE    ), &
       zd2valE(3,MAXEQNE,3,3), &
       zvalV(3,MAXEQNV),zdvalV( 3,MAXEQNV,3  ), &
       zd2valV(3,MAXEQNV,3,3), &
       zvalQ(  MAXEQNQ),zdvalQ(   MAXEQNQ,3  ), &
       zd2valQ(  MAXEQNQ,3,3)
  character(len=4) :: type
  integer :: &
       i, j, k, l, n, &
       ieq, ixi, ix, ivar, &
       nint, nrdofE, nrdofH, &
       iprint, iflag, icase
  real*8  :: wa, weight, rjac

  !---------------------------------------------------------------------
  iprint = 0 
  if (iprint.eq.1) then
     write(*,7001) Mdle,NODES(Mdle)%type
7001 format(' exact_error_h1: Mdle,type = ',i10,2x,a5)
  endif

  ! determine order of approximation
  call find_order(Mdle, norder)
  call find_orient(Mdle, nedge_orient,nface_orient)
  !
  call nodcor(Mdle,xnod)

  ! determine approximate solution dofs
  call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
  !
  Derr = 0.d0; Dnorm = 0.d0
  !
  !---------------------------------------------------------------------
  type  = NODES(Mdle)%type
  icase = NODES(Mdle)%case

  call set_3Dint(type,norder, nint,xiloc,wxi)
  do l=1,nint
     xi(1:3) = xiloc(1:3,l); wa = wxi(l)
     ! evaluate appropriate shape functions at the point
     call shape3H(type,xi,norder, &
          nedge_orient,nface_orient, &
          nrdofH,shapH,dshapH)

     x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0 
     do k=1,nrdofH
        x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
        do i=1,3
           dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
        enddo
     enddo
     
     ! evaluate the inverse derivatives and jacobian
     call geom(dxdxi, dxidx,rjac,iflag) 

     ! compute derivatives wrt physical coordinates
     do k=1,nrdofH
        dshapHx(1:3,k) = 0.d0
        do ixi=1,3
           dshapHx(1:3,k) = dshapHx(1:3,k) + dshapH(ixi,k)*dxidx(ixi,1:3)
        enddo
     enddo

     weight = wa*rjac

     ! evaluate the approximate solution
     zsolH(1:MAXEQNH) = ZERO ; zdsolH(1:MAXEQNH,1:3) = ZERO
     do k=1,nrdofH
        zsolH(1:MAXEQNH) = zsolH(1:MAXEQNH) + zdofH(1:MAXEQNH,k)*shapH(k)
        do ivar=1,3
           zdsolH(1:MAXEQNH,ivar) = zdsolH(1:MAXEQNH,ivar) + zdofH(1:MAXEQNH,k)*dshapHx(ivar,k)
        enddo
     enddo

     !  .....compute the exact solution 
     call exact(x,icase, zvalH,zdvalH,zd2valH, &
          zvalE,zdvalE,zd2valE, &
          zvalV,zdvalV,zd2valV, &
          zvalQ,zdvalQ,zd2valQ)

     !  .....accumulate for the error and norm
     do ieq=1,MAXEQNH
        Dnorm(ieq) = Dnorm(ieq) + abs(zvalH(ieq))**2*weight
        Derr(ieq) = Derr(ieq) + abs(zvalH(ieq)-zsolH(ieq))**2*weight
        if (H1_PROJ.eq.1) then
           do ix=1,3
              Dnorm(ieq) = Dnorm(ieq) + zdvalH(ieq,ix)**2*weight
              Derr(ieq) = Derr(ieq) + abs(zdvalH(ieq,ix)-zdsolH(ieq,ix))**2*weight
           enddo
        endif
     enddo
  enddo

  if (Derr(1).gt.1.0E-10) then
     write(*,7000)Mdle,Derr(1),Dnorm(1)
7000 format(' exact_error_h1: Mdle,Derr,Dnorm = ',i5,2x,3(e12.5,2x))
  endif

  if (iprint.eq.1) call pause

end subroutine exact_error_element
