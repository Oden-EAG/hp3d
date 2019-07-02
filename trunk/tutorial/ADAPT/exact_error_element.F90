!------------------------------------------------------------------------------------
!> Purpose : compute H1 error over element
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[out] Derr  - norm of          error^2 over element Mdle
!! @param[out] Dnorm - norm of exact solution^2 over element Mdle
!------------------------------------------------------------------------------------
subroutine exact_error_element(Mdle, Derr,Dnorm)
  !
  use control , only : EXGEOM
  use element_data
  use data_structure3D
  use lapl
  !------------------------------------------------------------------------------------
  implicit none
  integer,                  intent(in)  :: Mdle
  real*8,dimension(MAXEQNH),intent(out) :: Derr
  real*8,dimension(MAXEQNH),intent(out) :: Dnorm
  !------------------------------------------------------------------------------------
  !  ...element type
  character(len=4) :: etype
  !       
  !  ...element order of approximation
  integer,dimension(19) :: norder
  !
  !  ...node orientations
  integer,dimension(12) :: nedge_orient
  integer,dimension(6)  :: nface_orient
  !
  !  ...geometry    
  real*8,dimension(3,MAXbrickH) :: xnod
  real*8,dimension(3)           :: xi,x
  real*8,dimension(3,3)         :: dxdxi,dxidx
  !
  !  ...approximate solution dof's
  real*8,    dimension(MAXEQNH,MAXbrickH) :: zdofH
  real*8,    dimension(MAXEQNE,MAXbrickE) :: zdofE
  real*8,    dimension(MAXEQNV,MAXbrickV) :: zdofV
  real*8,    dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
  !
  !  ...approximate solution and its derivatives
  real*8,    dimension(MAXEQNH  ) :: zsolH
  real*8,    dimension(MAXEQNH,3) :: zdsolH
  !
  !  ...shape functions and their derivatives
  real*8,dimension(  MAXbrickH) :: shapH
  real*8,dimension(3,MAXbrickH) :: dshapH,dshapHx
  !
  !  ...3D quadrature data
  real*8,dimension(3,MAX_NINT3) :: xiloc
  real*8,dimension(MAX_NINT3)   :: wxi
  !
  !  ...exact solution
  real*8,    dimension(  MAXEQNH    ) ::   zvalH
  real*8,    dimension(  MAXEQNH,3  ) ::  zdvalH
  real*8,    dimension(  MAXEQNH,3,3) :: zd2valH
  real*8,    dimension(3,MAXEQNE    ) ::   zvalE
  real*8,    dimension(3,MAXEQNE,3  ) ::  zdvalE
  real*8,    dimension(3,MAXEQNE,3,3) :: zd2valE
  real*8,    dimension(3,MAXEQNV    ) ::   zvalV
  real*8,    dimension(3,MAXEQNV,3  ) ::  zdvalV
  real*8,    dimension(3,MAXEQNV,3,3) :: zd2valV
  real*8,    dimension(  MAXEQNQ    ) ::   zvalQ
  real*8,    dimension(  MAXEQNQ,3  ) ::  zdvalQ
  real*8,    dimension(  MAXEQNQ,3,3) :: zd2valQ

  !  ...number of H1, Hcurl, Hdiv, L2 dof's
  integer :: nrdofH,nrdofE,nrdofV,nrdofQ
  !
  !  ...miscellanea
  integer :: i,k,l,nint,icase,ivar,ix,iflag,ixi
  real*8  :: wa,rjac,weight
  !      
  integer :: iprint
  !------------------------------------------------------------------------------------
  !
  iprint = 0
  if (iprint.eq.1) then
     write(*,7001) Mdle,NODES(Mdle)%type
7001 format('exact_error_h1: Mdle,type = ',i10,2x,a5)
  endif
  !
  !  ...order of approximation, orientations, geometry dof's
  call find_order     (Mdle, norder)
  call find_orient(Mdle, nedge_orient,nface_orient)
  call nodcor     (Mdle,xnod)
  call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
  !
  !  ...initialize
  Derr = 0.d0; Dnorm = 0.d0
  !
  !------------------------------------------------------------------------------------
  !     E L E M E N T    I N T E G R A L                                              ! 
  !------------------------------------------------------------------------------------
  !
  !  ...set up the element quadrature
  etype = NODES(Mdle)%type
  icase = NODES(Mdle)%case
  !
  call set_3Dint(etype,norder, nint,xiloc,wxi)
  !
  !  ...loop through integration points
  do l=1,nint
     xi(1:3) = xiloc(1:3,l); wa = wxi(l)
     !
     !  .....evaluate appropriate shape functions at the point
     call shape3H(etype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,dshapH)
     !
     !  .....geometry map
     x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0
     do k=1,nrdofH
        x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
        do i=1,3
           dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
        enddo
     enddo
     !
     !  .....evaluate the inverse derivatives and jacobian
     call geom(dxdxi, dxidx,rjac,iflag) 
     if (iflag.ne.0) then
        write(*,*) 'elem_elem_error: NEGATIVE JACOBIAN FOR Mdle = ', &
             Mdle
        write(*,*) '        rjac = ',rjac
        stop 1
     endif
     !        
     !  .....compute derivatives wrt physical coordinates
     do k=1,nrdofH
        dshapHx(1:3,k) = 0.d0
        do ixi=1,3
           dshapHx(1:3,k) = dshapHx(1:3,k) + dshapH(ixi,k)*dxidx(ixi,1:3)
        enddo
     enddo
     !     
     !  .....total weight
     weight = wa*rjac
     !
     !  .....evaluate the approximate solution
     zsolH(1:MAXEQNH) = ZERO; zdsolH(1:MAXEQNH,1:3) = ZERO
     do k=1,nrdofH
        ZsolH(1:MAXEQNH) = ZsolH(1:MAXEQNH) + ZdofH(1:MAXEQNH,k)*shapH(k)
        do ivar=1,3
           ZdsolH(1:MAXEQNH,ivar) = ZdsolH(1:MAXEQNH,ivar) +  &
                ZdofH(1:MAXEQNH,k)*dshapHx(ivar,k)
        enddo
     enddo
     !
     !  .....compute the exact solution 
     call exact(x,icase, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE,  &
          zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
     if (iprint.eq.1) then
        write(*,7003) l,x(1:3),zsolH(1),zvalH(1)
7003    format('exact_error_h1: l,x = ',i4,2x,3f8.3,2x,  &
             ' APPROX SOLUTION = ',e12.5,' EXACT SOLUTION = ',e12.5)
     endif
     !
     !  .....loop over H1 equations
     do ivar=1,MAXEQNH

        Dnorm(ivar) = Dnorm(ivar) + abs(zvalH(ivar)              )**2*weight
        Derr(ivar)  = Derr(ivar)  + abs(zvalH(ivar) - zsolH(ivar))**2*weight

        if (IERROR_PROB.eq.IERROR_H1) then
           do ix=1,3
              Dnorm(ivar) = Dnorm(ivar) &
                   + abs(zdvalH(ivar,ix)                  )**2*weight
              Derr(ivar)  = Derr(ivar)  &
                   + abs(zdvalH(ivar,ix) - zdsolH(ivar,ix))**2*weight
           end do
        end if

        !  .....loop over H1 equations        
     end do
     !
     !  ...end of loop through integration points
  enddo
  !
  !
end subroutine exact_error_element
