!------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for elasticity problem 
!! @param[in]  Mdle      - middle node number
!! @param[out] Bloc      - elem load vector(s)
!! @param[in]  Nrow_Bloc - number of rows of load vector
!! @param[out] Aloc      - elem stiffness matrix
!! @param[in]  Ncol_Aloc - number of columns of stiffness matrix
!------------------------------------------------------------------------------------
!
subroutine elem_lapl(Mdle,Bloc,Nrow_Bloc,Aloc,Ncol_Aloc)
  use data_structure3D
  use element_data
  use lapl
  use assembly , only: NR_RHS
  !------------------------------------------------------------------------------------
  implicit none
  integer,                               intent(in)  :: Mdle
  integer,                               intent(in)  :: Nrow_Bloc,Ncol_Aloc
  real*8, dimension(Nrow_Bloc,NR_RHS),   intent(out) :: Bloc
  real*8, dimension(Nrow_Bloc,Ncol_Aloc),intent(out) :: Aloc
  !------------------------------------------------------------------------------------
  !
  !  ...element and face type
  character(len=4) :: etype,ftype
  !
  !  ...element order, face order, edge and face orientations
  integer,dimension(19) :: norder
  integer,dimension(5)  :: nordf
  integer,dimension(12) :: nedge_orient
  integer,dimension(6)  :: nface_orient
  !
  !  ...shape functions and their derivatives 
  real*8,dimension(  MAXbrickH) :: shapH
  real*8,dimension(3,MAXbrickH) :: dshapH,dshapHx
  !
  !  ...geometry 
  real*8,dimension(3,MAXbrickH) :: xnod
  real*8,dimension(3)           :: xi,x
  real*8,dimension(3,3)         :: dxdxi,dxidx
  real*8,dimension(3)           :: t
  real*8,dimension(3,2)         :: dxidt,dxdt
  !
  !  ...2D and 3D quadrature data
  real*8,dimension(3,MAXbrickH) :: xiloc
  real*8,dimension(  MAXbrickH) :: wxi
  !
  !  ...exact solution routine
  real*8 :: &
       zvalH(    MAXEQNH    ), &
       zdvalH(   MAXEQNH,3  ), &
       zd2valH(  MAXEQNH,3,3), &
       zvalE(  3,MAXEQNE    ), &
       zdvalE( 3,MAXEQNE,3  ), &
       zd2valE(3,MAXEQNE,3,3), &
       zvalV(  3,MAXEQNV    ), &
       zdvalV( 3,MAXEQNV,3  ), &
       zd2valV(3,MAXEQNV,3,3), &
       zvalQ(    MAXEQNQ    ), &
       zdvalQ(   MAXEQNQ,3  ), &
       zd2valQ(  MAXEQNQ,3,3)
  !
  !  ...miscellanea
  integer :: icase,iflag,ivar,  i,j,k,l,  nint,nrdofH,  ixi,  k1,k2, &
       ibeg,iend,jbeg,jend,l1,l2
  real*8 ::  rjac,wa,weight,zs
  !
  integer :: iprint
  !
  !------------------------------------------------------------------------------------
  !     
  !  ...order of approximation, element nodes, orientations
  call find_order     (Mdle, norder)
  call find_orient(Mdle, nedge_orient,nface_orient)
  call nodcor     (Mdle, xnod)
  !
  !  ...clear spaces for the element matrices                    
  Bloc = ZERO
  Aloc = ZERO
  !
  !------------------------------------------------------------------------------------
  !     E L E M E N T   I N T E G R A L S
  !------------------------------------------------------------------------------------
  !
  !  ...set up the element quadrature
  etype = NODES(Mdle)%type
  call set_3Dint(etype,norder, nint,xiloc,wxi)
  !
  !  ...loop through integration points
  do l=1,nint
     xi(1:3) = xiloc(1:3,l); wa = wxi(l)
     !
     !  .....derivatives and values of the shape functions at the point
     call shape3H(etype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,dshapH)
     !
     !  .....geometry map
     ! iso parameterization
     x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0 
     do k=1,nrdofH
        x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
        do i=1,3
           dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
        enddo
     enddo
     !
     !  .....evaluate the inverse derivatives and jacobian
     call geom(dxdxi, dxidx,rjac, iflag) 
     if (iflag.ne.0) then
        write(*,*) 'elem_proj: NEGATIVE JACOBIAN FOR Mdle = ',Mdle
        write(*,*) '        rjac = ',rjac
        stop 1
     endif
     !
     !  .....derivatives wrt physical coordinates
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
     !  .....get rhs terms
     icase=1
     call exact(x,icase, &
          zvalH,zdvalH,zd2valH, &
          zvalE,zdvalE,zd2valE, &
          zvalV,zdvalV,zd2valV, &
          zvalQ,zdvalQ,zd2valQ)
     !
     !  .....loop over TEST functions
     do k1=1,nrdofH
        !
        !         L O A D   V E C T O R
        Bloc(k1,1) = Bloc(k1,1) + zvalH(1)*shapH(k1)*weight
        if (IERROR_PROB.eq.IERROR_H1) then
           do ivar=1,3
              Bloc(k1,1) = Bloc(k1,1) + zdvalH(1,ivar)*dshapHx(ivar,k1)*weight
           enddo
        endif
        !
        !  .......loop over TRIAL functions
        do k2=1,nrdofH
           !
           !           S T I F F N E S S   M A T R I X
           zs = shapH(k1)*shapH(k2)
           if (IERROR_PROB.eq.IERROR_H1) then
              zs = zs &
                   + dshapHx(1,k1)*dshapHx(1,k2) &
                   + dshapHx(2,k1)*dshapHx(2,k2) &
                   + dshapHx(3,k1)*dshapHx(3,k2)
           endif
           Aloc(k1,k2) = Aloc(k1,k2) + zs*weight
           !
           !  .......TRIAL functions loop
        enddo
        !  .....TEST functions loop
     enddo
     !
     !  ...end of loop through integration points
  enddo
  !
  iprint=0
  if (iprint.eq.1) then
     !
     !
50   write(*,*) 'elem_h1: Bloc = '
     write(*,*) 'elem_h1: SET jbeg,jend,l1,l2'
     read(*,*) jbeg,jend,l1,l2
     if (jbeg.eq.0) go to 55
     do j=jbeg,jend
        write(*,7005) j,Bloc(j,l1:l2)
     enddo
     go to 50
55   write(*,*) 'elem_h1: Aloc = '
     write(*,*) 'elem_h1: SET jbeg,jend,ibeg,iend'
     read(*,*) jbeg,jend,ibeg,iend
     if (jbeg.eq.0) go to 60
     do j=jbeg,jend
        write(*,7005) j,Aloc(j,ibeg:iend)
#if C_MODE
7005    format(i3,2x,5(2e11.4,2x),10(/,5x,5(2e11.4,2x)))
#else
7005    format(i3,2x,10e12.5,10(/,5x,10e12.5))
#endif
     enddo
     go to 55
60   continue
  endif
  !
  !
end subroutine elem_lapl
