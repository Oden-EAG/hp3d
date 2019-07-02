!> Purpose : calcualte geometry error on element
!! @param[in]  Mdle - middle node number
!!
!! @param[inout] Derr    - error
!! @param[inout] Dnorm   - norm of solution
subroutine geometry_error_element(Mdle, Derr,Dnorm)
  use data_structure3D
  implicit none
  !
  ! ** Arguements
  integer, intent(in)    :: Mdle
  real*8,  intent(inout) :: Dnorm, Derr
  
  ! ** Locals
  ! order and soltuion
  integer, dimension(19)          :: norder
  real*8,  dimension(3,MAXbrickH) :: xnod

  integer :: nedge_orient(12), nface_orient(6)

  ! shape function
  real*8 :: &
       shapH(MAXbrickH),dshapH(3,MAXbrickH)

  ! geometry
  real*8, dimension(3)   :: xi,x,x_geom
  real*8, dimension(3,3) :: dxdxi,dxidx,dxidx_geom,dxdxi_geom

  ! quadrature
  real*8 :: xiloc(3,maxint),wxi(maxint)

  character(len=4) :: type
  integer :: i, j, k, l, nint, nrdofH, iprint, iflag
  real*8  :: wa, weight, rjac

  !---------------------------------------------------------------------
  iprint = 0
  if (iprint.eq.1) then
     write(*,7001) Mdle,NODES(Mdle)%type
7001 format('geometry_error_element: Mdle,type = ',i10,2x,a5)
  endif

  ! determine order of approximation
  call find_order(Mdle, norder)
  call find_orient(Mdle, nedge_orient,nface_orient)

  call nodcor(Mdle,xnod)

  Derr = 0.d0; Dnorm = 0.d0
  !-----------------------------------------------------------------------

  ! set up the element quadrature
  type  = NODES(Mdle)%type
  !
  call set_3Dint(type,norder, nint,xiloc,wxi)
  do l=1,nint
     xi(1:3) = xiloc(1:3,l); wa = wxi(l)
     !
     ! evaluate appropriate shape functions at the point
     call shape3H(type,xi,norder, &
          nedge_orient,nface_orient, &
          nrdofH,shapH,dshapH)
     
     ! mapping to master
     x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0 
     do k=1,nrdofH
        x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
        do i=1,3
           dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
        enddo
     enddo
     
     ! evaluate the inverse derivatives and jacobian
     call geom(dxdxi, dxidx,rjac,iflag) 

     !  .....total weight
     weight = wa*rjac
     
     call exact_geom(Mdle, xi, x_geom, dxdxi_geom) 
     
     if (iprint.eq.1) then
        write(*,7010) l,xi
7010    format('geometry_error: l,xi = ',i3,2x,3f8.3)
        write(*,7011) x
7011    format('approx = ',(3(e12.5,2x),3x))
        write(*,7012) x_geom
7012    format('exact  = ',(3(e12.5,2x),3x))
        call pause
     endif

     do i=1,3
        Dnorm = Dnorm + abs(x(i))**2*weight 
        Derr  = Derr  + abs(x_geom(i)-x(i))**2*weight 
     enddo
  enddo
  
  
end subroutine geometry_error_element
