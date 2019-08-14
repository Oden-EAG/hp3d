!------------------------------------------------------
subroutine volume_hp(Vol)
  !
  use data_structure3D
  use element_data
  !
  implicit none
  !
  !  ...dummy arguments
  real*8,intent(out) :: Vol
  real*8 :: vol_mdle
  integer :: iprint, iel, mdle, i
  !----------------------------------------------------
  !
  iprint = 0; Vol = 0.d0

  !  loop over active elements
  do iel=1,NRELES
     mdle = ELEM_ORDER(iel)
     call volume_hp_mdle(mdle, vol_mdle)
     Vol = Vol + vol_mdle
  enddo
endsubroutine volume_hp

!-----------------------------------------------------------------------      
subroutine volume_hp_mdle(Mdle, Vol)
  use data_structure3D
  use element_data
  implicit none

  integer, intent(in)  :: Mdle
  real*8,  intent(out) :: Vol
  !
  integer :: norder(19), nedge_orient(12), nface_orient(6)
  real*8 :: xnod(3,MAXbrickH), shapH(MAXbrickH),dshapH(3,MAXbrickH)
  !
  character(len=4) :: type
  real*8 :: xiloc(3,MAX_NINT3),wxi(MAX_NINT3)

  ! geometry
  real*8, dimension(3)   :: xi,x
  real*8, dimension(3,3) :: dxdxi,dxidx
  real*8  :: wa,rjac
  integer :: iprint,nint,i,k,l,iflag,nrdofH
  !
  iprint = 0;  
  Vol=0.d0

  !  element order, orientations, and gdofs
  call find_elem_nodes(mdle, norder, nedge_orient,nface_orient)
  call nodcor(mdle,xnod)

  !  set integration points        
  type = NODES(mdle)%type
  call set_3Dint(type,norder, nint,xiloc,wxi)

  !  loop over integration points        
  do l=1,nint
     xi(1:3) = xiloc(1:3,l); wa = wxi(l)

     !  evaluate appropriate shape functions at the point
     call shape3H(type,xi, &
          norder,nedge_orient,nface_orient, &
          nrdofH,shapH,dshapH)

     !  mapping to master
     x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0 
     do k=1,nrdofH
        x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
        do i=1,3
           dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
        enddo
     enddo

     !  evaluate the inverse derivatives and jacobian
     call geom(dxdxi, dxidx,rjac,iflag) 
     if (iflag.ne.0) then
        write(*,*)'volume_hp_mdle: NEGATIVE JACOBIAN Mdle = ',Mdle, rjac
        call pause
     endif

     !  accumulate
     Vol = Vol + wa*rjac
  enddo
  
endsubroutine volume_hp_mdle
