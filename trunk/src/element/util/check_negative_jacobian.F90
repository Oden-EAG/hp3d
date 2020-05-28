subroutine check_negative_jacobian(Nodesl, Nsize)

  !-------------------------------------------------------------------
  use control
  use data_structure3D
  implicit none
  integer, dimension(1000), intent(inout) :: Nodesl
  integer, intent(inout) :: Nsize
  !-------------------------------------------------------------------
  ! reference and physical coordinates
  real(8), dimension(3)   :: xi,x
  real(8), dimension(3,3) :: dxdxi,dxidx
  ! Gauss points and weights
  real(8), dimension(3,MAX_NINT3) :: xiloc
  real(8), dimension(  MAX_NINT3) :: wxi
  ! shape function
  real(8), dimension(  MAXbrickH) ::  shapH
  real(8), dimension(3,MAXbrickH) :: dshapH
  real(8), dimension(3,MAXbrickH) :: xnod
  integer, dimension(12) :: nedge_orient
  integer, dimension(6)  :: nface_orient
  integer, dimension(19) :: norder
  ! miscellanea
  real(8) :: rjac
  integer :: mdle, iel, i, k, nint, int_back, l, iflag, ndom, ic, nrdofH
  character(len=4) :: type
  !-------------------------------------------------------------------

  write(*,*) 'checking elements negative Jacobians...', NRELES

  int_back = INTEGRATION

  ! loop over active elements
  ic = 0
  do iel=1,NRELES
     mdle = ELEM_ORDER(iel)
     call find_elem_nodes(mdle, norder, nedge_orient,nface_orient)
     call nodcor(mdle, xnod)

     type = NODES(mdle)%type
     ! over integration to detect negative Jacobian when refined
     !
     INTEGRATION = 2
     call set_3Dint(type,norder, nint,xiloc,wxi)
     INTEGRATION = 0
     !
     ! loop over integration points
     do l=1,nint
        xi(1:3) = xiloc(1:3,l)
        select case(EXGEOM)
        case(0)
           call shape3DH(type,xi,norder, &
                         nedge_orient,nface_orient, &
                         nrdofH,shapH,dshapH)
           x    (1:3)    =0.d0
           dxdxi(1:3,1:3)=0.d0
           do k=1,nrdofH
              x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
              do i=1,3
                 dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
              end do
           end do
        case(1)
           call exact_geom(mdle,xi, x,dxdxi)
        end select
        call geom(dxdxi, dxidx,rjac,iflag)

        ! check Jacobian
        if (iflag.ne.0) then
           ic = ic + 1
           Nodesl(ic) = mdle
           write(*,*) 'mdle, type, rjac = ', mdle, NODES(mdle)%type, rjac
           exit
        endif
     enddo
  enddo
  Nsize = ic
  write(*,*) Nsize, ' elements Jacobians are negative '

  INTEGRATION = int_back
endsubroutine check_negative_jacobian
