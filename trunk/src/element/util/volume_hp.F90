!
!-----------------------------------------------------------------------
subroutine volume_hp(Vol)
!
   use data_structure3D
   use element_data
!
   implicit none
!
   real(8), intent(out) :: Vol
!
   real(8) :: vol_mdle
   integer :: iel, mdle
!
!------------------------------------------------------
!
   Vol = 0.d0
!
!..loop over active elements
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      call volume_hp_mdle(mdle, vol_mdle)
      Vol = Vol + vol_mdle
   enddo
!
end subroutine volume_hp
!
!-----------------------------------------------------------------------
!
subroutine volume_hp_mdle(Mdle, Vol)
!
   use data_structure3D
   use element_data
   implicit none
!
   integer, intent(in)  :: Mdle
   real(8), intent(out) :: Vol
!
   integer :: norder(19),nedge_orient(12),nface_orient(6)
   real(8) :: xnod(3,MAXbrickH),shapH(MAXbrickH),gradH(3,MAXbrickH)
!
   character(len=4) :: type
   real(8) :: xiloc(3,MAX_NINT3),wxi(MAX_NINT3)
!
!..geometry
   real(8), dimension(3)   :: xi,x
   real(8), dimension(3,3) :: dxdxi,dxidx
   real(8) :: wa,rjac
   integer :: nint,i,k,l,iflag,nrdofH
!
!------------------------------------------------------
!
   Vol = 0.d0
!
!..element order, orientations, and gdofs
   call find_elem_nodes(mdle, norder,nedge_orient,nface_orient)
   call nodcor(mdle,xnod)
!
!..set integration points
   type = NODES(mdle)%type
   call set_3Dint(type,norder, nint,xiloc,wxi)
!
!..loop over integration points
   do l=1,nint
      xi(1:3) = xiloc(1:3,l); wa = wxi(l)
!
!  ...evaluate appropriate shape functions at the point
      call shape3DH(type,xi,norder,nedge_orient,nface_orient, &
                    nrdofH,shapH,gradH)
!
!  ...geometry mapping
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, &
                  x,dxdxi,dxidx,rjac,iflag)
!
!  ...accumulate
      Vol = Vol + wa*rjac
!
   enddo
!
end subroutine volume_hp_mdle

