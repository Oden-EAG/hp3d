!--------------------------------------------------------------------------
subroutine flatten_tet(Ntet)
  !
  use GMP , only : TETRAS , TRIANGLES , CURVES
  !
  implicit none
  integer,intent(in) :: Ntet
  integer :: i,nt,nvoid,nc,j,lab
  !
  !--------------------------------------------------------------------------
  !
  ! loop over faces
  do i=1,4
     call decode(TETRAS(Ntet)%FigNo(i), nt,nvoid)
     if (TRIANGLES(nt)%Type.eq.'G1RecTri') then
        call flatten_trian(nt)
     end if
  enddo
  !
  ! loop over edges
  do i=1,6
     nc=iabs(TETRAS(Ntet)%EdgeNo(i))
     if (CURVES(nc)%Type.eq.'5Bezier') then
        call straighten_curve(nc)
        !  ...loop over connected figures
        do j=1,CURVES(nc)%NrFig
           call decode(iabs(CURVES(nc)%FigNo(i)), nt,lab)
           !  ...if connected figure is a reconstructed triangle
           if (lab.eq.1) then
              if (TRIANGLES(nt)%Type.eq.'G1RecTri') then
                 call flatten_trian(nt)
              end if
           endif
        enddo
     end if
  enddo
  !
  !
end subroutine flatten_tet
!
!
!
!--------------------------------------------------------------------------
subroutine flatten_trian(Nt)
  !
  use GMP , only : TRIANGLES
  !
  implicit none
  integer,intent(in)   :: Nt
  real(8),dimension(3) :: temp,v1,v2,v3
  integer              :: i,nc,ie,neig,j
  real(8)              :: x,y
  integer, parameter   :: deg=7
  integer, external    :: bijec
  !
  !--------------------------------------------------------------------------
  !  STEP 1 : straighten edges                                              |
  !--------------------------------------------------------------------------
  !
  do i=1,3
     nc=iabs(TRIANGLES(Nt)%EdgeNo(i))
     call straighten_curve(nc)
  enddo
  !
  !--------------------------------------------------------------------------
  !  STEP 2 : modify control points of adjacent triangles                   |
  !--------------------------------------------------------------------------
  !
  do i=1,3
     nc=iabs(TRIANGLES(Nt)%EdgeNo(i))
     call curve2trian(nc,Nt, neig,ie)
     call modify_G1trian(neig,ie)
  enddo
  !
  !--------------------------------------------------------------------------
  !  STEP 3 : modify control points of triangle                             |
  !--------------------------------------------------------------------------
  !
  TRIANGLES(Nt)%Rdata=0.d0
  call trian2verts(Nt, v1,v2,v3)
  do j=0,deg
     do i=0,(deg-j)
        x=float(i)/float(deg)
        y=float(j)/float(deg)
        temp(1:3)=(1.d0-x-y)*v1(1:3) + x*v2(1:3) + y*v3(1:3)
        TRIANGLES(Nt)%Rdata(bijec(i,j):bijec(i,j)+2)=temp(1:3)
     enddo
  enddo
  !
  !
end subroutine flatten_trian
!--------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------
subroutine straighten_curve(Nc)
  !
  use GMP , only : CURVES,POINTS
  !
  implicit none
  integer, intent(in)  :: Nc
  integer              :: deg,np,i
  real(8),dimension(3) :: v1,v2,temp
  real(8)              :: x
  !
  !--------------------------------------------------------------------------
  !  ...determine degree, check curve type
  select case(CURVES(Nc)%Type)
  case('5Bezier') ; deg=5
  case('6Bezier') ; deg=6
  case('7Bezier') ; deg=7
  case default
     write(*,7000) Nc,CURVES(Nc)%Type
7000 format(' straighten_curve: unsupported! Nc,Type = ',i7,2x,a10)
     return
  end select
  !
  !  ...determine curve endpoints
  np=CURVES(Nc)%EndPoNo(1)
  v1(1:3)=POINTS(np)%Rdata(1:3)
  np=CURVES(Nc)%EndPoNo(2)
  v2(1:3)=POINTS(np)%Rdata(1:3)
  !
  !  ...redefine control points
  CURVES(Nc)%Rdata=0.d0
  do i=0,deg
     x=float(i)/float(deg)
     temp(1:3)=(1.d0-x)*v1(1:3) + x*v2(1:3)
     CURVES(Nc)%Rdata(3*i:3*i+2)=temp(1:3)
  enddo
  !
  !
end subroutine straighten_curve
!--------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------
subroutine curve2trian(Nc,Nt, Neig,Ie)
  !
  use GMP , only : CURVES , TRIANGLES
  !
  implicit none
  integer, intent(in ) :: Nc,Nt
  integer, intent(out) :: Neig,Ie
  integer              :: icount,i,lab
  !
  !--------------------------------------------------------------------------
  !
  !  ...loop over connected figures
  icount=0
  do i=1,CURVES(Nc)%NrFig
     call decode(iabs(CURVES(Nc)%FigNo(i)), Neig,lab)
     !  ...if attached figure is a triangle
     if (lab.eq.1) then
        !  ...increment counter if a neighbor was found
        if ((Neig.ne.Nt).and.(TRIANGLES(Neig)%Type.eq.'G1RecTri')) then
           icount=1 ; exit
        endif
     endif
  enddo
  !
  !  ...check that exactly 1 neighbor was found
  if (icount.ne.1) then
     write(*,7000)Nc
7000 format(' curve2trian: no G1 neighbor found for Nc = ',i7)
     stop
  endif
  !
  !  ...determine local number for the edge
  icount=0
  do i=1,3
     if (iabs(TRIANGLES(Neig)%EdgeNo(i)).eq.Nc) then
        Ie=i ; icount=1 ; exit
     endif
  enddo
  !
  !  ...check
  if (icount.eq.0) then
     write(*,7001)Nc
7001 format(' curve2trian: did not find local number for Nc = ',i7)
     stop
  endif
  !
  !
end subroutine curve2trian
!--------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------
subroutine modify_G1trian(Nt,Ie)
  !
  use GMP , only : TRIANGLES
  !
  implicit none
  integer,intent(in)   :: Nt,Ie
  real(8),dimension(3) :: temp,v1,v2,v3
  real(8)              :: x,y
  integer              :: i,j,k,iflag
  integer, parameter   :: deg=7
  integer, external    :: bijec
  !
  !--------------------------------------------------------------------------
  !
  !  ...determine vertices
  call trian2verts(Nt, v1,v2,v3)
  !  ...loop over all control points
  do j=0,deg
     do i=0,(deg-j)
        k=deg-i-j
        x=float(i)/float(deg)
        y=float(j)/float(deg)
        temp(1:3)=(1.d0-x-y)*v1(1:3) + x*v2(1:3) + y*v3(1:3)
        select case(Ie)
        case(1) ; iflag=j
        case(2) ; iflag=k
        case(3) ; iflag=i
        end select

        !  .......if on the edge (btw, we are ALWAYS on the edge...)
        if (iflag.eq.0) then
           TRIANGLES(Nt)%Rdata(bijec(i,j):bijec(i,j)+2)=temp(1:3)
        endif
     enddo
  enddo

end subroutine modify_G1trian
!--------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------
subroutine trian2verts(Nt, V1,V2,V3)
  !
  use GMP , only : POINTS, TRIANGLES
  !
  implicit none
  integer,            intent(in )  :: Nt
  real(8),dimension(3),intent(out) :: V1,V2,V3
  integer                          :: i,np
  !
  !--------------------------------------------------------------------------
  !
  do i=1,3
     np=TRIANGLES(Nt)%VertNo(i)
     select case(i)
     case(1) ; V1(1:3)=POINTS(np)%Rdata(1:3)
     case(2) ; V2(1:3)=POINTS(np)%Rdata(1:3)
     case(3) ; V3(1:3)=POINTS(np)%Rdata(1:3)
     endselect
  enddo
  !
  !
end subroutine trian2verts
