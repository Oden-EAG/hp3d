subroutine poly_geom_util_verify

use data_structure3D_poly
use geometry_polydpg
use connectivity_poly
implicit none
integer :: nrv_f,nrv,jv,iexact_n,iexact_a,iexact_g,iexact_r,iexact_rx
integer, allocatable :: nverts_f(:),nverts(:),nedges(:),nfaces(:)
real*8, allocatable :: xverts_f(:,:),xverts(:,:),rotverts_f_ex(:,:),rotverts_f_calc(:,:)
real*8 :: aux,area_ex,area_calc,rad_ex,rad_calc
real*8, dimension(3) ::fn_ex,fn_calc,xg_ex,xg_calc
real*8, dimension(3,3) :: Qrot_f
real*8, parameter :: tolerance=1.d-12

open(unit=101,file='verify/face_verify', &
       form='formatted',access='sequential',status='unknown')
! read number of vertices of given face
read(101,*) nrv_f
allocate(xverts_f(3,nrv_f),rotverts_f_ex(3,nrv_f),rotverts_f_calc(3,nrv_f))
! read vertex coordinates
do jv=1,nrv_f
  read(101,*) xverts_f(1:3,jv)
enddo
! reads flag to know if exact normal is given in file
read(101,*) iexact_n
!if (iexact_n.eq.1) then
! read exact normal (or a dummy vector)
read(101,*) fn_ex(1:3)
!endif

! reads flag to know if exact area is given in file
read(101,*) iexact_a
!if (iexact_a.eq.1) then
! read exact area (or a dummy scalar)
read(101,*) area_ex
!endif

! reads flag to know if exact centroid is given in file
read(101,*) iexact_g
!if (iexact_g.eq.1) then
! read exact centroid (or a dummy vector)
read(101,*) xg_ex(1:3)
!endif

! reads flag to know if exact radius is given in file
read(101,*) iexact_r
!if (iexact_r.eq.1) then
! read exact radius (or a dummy scalar)
read(101,*) rad_ex
!endif
! reads flag to know if exact ROTATED vertex oordinates are given in file
read(101,*) iexact_rx
! read exact ROTATED vertex coordinates
do jv=1,nrv_f
  read(101,*) rotverts_f_ex(1:3,jv)
enddo
close(101)

 2000 format(3(e12.5))
 2001 format(1(e12.5))
 2002 format('PASSED with absolute tolerance ',1(e10.3))
 2003 format('FAILED with absolute tolerance ',1(e10.3))

write(*,*) 'Calling subroutine face_normal with input xverts_f (just first 3 vertices)'
write(*,*) '[output:normal unit vector]'
call face_normal(xverts_f(:,1:3),fn_calc)


if (iexact_n.eq.1) then
  call norm(fn_calc-fn_ex,aux)
  write(*,*) 'face_normal TEST :'
  if (aux.le.tolerance) then
    write(*,2002) tolerance
    write(*,*) ''
  else
    write(*,2003) tolerance
    call pause
    write (*,*) 'unit normal computed='
    write(*,2000) fn_calc
    write(*,*) ''
  endif
else
  write (*,*) 'NO EXACT VALUE TO COMPARE. Unit normal computed='
  write(*,2000) fn_calc
  write(*,*) ''
endif

write(*,*) 'Calling subroutine face_centroid with input xverts_f and nrv_f'
write(*,*) '[output: face centroid and area]'
call face_centroid(xverts_f,nrv_f,xg_calc,area_calc)


if (iexact_g.eq.1) then
  call norm(xg_calc-xg_ex,aux)
  write(*,*) 'face_centroid TEST 1:'
  if (aux.le.tolerance) then
    write(*,2002) tolerance
    write(*,*) ''
  else
    write(*,2003) tolerance
    write (*,*) 'centroid computed='
    write(*,2000) xg_calc
    write(*,*) ''
    call pause
  endif
else
  write (*,*) 'NO EXACT VALUE TO COMPARE. Centroid computed='
  write(*,2000) xg_calc
  write(*,*) ''
endif

if (iexact_a.eq.1) then
  aux=abs(area_calc-area_ex)
  write(*,*) 'face_centroid TEST 2:'
  if (aux.le.tolerance) then
    write(*,2002) tolerance
    write(*,*) ''
  else
    write(*,2003) tolerance
    write (*,*) 'area computed='
    write(*,2001) area_calc
    write(*,*) ''
    call pause
  endif
else
  write (*,*) 'NO EXACT VALUE TO COMPARE. Area computed='
  write(*,2001) area_calc
  write(*,*) ''
endif

write(*,*) 'Calling subroutine polyhedron_radius with input xverts_f, nrv_f, xg'
write(*,*) '[output: face radius w.r.t. centroid]'
call polyhedron_radius(Xverts_f,Nrv_f,xg_calc,Rad_calc)

if (iexact_r.eq.1) then
  aux=abs(rad_calc-rad_ex)
  write(*,*) 'polyhedron_radius TEST :'
  if (aux.le.tolerance) then
    write(*,2002) tolerance
    write(*,*) ''
  else
    write(*,2003) tolerance
    write (*,*) 'radius computed='
    write(*,2001) rad_calc
    write(*,*) ''
    call pause
  endif
else
  write (*,*) 'NO EXACT VALUE TO COMPARE. Radius computed='
  write(*,2001) rad_calc
  write(*,*) ''
endif


do jv=1,nrv_f
  xverts_f(1:3,jv)=xverts_f(1:3,jv)-xg_calc(1:3)
enddo

write(*,*) 'Calling subroutine face_rotation_to_2D with input xverts_f,nrv_f,fn_calc'
write(*,*) '[output:rotverts_f_calc,Qrot_f]'
call face_rotation_to_2D(Xverts_f,Nrv_f,Fn_calc,rotverts_f_calc,Qrot_f)

if (iexact_n.eq.1) then
  call norm(fn_calc-fn_ex,aux)
  write(*,*) 'face_rotation_to_2D TEST :'
  if (aux.le.tolerance) then
    write(*,2002) tolerance
    write(*,*) ''
  else
    write(*,2003) tolerance
    write (*,*) 'Rotated vertices computed='
    do jv=1,nrv_f
      write(*,2000) rotverts_f_calc(1:3,jv)
    enddo
    write(*,*) ''
    call pause
  endif
else
write (*,*) 'NO EXACT VALUE TO COMPARE. Rotated vertices computed='
do jv=1,nrv_f
  write(*,2000) rotverts_f_calc(1:3,jv)
enddo
write(*,*) ''
endif



end subroutine poly_geom_util_verify