!---------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!---------------------------------------------------------------------------

module geometry_transformations


contains

!----------------------------------------------------------------------------
!> Purpose : apply geometrical transformation to GMP points to get the
!            effect of translation, rotation, etc to a the geometry defined
!            in the input file.
!
!  IMPORTANT: Only makes sense with geometries defined entirely only
!  in terms of points. I.e., not valid when ellipse, sphere, etc., are
!  present.

!! @param[in] Tras_vector - Vector used to define a translation
!----------------------------------------------------------------------------
!
subroutine geometry_transformation_GMPpoints
  !
  use GMP
  use common_prob_data
  !----------------------------------------------------------------------------
  implicit none
  !----------------------------------------------------------------------------
  ! parameters
!!  double precision, dimension(3), intent(in) :: Trasvector

  ! local variables
  double precision, dimension (3,NRPOINT) :: coordinates
  double precision, dimension (3,3) :: rot_matrix
  double precision, dimension (3,3) :: rot_matrix_x
  double precision, dimension (3,3) :: rot_matrix_y
  double precision, dimension (3,3) :: rot_matrix_z
  double precision :: theta

  integer :: ii
  !----------------------------------------------------------------------------

  ! store point coordinates in local variable
  do ii=1,NRPOINT
     coordinates(1:3,ii) = POINTS(ii)%Rdata(1:3)
  end do

  ! translation
  if (COORD_TRANS_TRAS) then
     write(*,*) 'Translation vector: ', TRAS_VECTOR
     do ii=1,NRPOINT
        coordinates(1:3,ii) = coordinates(1:3,ii) + TRAS_VECTOR
     end do
  endif

  ! Rotation
  if (COORD_TRANS_ROT) then
     write(*,*) 'Rotation angles (degrees): ', ROT_ANGLES
     ! rotation around x axis
     rot_matrix_x=&
          RESHAPE((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))
     rot_matrix_y=rot_matrix_x
     rot_matrix_z=rot_matrix_x

     theta=ROT_ANGLES(1) ! theta in degrees
     theta=theta*PI/180.d0
     rot_matrix_x(1:2,1:2)=&
          RESHAPE((/cos(theta),sin(theta),-sin(theta),cos(theta)/),(/2,2/))

     theta=ROT_ANGLES(2) ! theta in degrees
     theta=theta*PI/180.d0
     rot_matrix_y((/1,3/),(/1,3/))=&
          RESHAPE((/cos(theta),sin(theta),-sin(theta),cos(theta)/),(/2,2/))

     theta=ROT_ANGLES(3) ! theta in degrees
     theta=theta*PI/180.d0
     rot_matrix_z(2:3,2:3)=&
          RESHAPE((/cos(theta),sin(theta),-sin(theta),cos(theta)/),(/2,2/))

     rot_matrix=matmul(matmul(rot_matrix_x,rot_matrix_y),rot_matrix_z)
     do ii=1,NRPOINT
        coordinates(1:3,ii) = matmul(rot_matrix,coordinates(1:3,ii))
     end do
  endif

  ! we return point coordinates after transformation to GMP data structure
  do ii=1,NRPOINT
      POINTS(ii)%Rdata(1:3) = coordinates(1:3,ii)
  end do

end subroutine geometry_transformation_GMPpoints



end module geometry_transformations
