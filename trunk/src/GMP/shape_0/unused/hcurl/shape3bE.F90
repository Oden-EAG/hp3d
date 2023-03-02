!----------------------------------------------------------------------------
!> Purpose : H(curl) shape functions for Hexahedron (Nedelec first family)
!!
!! @param[in] Xi        - master element coordinates
!! @param[in] Norder    - order of approximation for the nodes
!! @param[in] Ne_orient - edge orientations
!! @param[in] Nf_orient - face orientations

!! @param[out]            NrdofE    - number of element dof
!! @param[out]            ShapE     - values of shape functions
!! @param[out]            CurlE     - values of derivatives of the shape function
!----------------------------------------------------------------------------
!

subroutine shape3bE(Xi,Norder,Nedge_orient,Nface_orient, NrdofE,ShapE,CurlE)

  use element_data , only : NFAXES,IJKV,IXIEDGE,IBLENDE,NBLENDE, &
                                        IXIFACE,IBLENDF,NBLENDF
  use parameters   , only : MAXP,MAXbrickE

  implicit none
  !-------------------------------------------------------------------
  ! subroutine  arguments
  real(8),dimension(3)         ,intent(in ) :: Xi
  integer,dimension(19)        ,intent(in ) :: Norder
  integer,dimension(12)        ,intent(in ) :: Nedge_orient
  integer,dimension(6)         ,intent(in ) :: Nface_orient
  integer                      ,intent(out) :: NrdofE
  real(8),dimension(3,MAXbrickE),intent(out) :: ShapE
  real(8),dimension(3,MAXbrickE),intent(out) :: CurlE
  !-------------------------------------------------------------------
  ! local variables
  !-------------------------------------------------------------------

  logical, save :: initialized_curl_signs = .FALSE.
  integer, save, dimension(3,3) :: sign_curl_matrix = 0

  double precision :: shapE1D(3,MAXP,0:1) ! chequear dimension MAXP
  double precision :: shapH1D(3,MAXP,0:1),dshapH1D(3,MAXP,0:1) ! chequear dimension MAXP

  double precision :: shapE1Daux(MAXP), &
       shapH1D_1(MAXP), shapH1D_2(MAXP), dshapH1D_1(MAXP), dshapH1D_2(MAXP)

  double precision ::  vshap(2,3),dvshap(2,3) ! por limpiar ....

  integer :: direction, orientation, other_direction, &
       orient_direction, orient_other_direction, &
       other_direction_1,other_direction_2, &
       cross_direction, cross_direction_vec(1)
  integer :: nord_interior(3),nordh,nordv,nord1,nord2,nord3, &
       nord_face(2), nord_aux1, nord_aux2

  integer :: nvoid, ie, if, ibl1, ibl2, ixi, ixi1, ixi2, nv1, nv2, &
       ibl, nv, face_direction,  &
       k, kk, cont_diraux, cont_diraux1, cont_diraux2, &
       j

  integer :: nrdofEe, nrdofExy(2), nrdofExyz(3), nrdofH1xy(2), nrdofH1xyz(3)

  !-----------------------------------------------------------------------------


  integer, save :: iflag=0

  ! shape Hcurl functions and "derivatives" and others FOR JASON ROUTINE
  real(8),dimension(3,MAXbrickE) :: shapEjason
  real(8),dimension(3,MAXbrickE) :: curlEjason
  integer :: nrdofEjason, fakenrdofP
  real(8),dimension(MAXbrickE)   :: fakeShapP
  real(8),dimension(3,MAXbrickE) :: fakeGradP



  iflag=0
  if (iflag == 1) then
     write(*,*) 'JASON SHAPE FUNCTIONS'
     call shape3bEM2(xi,Norder, fakeNrdofP,fakeShapP,fakeGradP, &
          NrdofEjason,ShapEjason,CurlEjason)
     nrdofE=nrdofEjason
     ShapE(1:3,1:nrdofE)=ShapEjason(1:3,1:nrdofE)
     CurlE(1:3,1:nrdofE)=CurlEjason(1:3,1:nrdofE)
     return
!!$  else
!!$     write(*,*) 'LUISE SHAPE FUNCTIONS'
  endif


  !-----------------------------------------------------------------------------
  !     decode middle node order
  call decode(Norder(19), nordh, nord3)
  call decode(nordh, nord1,nord2)
  nord_interior=(/ nord1, nord2, nord3 /)
  if (ANY(nord_interior > MAXP) ) then
     write(*,*) 'ERROR shape3bE: interior order higher than MAXP'
     stop
  endif
  !
  ! evaluate 1D H1 shape functions and 1D E shape functions for three
  ! axis (directions) and orientations for maximums orders (minimum
  ! rule assures orders on edges and faces are lower than its
  ! respective interior directions)
  orientation_loop_0: do orientation=0,1

     direction_loop_0: do direction=1,3

        select case(orientation)
        case(0)
           call shape1_Jason (Xi(direction),nord_interior(direction), nvoid,&
             shapH1D(direction,:,0),dshapH1D(direction,:,0))
           call shape1E_Jason(Xi(direction),nord_interior(direction), nvoid, &
             shapE1D(direction,:,0))
        case(1)
           call shape1_Jason (1-Xi(direction),nord_interior(direction), nvoid,&
             shapH1D(direction,:,1),dshapH1D(direction,:,1))
           dshapH1D(direction,:,1) = - dshapH1D(direction,:,1)
           call shape1E_Jason(1-Xi(direction),nord_interior(direction), nvoid, &
                shapE1D(direction,:,1))
           shapE1D(direction,:,1) = - shapE1D(direction,:,1)
        end select

     end do direction_loop_0

  enddo orientation_loop_0

  !-----------------------------------------------------------------------------

  ShapE = 0.d0  ! better initialize only up to the required
                ! order: ShapE(:,1:"orden")
  CurlE = 0.d0  ! idem


  ! initialization (if not done in a previous call) of matrix with
  ! "curl_component_signs" used in computation of CurlE
  if (.not.initialized_curl_signs) then
     call sign_curl_component_initialization(3,Sign_curl_matrix)
     initialized_curl_signs = .TRUE.
  endif


  !-----------------------------------------------------------------------------


  !    GATTO:
  !  ...calculate 1D linear shape functions (used as blending functions)
  ! TO CLEAN
  do ixi=1,3 ! direction x,y,z

     vshap(1,ixi) = shapH1D(ixi,1,0);
     dvshap(1,ixi)=dshapH1D(ixi,1,0)
     vshap(2,ixi) = shapH1D(ixi,1,1);
     dvshap(2,ixi)=dshapH1D(ixi,1,1)

!!$     vshap(1,ixi) = 1.d0 - Xi(ixi); dvshap(1,ixi) = -1.d0
!!$     vshap(2,ixi) =        Xi(ixi); dvshap(2,ixi) =  1.d0
  enddo
  !
  k=0
  !
  !--------------------------------------------------------------
  ! calculate edge shape functions
  do ie=1,12
     ! determines if edge is along x,  y or z
     ixi = IXIEDGE(ie)
     !  linear blending functions along directions x, y, or z
     ibl1 = IBLENDE(1,ie);
     ibl2 = IBLENDE(2,ie)
     ! determines in what extreme the linear blending function is one
     ! or zero
     nv1 = NBLENDE(1,ie);
     nv2 = NBLENDE(2,ie)

     nrdofEe=Norder(ie)
     if (Norder(ie) > MAXP) then
        write(*,*) 'ERROR shape3bE: edge order higher than MAXP'
        stop
     endif
     shapE1Daux=shapE1D(ixi,:,Nedge_orient(ie))
     do j=1,nrdofEe ! nrdofEe = Norder(ie)
        k=k+1

        ShapE(ixi,k) = shapE1Daux(j)*vshap(nv1,ibl1)*vshap(nv2,ibl2)

        CurlE(ibl1,k) = shapE1Daux(j)*vshap(nv1,ibl1)*dvshap(nv2,ibl2) &
             * sign_curl_matrix(ixi,ibl2)
        CurlE(ibl2,k) = shapE1Daux(j)*dvshap(nv1,ibl1)*vshap(nv2,ibl2) &
             * sign_curl_matrix(ixi,ibl1)
     enddo
  enddo


  !--------------------------------------------------------------
  ! face node shape functions

  !  ...calculate face shape functions
  do if=1,6
     ! determines directions of the  two local axis, i.e., along x, y, or z
     ixi1 = IXIFACE(1,if)
     ixi2 = IXIFACE(2,if)
!     ixiface_aux=IXIFACE(1:2,if)  !(/ ixi1, ixi2 /) ! CLEAN
     !  linear blending function along directions x, y, or z
     ibl = IBLENDF(if);
     ! determines in what extreme the linear blending function is one
     ! or zero
     nv = NBLENDF(if)
     call decode(Norder(12+if), nordh,nordv)
     if (ANY( (/ nordh, nordv  /) > MAXP) ) then
        write(*,*) 'ERROR shape3bE: face order higher than MAXP'
        stop
     endif

!!$     if ( (if==1) .and. (Nface_orient(if) == 3)) then
!!$        write(*,*) 'YA ESTAMOS AQUI'
!!$     endif


!!!!! QUIZAS ESTO VAYA DENTRO DE BUCLE FACE_DIRECTION_LOOP
!!!! O NO HAGA FALTA!!!
     select case(NFAXES(3,Nface_orient(if)))
        !
        !  .....the face axes have NOT been reversed
     case(0)
        nord_face=(/ nordh, nordv /)
        !  ....the face axes HAVE been reversed
     case(1)
        nord_face=(/ nordv, nordh /)
     case default
        write(*,*) 'shape3bE. ERROR: orientation different from 0,1'
        stop
     end select

     face_direction_loop_pre: do face_direction=1,2
        nrdofExy(face_direction)=nord_face(face_direction)
        nrdofH1xy(face_direction)=nord_face(face_direction)-1
     enddo face_direction_loop_pre


     face_direction_loop: do face_direction=1,2
        direction=IXIFACE(face_direction,if)
        orient_direction=NFAXES(direction,Nface_orient(if))

        other_direction=IXIFACE(next_direction_number(face_direction,2),if)
        orient_other_direction=NFAXES(other_direction,Nface_orient(if))

        select case(NFAXES(3,Nface_orient(if)))
           !
           !  .....the face axes have NOT been reversed
        case(0)

           do cont_diraux=1,nrdofH1xy(next_direction_number(face_direction,2))
              do kk=1,nrdofExy(face_direction)

                 k=k+1

                 ShapE(direction,k) = shapE1D(direction,kk,orient_direction)* &
                      shapH1D(other_direction,cont_diraux+2,orient_other_direction)*vshap(nv,ibl)

                 CurlE(other_direction,k) = shapE1D(direction,kk,orient_direction)* &
                      shapH1D(other_direction,cont_diraux+2,orient_other_direction)*dvshap(nv,ibl) &
                      * sign_curl_matrix(direction,ibl)
                 CurlE(ibl,k) = shapE1D(direction,kk,orient_direction)* &
                      dshapH1D(other_direction,cont_diraux+2,orient_other_direction)*vshap(nv,ibl) &
                      * sign_curl_matrix(direction,other_direction)

              enddo
           enddo

        case(1)
           ! the face axes HAVE been reversed =>  order of loops gets REVERSED

           !¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡ CAMBIAR variable k N VEZ DE COPIAR Y PEGAR LOOP

           do cont_diraux=1,nrdofH1xy(next_direction_number(face_direction,2))
              do kk=1,nrdofExy(face_direction)

                 k=k+1

                 ShapE(other_direction,k) = shapE1D(other_direction,kk,orient_other_direction)* &
                      shapH1D(direction,cont_diraux+2,orient_direction)*vshap(nv,ibl)

                 CurlE(direction,k) = shapE1D(other_direction,kk,orient_other_direction)* &
                      shapH1D(direction,cont_diraux+2,orient_direction)*dvshap(nv,ibl) &
                      * sign_curl_matrix(other_direction,ibl)
                 CurlE(ibl,k) = shapE1D(other_direction,kk,orient_other_direction)* &
                      dshapH1D(direction,cont_diraux+2,orient_direction)*vshap(nv,ibl) &
                      * sign_curl_matrix(other_direction,direction)
              enddo
           enddo

        case default
           write(*,*) 'shape3bE. ERROR: orientation different from 0,1'
           stop
        end select

     enddo face_direction_loop

  enddo ! loop over faces

  !--------------------------------------------------------------
  ! calculate middle node shape functions

  direction_loop_pre: do direction=1,3
     nrdofExyz(direction)=nord_interior(direction)
     nrdofH1xyz(direction)=nord_interior(direction)-1
  end do direction_loop_pre

  direction_loop_interior: do direction=1,3

!!$     ! edges on directions x,y,z are 1,4,12, respectively
!!$     ibl_interior=IBLENDE(:,(/1, 4, 12/))
!!$     other_direction_1=ibl_interior(1,direction)
!!$     other_direction_2=ibl_interior(2,direction)

     ! Another way to compute the "other" two directions
     other_direction_1=next_direction_number(direction,3)
     other_direction_2=next_direction_number(direction+1,3)


     ! CAMBIAR FILAS POR COLUMNAS EN TODOS LOS ARRAYS DE RANK>=2
     nord_aux1=nord_interior(other_direction_1)
     shapH1D_1=shapH1D(other_direction_1,:,0)
     dshapH1D_1=dshapH1D(other_direction_1,:,0)

     nord_aux2=nord_interior(other_direction_2)
     shapH1D_2 = shapH1D(other_direction_2,:,0)
     dshapH1D_2=dshapH1D(other_direction_2,:,0)


     do kk=1,nrdofExyz(direction) ! = nord(direction)
        do cont_diraux1=1,nrdofH1xyz(other_direction_1)
           do cont_diraux2=1,nrdofH1xyz(other_direction_2)
              k=k+1

              ShapE(direction,k) = shapE1D(direction,kk,0)* &
                   shapH1D_1(cont_diraux1+2)*shapH1D_2(cont_diraux2+2)

              CurlE(other_direction_1,k) = shapE1D(direction,kk,0)* &
                   shapH1D_1(cont_diraux1+2)*dshapH1D_2(cont_diraux2+2) &
                   * sign_curl_matrix(direction,other_direction_2)

              CurlE(other_direction_2,k) = shapE1D(direction,kk,0)* &
                   dshapH1D_1(cont_diraux1+2)*shapH1D_2(cont_diraux2+2) &
                   * sign_curl_matrix(direction,other_direction_1)

           enddo
        enddo
     enddo

  enddo direction_loop_interior

  !--------------------------------------------------------------

  ! count of total number of dof
  NrdofE=k

  !--------------------------------------------------------------

CONTAINS

  !----------------------------------------------------------------------------
  !> Purpose : compute the sign (1 or -1) of terms dE_i/dx_j that
  !appear on calculation of 3D curl
  !
  ! Example: Assuming 3D Cartesian axis numerated such as: x->x_1,
  ! y->x_2, z-> x_3
  !
  !  The curl produced by the first component E_1 involves terms
  !  dE_1/dx_2 and dE_1/dx3; actually,
  !  curlE(1)=0
  !  curlE(2) = + dE_1/dx_3
  !  curlE(3) = - dE_1/dx_2
  !
  ! Then, this function computes the sign (+1,-1) above, i.e.,
  ! curlE(2)=sign_curl_component(1,3) dE_1/dx_3
  ! curlE(3)=sign_curl_component(1,2) dE_1/dx_2
  !
  ! Note: Same is valid for E_2 and E_3 components, for instance,
  !       curlE(1)=sign_curl_component(2,3) dE_2/dx_3
  !
  !! @param[in] direction1              - intepreted as 1->x, 2->y, 3->z
  !! @param[in] direction2              - idem

  !! @param[out] sign_curl_component    - value  +1 or -1
  !----------------------------------------------------------------------------
  !
  function sign_curl_component(Direction1,Direction2)

    use cross_product_module
    implicit none
    integer, intent(in) :: Direction1,Direction2
    integer             :: sign_curl_component
    !----------------------------------------------

    integer, dimension(3) :: vec1,vec2,vec_result
    logical, dimension(3) :: mask

    !----------------------------------------------

    if ((Direction1 > 3) .or. (Direction1 < 0) ) then
       write(*,*) 'function sign_curl_component. ERROR: ', &
            'Input parameter direction1 out of limits'
       stop
    endif
    if ((Direction2 > 3) .or. (Direction2 < 0) ) then
       write(*,*) 'function sign_curl_component. ERROR: ', &
            'Input parameter direction2 out of limits'
       stop
    endif

    vec1=0
    vec1(Direction1)=1

    vec2=0
    vec2(Direction2)=1

    call cross_product3D_int(vec1,vec2,vec_result)
    mask = vec_result.eq.(/0,0,0/)

    ! check that vec_result is a vector with only one non-null component
    ! and its value is 1 or -1
    if ( (count(mask) /= 2).or.(abs(sum(vec_result)) /= 1) ) then
       write(*,*) 'function sign_curl_component. ERROR: ', &
            'wrong output from cross_product3D_int'
       stop
    endif

    ! we detect which component != 0
    cross_direction_vec=maxloc(abs(vec_result))
    cross_direction=cross_direction_vec(1)

    ! the sign is the negative sign of the non-null component of vec_result
    sign_curl_component = - sign(1,vec_result(cross_direction))

  end  function sign_curl_component


  !----------------------------------------------------------------------------
  !> Purpose : Precompute matrix containing the "sign_curl_components"
  !(see function "sign_curl_component") for a Cartesian system of
  !coordinates of dimension Ndim
  !
  ! Note: See  function "sign_curl_component"
  !
  !! @param[in] Ndim                - Dimension of the space, e.g., 3->3D)

  !! @param[out] Sign_curl_matrix        - Matrix with ij component =
  !!                                          sign_curl_component(i,j)
  !----------------------------------------------------------------------------
  !
  subroutine sign_curl_component_initialization(Ndim,Sign_curl_matrix)

    implicit none
    integer :: Ndim
    integer :: Sign_curl_matrix(Ndim,Ndim)

    integer                  :: ii, jj

    do ii=1,Ndim
       do jj=ii+1,Ndim
           Sign_curl_matrix(ii,jj) = &
                sign_curl_component(ii,jj)
           Sign_curl_matrix(jj,ii) = - Sign_curl_matrix(ii,jj)
        enddo
     enddo

  end subroutine sign_curl_component_initialization

  !----------------------------------------------------------------------------
  !> Purpose : compute the next "Direction number" in a Ndim dimensional space
  !
  ! Example: Assuming 3D Cartesian axis numerated such as: x->1, y->2, z-> 3
  !   next_direction_number(1->"x") = 2->"y"
  !   next_direction_number(2->"y") = 3->"z"
  !   next_direction_number(3->"z") = 1->"x"
  !   next_direction_number(4->"x") = 2->"y"
  !!
  !! @param[in] Direction_number    - interpreted as 1->x, 2->y, 3->z
  !! @param[in] Ndim                - dimensional of the space (=3 for 3D)

  !! @param[out]                    - next_direction_number
  !----------------------------------------------------------------------------
  !
  function next_direction_number(Direction_number,Ndim)

    implicit none
    integer, intent(in) :: Direction_number,Ndim
    integer             :: Next_direction_number

    if ( (Direction_number <= 0) .or. (Ndim <= 0) ) then
       write(*,*) 'function next_direction. ERROR: input ', &
            'parameter must be positive'
       stop
    endif
    Next_direction_number=mod(Direction_number+1,Ndim)
    if (Next_direction_number == 0) then
       Next_direction_number=Ndim
    endif

  end function next_direction_number


end subroutine shape3bE


