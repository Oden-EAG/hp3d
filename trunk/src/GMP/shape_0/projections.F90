!----------------------------------------------------------------------
!> Purpose : triangle to edge projection
!
!  @param[in ] Iv1,2   - vertices defining the edge
!  @param[in ] Vshape  - master triangle vertex shape functions,
!                        computed at the point of interest
!  @param[in ] Dvshape - derivatives of vertex shape functions
!  @param[out] T       - edge coordinate
!  @param[out] Dtdxi   - derivative of edge coordinate
!
!  @revision Nov 12
!----------------------------------------------------------------------
   subroutine proj_t2e(Iv1,Iv2,Vshape,Dvshape, T,Dtdxi)
!
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2
      real(8),dimension(  3),intent(in ) :: Vshape
      real(8),dimension(2,3),intent(in ) :: Dvshape
      real(8)               ,intent(out) :: T
      real(8),dimension(2)  ,intent(out) :: Dtdxi
!----------------------------------------------------------------------
!
      T = (Vshape(Iv2)-Vshape(Iv1)+1.d0)*.5d0
      Dtdxi(1:2) = (Dvshape(1:2,Iv2)-Dvshape(1:2,Iv1))*.5d0
!
   end subroutine proj_t2e
!
!
!----------------------------------------------------------------------
!> Purpose : quad to edge projection
!
!  @param[in ] Iv1,2   - vertices defining the edge
!  @param[in ] Vshape  - master quad vertex shape functions, computed
!                        at the point of interest
!  @param[in ] Dvshape - derivatives of vertex shape functions
!  @param[out] T       - edge coordinate
!  @param[out] Dtdxi   - derivative of edge coordinate
!
!  @revision Nov 12
!----------------------------------------------------------------------
   subroutine proj_r2e(Iv1,Iv2,Vshape,Dvshape, T,Dtdxi)
!
      use element_data , only : QUADR_EDGE_TO_VERT
!
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2
      real(8),dimension(  4),intent(in ) :: Vshape
      real(8),dimension(2,4),intent(in ) :: Dvshape
      real(8)               ,intent(out) :: T
      real(8),dimension(2)  ,intent(out) :: Dtdxi
!
      integer :: iedge,i,iv
!----------------------------------------------------------------------
!
      call quad_aux(Iv1,Iv2, iedge)
!
!  ...loop over edge vertices
      T=0.d0 ; Dtdxi(1:2)=0.d0
      do i=1,2
        iv=QUADR_EDGE_TO_VERT(i,iedge)
        T          = T          +  Vshape(    iv)
        Dtdxi(1:2) = Dtdxi(1:2) + Dvshape(1:2,iv)
      enddo
!
   end subroutine proj_r2e
!
!
!----------------------------------------------------------------------
!> Purpose : tetrahedroN to edge projection
!
!  @param[in ] Iv1,2   - vertices defining the edge
!  @param[in ] Vshape  - master tet vertex shape functions, computed
!                        at the point of interest
!  @param[in ] Dvshape - derivatives of vertex shape functions
!  @param[out] T       - edge coordinate
!  @param[out] Dtdxi   - derivative of edge coordinate
!
!  @revision Nov 12
!----------------------------------------------------------------------
   subroutine proj_n2e(Iv1,Iv2,Vshape,Dvshape, T,Dtdxi)
!
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2
      real(8),dimension(  4),intent(in ) :: Vshape
      real(8),dimension(3,4),intent(in ) :: Dvshape
      real(8)               ,intent(out) :: T
      real(8),dimension(3)  ,intent(out) :: Dtdxi
!----------------------------------------------------------------------
!
      T         =( Vshape(    Iv2)- Vshape(    Iv1)+1.d0)*0.5d0
      Dtdxi(1:3)=(Dvshape(1:3,Iv2)-Dvshape(1:3,Iv1)     )*0.5d0
!
   end subroutine proj_n2e
!
!
!----------------------------------------------------------------------
!> Purpose : tetrahedroN to face projection
!
!  @param[in ] Iv1,2,3   - vertices defining the face
!  @param[in ] Vshape    - master tet vertex shape functions, computed
!                          at the point of interest
!  @param[in ] Dvshape   - derivatives of vertex shape functions
!  @param[out] T         - face coordinate
!  @param[out] Dtdxi     - derivative of face coordinate
!
!  @revision Nov 12
!----------------------------------------------------------------------
   subroutine proj_n2f(Iv1,Iv2,Iv3,Vshape,Dvshape, T,Dtdxi)
!
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2,Iv3
      real(8),dimension(  4),intent(in ) :: Vshape
      real(8),dimension(3,4),intent(in ) :: Dvshape
      real(8),dimension(2)  ,intent(out) :: T
      real(8),dimension(2,3),intent(out) :: Dtdxi
!
      real(8)              :: tsum
      real(8),dimension(3) :: dtsum
!----------------------------------------------------------------------
!
      tsum  =  1.d0 - Vshape(Iv1) - Vshape(Iv2) - Vshape(Iv3)
      dtsum(1:3) = -Dvshape(1:3,Iv1)-Dvshape(1:3,Iv2)-Dvshape(1:3,Iv3)
      T(1) = Vshape(Iv2) + tsum/3.d0
      Dtdxi(1,1:3) = Dvshape(1:3,Iv2) + dtsum(1:3)/3.d0
      T(2) = Vshape(Iv3) + tsum/3.d0
      Dtdxi(2,1:3) = Dvshape(1:3,Iv3) + dtsum(1:3)/3.d0
!
!!      T(1) = Vshap(Iv2); Dtdxi(1,1:3) = Dvshap(1:3,Iv2)
!!      T(2) = Vshap(Iv3); Dtdxi(2,1:3) = Dvshap(1:3,Iv3)
!!      do iv=1,4
!!        if ((iv.ne.Iv1).and.(iv.ne.Iv2).and.(iv.ne.Iv3)) then
!!        T(1) = T(1) + Vshap(iv)/3.d0
!!          Dtdxi(1,1:3) = Dtdxi(1,1:3) + Dvshap(1:3,iv)/3.d0
!!          T(2) = T(2) + Vshap(iv)/3.d0
!!          Dtdxi(2,1:3) = Dtdxi(2,1:3) + Dvshap(1:3,iv)/3.d0
!!        endif
!!      enddo
!!      if (iprint.eq.1) then
!!        write(*,7001) Xi(1:3),T(1:2)
!! 7001   format('proj_n2f: Xi = ',3e12.5,' T = ',2e12.5)
!!        call pause
!!      endif
!
   end subroutine proj_n2f
!
!
!----------------------------------------------------------------------
!
!   routine name       - proj_d2e
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - define pyramiD to edge projection used
!                        for edge to pyramiD extensions
!
!   arguments :
!     in:
!             Xi       - master element coordinates
!             Ie       - edge number
!             Vshap    - vertex shape functions
!             Dvshap   - derivatives of vertex shape functions wrt
!                        master pyramid coordinates
!     out:
!             T        - local edge coordinate
!             Dtdxi    - derivatives of the local edge coordinate
!                        wrt master tetrahedron coordinates
!
!----------------------------------------------------------------------
!
   subroutine proj_d2e(Xi,Ie,Vshap,Dvshap, T,Dtdxi)
!
      implicit none
!
      integer, intent(in)  :: Ie
      real(8), intent(in)  :: Xi(3),Vshap(5),Dvshap(3,5)
      real(8), intent(out) :: T,Dtdxi(3)
!
      integer :: iv
!
      select case(Ie)
      case(1,3)
        T = Xi(1) + Xi(3)*.5d0
        Dtdxi(1) = 1.d0; Dtdxi(2) = 0.d0; Dtdxi(3) = .5d0
      case(2,4)
        T = Xi(2) + Xi(3)*.5d0
        Dtdxi(1) = 0.d0; Dtdxi(2) = 1.d0; Dtdxi(3) = .5d0
      case(5,6,7,8)
        iv = Ie-4
        T = (Vshap(5) - Vshap(iv) + 1.d0)/2.d0
        Dtdxi(1:3) = (Dvshap(1:3,5) - Dvshap(1:3,iv))/2.d0
      end select
!
!
   end subroutine proj_d2e
!
!----------------------------------------------------------------------
!
!   routine name       - proj_d2f
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - define pyramiD to face projection used
!                        for face to pyramiD extensions
!
!   arguments :
!     in:
!         Xi           - master element coordinates
!         If           - face number
!                        master tetrahedron coordinates
!     out:
!         T            - local face coordinates
!         Dtdxi        - derivatives of the local face coordinates
!                        wrt master tetrahedron coordinates
!
!----------------------------------------------------------------------
!
   subroutine proj_d2f(Xi,If, T,Dtdxi)
!
      implicit none
!
      integer, intent(in)  :: If
      real(8), intent(in)  :: Xi(3)
      real(8), intent(out) :: T(2),Dtdxi(2,3)
!
!
      T(2) = Xi(3)
      Dtdxi(2,1:2) = 0.d0; Dtdxi(2,3) = 1.d0; Dtdxi(1,3) = 0.d0
      select case(If)
      case(2,4)
        T(1) = Xi(1); Dtdxi(1,1) = 1.d0; Dtdxi(1,2) = 0.d0
      case(3,5)
        T(1) = Xi(2); Dtdxi(1,1) = 0.d0; Dtdxi(1,2) = 1.d0
      case default
        write(*,7001) If
 7001   format('proj_d2f: WRONG If = ',i4)
        stop 1
      end select
!
!
   end subroutine proj_d2f
!
!
!----------------------------------------------------------------------
!> Purpose : brick to edge projection
!
!  @param[in ] Iv1,2   - vertices defining the edge
!  @param[in ] Vshape  - master brick vertex shape functions, computed
!                        at the point of interest
!  @param[in ] Dvshape - derivatives of vertex shape functions
!  @param[out] T       - edge coordinate
!  @param[out] Dtdxi   - derivative of edge coordinate
!
!  @revision Nov 12
!----------------------------------------------------------------------
   subroutine proj_b2e(Iv1,Iv2,Vshape,Dvshape, T,Dtdxi)
!
      use element_data , only : BRICK_FACE_TO_VERT
!
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2
      real(8),dimension(  8),intent(in ) :: Vshape
      real(8),dimension(3,8),intent(in ) :: Dvshape
      real(8)               ,intent(out) :: T
      real(8),dimension(3  ),intent(out) :: Dtdxi
!
      integer :: iface,i,iv
!----------------------------------------------------------------------
!
      call hexa_aux(Iv1,Iv2,iface)
!
      T=0.d0 ; Dtdxi(1:3)=0.d0
!  ...loop over face vertices
      do i=1,4
        iv=BRICK_FACE_TO_VERT(i,iface)
        T          = T          +  Vshape(    iv)
        Dtdxi(1:3) = Dtdxi(1:3) + Dvshape(1:3,iv)
      enddo
!
!
   end subroutine proj_b2e
!
!
!----------------------------------------------------------------------
!> Purpose : brick to face projection
!
!  @param[in ] Iv1,2,4 - vertices defining the face orientation
!  @param[in ] Vshape  - master brick vertex shape functions, computed
!                        at the point of interest
!  @param[in ] Dvshape - derivatives of vertex shape functions
!  @param[out] T       - face coordinates
!  @param[out] Dtdxi   - derivative of face coordinates
!
!  @revision Nov 12
!----------------------------------------------------------------------
   subroutine proj_b2f(Iv1,Iv2,Iv4,Vshape,Dvshape, T,Dtdxi)
!
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2,Iv4
      real(8),dimension(  8),intent(in ) :: Vshape
      real(8),dimension(3,8),intent(in ) :: Dvshape
      real(8),dimension(2  ),intent(out) :: T
      real(8),dimension(2,3),intent(out) :: Dtdxi
!----------------------------------------------------------------------
!
!  ...project on 1st axis
      call proj_b2e(Iv1,Iv2,Vshape,Dvshape, T(1),Dtdxi(1,1:3))
!  ...project on 2nd axis
      call proj_b2e(Iv1,Iv4,Vshape,Dvshape, T(2),Dtdxi(2,1:3))
!
!
   end subroutine proj_b2f
