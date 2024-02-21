c----------------------------------------------------------------------
c> Purpose : triangle to edge projection
c
c  @param[in ] Iv1,2   - vertices defining the edge
c  @param[in ] Vshape  - master triangle vertex shape functions,
c                        computed at the point of interest
c  @param[in ] Dvshape - derivatives of vertex shape functions
c  @param[out] T       - edge coordinate
c  @param[out] Dtdxi   - derivative of edge coordinate
c
c  @revision Nov 12
c----------------------------------------------------------------------
      subroutine proj_t2e(Iv1,Iv2,Vshape,Dvshape, T,Dtdxi)
c
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2
      real(8),dimension(  3),intent(in ) :: Vshape
      real(8),dimension(2,3),intent(in ) :: Dvshape
      real(8)               ,intent(out) :: T
      real(8),dimension(2)  ,intent(out) :: Dtdxi
c----------------------------------------------------------------------
c
      T = (Vshape(Iv2)-Vshape(Iv1)+1.d0)*.5d0
      Dtdxi(1:2) = (Dvshape(1:2,Iv2)-Dvshape(1:2,Iv1))*.5d0
c
      end subroutine proj_t2e
c
c
c----------------------------------------------------------------------
c> Purpose : quad to edge projection
c
c  @param[in ] Iv1,2   - vertices defining the edge
c  @param[in ] Vshape  - master quad vertex shape functions, computed
c                        at the point of interest
c  @param[in ] Dvshape - derivatives of vertex shape functions
c  @param[out] T       - edge coordinate
c  @param[out] Dtdxi   - derivative of edge coordinate
c
c  @revision Nov 12
c----------------------------------------------------------------------
      subroutine proj_r2e(Iv1,Iv2,Vshape,Dvshape, T,Dtdxi)
c
      use element_data , only : QUADR_EDGE_TO_VERT
c
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2
      real(8),dimension(  4),intent(in ) :: Vshape
      real(8),dimension(2,4),intent(in ) :: Dvshape
      real(8)               ,intent(out) :: T
      real(8),dimension(2)  ,intent(out) :: Dtdxi
c
      integer :: iedge,i,iv
c----------------------------------------------------------------------
c
      call quad_aux(Iv1,Iv2, iedge)
c
c  ...loop over edge vertices
      T=0.d0 ; Dtdxi(1:2)=0.d0
      do i=1,2
        iv=QUADR_EDGE_TO_VERT(i,iedge)
        T          = T          +  Vshape(    iv)
        Dtdxi(1:2) = Dtdxi(1:2) + Dvshape(1:2,iv)
      enddo
c
      end subroutine proj_r2e
c
c
c----------------------------------------------------------------------
c> Purpose : tetrahedroN to edge projection
c
c  @param[in ] Iv1,2   - vertices defining the edge
c  @param[in ] Vshape  - master tet vertex shape functions, computed
c                        at the point of interest
c  @param[in ] Dvshape - derivatives of vertex shape functions
c  @param[out] T       - edge coordinate
c  @param[out] Dtdxi   - derivative of edge coordinate
c
c  @revision Nov 12
c----------------------------------------------------------------------
      subroutine proj_n2e(Iv1,Iv2,Vshape,Dvshape, T,Dtdxi)
c
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2
      real(8),dimension(  4),intent(in ) :: Vshape
      real(8),dimension(3,4),intent(in ) :: Dvshape
      real(8)               ,intent(out) :: T
      real(8),dimension(3)  ,intent(out) :: Dtdxi
c----------------------------------------------------------------------
c
      T         =( Vshape(    Iv2)- Vshape(    Iv1)+1.d0)*0.5d0
      Dtdxi(1:3)=(Dvshape(1:3,Iv2)-Dvshape(1:3,Iv1)     )*0.5d0
c
      end subroutine proj_n2e
c
c
c----------------------------------------------------------------------
c> Purpose : tetrahedroN to face projection
c
c  @param[in ] Iv1,2,3   - vertices defining the face
c  @param[in ] Vshape    - master tet vertex shape functions, computed
c                          at the point of interest
c  @param[in ] Dvshape   - derivatives of vertex shape functions
c  @param[out] T         - face coordinate
c  @param[out] Dtdxi     - derivative of face coordinate
c
c  @revision Nov 12
c----------------------------------------------------------------------
      subroutine proj_n2f(Iv1,Iv2,Iv3,Vshape,Dvshape, T,Dtdxi)
c
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2,Iv3
      real(8),dimension(  4),intent(in ) :: Vshape
      real(8),dimension(3,4),intent(in ) :: Dvshape
      real(8),dimension(2)  ,intent(out) :: T
      real(8),dimension(2,3),intent(out) :: Dtdxi
c
      real(8)              :: tsum
      real(8),dimension(3) :: dtsum
c----------------------------------------------------------------------
c
      tsum  =  1.d0 - Vshape(Iv1) - Vshape(Iv2) - Vshape(Iv3)
      dtsum(1:3) = -Dvshape(1:3,Iv1)-Dvshape(1:3,Iv2)-Dvshape(1:3,Iv3)
      T(1) = Vshape(Iv2) + tsum/3.d0
      Dtdxi(1,1:3) = Dvshape(1:3,Iv2) + dtsum(1:3)/3.d0
      T(2) = Vshape(Iv3) + tsum/3.d0
      Dtdxi(2,1:3) = Dvshape(1:3,Iv3) + dtsum(1:3)/3.d0
c
cc      T(1) = Vshap(Iv2); Dtdxi(1,1:3) = Dvshap(1:3,Iv2)
cc      T(2) = Vshap(Iv3); Dtdxi(2,1:3) = Dvshap(1:3,Iv3)
cc      do iv=1,4
cc        if ((iv.ne.Iv1).and.(iv.ne.Iv2).and.(iv.ne.Iv3)) then
cc        T(1) = T(1) + Vshap(iv)/3.d0
cc          Dtdxi(1,1:3) = Dtdxi(1,1:3) + Dvshap(1:3,iv)/3.d0
cc          T(2) = T(2) + Vshap(iv)/3.d0
cc          Dtdxi(2,1:3) = Dtdxi(2,1:3) + Dvshap(1:3,iv)/3.d0
cc        endif
cc      enddo
cc      if (iprint.eq.1) then
cc        write(*,7001) Xi(1:3),T(1:2)
cc 7001   format('proj_n2f: Xi = ',3e12.5,' T = ',2e12.5)
cc        call pause
cc      endif
c
      end subroutine proj_n2f
c
c
c----------------------------------------------------------------------
c
c   routine name       - proj_d2e
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - define pyramiD to edge projection used
c                        for edge to pyramiD extensions
c
c   arguments :
c     in:
c             Xi       - master element coordinates
c             Ie       - edge number
c             Vshap    - vertex shape functions
c             Dvshap   - derivatives of vertex shape functions wrt
c                        master pyramid coordinates
c     out:
c             T        - local edge coordinate
c             Dtdxi    - derivatives of the local edge coordinate
c                        wrt master tetrahedron coordinates
c
c----------------------------------------------------------------------
c
      subroutine proj_d2e(Xi,Ie,Vshap,Dvshap, T,Dtdxi)
c
      implicit none
c
      integer, intent(in)  :: Ie
      real(8), intent(in)  :: Xi(3),Vshap(5),Dvshap(3,5)
      real(8), intent(out) :: T,Dtdxi(3)
c
      integer :: iv
c
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
c
c
      end subroutine proj_d2e
c
c----------------------------------------------------------------------
c
c   routine name       - proj_d2f
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - define pyramiD to face projection used
c                        for face to pyramiD extensions
c
c   arguments :
c     in:
c         Xi           - master element coordinates
c         If           - face number
c                        master tetrahedron coordinates
c     out:
c         T            - local face coordinates
c         Dtdxi        - derivatives of the local face coordinates
c                        wrt master tetrahedron coordinates
c
c----------------------------------------------------------------------
c
      subroutine proj_d2f(Xi,If, T,Dtdxi)
c
      implicit none
c
      integer, intent(in)  :: If
      real(8), intent(in)  :: Xi(3)
      real(8), intent(out) :: T(2),Dtdxi(2,3)
c
c
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
c
c
      end subroutine proj_d2f
c
c
c----------------------------------------------------------------------
c> Purpose : brick to edge projection
c
c  @param[in ] Iv1,2   - vertices defining the edge
c  @param[in ] Vshape  - master brick vertex shape functions, computed
c                        at the point of interest
c  @param[in ] Dvshape - derivatives of vertex shape functions
c  @param[out] T       - edge coordinate
c  @param[out] Dtdxi   - derivative of edge coordinate
c
c  @revision Nov 12
c----------------------------------------------------------------------
      subroutine proj_b2e(Iv1,Iv2,Vshape,Dvshape, T,Dtdxi)
c
      use element_data , only : BRICK_FACE_TO_VERT
c
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2
      real(8),dimension(  8),intent(in ) :: Vshape
      real(8),dimension(3,8),intent(in ) :: Dvshape
      real(8)               ,intent(out) :: T
      real(8),dimension(3  ),intent(out) :: Dtdxi
c
      integer :: iface,i,iv
c----------------------------------------------------------------------
c
      call hexa_aux(Iv1,Iv2,iface)
c
      T=0.d0 ; Dtdxi(1:3)=0.d0
c  ...loop over face vertices
      do i=1,4
        iv=BRICK_FACE_TO_VERT(i,iface)
        T          = T          +  Vshape(    iv)
        Dtdxi(1:3) = Dtdxi(1:3) + Dvshape(1:3,iv)
      enddo
c
c
      end subroutine proj_b2e
c
c
c----------------------------------------------------------------------
c> Purpose : brick to face projection
c
c  @param[in ] Iv1,2,4 - vertices defining the face orientation
c  @param[in ] Vshape  - master brick vertex shape functions, computed
c                        at the point of interest
c  @param[in ] Dvshape - derivatives of vertex shape functions
c  @param[out] T       - face coordinates
c  @param[out] Dtdxi   - derivative of face coordinates
c
c  @revision Nov 12
c----------------------------------------------------------------------
      subroutine proj_b2f(Iv1,Iv2,Iv4,Vshape,Dvshape, T,Dtdxi)
c
      implicit none
      integer               ,intent(in ) :: Iv1,Iv2,Iv4
      real(8),dimension(  8),intent(in ) :: Vshape
      real(8),dimension(3,8),intent(in ) :: Dvshape
      real(8),dimension(2  ),intent(out) :: T
      real(8),dimension(2,3),intent(out) :: Dtdxi
c----------------------------------------------------------------------
c
c  ...project on 1st axis
      call proj_b2e(Iv1,Iv2,Vshape,Dvshape, T(1),Dtdxi(1,1:3))
c  ...project on 2nd axis
      call proj_b2e(Iv1,Iv4,Vshape,Dvshape, T(2),Dtdxi(2,1:3))
c
c
      end subroutine proj_b2f
