#if HP3D_USE_X11

!-------------------------------------------------------------------------------------
!> @brief      routine computes physical coordinates and quantity to display
!!
!> @param[in]  Numlev       - 0:coordinate only, >0 compute quantity to display
!> @param[in]  Mdle         - element middle node
!> @param[in]  Iface        - face number
!> @param[in]  Nedge_orient - edge orientation
!> @param[in]  Nface_orient - face orientation
!> @param[in]  Norder       - order of approximation for the element
!> @param[in]  Xnod         - geometry dof for the element
!> @param[in]  ZdofH        - H1 dofs
!> @param[in]  ZdofE        - H(curl) dofs
!> @param[in]  ZdofV        - H(div) dofs
!> @param[in]  ZdofQ        - L^2 dofs
!> @param[in]  T            - master face coordinates
!> @param[out] X            - physical coordinates
!> @param[out] Val          - value of the quantity to display
!!
!> @date       Feb 2023
!-------------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine compute_face(Numlev,Mdle,Iface,Nedge_orient,Nface_orient,Norder, &
                        Xnod,ZdofH,ZdofE,ZdofV,ZdofQ,T, X,Val)
!
      use data_structure3D
!
      implicit none
      integer,                             intent(in)  :: Numlev,Mdle,Iface
      integer, dimension(12),              intent(in)  :: Nedge_orient
      integer, dimension(6),               intent(in)  :: Nface_orient
      integer, dimension(19),              intent(in)  :: Norder
      real(8), dimension(3,MAXbrickH),     intent(in)  :: Xnod
      VTYPE, dimension(MAXEQNH,MAXbrickH), intent(in)  :: ZdofH
      VTYPE, dimension(MAXEQNE,MAXbrickE), intent(in)  :: ZdofE
      VTYPE, dimension(MAXEQNV,MAXbrickV), intent(in)  :: ZdofV
      VTYPE, dimension(MAXEQNQ,MAXbrickQ), intent(in)  :: ZdofQ
      real(8), dimension(2),               intent(in)  :: T
      real(8), dimension(3),               intent(out) :: X
      real(8),                             intent(out) :: Val
!
      VTYPE, dimension(  MAXEQNH  ) :: zsolH
      VTYPE, dimension(  MAXEQNH,3) :: zgradH
      VTYPE, dimension(3,MAXEQNE  ) :: zsolE
      VTYPE, dimension(3,MAXEQNE  ) :: zcurlE
      VTYPE, dimension(3,MAXEQNV  ) :: zsolV
      VTYPE, dimension(  MAXEQNV  ) :: zdivV
      VTYPE, dimension(  MAXEQNQ  ) :: zsolQ
!
      real(8), dimension(3)   :: xi
      real(8), dimension(3,2) :: dxidt
      real(8), dimension(3,3) :: dxdxi
      real(8), dimension(3,2) :: dxdt
      real(8), dimension(3) :: rn
!
      integer :: i, j, ntype
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!--------------------------------------------------------------------------------------------------
!
      ntype=NODES(Mdle)%ntype
      call face_param(ntype,Iface,T, xi,dxidt)
!
!     evaluate the physical coordinates of the point and the value of the solution
      call soleval(Mdle,xi,Nedge_orient,Nface_orient,Norder,Xnod,ZdofH,ZdofE,ZdofV,ZdofQ, &
                   Numlev, X,dxdxi,zsolH,zgradH,zsolE,zcurlE,zsolV,zdivV,zsolQ            )
!
#if HP3D_DEBUG
!     printing
      if (iprint.eq.1) then
        write(*,*)'compute_face:'
        do i=1,3
          write(*,1112)i,dxdxi(i,1:3)
1112      format(' i,dxdxi(i,:) = ',i1,2x,3(e12.5,2x))
        enddo
        do i=1,2
          write(*,1113)i,dxidt(1:3,i)
1113      format(' i,dxidt(:,i) = ',i1,2x,3(e12.5,2x))
        enddo
      endif
#endif
!
!     face parameterization
      dxdt(1:3,1:2)=0.d0
      do i=1,2
        do j=1,3
          dxdt(1:3,i) = dxdt(1:3,i) + dxdxi(1:3,j)*dxidt(j,i)
        enddo
      enddo
!
#if HP3D_DEBUG
!     printing
      if (iprint.eq.1) then
        do i=1,2
          write(*,1111)i,dxdt(1:3,i)
1111      format(' i, dxdt(:,i) = ',i1,2x,3(e12.5,2x))
        enddo
      endif
#endif
!
!     normal vector
      call cross_product(dxdt(1:3,1),dxdt(1:3,2), rn)
      call normalize(rn)
      rn(1:3) = rn(1:3)*Nsign_param(ntype,Iface)
!
!     select the quantity to display
      Val=0.d0
      if (Numlev.gt.0) then
        call soldis(Mdle,xi,X,rn,zsolH,zgradH,zsolE,zcurlE,zsolV,zdivV,zsolQ, Val)
      endif
!
!
end subroutine compute_face

#endif
