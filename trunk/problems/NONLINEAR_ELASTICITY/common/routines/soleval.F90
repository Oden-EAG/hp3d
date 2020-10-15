!-------------------------------------------------------------------------------
!> Purpose : routine evaluates physical coordinates and solution at a point
!!
!> @param[in]  Mdle         - element (middle node) number
!> @param[in]  Xi           - master element coordinates of a point
!> @param[in]  Nedge_orient - edge orientation
!> @param[in]  Nface_orient - face orientation
!> @param[in]  Norder       - order of approximation 
!> @param[in]  Xnod         - geometry dof
!> @param[in]  ZdofH        - H1      dofs
!> @param[in]  ZdofE        - H(curl) dofs
!> @param[in]  ZdofV        - H(div)  dofs
!> @param[in]  ZdofQ        - L^2     dofs
!> @param[in]  Nflag        - 0 : geometry map only  
!!
!> @param[out] X            - physical coordinates
!> @param[out] Dxdxi        - derivatives of physical coordinates wrt
!!                            element master coordinates
!> @param[out] ZsolH,ZgradH - H1 solution values
!> @param[out] ZsolE,ZcurlE - H(curl) solution
!> @param[out] ZsolV,ZdiV   - H(div) solution
!> @param[out] ZsolQ        - L2 solution
!!
!> @data Nov 14
!-------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine soleval(Mdle,Xi,Nedge_orient,Nface_orient,Norder,Xnod,ZdofH,ZdofE,ZdofV,ZdofQ,Nflag, &
                   X,Dxdxi,ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ                          )
!     
      use control
      use data_structure3D_poly
!  
      implicit none
      integer,                             intent(in)  :: Mdle, Nflag
      real*8,  dimension(3),               intent(in)  :: Xi
      integer, dimension(12),              intent(in)  :: Nedge_orient
      integer, dimension(6),               intent(in)  :: Nface_orient
      integer, dimension(19),              intent(in)  :: Norder
      real*8,  dimension(3,MAXbrickH),     intent(in)  :: Xnod
      VTYPE, dimension(MAXEQNH,MAXbrickH), intent(in)  :: ZdofH
      VTYPE, dimension(MAXEQNE,MAXbrickE), intent(in)  :: ZdofE
      VTYPE, dimension(MAXEQNV,MAXbrickV), intent(in)  :: ZdofV
      VTYPE, dimension(MAXEQNQ,MAXbrickQ), intent(in)  :: ZdofQ
      real*8, dimension(3),                intent(out) :: X
      real*8, dimension(3,3),              intent(out) :: Dxdxi
      VTYPE, dimension(  MAXEQNH  ),       intent(out) :: ZsolH
      VTYPE, dimension(  MAXEQNH,3),       intent(out) :: ZgradH
      VTYPE, dimension(3,MAXEQNE  ),       intent(out) :: ZsolE
      VTYPE, dimension(3,MAXEQNE  ),       intent(out) :: ZcurlE
      VTYPE, dimension(3,MAXEQNV  ),       intent(out) :: ZsolV
      VTYPE, dimension(  MAXEQNV  ),       intent(out) :: ZdivV
      VTYPE, dimension(  MAXEQNQ  ),       intent(out) :: ZsolQ
!
#ifdef _PYHP3D
      integer :: etype
#else
      character(len=4) :: etype
#endif
!
      real*8, dimension(3,3) :: dxidx
      real*8,dimension(  MAXbrickH) :: shapH
      real*8,dimension(3,MAXbrickH) :: gradH,gradHx
      real*8,dimension(3,MAXbrickE) :: shapE,shapEx,curlE,curlEx
      real*8,dimension(3,MAXbrickV) :: shapV,shapVx
      real*8,dimension(  MAXbrickV) :: divV,divVx
      real*8,dimension(  MAXbrickQ) :: shapQ
!
      integer :: iprint,iflag,i,j,k,n, ivar,nrdofH,nrdofE,nrdofV,nrdofQ
      real*8  :: s, rjac
!-------------------------------------------------------------------------------
!  
      iprint=0
!
!     evaluate H1 shape functions
      etype=NODES(Mdle)%type
      call shape3H(etype,Xi,Norder,Nedge_orient,Nface_orient, nrdofH,shapH,gradH)
!      
!     geometry map
      select case(EXGEOM)
!      
!     -- ISOPARAMETRIC GEOMETRY MAP --
      case(0)
        X(1:3)=0.d0 ; Dxdxi(1:3,1:3)=0.d0
        do k=1,nrdofH
          X(1:3) = X(1:3) + Xnod(1:3,k)*shapH(k)
          do i=1,3
            Dxdxi(1:3,i) = Dxdxi(1:3,i) + Xnod(1:3,k)*gradH(i,k)
          enddo
        enddo
!
!     -- EXACT GEOMETRY MAP --
      case(1) ; call exact_geom(Mdle,Xi, X,Dxdxi)
      endselect
!
!     check Jacobian
      call geom(Dxdxi, dxidx,rjac,iflag)
      if (iflag /= 0) then
         write(*,8000) Mdle, iflag, Xi, rjac
 8000    format(' soleval: negative jacobian, Mdle, iflag, Xi, rjac = ', &
         i4,2x,i1,2x,4(e12.5,2x))
         write(*,*)'Exiting routine'
         call pause
         return
      endif
!      
!     if only geometry map is needed, return
      if (Nflag.eq.0) return
!
!     initialize
      ZsolH(    1:MAXEQNH)=ZERO ; ZgradH(    1:MAXEQNH,1:3)=ZERO
      ZsolE(1:3,1:MAXEQNE)=ZERO ; ZcurlE(1:3,1:MAXEQNE    )=ZERO
      ZsolV(1:3,1:MAXEQNV)=ZERO ; ZdivV(     1:MAXEQNV    )=ZERO
      ZsolQ(    1:MAXEQNQ)=ZERO
! !
! !===============================================================================
! !  H1 SOLUTION                                                                 | 
! !===============================================================================
! !
! !     compute derivatives wrt physical coordinates
!       do k=1,nrdofH
!         gradHx(1:3,k) = 0.d0
!         do i=1,3
!           gradHx(1:3,k) = gradHx(1:3,k) + gradH(i,k)*dxidx(i,1:3)
!         enddo
!       enddo

! !     loop through components
!       do n=1,MAXEQNH
! !
! !       loop through dof's      
!         do k=1,nrdofH
!           ZsolH( n    ) = ZsolH( n    ) + ZdofH(n,k)*shapH(     k)
!           ZgradH(n,1:3) = ZgradH(n,1:3) + ZdofH(n,k)*gradHx(1:3,k)
!         enddo
!       enddo
! !
! !     printing
!       if (iprint.eq.1) then
!         write(*,7001) Mdle,Xi(1:3)
! 7001    format('soleval: Mdle,Xi = ',i8,2x,3f8.3)
!         write(*,7002) X(1:3)
! 7002    format('         X = ',3f8.3)
!         do ivar=1,MAXEQNH
!           write(*,7003) ivar,ZsolH(ivar),ZgradH(ivar,1:3)
! #if C_MODE
! 7003     format('         ivar, ZsolH(ivar) = ',i2,2x,2e12.5,2x,3(2e12.5,2x))
! #else
! 7003     format('         ivar, ZsolH(ivar) = ',i2,2x,e12.5,3(e12.5,2x))
! #endif
!         enddo
!         call pause
!       endif
! !
! !===============================================================================
! !  H(curl) SOLUTION                                                            |
! !===============================================================================
! !
! !     H(curl) shape functions
!       call shape3E(etype,Xi,Norder,Nedge_orient,Nface_orient, nrdofE,shapE,curlE)
! !
! !     Piola transform
!       shapEx=ZERO ; curlEx=ZERO
!       do k=1,nrdofE
!          do i=1,3
!             do j=1,3
!                shapEx(i,k) = shapEx(i,k)+ (dxidx(j,i)*shapE(j,k))
!                curlEx(i,k) = curlEx(i,k)+(dxdxi(i,j)*curlE(j,k))/rjac
!             enddo
!          enddo
!       enddo
! !    
! !     H(curl) solution and its curl
!       do n=1,MAXEQNE
!          do k=1,nrdofE
!             ZsolE (1:3,n) = ZsolE (1:3,n) + ZdofE(n,k)*shapEx(1:3,k)
!             ZcurlE(1:3,n) = ZcurlE(1:3,n) + ZdofE(n,k)*curlEx(1:3,k)
!          enddo
!       enddo
! !
! !     printing
!       if (iprint.eq.2) then
!          write(*,7001) Mdle,Xi(1:3)
!          write(*,7002) X(1:3)
!          do ivar=1,MAXEQNE
!             write(*,7004) ivar,ZsolE(1:3,ivar),ZcurlE(1:3,ivar)
! #if C_MODE
! 7004        format('         ivar, ZsolE(1:3,ivar) = ',i2,2x,3(2e12.5,2x), &
!                    ' ZcurlE(1:3,ivar) = ',2x,3(2e12.5,2x))
! #else
! 7004        format('         ivar, ZsolE(1:3,ivar) = ',i2,2x,3(e12.5,2x), &
!                    ' ZcurlE(1:3,ivar) = ',2x,3(e12.5,2x))
! #endif
!          enddo
!          call pause
!       endif
! !
! !===============================================================================
! !  H(div) SOLUTION                                                             |
! !===============================================================================
! !
! !     H(div) shape functions
!       call shape3V(etype,Xi,Norder,Nface_orient, nrdofV,shapV,divV)
! !      
! !     Piola transform (H(div) shape functions are transformed as H(curl) curls)
!       shapVx=ZERO ; divVx=ZERO
!       do k=1,nrdofV
!         do i=1,3
!           do j=1,3
!             shapVx(i,k) = shapVx(i,k)+ dxdxi(i,j)*shapV(j,k)/rjac
!           enddo
!         enddo
!         divVx(k) = divV(k)/rjac
!       enddo
! !       
! !     H(div) solution and its divergence 
!       ZsolV(1:3,1:MAXEQNV)=ZERO ; ZdivV(1:MAXEQNV)=ZERO
!       do k=1,nrdofV
!         do ivar=1,MAXEQNV
!           ZsolV(1:3,ivar) = ZsolV(1:3,ivar) + zdofV(ivar,k)*shapVx(1:3,k)
!           ZdivV(    ivar) = ZdivV(    ivar) + zdofV(ivar,k)* divVx(    k)
!         enddo
!       enddo
! !
! !     printing
!       if (iprint.eq.2) then
!         write(*,7001) Mdle,Xi(1:3)
!         write(*,7002) X(1:3)
!         do ivar=1,MAXEQNV
!           write(*,7005) ivar,ZsolV(1:3,ivar),ZdivV(ivar)
! #if C_MODE
! 7005      format('         ivar, ZsolV(1:3,ivar) = ',i2,2x,3(2e12.5,2x), &
!                ' ZdivV(ivar) = ',2x,2(e12.5))
! #else
! 7005      format('         ivar, ZsolV(1:3,ivar) = ',i2,2x,3(e12.5,2x), &
!                ' ZdivV(ivar) = ',2x,e12.5)
! #endif
!         enddo
!         call pause
!       endif
!
!===============================================================================
!  L2 SOLUTION                                                                 |
!===============================================================================
!
!     L2 shape functions
      call shape3Q(etype,Xi,Norder, nrdofQ,shapQ)
!      
!     Piola transform
      shapQ(1:nrdofQ)=shapQ(1:nrdofQ)/rjac
!      
!     evaluate the approximate solution
      ZsolQ(1:MAXEQNQ)=ZERO
      do k=1,nrdofQ
        do ivar=1,MAXEQNQ
          ZsolQ(ivar) = ZsolQ(ivar) + zdofQ(ivar,k)*shapQ(k)
        enddo
      enddo
!
!     printing
      if (iprint.eq.2) then
        write(*,7001) Mdle,Xi(1:3)
        write(*,7002) X(1:3)
        do ivar=1,MAXEQNQ
          write(*,7006) ivar,ZsolQ(ivar)
#if C_MODE
7006      format('         ivar, ZsolQ(ivar) = ',i2,2x,2e12.5,2x)
#else
7006      format('         ivar, ZsolV(ivar) = ',i2,2x,e12.5,2x)
#endif
        enddo
        call pause
      endif
!
!
end subroutine soleval
