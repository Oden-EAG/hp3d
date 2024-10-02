!-------------------------------------------------------------------------------
!> @brief      routine evaluates physical coordinates and solution at a point
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
!!                            1 : all variables
!!                           <0 : 4-digit encoding specifying variables
!!                                e.g.: -1010: H1, Hdiv
!!
!> @param[out] X            - physical coordinates
!> @param[out] Dxdxi        - derivatives of physical coordinates wrt
!!                            element master coordinates
!> @param[out] ZsolH,ZgradH - H1 solution values
!> @param[out] ZsolE,ZcurlE - H(curl) solution
!> @param[out] ZsolV,ZdiV   - H(div) solution
!> @param[out] ZsolQ        - L2 solution
!!
!> @date       Apr 2024
!-------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine soleval(Mdle,Xi,Nedge_orient,Nface_orient,Norder,Xnod,ZdofH,ZdofE,ZdofV,ZdofQ,Nflag, &
                   X,Dxdxi,ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ                          )
!
   use control
   use data_structure3D
!
   implicit none
   integer,                             intent(in)  :: Mdle
   real(8), dimension(3),               intent(in)  :: Xi
   integer, dimension(12),              intent(in)  :: Nedge_orient
   integer, dimension(6),               intent(in)  :: Nface_orient
   integer, dimension(19),              intent(in)  :: Norder
   real(8), dimension(3,MAXbrickH),     intent(in)  :: Xnod
   VTYPE, dimension(MAXEQNH,MAXbrickH), intent(in)  :: ZdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE), intent(in)  :: ZdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV), intent(in)  :: ZdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ), intent(in)  :: ZdofQ
   integer,                             intent(in)  :: Nflag
   real(8), dimension(3),               intent(out) :: X
   real(8), dimension(3,3),             intent(out) :: Dxdxi
   VTYPE, dimension(  MAXEQNH  ),       intent(out) :: ZsolH
   VTYPE, dimension(  MAXEQNH,3),       intent(out) :: ZgradH
   VTYPE, dimension(3,MAXEQNE  ),       intent(out) :: ZsolE
   VTYPE, dimension(3,MAXEQNE  ),       intent(out) :: ZcurlE
   VTYPE, dimension(3,MAXEQNV  ),       intent(out) :: ZsolV
   VTYPE, dimension(  MAXEQNV  ),       intent(out) :: ZdivV
   VTYPE, dimension(  MAXEQNQ  ),       intent(out) :: ZsolQ
!
   real(8),dimension(3,3)         :: dxidx
   real(8),dimension(  MAXbrickH) :: shapH
   real(8),dimension(3,MAXbrickH) :: gradH,gradHx
   real(8),dimension(3,MAXbrickE) :: shapE,shapEx,curlE,curlEx
   real(8),dimension(3,MAXbrickV) :: shapV,shapVx
   real(8),dimension(  MAXbrickV) :: divV
   real(8),dimension(  MAXbrickQ) :: shapQ
!
   integer :: iflag,k,n,nrdofH,nrdofE,nrdofV,nrdofQ,ntype
   integer :: mod,len,vflag(4)
   real(8) :: rjac
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
#endif
!
!-------------------------------------------------------------------------------
!
!..evaluate H1 shape functions
   ntype=NODES(Mdle)%ntype
   call shape3DH(ntype,Xi,Norder,Nedge_orient,Nface_orient, nrdofH,shapH,gradH)
!
!..geometry map
   call geom3D(Mdle,Xi,Xnod,shapH,gradH,nrdofH, &
               X,Dxdxi,dxidx,rjac,iflag)
!
!..if only geometry map is needed, return
   if (Nflag.eq.0) return
!
   if (Nflag < 0) then
      mod = 10; len = 4
      call decod(-Nflag,mod,len, vflag)
   else
      vflag(1:4) = 1
   endif
!
!..initialize
   ZsolH(    1:MAXEQNH)=ZERO ; ZgradH(    1:MAXEQNH,1:3)=ZERO
   ZsolE(1:3,1:MAXEQNE)=ZERO ; ZcurlE(1:3,1:MAXEQNE    )=ZERO
   ZsolV(1:3,1:MAXEQNV)=ZERO ; ZdivV(     1:MAXEQNV    )=ZERO
   ZsolQ(    1:MAXEQNQ)=ZERO
!
!===============================================================================
!  H1 SOLUTION                                                                 |
!===============================================================================
!
   if (MAXEQNH.eq.0 .or. vflag(1).eq.0) goto 10
!
!..Piola transform
!   gradHx = 0.d0
!   do k=1,nrdofH
!      do i=1,3
!         gradHx(1:3,k) = gradHx(1:3,k) + gradH(i,k)*dxidx(i,1:3)
!      enddo
!   enddo
   call piola_hcurl(dxidx,nrdofH,gradH(1:3,1:nrdofH), gradHx(1:3,1:nrdofH))
!
!..H1 solution and its grad
   do k=1,nrdofH
      do n=1,MAXEQNH
         ZsolH( n    ) = ZsolH( n    ) + ZdofH(n,k)*shapH(     k)
         ZgradH(n,1:3) = ZgradH(n,1:3) + ZdofH(n,k)*gradHx(1:3,k)
      enddo
   enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) Mdle,Xi(1:3)
7001    format('soleval: Mdle,Xi = ',i8,2x,3f8.3)
        write(*,7002) X(1:3)
7002    format('         X = ',3f8.3)
        do n=1,MAXEQNH
          write(*,7003) n,ZsolH(n),ZgradH(n,1:3)
#if HP3D_COMPLEX
7003     format('         n, ZsolH(n) = ',i2,2x,2e12.5,2x,3(2e12.5,2x))
#else
7003     format('         n, ZsolH(n) = ',i2,2x,e12.5,3(e12.5,2x))
#endif
        enddo
        call pause
      endif
#endif
!
   10 continue
!
!===============================================================================
!  H(curl) SOLUTION                                                            |
!===============================================================================
!
   if (MAXEQNE.eq.0 .or. vflag(2).eq.0) goto 20
!
!..H(curl) shape functions
   call shape3DE(ntype,Xi,Norder,Nedge_orient,Nface_orient, nrdofE,shapE,curlE)
!
!..Piola transform
!   shapEx=0.d0 ; curlEx=0.d0
!   do k=1,nrdofE
!      do i=1,3
!         do j=1,3
!            shapEx(i,k) = shapEx(i,k)+(dxidx(j,i)*shapE(j,k))
!            curlEx(i,k) = curlEx(i,k)+(Dxdxi(i,j)*curlE(j,k))/rjac
!         enddo
!      enddo
!   enddo
   call piola_hcurl(dxidx,     nrdofE,shapE(1:3,1:nrdofE), shapEx(1:3,1:nrdofE))
   call piola_hdiv (Dxdxi,rjac,nrdofE,curlE(1:3,1:nrdofE), curlEx(1:3,1:nrdofE))
!
!..H(curl) solution and its curl
   do k=1,nrdofE
      do n=1,MAXEQNE
         ZsolE (1:3,n) = ZsolE (1:3,n) + ZdofE(n,k)*shapEx(1:3,k)
         ZcurlE(1:3,n) = ZcurlE(1:3,n) + ZdofE(n,k)*curlEx(1:3,k)
      enddo
   enddo
!
#if HP3D_DEBUG
      if (iprint.eq.2) then
         write(*,7001) Mdle,Xi(1:3)
         write(*,7002) X(1:3)
         do n=1,MAXEQNE
            write(*,7004) n,ZsolE(1:3,n),ZcurlE(1:3,n)
#if HP3D_COMPLEX
7004        format('         n, ZsolE(1:3,n) = ',i2,2x,3(2e12.5,2x), &
                   ' ZcurlE(1:3,n) = ',2x,3(2e12.5,2x))
#else
7004        format('         n, ZsolE(1:3,n) = ',i2,2x,3(e12.5,2x), &
                   ' ZcurlE(1:3,n) = ',2x,3(e12.5,2x))
#endif
         enddo
         call pause
      endif
#endif
!
   20 continue
!
!===============================================================================
!  H(div) SOLUTION                                                             |
!===============================================================================
!
   if (MAXEQNV.eq.0 .or. vflag(3).eq.0) goto 30
!
!..H(div) shape functions
   call shape3DV(ntype,Xi,Norder,Nface_orient, nrdofV,shapV,divV)
!
!..Piola transform
!   shapVx=0.d0 ; divVx=0.d0
!   do k=1,nrdofV
!      do j=1,3
!         do i=1,3
!            shapVx(i,k) = shapVx(i,k) + Dxdxi(i,j)*shapV(j,k)/rjac
!         enddo
!      enddo
!      divVx(k) = divV(k)/rjac
!   enddo
   call piola_hdiv(Dxdxi,rjac,nrdofV,shapV(1:3,1:nrdofV), shapVx(1:3,1:nrdofV))
   divV(1:nrdofV) = divV(1:nrdofV)/rjac
!
!..H(div) solution and its divergence
   do k=1,nrdofV
      do n=1,MAXEQNV
         ZsolV(1:3,n) = ZsolV(1:3,n) + zdofV(n,k)*shapVx(1:3,k)
         ZdivV(    n) = ZdivV(    n) + zdofV(n,k)* divV (    k)
      enddo
   enddo
!
#if HP3D_DEBUG
      if (iprint.eq.2) then
        write(*,7001) Mdle,Xi(1:3)
        write(*,7002) X(1:3)
        do n=1,MAXEQNV
          write(*,7005) n,ZsolV(1:3,n),ZdivV(n)
#if HP3D_COMPLEX
7005      format('         n, ZsolV(1:3,n) = ',i2,2x,3(2e12.5,2x), &
               ' ZdivV(n) = ',2x,2(e12.5))
#else
7005      format('         n, ZsolV(1:3,n) = ',i2,2x,3(e12.5,2x), &
               ' ZdivV(n) = ',2x,e12.5)
#endif
        enddo
        call pause
      endif
#endif
!
   30 continue
!
!===============================================================================
!  L2 SOLUTION                                                                 |
!===============================================================================
!
   if (MAXEQNQ.eq.0 .or. vflag(4).eq.0) goto 40
!
!..L2 shape functions
   call shape3DQ(ntype,Xi,Norder, nrdofQ,shapQ)
!
!..Piola transform
   shapQ(1:nrdofQ) = shapQ(1:nrdofQ)/rjac
!
!..evaluate the approximate solution
   do k=1,nrdofQ
      do n=1,MAXEQNQ
         ZsolQ(n) = ZsolQ(n) + zdofQ(n,k)*shapQ(k)
      enddo
   enddo
!
#if HP3D_DEBUG
      if (iprint.eq.2) then
        write(*,7001) Mdle,Xi(1:3)
        write(*,7002) X(1:3)
        do n=1,MAXEQNQ
          write(*,7006) n,ZsolQ(n)
#if HP3D_COMPLEX
7006      format('         n, ZsolQ(n) = ',i2,2x,2e12.5,2x)
#else
7006      format('         n, ZsolV(n) = ',i2,2x,e12.5,2x)
#endif
        enddo
        call pause
      endif
#endif
!
   40 continue
!
end subroutine soleval
