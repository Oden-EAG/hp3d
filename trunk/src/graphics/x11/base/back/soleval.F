c----------------------------------------------------------------------
c
c   routine name       - soleval
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 11
c
c   purpose            - routine evaluates physical coordinates and
c                        solution at a point in a 3D element
c
c   arguments :
c     in:
c         Mdle         - element (middle node) number
c         Xi           - master element coordinates of a point
c         Nedge_orient - edge orientation
c         Nface_orient - face orientation
c         Norder       - order of approximation 
c         Xnod         - geometry dof
c         ZdofH,ZdofE,ZdofV,ZdofQ - H1,H(curl),H(div) and L2 dof
c         Nflag        = 0  coordinates only
c                        else  both coordinates and the solution
c     out:
c         X            - physical coordinates
c         Dxdxi        - derivatives of physical coordinates wrt
c                        element master coordinates
c         ZsolH        - H1 solution values
c         ZgradH       - H1 solution gradient
c         ZsolE        - H(curl) solution
c         ZcurlE       - H(curl) solution curl
c         ZsolV        - H(div) solution
c         ZdivV        - H(div) solution divergence
c         ZsolQ        - L2 solution
c
c----------------------------------------------------------------------
c
      subroutine soleval(Mdle,Xi,
     .                   Nedge_orient,Nface_orient,Norder,
     .                   Xnod,ZdofH,ZdofE,ZdofV,ZdofQ,Nflag, 
     .                   X,Dxdxi,
     .                   ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ)
c
      use control
      use data_structure3D
#include "syscom.blk"
c
      character(len=4) type
      dimension Xi(NDIMEN),Xnod(NDIMEN,MAXbrickH),
     .          ZdofH(MAXEQNH,MAXbrickH),ZdofE(MAXEQNE,MAXbrickE),
     .          ZdofV(MAXEQNV,MAXbrickV),ZdofQ(MAXEQNQ,MAXbrickQ)
      dimension Nedge_orient(12),Nface_orient(6),Norder(19)
      dimension X(NDIMEN),Dxdxi(NDIMEN,NDIMEN),
     .          ZsolH(MAXEQNH),ZgradH(MAXEQNH,3),
     .          ZsolE(3,MAXEQNE),ZcurlE(3,MAXEQNE),
     .          ZsolV(3,MAXEQNV),ZdivV(MAXEQNV),
     .          ZsolQ(MAXEQNQ)
c
      dimension shapH(MAXbrickH),gradH(3,MAXbrickH),
     .          shapE(3,MAXbrickE),curlE(3,MAXbrickE),
     .          shapV(3,MAXbrickV),divV(MAXbrickV),
     .          shapQ(MAXbrickQ)
      dimension zgradH_dxi(3),dxidx(3,3),zsolExi(3),zcurlExi(3)
c      
c----------------------------------------------------------------------
c
      iprint=0
      type=NODES(Mdle)%type
c     
c  ...evaluate shape functions at the point
      call shape3H(Type,Xi,Norder,Nedge_orient,Nface_orient,
     .             nrdofH,shapH,gradH)
c
c  ...printing      
      if (iprint.eq.1) then
        write(*,*) 'soleval: Xnod = '
        do k=1,nrdofH
          write(*,6001) k,Xnod(1:3,k)
 6001     format('k,Xnod(1:3,k) = ',i3,2x,3f8.3)
        enddo
        call pause
        write(*,*) 'soleval: ZdofH = '
        do k=1,nrdofH
          if (iprint.eq.1) then
            write(*,6002) k,shapH(k),ZdofH(1:MAXEQNH,k)
 6002       format('soleval: k,shapH,ZdofH = ',i3,2x,f8.3,
#if C_MODE
     .            /,5('(',2e12.5,')',2x))
#else
     .            /,10e12.5)
#endif
          endif
        enddo
      endif
c
      select case(EXGEOM)
c  ...parametric element      
      case(0)
        X(1:3)=0.d0 ; Dxdxi(1:3,1:3)=0.d0
        do k=1,nrdofH
          X(1:3) = X(1:3) + Xnod(1:3,k)*shapH(k)
          do i=1,3
            Dxdxi(1:3,i) = Dxdxi(1:3,i) + Xnod(1:3,k)*gradH(i,k)
          enddo
        enddo
c  ...exact geometry element        
      case(1)
        call exact_geom(Mdle,Xi, X,Dxdxi)
      end select
c
c  ...if geometry only, return
      if (Nflag.eq.0) return
c      
      call geom(Dxdxi, dxidx,rjac,iflag)
      if (iflag.ne.0) then
        write(*,*) 'soleval: iflag = ',iflag
        write(*,*) 'soleval: Xi = ',Xi
        call pause
      endif
c-----------------------------------------------------------------------
c
c  ...evaluate the solution and its derivatives
      ZsolH(1:MAXEQNH)      = ZERO
      ZgradH(1:MAXEQNH,1:3) = ZERO
c
c  ...loop through H1 solution components
      do i=1,MAXEQNH
c
        ZsolH(i)        = ZERO
        zgradH_dxi(1:3) = ZERO
        do k=1,nrdofH
          ZsolH(i) = ZsolH(i) 
     .             + ZdofH(i,k)*shapH(k)
          zgradH_dxi(1:3) = zgradH_dxi(1:3) 
     .                    + ZdofH(i,k)*gradH(1:3,k)
        enddo
c
        do ivar=1,3
          ZgradH(i,ivar) = zgradH_dxi(1)*dxidx(1,ivar) 
     .                   + zgradH_dxi(2)*dxidx(2,ivar)  
     .                   + zgradH_dxi(3)*dxidx(3,ivar) 
        enddo
c
      enddo

      if (iprint.eq.1) then
        write(*,7001) Mdle,Xi(1:3)
 7001   format('soleval: Mdle,Xi = ',i8,2x,3f8.3)
        write(*,7002) X(1:3)
 7002   format('         X = ',3f8.3)
        do ivar=1,MAXEQNH
          write(*,7003) ivar,ZsolH(ivar),ZgradH(ivar,1:3)
#if C_MODE
 7003     format('         ivar, ZsolH(ivar) = ',i2,2x,2e12.5,
     .           2x,3(2e12.5,2x))
#else
 7003     format('         ivar, ZsolH(ivar) = ',i2,2x,e12.5,
     .           3(e12.5,2x))
#endif
        enddo
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...sanity check for jacobians
      if (iprint.eq.2) then
        do i=1,3
        do j=1,3
          s = 0.d0
          do k=1,3
            s = s + dxdxi(i,k)*dxidx(k,j)
          enddo
          write(*,6034) i,j,s
 6034     format('soleval: i,j,s = ',2i2,2x,e12.5)
        enddo
        enddo
        call pause
      endif
c
c  ...evaluate shape functions at the point
      call shape3E(type,Xi,Norder,Nedge_orient,Nface_orient,
     .             nrdofE,shapE,curlE)
c 
c  ...evaluate the H(curl) solution and its curl
      ZsolE(1:3,1:MAXEQNE)  = ZERO
      ZcurlE(1:3,1:MAXEQNE) = ZERO
c
c  ...loop through H(curl) solution components
      do i=1,MAXEQNE
c
        ZsolExi(1:3)  = ZERO
        zcurlExi(1:3) = ZERO
        do k=1,nrdofE
          zsolExi(1:3)  = ZsolExi(1:3) 
     .                  + ZdofE(i,k)*shapE(1:3,k)
          zcurlExi(1:3) = zcurlExi(1:3) 
     .                  + ZdofE(i,k)*curlE(1:3,k)
        enddo
c
c  .....Piola transform
        ZsolE(1:3,i)  = ( zsolExi(1)*dxidx(1,1:3) 
     .                +   zsolExi(2)*dxidx(2,1:3)  
     .                +   zsolExi(3)*dxidx(3,1:3) )
        ZcurlE(1:3,i) = ( dxdxi(1:3,1)*zcurlExi(1)
     .                +   dxdxi(1:3,2)*zcurlExi(2)
     .                +   dxdxi(1:3,3)*zcurlExi(3) )/rjac
c
      enddo

      iprint =0 
      
      if (iprint.eq.2) then
        write(*,7001) Mdle,Xi(1:3)
        write(*,7002) X(1:3)
        do ivar=1,MAXEQNE
          write(*,7004) ivar,ZsolE(1:3,ivar),ZcurlE(1:3,ivar)
#if C_MODE
 7004     format('         ivar, ZsolE(1:3,ivar) = ',i2,2x,3(2e12.5,2x),
     .           ' ZcurlE(1:3,ivar) = ',2x,3(2e12.5,2x))
#else
 7004     format('         ivar, ZsolE(1:3,ivar) = ',i2,2x,3(e12.5,2x),
     .           ' ZcurlE(1:3,ivar) = ',2x,3(e12.5,2x))
#endif
        enddo
        call pause
      endif
c
      end subroutine
