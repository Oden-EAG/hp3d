subroutine exact_error_hcurl(Mdle, Derr,Dnorm)
!
      use control
      use element_data
      use data_structure3D
      use GMP
      use PROJ
#include "syscom.blk"
!
!  ...norm of exact solution and relative error
      dimension Dnorm(MAXEQNE)
      dimension Derr (MAXEQNE)

!  ...element order (use brick)
      dimension norder(19)
!  ...element goemetry dof's      
      dimension xnod(3,MAXbrickH)
!  ...element solution dof's      
      dimension zdofH(MAXEQNH,MAXbrickH), &
                zdofE(MAXEQNE,MAXbrickE), &
                zdofV(MAXEQNV,MAXbrickV), &
                zdofQ(MAXEQNQ,MAXbrickQ)
!
!  ...node orientations (use brick)
      dimension nedge_orient(12), nface_orient(6)
!
!  ...shape functions and their derivatives 
      dimension shapH(MAXbrickH),dshapH(3,MAXbrickH), &
                dshapHx(3,MAXbrickH), &
                shapE(3,MAXbrickE),curlE(3,MAXbrickE), &
                shapEx(3,MAXbrickE),curlEx(3,MAXbrickE)
      
!
!  ...geometry 
      dimension xi(3),x(3),dxdxi(3,3),dxidx(3,3)
!
!  ...3D quadrature data
      dimension xiloc(3,MAX_NINT3),wxi(MAX_NINT3)
!
!  ...approximate solution before and after Piola transform
      dimension zsolH(  MAXEQNH),zdsolH(  MAXEQNH,3), &
                zsolE(3,MAXEQNE),zcurlE(3,MAXEQNE  )
!
      dimension zaux1(3),zaux2(3)      
!     
! ...exact solution
      dimension & 
                zvalH(    MAXEQNH    ),& 
                zdvalH(   MAXEQNH,3  ),&
                zd2valH(  MAXEQNH,3,3),&
                zvalE(  3,MAXEQNE    ),&
                zcvalE( 3,MAXEQNE    ),& 
                zdvalE( 3,MAXEQNE,3  ),&
                zd2valE(3,MAXEQNE,3,3),&
                zvalV(  3,MAXEQNV    ),&
                zdvalV( 3,MAXEQNV,3  ),&
                zd2valV(3,MAXEQNV,3,3),&
                zvalQ(    MAXEQNQ    ),&
                zdvalQ(   MAXEQNQ,3  ),&
                zd2valQ(  MAXEQNQ,3,3)
!
!  ...for computing L2 error only
      logical :: L2PROJ=.FALSE.
!
!---------------------------------------------------------------------
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  ...HACK FOR HEXA
      if ( (NODES(Mdle)%type.eq.'bric').or.         &
           (NODES(Mdle)%type.eq.'mdlb')     ) then
        Derr=0.d0 ; Dnorm=0.d0
        return
      endif
!      
!  ...HACK FOR PYRAMID
      if ( (NODES(Mdle)%type.eq.'pyra').or.         &
           (NODES(Mdle)%type.eq.'mdld')     ) then
        Derr=0.d0 ; Dnorm=0.d0
        return
      endif
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
      select case(Mdle)
      case(1)      ; iprint=0
      case default ; iprint=0
      endselect
!
      if (iprint.eq.1) then
        write(*,7001) Mdle,NODES(Mdle)%type
 7001   format(' exact_error_hcurl: Mdle,type = ',i10,2x,a5)
      endif
!
!  ...determine order of approximation
      call find_order(Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
!
      call nodcor(Mdle,xnod)
      if (iprint.eq.1) then
        call celndof(NODES(Mdle)%type,Norder,  &
                     nrdofH,nrdofE,nrdofV,nrdofQ)
        do k=1,nrdofH
          write(*,7002) k,xnod(1:3,k)
 7002     format('exact_error_hcurl: k,xnod(1:3,k) = ',i3,3e12.5)
        enddo
      endif
!
!  ...determine solution dof
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
      if (iprint.eq.1) then
        write(*,*) 'exact_error_hcurl: zdofE = '
        do k=1,nrdofE
          write(*,7008) k,zdofE(1,k)
 7008     format('k = ',i4,' zdofE(k,1) = ',e12.5)
        enddo
        call pause
      endif
!
      Derr = 0.d0; Dnorm = 0.d0
!
!-----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L      
!-----------------------------------------------------------------------
!
!  ...set up the element quadrature
      call set_3Dint(NODES(Mdle)%type,norder, nint,xiloc,wxi)
!
!  ...loop through integration points
      do l=1,nint
        xi(1:3) = xiloc(1:3,l); wa = wxi(l)
!
!  .....evaluate appropriate shape functions at the point
        call shape3H(NODES(Mdle)%type,xi,norder, &
                     nedge_orient,nface_orient, &
                     nrdofH,shapH,dshapH)
        call shape3E(NODES(Mdle)%type,xi,norder, &
                     nedge_orient,nface_orient, &
                     nrdofE,shapE,curlE)
!
!  .....geometry map
        select case(EXGEOM)
        case(0)
          x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0 
          do k=1,nrdofH
            x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
            do i=1,3
              dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
            enddo
          enddo
        case(1)
          call exact_geom(Mdlem,xi, x,dxdxi)
        endselect
!
!  .....evaluate the inverse derivatives and jacobian
        call geom(dxdxi, dxidx,rjac,iflag) 
        if (iflag.ne.0) then
          write(*,*) 'elem_elem_error: NEGATIVE JACOBIAN FOR Mdle = ', Mdle
          write(*,*) '        rjac = ',rjac
          stop 1
        endif
!        
!   ....shapE and curlE wrt physical coordinate (Piola transform)
        shapEx = ZERO; curlEx = ZERO;
        do k=1,nrdofE
          do i=1,3
            do j=1,3
              shapEx(i,k) = shapEx(i,k)+ &
                            (dxidx(j,i)*shapE(j,k))
              curlEx(i,k) = curlEx(i,k)+ &
                            (dxdxi(i,j)*curlE(j,k))/rjac
            enddo
          enddo
        enddo
        
!
!  .....total weight
        weight = wa*rjac

        zsolE(1:3,1:MAXEQNE) = ZERO; zcurlE(1:3,1:MAXEQNE) = ZERO
!
!  .....compute the approxiamte solution
        do n=1,MAXEQNE
          do k=1,nrdofE
            zsolE (1:3,n) = zsolE (1:3,n) + zdofE(n,k)*shapEx(1:3,k)
            zcurlE(1:3,n) = zcurlE(1:3,n) + zdofE(n,k)*curlEx(1:3,k)
          enddo 
!          
!  .....end of loop over equations            
        enddo
!
!  .....compute the exact solution 
        icase=1
        call exact(x,icase, zvalH,zdvalH,zd2valH, &
                            zvalE,zdvalE,zd2valE, &
                            zvalV,zdvalV,zd2valV, &
                            zvalQ,zdvalQ,zd2valQ)
!  .....curl of exact solution
        zcvalE(1,1:MAXEQNE) = zdvalE(3,1:MAXEQNE,2) - &
                              zdvalE(2,1:MAXEQNE,3)
        zcvalE(2,1:MAXEQNE) = zdvalE(1,1:MAXEQNE,3) - &
                              zdvalE(3,1:MAXEQNE,1)
        zcvalE(3,1:MAXEQNE) = zdvalE(2,1:MAXEQNE,1) - &
                              zdvalE(1,1:MAXEQNE,2)
!
!  .....loop over equations                
        do n=1,MAXEQNE
          if (iprint.eq.1) then
            write(*,7010) l,n,x
 7010       format('exact_error_hcurl: l,n,x = ',2i3,2x,3f8.3)
            write(*,7011) zsolE(1:3,n),zcurlE(1:3,n)
#if(C_MODE)
 7011       format('approx(grad,curl) = ',2(3(2e12.5,2x),3x))
#else
 7011       format('approx(grad,curl) = ',2(3e12.5,2x))
#endif
            write(*,7012) zvalE(1:3,n),zcvalE(1:3,n)
#if(C_MODE)
 7012       format('exact(grad,curl)  = ',2(3(2e12.5,2x),3x))
#else
 7012       format('exact(grad,curl)  = ',2(3e12.5,2x))
#endif
            call pause
          endif
!
!  .......compute norm of the hp solution, and norm of the error          
          do i=1,3
            Dnorm(n) = Dnorm(n) + abs( zvalE(i,n)            )**2*weight
            Derr(n)  = Derr(n)  + abs( zvalE(i,n)- zsolE(i,n))**2*weight
!  .........L2 norm of the gradient
            if (.NOT. L2PROJ) then
            Dnorm(n) = Dnorm(n) + abs(zcvalE(i,n)            )**2*weight
            Derr(n)  = Derr(n)  + abs(zcvalE(i,n)-zcurlE(i,n))**2*weight
            endif
          enddo    
!
!  .......end of loop over equations
          enddo
!        
!  ...end of loop through integration points
      enddo
!
      if (iprint.eq.1) call pause
!
!
endsubroutine exact_error_hcurl
