subroutine exact_error_h1(Mdle, Derr,Dnorm)
!
      use control
      use element_data
      use data_structure3D
      use GMP
      use PROJ
#include "syscom.blk"
!
!  ...norm of exact solution and relative error
      dimension Dnorm(MAXEQNH)
      dimension Derr (MAXEQNH)

!  ...element order, geometry and solution dof
      dimension norder(19)
      dimension xnod(3,MAXbrickH), &
                zdofH(MAXEQNH,MAXbrickH),zdofE(MAXEQNE,MAXbrickE), &
                zdofV(MAXEQNV,MAXbrickV),zdofQ(MAXEQNQ,MAXbrickQ)
!
!  ...node orientations
      dimension nedge_orient(12), nface_orient(6)
!
!  ...shape functions and their derivatives
      dimension shapH(MAXbrickH),dshapH(3,MAXbrickH), &
                dshapHx(3,MAXbrickH)
!
!
!  ...geometry
      dimension xi(3),x(3),dxdxi(3,3),dxidx(3,3)
!
!  ...3D quadrature data
      dimension xiloc(3,MAX_NINT3),wxi(MAX_NINT3)
!
!  ...approximate solution and its derivatives
      dimension zsolH(MAXEQNH),zdsolH(MAXEQNH,3)
!
!  ...exact solution routine
      dimension &
                zvalH(    MAXEQNH    ),& 
                zdvalH(   MAXEQNH,3  ),&
                zd2valH(  MAXEQNH,3,3),&
                zvalE(  3,MAXEQNE    ),&
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
      select case(Mdle)
      case(1,40)
        iprint=0
      case default
        iprint=0
      end select
!
!
      if (iprint.eq.1) then
        write(*,7001) Mdle,NODES(Mdle)%type
 7001   format(' exact_error_h1: Mdle,type = ',i10,2x,a5)
      endif
!
!  ...determine order of approximation
      call find_order(Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
!
      call nodcor(Mdle,xnod)
      if (iprint.eq.1) then
        call celndof(NODES(Mdle)%type,Norder, &
                     nrdofH,nrdofE,nrdofV,nrdofQ)
        do k=1,nrdofH
          write(*,7002) k,xnod(1:3,k)
 7002     format(' exact_error_h1: k,xnod(1:3,k) = ',i3,3e12.5)
        enddo
      endif
!
!  ...determine approximate solution dofs
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
      if (iprint.eq.1) then
        write(*,*) 'exact_error_h1: zdofH = '
        do k=1,nrdofH
          write(*,7005) k,zdofH(1,k)
 7005     format(' k = ',i4,' zdofH(k,1) = ',e12.5)
        enddo
        call pause
      endif
!
      Derr=0.d0 ; Dnorm=0.d0
!
!-----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L
!-----------------------------------------------------------------------
!      
      INTEGRATION=2
      call set_3Dint(NODES(Mdle)%type,norder, nint,xiloc,wxi)
      INTEGRATION=0
!
!  ...loop through integration points
      do l=1,nint
        xi(1:3) = xiloc(1:3,l); wa = wxi(l)
!
!  .....evaluate appropriate shape functions at the point
        call shape3H(NODES(Mdle)%type,xi,norder,nedge_orient,nface_orient, &
                    nrdofH,shapH,dshapH)
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
          call exact_geom(Mdle,xi, x,dxdxi)
        endselect
!
!  .....evaluate the inverse derivatives and jacobian
        call geom(dxdxi, dxidx,rjac,iflag) 
        if (iflag.ne.0) then
          write(*,*) 'elem_elem_error: NEGATIVE JACOBIAN FOR Mdle = ',Mdle
          write(*,*) '        rjac = ',rjac
          stop 1
        endif
!        
!  .....compute derivatives wrt physical coordinates
        do k=1,nrdofH
          dshapHx(1:3,k) = 0.d0
          do ixi=1,3
            dshapHx(1:3,k) = dshapHx(1:3,k) &
                           + dshapH(ixi,k)*dxidx(ixi,1:3)
          enddo
        enddo
!     
!  .....total weight
        weight = wa*rjac
!
!  .....evaluate the approximate solution
        zsolH(1:MAXEQNH) = ZERO ; zdsolH(1:MAXEQNH,1:3) = ZERO
        do k=1,nrdofH
          zsolH(1:MAXEQNH) = zsolH(1:MAXEQNH) &
                           + zdofH(1:MAXEQNH,k)*shapH(k)
          do ivar=1,3
            zdsolH(1:MAXEQNH,ivar) = zdsolH(1:MAXEQNH,ivar) &
                           + zdofH(1:MAXEQNH,k)*dshapHx(ivar,k)
          enddo
        enddo
!
!  .....compute the exact solution 
        icase=1
        call exact(x,icase, zvalH,zdvalH,zd2valH, &
                            zvalE,zdvalE,zd2valE, &
                            zvalV,zdvalV,zd2valV, &
                            zvalQ,zdvalQ,zd2valQ)
        if (iprint.eq.1) then
          write(*,7003) l,x(1:3),zsolH(1),zvalH(1)
 7003     format('exact_error_h1: l,x = ',i4,2x,3f8.3,2x, &
                 ' APPROX SOLUTION = ',e12.5,' EXACT SOLUTION = ',e12.5)
        endif
!
!  .....accumulate for the error and norm
        do ieq=1,MAXEQNH
          Dnorm(ieq) = Dnorm(ieq) + abs(zvalH(ieq))**2*weight
          Derr(ieq) = Derr(ieq) + abs(zvalH(ieq)-zsolH(ieq))**2*weight
          if (.NOT. L2PROJ) then
            do ix=1,3
              Dnorm(ieq) = Dnorm(ieq) + abs(zdvalH(ieq,ix))**2*weight
              Derr(ieq) = Derr(ieq) &
                + abs(zdvalH(ieq,ix)-zdsolH(ieq,ix))**2*weight
            enddo
          endif 
        enddo
!
!  ...end of loop through integration points
      enddo

      if (Derr(1).gt.1.0E-10) then
!cc        write(*,7000)Mdle,Derr(1),Dnorm(1)
7000    format(' exact_error_h1: Mdle,Derr,Dnorm = ',i5,2x,3(e12.5,2x))
!cc        call pause
      endif
!
!  ...store in data structure
      NODES(Mdle)%error(1,1)=Derr(1)      
!
!
      if (iprint.eq.1) call pause
!
!
endsubroutine exact_error_h1
