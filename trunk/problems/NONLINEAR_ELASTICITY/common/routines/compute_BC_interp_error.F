c--------------------------------------------------------------------
c
c     routine name      - elem_interp_error
c
c--------------------------------------------------------------------
c
c     latest revision:  - Mar 15
c
c     purpose:          - routine returns element interpolation
c                         error for H1,H(curl and H(div) BC data
c                         (squared); a debugging routine verifying
c                         PB interpolation of Dirichlet data;
c                         at the moment, the routine is hardwired
c                         to the specific problem, it should
c                         be generalized to a system routine
c
c     arguments:
c
c     in:
c             Mdle      - an element middle node number, identified
c                         with the element
c     out:
c             ErrorH, ErrorE, ErrorV - element errors
c
c---------------------------------------------------------------------
c
      subroutine elem_interp_error(Mdle, ErrorH,ErrorE,ErrorV)
c
      use control
      use parametersDPG
      use element_data
      use data_structure3D
      ! use DPGHcurl
#include "syscom.blk"
c
      character(len=4) :: etype,ftype
c
c  ...element order, orientation for edges and faces
      dimension norder(19),norient_edge(12),norient_face(6)
c
c  ...face order
      dimension norderf(5)
c
c  ...geometry dof
      dimension xnod(3,MAXbrickH)
c
c  ...solution dof
      dimension zdofH(MAXEQNH,MAXbrickH),
     .          zdofE(MAXEQNE,MAXbrickE),
     .          zdofV(MAXEQNV,MAXbrickV),
     .          zdofQ(MAXEQNQ,MAXbrickQ)
c
c  ...geometry
      dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3),
     .          dxdt(3,2),rt(3,2),rn(3),t(2)
c
c  ...H1 shape functions
      dimension shapH(MAXbrickH),gradH(3,MAXbrickH)
c
c  ...Hcurl shape functions
      dimension shapE(3,MAXbrickE),curlE(3,MAXbrickE)
c
c  ...Hdiv shape functions
      dimension shapV(3,MAXbrickV),divV(MAXbrickV)
c
c  ...2D quadrature data
      dimension tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
c
c  ...derivatives wrt physical coordinates
      dimension dv1(3),dv2(3)
c
c  ...approximate solution
      dimension zgradHxi(3),zgradH(3),
     .          zsolExi(3), zsolE(3),
     .          zsolVxi(3), zsolV(3)
c
c  ...tangential trace
      dimension zsolE_trace(3),zvalE_trace(3)
c
c  ...Dirichlet data
      dimension  zvalH(  MAXEQNH    ),
     .          zdvalH(  MAXEQNH,3  ),
     .           zvalE(3,MAXEQNE    ),
     .          zdvalE(3,MAXEQNE,3  ),
     .           zvalV(3,MAXEQNV    ),
     .          zdvalV(3,MAXEQNV,3  )
c
c  ...BC's flags
      dimension ibc(6,NR_PHYSA)
c
c---------------------------------------------------------------------
c
      select case(Mdle)
      case(1)
        iprint=0
      case default
        iprint=0
      end select
c
c  ...element type
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
c
c  ...case number
      icase = NODES(Mdle)%case
c
c  ...determine order of approximation
      call find_order(Mdle, norder)
c
c  ...determine edge and face orientations
      call find_orient( Mdle, norient_edge,norient_face)
c
c  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
c
c  ...determine solution dof
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
      if (iprint.eq.1) then
c
c  .....print out the lowest order dof
        write(*,7020) xnod(1,1:nrv),xnod(2,1:nrv),xnod(3,1:nrv)
 7020   format('elem_interp_error: xnod  = ',8(f8.3,2x),
     . 2(    /,'                           ',8(f8.3,2x)))
        if (NREQNH(icase).gt.0)
     .  write(*,7030) zdofH(1,1:nrv)
#if C_MODE
 7030   format('elem_interp_error: zdofH = ',10(/,6(2e11.4,2x)))
#else
 7030   format('elem_interp_error: zdofH = ',10(/,6( e12.5,2x)))
#endif
        if (NREQNE(icase).gt.0)
     .  write(*,7040) zdofE(1,1:nre)
#if C_MODE
 7040   format('elem_interp_error: zdofE = ',10(/,6(2e11.4,2x)))
#else
 7040   format('elem_interp_error: zdofE = ',10(/,6( e12.5,2x)))
#endif
        if (NREQNV(icase).gt.0)
     .  write(*,7050) zdofV(1,1:nrf)
#if C_MODE
 7050   format('elem_interp_error: zdofV = ',10(/,6(2e11.4,2x)))
#else
 7050   format('elem_interp_error: zdofV = ',10(/,6( e12.5,2x)))
#endif
      endif
c
c  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
      if (iprint.ge.1) then
        write(*,7001) Mdle
 7001   format('elem_inter_error: BC FLAGS FOR Mdle = ',i5)
        do i=1,NR_PHYSA
          write(*,7002) PHYSA(i), ibc(1:nrf,i)
 7002     format('                        ATTRIBUTE = ',a6,
     .           ' FLAGS = ',6i2)
        enddo
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...boundary integrals
      ErrorH = 0.d0; ErrorE = 0.d0; ErrorV = 0.d0
c
c  ...loop through element faces
      do if=1,nrf
        iflag=0
        do i=1,NR_PHYSA
          if (ibc(if,i).eq.1) iflag=1
          call locate(ibc(if,i),DIRICHLET_LIST,NR_DIRICHLET_LIST, loc)
          if (loc.ne.0) iflag=1
        enddo
c
c  .....skip if no Dirichlet BC have been found
        if (iflag.eq.0) cycle
c
c  .....sign factor to determine the OUTWARD normal unit vector
        nsign = nsign_param(etype,if)
c
c  .....face type
        ftype = face_type(etype,if)
c
c  .....face order of approximation
        call face_order(etype,if,norder, norderf)
c
c  .....set 2D quadrature; overintegrate
        INTEGRATION = 2
        call set_2Dint(ftype,norderf, nint,tloc,wtloc)
        INTEGRATION = 0
c
c  .....loop through integration points
        do l=1,nint
c
c  .......face coordinates
          t(1:2) = tloc(1:2,l)
c
c  .......face parametrization
          call face_param(etype,if,t, xi,dxidt)
c
c  .......determine element H1 shape functions
          call shape3H(etype,xi,norder,norient_edge,norient_face,
     .                 nrdofH,shapH,gradH)
c
c  .......determine element Hcurl shape functions
          call shape3E(etype,xi,norder,norient_edge,norient_face,
     .                 nrdofE,shapE,curlE)
c
c  .......determine element Hdiv shape functions
          call shape3V(etype,xi,norder,norient_face,
     .                 nrdofV,shapV,divV)
c
c  .......geometry
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,
     .                 x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
          weight = bjac*wtloc(l)
c
c  .......compute approximate trace at the point
          if (NREQNH(icase).gt.0) then
            zsolH = ZERO
            do k=1,nrdofH
              zsolH = zsolH + zdofH(1,k)*shapH(k)
            enddo
          endif
c
c  .......compute approximate tangential trace at the point
          if (NREQNE(icase).gt.0) then
            zsolExi(1:3) = ZERO
            do k=1,nrdofE
              zsolExi(1:3) = zsolExi(1:3) + zdofE(1,k)*shapE(1:3,k)
            enddo
            zsolE(1:3) = zsolExi(1)*dxidx(1,1:3)
     .                 + zsolExi(2)*dxidx(2,1:3)
     .                 + zsolExi(3)*dxidx(3,1:3)
            call zcross_product(rn,zsolE, zsolE_trace)
          endif
c
c  .......compute approximate normal trace at the point
          if (NREQNV(icase).gt.0) then
            zsolVxi(1:3) = ZERO
            do k=1,nrdofV
              zsolVxi(1:3) = zsolVxi(1:3) + zdofV(1,k)*shapV(1:3,k)
            enddo
            zsolV(1:3) = (dxdxi(1:3,1)*zsolVxi(1)
     .                   +dxdxi(1:3,2)*zsolVxi(2)
     .                   +dxdxi(1:3,3)*zsolVxi(3))/rjac
            zsolVn = ZsolV(1)*rn(1)+ZsolV(2)*rn(2)+ZsolV(3)*rn(3)
          endif
c
          call dirichlet(mdle,x,icase,
     .                   zvalH,zdvalH, zvalE,zdvalE, zvalV,zdvalV)
          call zcross_product(rn,zvalE(1:3,1), zvalE_trace)
          zvalVn = zvalV(1,1)*rn(1)+zvalV(2,1)*rn(2)+zvalV(3,1)*rn(3)
c
          if (iprint.eq.1) then
            write(*,8040) if,x(1:3)
 8040       format('elem_interp_error: if = ',i6,' x = ',3f8.3)
            if (NREQNH(icase).gt.0) then
              write(*,8050) zsolH
 8050         format('  approximate H trace = ',2e12.5)
              write(*,8055) zvalH(1)
 8055         format('  exact H trace       = ',2e12.5)
            endif
            if (NREQNE(icase).gt.0) then
              write(*,8060) zsolE_trace(1:3)
              write(*,8065) zvalE_trace(1:3)
#if C_MODE
 8060         format('  approximate E trace = ',3(2e12.5,2x))
 8065         format('  exact E trace       = ',3(2e12.5,2x))
#else
 8060         format('  approximate E trace = ',3( e12.5,2x))
 8065         format('  exact E trace       = ',3( e12.5,2x))
#endif
            endif
            if (NREQNV(icase).gt.0) then
              write(*,8070) zsolVn
 8070         format('  approximate flux  = ',2e12.5)
              write(*,8075) zvalVn
 8075         format('  exact flux        = ',2e12.5)
            endif
          endif
c
c  .......accumulate for the interpolation errors
          if (ibc(if,1).ne.0) then
            ErrorE = ErrorE
     .             + (abs(zsolE_trace(1)-zvalE_trace(1))**2
     .              + abs(zsolE_trace(2)-zvalE_trace(2))**2
     .              + abs(zsolE_trace(3)-zvalE_trace(3))**2
     .               )*weight
          endif
        enddo
        if (iprint.eq.1) call pause
      enddo
c
c
      if (iprint.ge.1) then
        write(*,7010) Mdle, ErrorH,ErrorE,ErrorV
 7010   format('elem_interp_error: Mdle, ErrorH,ErrorE,ErrorV = ',
     .         i5,3x,3e12.5)
        call pause
      endif
c
      end subroutine elem_interp_error

c--------------------------------------------------------------------
c
c     routine name      - compute_BC_interp_error
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 14
c
c     purpose:          - routine returns global BC interpolation
c                         error for trace and flux BC data
c                         (a verification for interpolating
c                          and exact solution routines)
c
c---------------------------------------------------------------------
c
      subroutine compute_BC_interp_error
c
      use data_structure3D
      use environment      , only : QUIET_MODE
#include "syscom.blk"
c
c  ...visitation flag
      integer, save :: ivis = 0
c
c  ...errors for the old mesh
      real*8,  save :: errorHold, errorEold, errorVold
c
c  ...residuals and rates to display
      real*8 , dimension(6,10), save :: rwork
c
      iprint=0
c
c  ...compute total residual and number of dof
      errorH =  0.d0;errorE =  0.d0; errorV = 0.d0
      mdle = 0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_interp_error(mdle, derrorH,derrorE,derrorV)
        if ((derrorH.eq.0.d0).and.(derrorE.eq.0.d0).and.
     .      (derrorV.eq.0.d0)) cycle
        if (iprint.eq.1) then
          write(*,7010) iel, mdle, derrorH,derrorE,derrorV
 7010     format('compute_interp_error: iel, mdle = ',2i5,
     .           ' element trace errors = ',3e12.5)
        endif
        errorH = errorH + derrorH
        errorE = errorE + derrorE
        errorV = errorV + derrorV
      enddo
      errorH = sqrt(errorH)
      errorE = sqrt(errorE)
      errorV = sqrt(errorV)
c
c  ...compute rate
      rateH = 0.d0; rateE = 0.d0; rateV = 0.d0
      if (ivis.ne.0) then
        rateH = log(errorHold/errorH)/log(2.d0)
        rateE = log(errorEold/errorE)/log(2.d0)
        rateV = log(errorVold/errorV)/log(2.d0)
      endif
c
c  ...save current data
      ivis = ivis+1
      errorHold = errorH
      errorEold = errorE
      errorVold = errorV
c
c  ...store data to display
      rwork(1,ivis) = errorH; rwork(2,ivis) = rateH
      rwork(3,ivis) = errorE; rwork(4,ivis) = rateE
      rwork(5,ivis) = errorV; rwork(6,ivis) = rateV
c
c  ...display the convergence history
      if (.NOT. QUIET_MODE) then
        write(*,*)''
        write(*,*)'         -- Error Report --'
        write(*,7100)
 7100   format(' Mesh,      errorH,  rateH,'
     .         '       errorE,  rateE, ',
     .         '       errorV,  rateV ')
        do i=1,ivis
          write(*,7110)i,rwork(1:6,i)
 7110     format(i3,3x,3(e12.5,f8.2,2x))
        enddo
        write(*,*)''
      endif
c
c
      end subroutine compute_BC_interp_error
c
c--------------------------------------------------------------------
c
c     routine name      - zcross_product
c
c--------------------------------------------------------------------
c
c     latest revision:  - May 14
c
c     purpose:          - compute cross product of real and
c                         (possibly) complex valued vectors
c
c---------------------------------------------------------------------
c
      subroutine zcross_product(A, Zb, Zc)
#include "syscom.blk"
      dimension A(3),Zb(3), Zc(3)
c
      Zc(1) =   A(2)*Zb(3) - A(3)*Zb(2)
      Zc(2) = - A(1)*Zb(3) + A(3)*Zb(1)
      Zc(3) =   A(1)*Zb(2) - A(2)*Zb(1)
c
      end subroutine zcross_product







