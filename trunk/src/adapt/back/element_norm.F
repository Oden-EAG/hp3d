c///////////////////////////////////////////////////////////////////////
c P. Gatto, Dec 2014 : prehistroric routine, probably copied from old
c                      code! Do not use it, unless you revise it!
c///////////////////////////////////////////////////////////////////////
c-----------------------------------------------------------------------
c
c     routine name      - element_norm
c
c-----------------------------------------------------------------------
c
c     latest revision:  - Jan 09
c
c     purpose:          - computes contributions to the norm (squared)
c                         of FE solution over an element
c
c     arguments:
c     in:
c              Mdle     - element middle node number
c              Nr_RHS   - number of loads (RHS)
c              Nshift   = 0 there is one solution only
c                       > 0 both fine grid solution u_hP and coarse
c                           grid solution u_Hp are available,
c                           Nshift indicates the location of the coarse
c                           grid solution
c     out:
c              Dnorm    - norms listed below for all physical attributes
c                         (summed up for all load components)
c                         1 - |u_hP|^2
c                         2 - |u_ex|^2
c                         3 - |u_ex - u_hP|^2
c                         4 - |u_ex - u_Hp|^2
c                         5 - |u_hP - u_Hp|^2
c
c-----------------------------------------------------------------------
c
      subroutine element_norm(Mdle,Nr_RHS,Nshift,Idec, Dnorm)
c
      use control
      use data_structure3D
      use element_data
#include "syscom.blk"
c
      dimension Dnorm(5,NR_PHYSA)
c
c  ...element order
      dimension norder(19)
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
c  ...node orientations
      dimension nedge_orient(12), nface_orient(6)
c
c  ...shape functions and their derivatives wrt master coordinates
      dimension shapH(MAXbrickH),gradH(3,MAXbrickH)
c
c  ...geometry
      dimension xi(3),x(3),dxdxi(3,3),dxidx(3,3)
c
c  ...3D quadrature data
      dimension xiloc(3,MAX_NINT3),wxi(MAX_NINT3)
c
c  ...approximate solution
      dimension zsolH(MAXEQNH),zgradH(MAXEQNH,3),
     .          zsolE(MAXEQNE,3),zcurlE(MAXEQNE,3),
     .          zsolV(MAXEQNV,3),zdivV(MAXEQNV),
     .          zsolQ(MAXEQNQ)
c
c  ...exact solution
      dimension ZvalH(MAXEQNH),
     .          ZdvalH(MAXEQNH,3),Zd2valH(MAXEQNH,3,3),
     .          ZvalE(3,MAXEQNE),
     .          ZdvalE(3,MAXEQNE,3),Zd2valE(3,MAXEQNE,3,3),
     .          ZvalV(3,MAXEQNV),
     .          ZdvalV(3,MAXEQNV,3),Zd2valV(3,MAXEQNV,3,3),
     .          ZvalQ(MAXEQNQ),
     .          ZdvalQ(MAXEQNQ,3),Zd2valQ(MAXEQNQ,3,3)
c
c  ...middle node index
      integer index(NRINDEX)
c
c  ...norms for various components
      dimension dnormH(5,MAXEQNH),dnormE(5,MAXEQNE),
     .          dnormV(5,MAXEQNV),dnormQ(5,MAXEQNQ)
c
c  ...work space
      dimension zaux(3)
c
c-----------------------------------------------------------------------
c
      if (Mdle.eq.3) then
        iprint=0
      else
        iprint=0
      endif
c
      dnormH(1:5,1:MAXEQNH) = 0.d0
      dnormE(1:5,1:MAXEQNE) = 0.d0
      dnormV(1:5,1:MAXEQNV) = 0.d0
      dnormQ(1:5,1:MAXEQNQ) = 0.d0
c
c  ...number of active variables for a mesh
      nvarH = NRHVAR*Nr_RHS; nsH = NRHVAR*Nshift
      nvarE = NREVAR*Nr_RHS; nsE = NREVAR*Nshift
      nvarV = NRVVAR*Nr_RHS; nsV = NRVVAR*Nshift
      nvarQ = NRQVAR*Nr_RHS; nsQ = NRQVAR*Nshift
c
c  ...determine order of approximation
      call find_order(Mdle, norder)
c
c  ...determine the node orientation
      call find_orient(Mdle, nedge_orient,nface_orient)
c
c  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
c
c  ...determine solution dof
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
c
c  ...set up the element quadrature
      call set_3Dint(NODES(Mdle)%type,norder, nint,xiloc,wxi)
c
c  ...loop through integration points
      do l=1,nint
        xi(1:3) = xiloc(1:3,l); wa = wxi(l)
c
c  .....evaluate derivatives and values of the shape functions
c       at the point
        call shape3H(NODES(Mdle)%type,xi,norder,
     .               nedge_orient,nface_orient,
     .               nrdofH,shapH,gradH)
c
c  .....determine physical coordinates and the derivatives of
c       the physical coordinates wrt master element coordinates
        x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0
        do k=1,nrdofH
          x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
          do i=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*gradH(i,k)
          enddo
        enddo
c
c  .....evaluate the inverse derivatives and jacobian
        call geom(dxdxi, dxidx,rjac,iflag)
        weight = wa*rjac
c
c  .....evaluate the exact solution
        if (NEXACT.ne.0)
     .    call exact(Xp,NODES(Mdle)%case,
     .               ZvalH,ZdvalH,Zd2valH,
     .               ZvalE,ZdvalE,Zd2valE,
     .               ZvalV,ZdvalV,Zd2valV,
     .               ZvalQ,ZdvalQ,Zd2valQ)
c
c  .....if there are H1 variables
        if (NRHVAR.gt.0) then
c
c  .......evaluate the approximate H1 solution
          zsolH(1:MAXEQNH) = ZERO
          zgradH(1:MAXEQNH,1:3) = ZERO
          do k=1,nrdofH
            zsolH(1:MAXEQNH) = zsolH(1:MAXEQNH)
     .                       + zdofH(MAXEQNH,k)*shapH(k)
            do ivar=1,3
              zgradH(1:MAXEQNH,ivar) = zgradH(1:MAXEQNH,ivar)
     .                               + zdofH(1:MAXEQNH,k)*gradH(ivar,k)
            enddo
          enddo
c
c  .......use chain formula to compute the gradient wrt physical
c         coordinates
          do ivar=1,MAXEQNH
            zaux(1:3) = zgradH(ivar,1:3)
            zgradH(ivar,1:3) = ZERO
            do ixi=1,3
              zgradH(ivar,1:3) = zgradH(ivar,1:3)
     .                         + zaux(ixi)*dxidx(ixi,1:3)
            enddo
          enddo
          if (iprint.eq.1) then
            do ivar=1,MAXEQNH
              write(*,7003) ivar,zsolH(ivar),zgradH(ivar,1:3)
 7003         format('element_norm: ivar = ',i3,' zsolH = ',
     .                2e12.5,' zgradH = ',3(2e12.5,2x))
            enddo
            call pause
          endif
c
c  .......norm of the fine (first) grid solution
          dnormH(1,1:nvarH) = dnormH(1,1:nvarH)
     .      + (abs(zsolH(1:nvarH))**2
     .      +  abs(zgradH(1:nvarH,1))**2
     .      +  abs(zgradH(1:nvarH,2))**2
     .      +  abs(zgradH(1:nvarH,3))**2)*weight
c
          if (NEXACT.ne.0) then
c
c  .........norm of the exact solution
            dnormH(2,1:nvarH) = dnormH(2,1:nvarH)
     .      + (abs(ZvalH(1:nvarH))**2
     .      +  abs(ZdvalH(1:nvarH,1))**2
     .      +  abs(ZdvalH(1:nvarH,2))**2
     .      +  abs(ZdvalH(1:nvarH,3))**2)*weight
c
c  .........norm of the difference between exact and fine(first) grid
c           solution
            dnormH(3,1:nvarH) = dnormH(3,1:nvarH)
     .      + (abs(ZvalH(1:nvarH)-zsolH(1:nvarH))**2
     .      +  abs(ZdvalH(1:nvarH,1)-zgradH(1:nvarH,1))**2
     .      +  abs(ZdvalH(1:nvarH,2)-zgradH(1:nvarH,2))**2
     .      +  abs(ZdvalH(1:nvarH,3)-zgradH(1:nvarH,3))**2)*weight
c
c  .........norm of the difference between exact and coarse(second) grid
c           solution
            if ((Nshift.gt.0).and.(Idec.eq.1))
     .      dnormH(4,1:nvarH) = dnormH(4,1:nvarH)
     .      + (abs(ZvalH(1:nvarH)-zsolH(nsh+1:nsh+nvarH))**2
     .      +  abs(ZdvalH(1:nvarH,1)-zgradH(nsh+1:nsh+nvarH,1))**2
     .      +  abs(ZdvalH(1:nvarH,2)-zgradH(nsh+1:nsh+nvarH,2))**2
     .      +  abs(ZdvalH(1:nvarH,3)-zgradH(nsh+1:nsh+nvarH,3))**2)
     .         *weight
          endif
c
c  .......norm of the difference between first and second grid solutions
          if ((Nshift.gt.0).and.(Idec.eq.1)) then
            dnormH(5,1:nvarH) = dnormH(5,1:nvarH)
     .      + (abs(zsolH(1:nvarH)-zsolH(nsh+1:nsh+nvarH))**2
     .      +  abs(zgradH(1:nvarH,1)-zgradH(nsh+1:nsh+nvarH,1))**2
     .      +  abs(zgradH(1:nvarH,2)-zgradH(nsh+1:nsh+nvarH,2))**2
     .      +  abs(zgradH(1:nvarH,3)-zgradH(nsh+1:nsh+nvarH,3))**2)
     .         *weight
          endif
          if (iprint.eq.1) then
            write(*,7004) dnormH(5,1:nvarH)
 7004       format('element_norm: dnormH = ',10(/,10e12.5))
            call pause
          endif
c
c  .....if there are H1 unknowns
        endif
c
c  .....if there are H(curl) variables
        if (NREVAR.gt.0) then
          write(*,*) 'element_norm: UNFINISHED'; stop 1
c
c  .....if there are H(curl) unknowns
        endif
c
c  .....if there are H(div) variables
        if (NRVVAR.gt.0) then
          write(*,*) 'element_norm: UNFINISHED'; stop 1
c
c  .....if there are H(div) unknowns
        endif
c
c  .....if there are L2 variables
        if (NRQVAR.gt.0) then
          write(*,*) 'element_norm: UNFINISHED'; stop 1
c
c  .....if there are L2 unknowns
        endif
c
c  ...end of loop through integration points
      enddo
c
c-----------------------------------------------------------------------
c
      Dnorm(1:5,1:NR_PHYSA) = 0.d0
c
c  ...get the middle node index
      call get_index(mdle, index)
c
c  ...initiate component number
      ivar=0; ivarH=0; ivarE=0; ivarV=0; ivarQ=0
c
c  ...loop through the physical attributes
      do iphys=1,NR_PHYSA
c
c  .....loop through components
        do icomp=1,NR_COMP(iphys)
c
c  .......update the global component
          ivar=ivar+1
c
          select case(DTYPE(iphys))
          case('contin')
            ivarH=ivarH+1
            if (index(ivar).ne.0) then
c
c  ...........loop through the loads
              do l=1,Nr_RHS
                Dnorm(1,iphys) = Dnorm(1,iphys)
     .                         + dnormH(1,(l-1)*NRHVAR+ivarH)
                if (NEXACT.ne.0) then
                  Dnorm(2,iphys) = Dnorm(2,iphys)
     .                           + dnormH(2,(l-1)*NRHVAR+ivarH)
                  Dnorm(3,iphys) = Dnorm(3,iphys)
     .                           + dnormH(3,(l-1)*NRHVAR+ivarH)
                  if ((Nshift.gt.0).and.(Idec.eq.1)) then
                    Dnorm(4,iphys) = Dnorm(4,iphys)
     .                             + dnormH(4,(l-1)*NRHVAR+ivarH)
                  endif
                endif
                if ((Nshift.gt.0).and.(Idec.eq.1)) then
                  Dnorm(5,iphys) = Dnorm(5,iphys)
     .                           + dnormH(5,(l-1)*NRHVAR+ivarH)
                endif
              enddo
            endif
          case('tangen')
          case('normal')
          case('discon')
          end select
        enddo
      enddo
c
      if (iprint.ge.1) then
        write(*,7010) Mdle,Index(1:NRINDEX)
 7010   format('element_norm: Mdle, Index = ',i6,3x,10i2)
        do iphys=1,NR_PHYSA
          write(*,7011) iphys,Dnorm(1:5,iphys)
 7011     format('              iphys = ',i2,' Dnorm = ',5e12.5)
        enddo
        call pause
      endif
c
c
      end
