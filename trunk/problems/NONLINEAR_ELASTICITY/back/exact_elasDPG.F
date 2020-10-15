c-----------------------------------------------------------------------
c
c   routine name       - exact
c
c-----------------------------------------------------------------------
c
c   latest revision    - Feb 11
c
c   purpose            - routine evaluates an exact solution for 
c                        the elasticity problem
c
c   arguments :
c     in:
c              Xp      - coordinates of a point
c     out:
c              ZvalH   - value of the H1 solution
c              ZdvalH  - corresponding first derivatives
c              Zd2valH - corresponding second derivatives
c              ZvalE   - value of the H(curl) solution
c              ZdvalE  - corresponding first derivatives
c              Zd2valE - corresponding second derivatives
c              ZvalV   - value of the H(div) solution
c              ZdvalV  - corresponding first derivatives
c              Zd2valV - corresponding second derivatives
c              ZvalQ   - value of the L2 solution
c              ZdvalQ  - corresponding first derivatives
c              Zd2valQ - corresponding second derivatives
c
c-----------------------------------------------------------------------
c
      subroutine exact(Xp,
     .                 ZvalH,ZdvalH,Zd2valH,
     .                 ZvalE,ZdvalE,Zd2valE,
     .                 ZvalV,ZdvalV,Zd2valV,
     .                 ZvalQ,ZdvalQ,Zd2valQ)
c
      use data_structure2D
      use control
      use problem
c
#include "syscom.blk"
      common /cexact/ iprint
c
      dimension Xp(2),
     .          ZvalH(MAXEQNH),
     .          ZdvalH(MAXEQNH,2),Zd2valH(MAXEQNH,2,2),
     .          ZvalE(2,MAXEQNE),
     .          ZdvalE(2,MAXEQNE,2),Zd2valE(2,MAXEQNE,2,2),
     .          ZvalV(2,MAXEQNV),
     .          ZdvalV(2,MAXEQNV,2),Zd2valV(2,MAXEQNV,2,2),
     .          ZvalQ(MAXEQNQ),
     .          ZdvalQ(MAXEQNQ,2),Zd2valQ(MAXEQNQ,2,2)
c
c  ...physics flags decoded
      dimension ncase(NR_PHYSA)
      dimension rstress1(2), rstress2(2)
      REAL, PARAMETER :: pi = 3.1415927
c
c
      iprint=0
c
      ZvalH=ZERO; ZdvalH=ZERO; Zd2valH=ZERO
      ZvalE=ZERO; ZdvalE=ZERO; Zd2valE=ZERO
      ZvalV=ZERO; ZdvalV=ZERO; Zd2valV=ZERO
      ZvalQ=ZERO; ZdvalQ=ZERO; Zd2valQ=ZERO
c
      x=Xp(1);
      y=Xp(2);
c-----------------------------------------------------------------------
c  
c  ...select form of displacement
c
c     1 - constant
c     2 - linear
c     3 - bilinear
c     4 - quadratic
c     5 - L-shaped domain
c     6 - sinusoid
c
      iselect = 6
c
c
c
      select case (iselect)
c
      case(1)
        ZvalQ(4) = 1.d0
        ZvalQ(5) = 2.d0
c
      case(2)
c
        ZvalQ(4) = 2.d0*x+y
        ZvalQ(5) = y
        ZdvalQ(4,1) = 2.d0
        ZdvalQ(4,2) = 1.d0
        ZdvalQ(5,2) = 1.d0
c
      case(3)
        ZvalQ(4) = x*y
        ZvalQ(5) = 0.d0
        ZdvalQ(4,1) = y
        ZdvalQ(4,2) = x
        Zd2valQ(4,1,2) = 1.d0
        Zd2valQ(4,2,1) = 1.d0
c
      case(4)
        ZvalQ(4) = x**2+y**2
        ZvalQ(5) = x**2
        ZdvalQ(4,1) = 2.d0*x
        ZdvalQ(4,2) = 2.d0*y
        ZdvalQ(5,1) = 2.d0*x
        Zd2valQ(4,1,1) = 2.d0
        Zd2valQ(4,2,2) = 2.d0
        Zd2valQ(5,1,1) = 2.d0
c       
      case(5)
        call singular_sol_df(x,y, rstress1,rdiv1,rstress2,rdiv2,
     .        solQ1,solQ2,solQ3,solu11,solu12,solu21,solu22)
        ZvalQ(4) = solQ1
        ZvalQ(5) = solQ2
        ZvalQ(6) = solQ3
        ZdvalQ(4,1) = solu11
        ZdvalQ(4,2) = solu12
        ZdvalQ(5,1) = solu21
        ZdvalQ(5,2) = solu22
c   
      case(6)
        ZvalQ(4) = sin(pi*x)*sin(pi*y)
        ZdvalQ(4,1) = pi*cos(pi*x)*sin(pi*y)
        ZdvalQ(4,2) = pi*sin(pi*x)*cos(pi*y)
        Zd2valQ(4,1,1) = -1.d0*(pi**2)*sin(pi*x)*sin(pi*y)
        Zd2valQ(4,2,2) = -1.d0*(pi**2)*sin(pi*x)*sin(pi*y)
        Zd2valQ(4,1,2) = (pi**2)*cos(pi*x)*cos(pi*y)
        Zd2valQ(4,2,1) =(pi**2)*cos(pi*x)*cos(pi*y)
        ZvalQ(5) = sin(pi*x)*sin(pi*y)
        ZdvalQ(5,1) = pi*cos(pi*x)*sin(pi*y)
        ZdvalQ(5,2) = pi*sin(pi*x)*cos(pi*y)
        Zd2valQ(5,1,1) = -1.d0*(pi**2)*sin(pi*x)*sin(pi*y)
        Zd2valQ(5,2,2) = -1.d0*(pi**2)*sin(pi*x)*sin(pi*y)
        Zd2valQ(5,1,2) = (pi**2)*cos(pi*x)*cos(pi*y)
        Zd2valQ(5,2,1) =(pi**2)*cos(pi*x)*cos(pi*y)
      end select
c
      if (iselect.eq.5) then
        ZvalQ(1) = rstress1(1)
        ZvalQ(2) = rstress1(2)
        ZvalQ(3) = rstress2(2)
        if (iprint.eq.1) then
          write(*,7010) 1,ZvalQ(1:3)
 7010     format('exact: STRESSES = ',i2,2x,3e12.5)
        endif
c
      else 
c  .....sigma_x
        ZvalQ(1) = 2.d0*RMU*ZdvalQ(4,1)
     .              +RLAMBDA*(ZdvalQ(4,1)+ZdvalQ(5,2))
        ZdvalQ(1,1) = 2.d0*RMU*Zd2valQ(4,1,1)
     .              +RLAMBDA*(Zd2valQ(4,1,1)+Zd2valQ(5,2,1))
        ZdvalQ(1,2) = 2.d0*RMU*Zd2valQ(4,1,2)
     .              +RLAMBDA*(Zd2valQ(4,1,2)+Zd2valQ(5,2,2))
c
c  .....sigma_xy
c
        ZvalQ(2) = RMU*(ZdvalQ(4,2)+ZdvalQ(5,1))
        ZdvalQ(2,1) = RMU*(Zd2valQ(4,2,1)+Zd2valQ(5,1,1))
        ZdvalQ(2,2) = RMU*(Zd2valQ(4,2,2)+Zd2valQ(5,1,2))
c
c  .....sigma_y
c
        ZvalQ(3) = 2.d0*RMU*ZdvalQ(5,2)
     .              +RLAMBDA*(ZdvalQ(4,1)+ZdvalQ(5,2))
        ZdvalQ(3,1) = 2.d0*RMU*Zd2valQ(5,2,1)
     .              +RLAMBDA*(Zd2valQ(4,1,1)+Zd2valQ(5,2,1))
        ZdvalQ(3,2) = 2.d0*RMU*Zd2valQ(5,2,2)
     .              +RLAMBDA*(Zd2valQ(4,1,2)+Zd2valQ(5,2,2))
c
c  .....p
c
        ZvalQ(6) = 0.5*(ZdvalQ(4,2)-ZdvalQ(5,1))
        ZdvalQ(6,1) = 0.5*(Zd2valQ(4,2,1)-Zd2valQ(5,1,1))
        ZdvalQ(6,2) = 0.5*(Zd2valQ(4,2,2)-Zd2valQ(5,1,2))
c
        if (iprint.eq.1) then
          write(*,7010) 2,ZvalQ(1:3)
          call pause
        endif
      endif
c
c
      if (NEXACT.eq.0) then
        write(*,*) 'exact: EXACT SOLUTION UNKNOWN'
        stop1
      endif
c
c
c
c
      if (iprint.eq.1) then
        write(*,7001) ncase
 7001   format('exact: ncase = ',2i2)
        write(*,7002) Xp
 7002   format('exact: Xp = ',2f8.3)
        do ivar=1,1
          write(*,7003) ivar, ZvalQ(ivar)
 7003     format('exact: ivar,ZvalQ(ivar) = ',i2,2x,2e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine

c-----------------------------------------------------------------------
c
c   routine name       - singular_sol_df
c
c-----------------------------------------------------------------------
c
c   latest revision    - Jun 10
c
c   purpose            - routine evaluates a manufactured (exact)
c                        solution to the L-shape linear elasticity
c                        problem
c   arguments :
c     in:
c              tnx,tny  - coordinates of a point 
c     out:
c              solV1    = (\sigma_xx, \sigma_xy)
c              Soldiv1  = \sigma_xx,x + \sigma_xy,y
c              SolV2    = (\sigma_yx, \sigma_yy)
c              Soldiv2  = \sigma_yx,x + \sigma_yy,y
c              SolQ1    = displacement u
c              SolQ2    = displacement v
c              SolQ3    = infinitezimal rotation 
c                         \omega_xy = 1/2(u_y-v_x)
c              Solu11   = u,x
c              Solu12   = u,y
c              Solu21   = v,x
c              Solu22   = v,y
c
c-----------------------------------------------------------------------
c
      subroutine singular_sol_df(Tnx,Tny, SolV1,Soldiv1,SolV2,Soldiv2,
     .              SolQ1,SolQ2,SolQ3,Solu11,Solu12,Solu21,Solu22)
c
      use problem
#include "syscom.blk"
      common /csingular/ iprint
c
      dimension SolV1(2),SolV2(2)
      dimension Zdval1(2),Zdval2(2,2)

      SolV1(1:2)=0.d0; Soldiv1=0.d0; SolV2(1:2)=0.d0; Soldiv2=0.d0
      SolQ1=0.d0; SolQ2=0.d0; SolQ3=0.d0
c
c     rotate coordinates 45 degrees counter-clockwise
      tmp_a = sqrt(2.d0)/2.d0
      x = tmp_a*(Tnx+Tny)
      y = tmp_a*(-Tnx+Tny) 
c
      r = sqrt(x**2+y**2)
      eps = 1.d-10
      if (r.lt.eps) then
      else
c
c  .....these material constants CANNOT be changed
        zmu = 66.16380367543577d0
        zlamda = 99.3d0
        if (iprint.eq.1) then
          write(*,*) 'zmu,zlamda = ',zmu,zlamda
        endif
        zpoisson = zlamda/(2.d0*(zlamda+zmu))
        zsigma = zpoisson/(1.d0+zpoisson)        
        Gc = zmu 
        pi = acos(-1.d0)
        if (y.eq.0.d0) then
          if (x.lt.0.d0) then
            theta = pi
          elseif (x.gt.0.d0) then
            theta = 0.d0
          endif
        elseif (y.gt.0.d0) then
          theta = acos(x/r)
        else
          theta = -acos(x/r)
        endif
c
        drdx = x/r ; dthetadx = -y/r**2 
        drdy = y/r ; dthetady = x/r**2
c
        a = 0.60404d0 
c       
        alpha = 3.d0*pi/4.d0  
        ak = (4.d0*(1-zsigma)-(a+1.d0))*sin((a-1.d0)*alpha)
     .                    /((a+1.d0)*sin((a+1.d0)*alpha))

        F = ak*sin((a+1.d0)*theta)+sin((a-1.d0)*theta)
        Fd = (a+1.d0)*ak*cos((a+1.d0)*theta)
     .           +(a-1.d0)*cos((a-1.d0)*theta)
        Fdd = -(a+1.d0)**2*ak*sin((a+1.d0)*theta)
     .           -(a-1.d0)**2*sin((a-1.d0)*theta)
        G = -4.d0/(a-1.d0)*cos((a-1.d0)*theta)
        Gd = 4.d0*sin((a-1.d0)*theta)
        Gdd = 4.d0*(a-1.d0)*cos((a-1.d0)*theta)  

        b =  1.d0/(2.d0*Gc)
        ur = b*r**a*(-(a+1.d0)*F+(1-zsigma)*Gd)
        utheta = b*r**a*(-Fd+(1.d0-zsigma)*(a-1.d0)*G)
        urdr = a* b*r**(a-1.d0)*(-(a+1.d0)*F+(1-zsigma)*Gd)
        urdtheta = b*r**a*(-(a+1.d0)*Fd+(1-zsigma)*Gdd)
        uthetadr = a*b*r**(a-1.d0)*(-Fd+(1.d0-zsigma)*(a-1.d0)*G)
        uthetadtheta = b*r**a*(-Fdd+(1.d0-zsigma)*(a-1.d0)*Gd) 
        if (iprint.eq.1) then
          write(*,7001) ur,utheta,theta
 7001     format('singular: ur,utheta,theta = ',3e12.5)
        endif
c
c  .....trr,tthetatheat,trtheta are the cylindrical stress components
        trr = r**(a-1.d0)*(Fdd+(a+1.d0)*F)
        tthetatheta = a*(a+1.d0)*r**(a-1.d0)*F
        trtheta = -a*r**(a-1.d0)*Fd
c
c  .....displacements in the rotated system of coordinates
        SolQ1 = ur*cos(theta)-utheta*sin(theta)
        SolQ2 = ur*sin(theta)+utheta*cos(theta)
c
c  .....tranform displacements to the rotated system of coordinates
        tmp_solQ1 = tmp_a*(SolQ1-SolQ2)
        tmp_solQ2 = tmp_a*(SolQ1+SolQ2)
        SolQ1 = tmp_solQ1
        SolQ2 = tmp_solQ2 
c
c  .....and their derivatives in the rotated system of coordinates
        Solu11 = (urdr*cos(theta)-uthetadr*sin(theta))*drdx
     .   +(urdtheta*cos(theta)-ur*sin(theta)-uthetadtheta*sin(theta)
     .     -utheta*cos(theta))*dthetadx
        Solu12 = (urdr*cos(theta)-uthetadr*sin(theta))*drdy
     .   +(urdtheta*cos(theta)-ur*sin(theta)-uthetadtheta*sin(theta)
     .     -utheta*cos(theta))*dthetady
        Solu21 = (urdr*sin(theta)+uthetadr*cos(theta))*drdx
     .   +(urdtheta*sin(theta)+ur*cos(theta)+uthetadtheta*cos(theta)
     .     -utheta*sin(theta))*dthetadx 
        Solu22 = (urdr*sin(theta)+uthetadr*cos(theta))*drdy
     .   +(urdtheta*sin(theta)+ur*cos(theta)+uthetadtheta*cos(theta)
     .     -utheta*sin(theta))*dthetady         
c
c  .....infinitesimal rotation component
        tmp1 = (urdr*cos(theta)-uthetadr*sin(theta))*drdy
     .     +(urdtheta*cos(theta)-uthetadtheta*sin(theta))*dthetady
     .     +(-ur*sin(theta)-utheta*cos(theta))*dthetady
        tmp2 =  (urdr*sin(theta)+uthetadr*cos(theta))*drdx
     .     +(urdtheta*sin(theta)+uthetadtheta*cos(theta))*dthetadx
     .     +(ur*cos(theta)-utheta*sin(theta))*dthetadx        
        SolQ3 = (tmp1-tmp2)/2.d0
c
c  .....transform displacement derivatives to the original system
        tmp_u11 = 0.5d0*(Solu11-Solu21-Solu12+Solu22)
        tmp_u12 = 0.5d0*(Solu11-Solu21+Solu12-Solu22)
        tmp_u21 = 0.5d0*(Solu11+Solu21-Solu12-Solu22)
        tmp_u22 = 0.5d0*(Solu11+Solu21+Solu12+Solu22)
c
        Solu11 = tmp_u11
        Solu12 = tmp_u12
        Solu21 = tmp_u21
        Solu22 = tmp_u22
c
c  .....transform stresses from the cylindrical to the rotated
c       Cartesian system of coordinates 
        SolV1(1) = trr*(cos(theta)**2)
     .      -2.d0*sin(theta)*cos(theta)*trtheta
     .            +tthetatheta*(sin(theta)**2)
        SolV2(2) = trr*(sin(theta)**2)
     .      +2.d0*sin(theta)*cos(theta)*trtheta
     .            +tthetatheta*(cos(theta)**2)
        SolV1(2) = trr*sin(theta)*cos(theta)
     .      +(cos(theta)**2-sin(theta)**2)*trtheta
     .            -tthetatheta*sin(theta)*cos(theta)
        SolV2(1) = SolV1(2)            
c
c  .....transform stresses to the original system of coordinates
        tmp_v11 = SolV1(1); tmp_v12 = SolV1(2)
        tmp_v21 = SolV2(1); tmp_v22 = SolV2(2)
c
        SolV1(1) = 0.5d0*(tmp_v11-tmp_v21-tmp_v12+tmp_v22)
        SolV1(2) = 0.5d0*(tmp_v11-tmp_v21+tmp_v12-tmp_v22)
        SolV2(1) = 0.5d0*(tmp_v11+tmp_v21-tmp_v12-tmp_v22)
        SolV2(2) = 0.5d0*(tmp_v11+tmp_v21+tmp_v12+tmp_v22)
c              
c 
      endif
c

      end 


      subroutine exact_flux(x,rt,rn, zvalT,zdvalT,zvalF)
c
      use data_structure2D
      use control
      use problem
c
#include "syscom.blk"
c
c  ...normal to the boundary implied by the edge orientation
      dimension rn(2),rt(2),x(2)
      dimension zvalT(MAXEQNH), zdvalT(MAXEQNH,1), zvalF(MAXEQNV)
c
c  ...the exact solution and its derivatives
      dimension zvalH(MAXEQNH),
     .          zdvalH(MAXEQNH,2),zd2valH(MAXEQNH,2,2),
     .          zvalE(2,MAXEQNE),
     .          zdvalE(2,MAXEQNE,2),zd2valE(2,MAXEQNE,2,2),
     .          zvalV(2,MAXEQNV),
     .          zdvalV(2,MAXEQNV,2),zd2valV(2,MAXEQNV,2,2),
     .          zvalQ(MAXEQNQ),
     .          zdvalQ(MAXEQNQ,2),zd2valQ(MAXEQNQ,2,2)
c
c
        call exact(x,
     .             zvalH,zdvalH,zd2valH,
     .             zvalE,zdvalE,zd2valE,
     .             zvalV,zdvalV,zd2valV,
     .             zvalQ,zdvalQ,zd2valQ)

        zvalT(1) = -zvalQ(4)
        zvalT(2) = -zvalQ(5)
        zvalF(1) = -(zvalQ(1)*rn(1) + zvalQ(2)*rn(2))
        zvalF(2) = -(zvalQ(2)*rn(1) + zvalQ(3)*rn(2))

      end subroutine

      
