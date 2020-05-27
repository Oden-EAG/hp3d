c----------------------------------------------------------------------
c
c   routine name       - disp_vtk
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Aug 07
c
c   purpose            - routine selects a quantity to be displayed
c                        on the screen; this is a sample user routine
c                        that should be customized for each particular
c                        application
c
c   arguments :
c     in:
c               Mdle   - element (middle node) number
c               Xp     - coordinates of a point
c               Icase  - case number
c               ZsolH  - approximate H1 solution vector at the point
c               ZdsolH - derivatives of the approximate solution
c                        wrt physical coordinates
c               ZsolE  - approximate H(curl) solution vector at the point
c               ZcurlE - the corresponding curl
c               ZsolV  - approximate H(div) solution vector at the point
c               ZdivE  - the corresponding div
c               ZsolQ  - approximate L2 solution vector at the point
c     out:
c               Val    - the quantity to display
c
c
c----------------------------------------------------------------------
c
      subroutine disp_vtk(Mdle,X,Iphys,Nrcomp,
     .                    ZH,ZgradH,ZE,ZcurlE,ZV,ZdivV,ZQ,
     .                    Zval)
c
      use data_structure3D
c
#include "syscom.blk"
c
      dimension X(2),
     .          ZH(NRHVAR),  ZgradH(NRHVAR,3),
     .          ZE(NREVAR,3),ZcurlE(NREVAR,3),
     .          ZV(NRVVAR,3),ZdivV(NRVVAR),
     .          ZQ(NRQVAR),
     .          Zval(Nrcomp)
c
c-----------------------------------------------------------------------
c
      if (Iphys>=1 .and. Iphys<=NR_PHYSA) then
c
        select case(DTYPE(Iphys))
        case('contin')
          do i=1,NR_COMP(Iphys)
            Zval(i) = ZH(ADRES(Iphys)+i)
          enddo
        case('tangen')
          Zval(:) = ZE(ADRES(Iphys)+1,:)
        case('normal')
          Zval(:) = ZV(ADRES(Iphys)+1,:)
        case('discon')
          do i=1,NR_COMP(Iphys)
            Zval(i) = ZQ(ADRES(Iphys)+i)
          enddo
        case default
          write(*,*)'disp_vtk: DTYPE(Iphys) =',DTYPE(Iphys)
          stop
        end select
c
c       for elasticity, add pressure as last component
        if (PHYSA(Iphys)=='elast' .and. NODES(Mdle)%case==2) then
c
          call getmat_elast(mdle, omega,ZEYoung,rnu,rho,ufactor)
c
          rmu     = dble(zEYoung)/(2.d0*(1.d0+rnu))
          lambda  = dble(zEYoung)*rnn/((1.d0+rnu)*(1.d0-2.d0*rnu))
          rk      = lambda + 2.d0*rmu/3.d0
c
          Zval(NR_COMP(Iphys)+1) = -rk*(ZgradH(1,1)
     .                           +      ZgradH(2,2)
     .                           +      ZgradH(3,3))
c
        endif
c
      else
c
c       GLOBAL PRESSURE
c
        select case(NODES(Mdle)%case)
c
c       acoustics
        case(1)

          Zval(1) = ZH(4)
c
c       elasticity
        case(2)
c
          call getmat_elast(mdle, omega,ZEYoung,rnu,rho,ufactor)
c
          rmu     = dble(zEYoung)/(2.d0*(1.d0+rnu))
          lambda  = dble(zEYoung)*rnn/((1.d0+rnu)*(1.d0-2.d0*rnu))
          rk      = lambda + 2.d0*rmu/3.d0
c
          Zval(1) = -rk*(ZgradH(1,1)+ZgradH(2,2)+ZgradH(3,3))
c
        case default
          write(*,*)'disp_vtk: UNKNOWN case =',NODES(Mdle)%case
          stop
        end select
c
      endif
c
      end subroutine disp_vtk
