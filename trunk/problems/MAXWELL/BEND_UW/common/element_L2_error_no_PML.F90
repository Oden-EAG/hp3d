#include "typedefs.h"
!
subroutine element_L2_error_no_PML(Mdle,Flag,errorQ,rnormQ)
!
      use control          , only : INTEGRATION
      use data_structure3D
      use environment      , only : L2PROJ
      use physics
!
      implicit none
      integer, intent(in)  :: Flag(NR_PHYSA)
      integer, intent(in)  :: Mdle
      real(8), intent(out) :: errorQ
      real(8), intent(out) :: rnormQ
!
!     node case (decimal form)
      integer,dimension(NR_PHYSA) :: icased
!
!     element, face order, geometry dof
      integer,dimension(19)          :: norder
      real(8),dimension(3,MAXbrickH) :: xnod
      integer,dimension(12)          :: nedge_orient
      integer,dimension(6)           :: nface_orient
!
!     geometry
      real(8),dimension(3)   :: xi,x
      real(8),dimension(3,3) :: dxidx,dxdxi
      real(8)                :: rjac
!
!     3D quadrature data
      real(8),dimension(3,MAX_NINT3) :: xiloc
      real(8),dimension(  MAX_NINT3) :: wxi
!
!     approximate solution dof's
      VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
      VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
      VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
      VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!
!     approximate solution
      VTYPE, dimension(  MAXEQNH  ) ::  zsolH
      VTYPE, dimension(  MAXEQNH,3) :: zdsolH
      VTYPE, dimension(3,MAXEQNE  ) ::  zsolE
      VTYPE, dimension(3,MAXEQNE  ) :: zcurlE
      VTYPE, dimension(3,MAXEQNV  ) ::  zsolV
      VTYPE, dimension(  MAXEQNV  ) ::  zdivV
      VTYPE, dimension(  MAXEQNQ  ) ::  zsolQ
!
!     exact solution
      VTYPE, dimension(  MAXEQNH    ) ::   zvalH
      VTYPE, dimension(  MAXEQNH,3  ) ::  zdvalH
      VTYPE, dimension(  MAXEQNH,3,3) :: zd2valH
      VTYPE, dimension(3,MAXEQNE    ) ::   zvalE
      VTYPE, dimension(3,MAXEQNE,3  ) ::  zdvalE
      VTYPE, dimension(3,MAXEQNE,3,3) :: zd2valE
      VTYPE, dimension(3,MAXEQNV    ) ::   zvalV
      VTYPE, dimension(3,MAXEQNV,3  ) ::  zdvalV
      VTYPE, dimension(3,MAXEQNV,3,3) :: zd2valV
      VTYPE, dimension(  MAXEQNQ    ) ::   zvalQ
      VTYPE, dimension(  MAXEQNQ,3  ) ::  zdvalQ
      VTYPE, dimension(  MAXEQNQ,3,3) :: zd2valQ
!
!     miscellanea
      integer :: nint,icase,iattr,l,i,j,ibeg,iflag,iload,icomp,ndom,ivar,nflag,activePML
      real(8) :: weight,wa
!
!---------------------------------------------------------------------------------------
!
!     initialize global quantities
      errorQ=0.d0 ; rnormQ=0.d0
!
!     order of approx, orientations, geometry dof's, solution dof's
      call find_order( Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor(     Mdle, xnod)
      call solelm(     Mdle, zdofH,zdofE,zdofV,zdofQ)
!
!     set up the element quadrature
      INTEGRATION=2
      call set_3Dint(NODES(Mdle)%ntype,norder, nint,xiloc,wxi)
      INTEGRATION=0
!
!     supported physical attributes
      icase=NODES(Mdle)%case
      call decod(icase,2,NR_PHYSA, icased)
!
!     loop over physical attributes
      do iattr=1,NR_PHYSA
!
!       if the error not needed, skip
        if (Flag(iattr) == 0) cycle
!
!       if attribute is absent, skip
        if (icased(iattr) == 0) cycle
!
!       address of the 1st component for the attribute
        ibeg=ADRES(iattr)
!
!===================================================================================
!  L2 ATTRIBUTE                                                                    |
!===================================================================================
!
!         loop through integration points
        do l=1,nint
!
!           Gauss point and weight
          xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
          nflag=1
          call soleval(Mdle,xi,nedge_orient,nface_orient,norder,   &
                        xnod,zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi, &
                        zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)

          ! We check whether the integration point lies in a PML zone. 
          ! This assumes all quadrature is placed in the interior of the elements, so it may fail if the point belongs to an element boundary
          call is_pml(Mdle,x,activePML)
          ! we proceed with the integration only if activePML==0, otherwise we leave
          if (activePML.gt.0) goto 999
!
!           -- EXACT SOLUTION --
          call exact(x,icase, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, &
                              zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
!
!           Jacobian
          call geom(dxdxi, dxidx,rjac,iflag)
          if (iflag /= 0) then
            call find_domain(Mdle, ndom)
            write(*,9997) Mdle,ndom,rjac
9997        format(' element_L2_error_no_PML: mdle,ndom,rjac = ',i8,2x,i2,2x,e12.5)
          endif

!           total weight
          weight=wa*rjac
!
!           loop over rhs's
          do iload=1,NRRHS
!
!             loop over components of the physical attribute
            do icomp=1,NR_COMP(iattr)
!
              i=(iload-1)*NRQVAR+ibeg+icomp
!
!               accumulate L2 norm
              ErrorQ = ErrorQ + abs(zvalQ(i) - zsolQ(i))**2 * weight
              RnormQ = RnormQ + abs(zvalQ(i)           )**2 * weight
!
            enddo
          enddo
!
!         loop over integration points
        enddo
!
!       loop over physical attributes
        enddo
!
 999 continue
!
end subroutine element_L2_error_no_pml
