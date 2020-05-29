!
#include "typedefs.h"
!
!-------------------------------------------------------------------------------------
!> Purpose - Selection for quantities to be displayed by graphics
!!
!> @data Nov 14
!-------------------------------------------------------------------------------------
!
subroutine soldis_select_system
!
      use physics          , only : NR_PHYSA, NR_COMP, DTYPE, PHYSA
      use data_structure3D , only : NRCOMS
      use graphmod         , only : ISELECT
!
      integer :: iattr,icomp,ireal,iload
!
!-------------------------------------------------------------------------------------
!
!     print attributes info
      write(*,*)'========================================='
      write(*,*)'     ATTRIBUTE | DISC. SPACE | COMPONENTS'
      do iattr=1,NR_PHYSA
!
        write(*,1000) iattr, PHYSA(iattr), DTYPE(iattr), NR_COMP(iattr)
 1000   format(i2,' - ',a6,7x,a6,8x,i10)
!
      enddo
      write(*,*)'========================================='
      write(*,*)' '
!
!     attribute
 11   continue
      write(*,*) 'Set attribute'
      read( *,*) iattr
      if ((iattr < 1) .OR. (iattr > NR_PHYSA))  goto 11
!
!     component
 12   continue
      write(*,*) 'Set component (use negative for imaginary part)'
      read( *,*) icomp
      if ((abs(icomp) < NR_COMP(iattr)) .OR. (abs(icomp) > NR_COMP(iattr)))  goto 12
!
!     remap to 1 - real ; 2 - imaginary
      if (icomp > 0) then ; ireal=1
      else                ; ireal=2
      endif
      icomp = abs(icomp)
!
!     load
 13   continue
      write(*,1001) NRCOMS
 1001 format(' Set load (NRCOMS  = ',i2,')')
      read(*,*) iload
      if ((iload < 1) .OR. (iload > NRCOMS))  goto 13
!
!     trace
      select case(DTYPE(iattr))
      case('tangen') ; write(*,*) 'Displaying magnitude of tangential trace...'
      case('normal') ; write(*,*) 'Displaying magnitude of normal trace...'
      endselect
!
!     store selection
      ISELECT = ireal*1000 + iattr*100 + icomp*10 + iload*1
!
!
end subroutine soldis_select_system
!
!
!-------------------------------------------------------------------------------------
!> Purpose - Compute quantity to be displayed by graphics
!!
!> @param[in]  Mdle   - middle node
!> @param[in]  Xi     - master element coordinates
!> @param[in]  X      - coordinates of a point
!> @param[in]  Rn     - outward normal unit vector
!> @param[in]  ZsolH  - value of H1      solution
!> @param[in]  ZgradH - grad  of H1      solution
!> @param[in]  ZsolE  - value of H(curl) solution
!> @param[in]  ZcurlE - curl  of H(curl) solution
!> @param[in]  ZsolV  - value of H(div)  solution
!> @param[in]  ZdivV  - div   of H(div)  solution
!> @param[in]  ZsolQ  - value of L^2     solution
!> @param[out] Val    - quantity to display
!!
!> @data Nov 14
!-------------------------------------------------------------------------------------
subroutine soldis_system(Mdle,Xi,X,Rn,SolH,GradH,SolE,CurlE,SolV,DivV,SolQ, Val)
!
      use data_structure3D , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,  &
                                    NRHVAR, NREVAR, NRVVAR, NRQVAR
      use graphmod         , only : ISELECT
      use physics
!
      implicit none
      integer,                     intent(in)  :: Mdle
      real(8),dimension(3),        intent(in)  :: Xi,X,Rn
      VTYPE,dimension(  MAXEQNH  ),intent(in)  :: SolH
      VTYPE,dimension(  MAXEQNH,3),intent(in)  :: GradH
      VTYPE,dimension(3,MAXEQNE  ),intent(in)  :: SolE
      VTYPE,dimension(3,MAXEQNE  ),intent(in)  :: CurlE
      VTYPE,dimension(3,MAXEQNV  ),intent(in)  :: SolV
      VTYPE,dimension(  MAXEQNV  ),intent(in)  :: DivV
      VTYPE,dimension(  MAXEQNQ  ),intent(in)  :: SolQ
      real(8),                     intent(out) :: Val
!
!!!!     exact solution
!!!      VTYPE,dimension(  MAXEQNH    ) ::   valH
!!!      VTYPE,dimension(  MAXEQNH,3  ) ::  dvalH
!!!      VTYPE,dimension(  MAXEQNH,3,3) :: d2valH
!!!      VTYPE,dimension(3,MAXEQNE    ) ::   valE
!!!      VTYPE,dimension(3,MAXEQNE,3  ) ::  dvalE
!!!      VTYPE,dimension(3,MAXEQNE,3,3) :: d2valE
!!!      VTYPE,dimension(3,MAXEQNV    ) ::   valV
!!!      VTYPE,dimension(3,MAXEQNV,3  ) ::  dvalV
!!!      VTYPE,dimension(3,MAXEQNV,3,3) :: d2valV
!!!      VTYPE,dimension(  MAXEQNQ    ) ::   valQ
!!!      VTYPE,dimension(  MAXEQNQ,3  ) ::  dvalQ
!!!      VTYPE,dimension(  MAXEQNQ,3,3) :: d2valQ
!!!      integer :: ivoid
!
      real(8) :: aux(3),aux_n(3)
      real(8) :: s
      integer :: iattr,icomp,ireal,ibeg,iload,isol
!
      real(8), external :: dreal_part,dimag_part
!
!-------------------------------------------------------------------------------------
!
!
!!!!     compute exact solution
!!!      ivoid=0
!!!      call exact(X,ivoid, valH,dvalH,d2valH,valE,dvalE,d2valE, &
!!!                          valV,dvalV,d2valV,valQ,dvalQ,d2valQ )
!
!     decode
      ireal = ISELECT/1000 ; iattr = ISELECT - ireal*1000
      iattr =   iattr/ 100 ; icomp = ISELECT - ireal*1000 - iattr*100
      icomp =   icomp/  10 ; iload = ISELECT - ireal*1000 - iattr*100 - icomp*10
!
!     address of 1st component for the attribute
      ibeg=ADRES(iattr)
!
!     discretization type
      select case(DTYPE(iattr))
!
!     -- H1 --
      case('contin')
!
        isol = (iload-1)*NRHVAR + ibeg + icomp
!
        if (ireal == 1) then ; Val = dreal_part(SolH(isol))
        else                 ; Val = dimag_part(SolH(isol))
        endif
!
!     -- H(curl) --
      case('tangen')
!
        isol = (iload-1)*NREVAR + ibeg + icomp
!
        if (ireal == 1) then ; aux(1) = dreal_part(SolE(1,isol))
                               aux(2) = dreal_part(SolE(2,isol))
                               aux(3) = dreal_part(SolE(3,isol))
        else                 ; aux(1) = dimag_part(SolE(1,isol))
                               aux(2) = dimag_part(SolE(2,isol))
                               aux(3) = dimag_part(SolE(3,isol))
        endif
!
!       normal component
        call scalar_product(aux,Rn, s)
        aux_n = s*Rn
!
!       tangential component
        aux = aux - aux_n
        call norm(aux, Val)
!
!     -- H(div) --
      case('normal')
!
        isol = (iload-1)*NRVVAR + ibeg + icomp
!
        if (ireal == 1) then ; aux(1) = dreal_part(SolV(1,isol))
                               aux(2) = dreal_part(SolV(2,isol))
                               aux(3) = dreal_part(SolV(3,isol))
        else                 ; aux(1) = dimag_part(SolV(1,isol))
                               aux(2) = dimag_part(SolV(2,isol))
                               aux(3) = dimag_part(SolV(3,isol))
        endif

        call scalar_product(aux,Rn, s)
        Val = abs(s)
!
!     -- L2 --
      case('discon')
!
        isol = (iload-1)*NRQVAR + ibeg + icomp
!
        if (ireal == 1) then ; Val = dreal_part(SolQ(isol))
        else                 ; Val = dimag_part(SolQ(isol))
        endif
      endselect
!
!
end subroutine soldis_system
