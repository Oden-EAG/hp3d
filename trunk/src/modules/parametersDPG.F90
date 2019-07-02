!>  Purpose :
!!    module defines various parameters used in the DPG method
!!
module parametersDPG
  use parameters
  implicit none
  !
  !  ...increment in order to define the enriched spaces
  integer, parameter :: MAXNORD_ADD=1
  integer :: NORD_ADD
  !
  !  ...maximum enriched order
  integer, parameter :: MAXPP = MAXP+MAXNORD_ADD
  !
  !  ...maximum number of quadrature points for an enriched  2D element
  integer, parameter :: MAXNINT2ADD = (MAXPP+1)**2
  !
  !  ...maximum number of quadrature points for an enriched  3D element
  integer, parameter :: MAXNINT3ADD = (MAXPP+1)**2*(MAXPP+2)
  !
  !----------------------------------------------------------------------
  !  === ELEMENT ====
  !  ...max number of local dof for a 2D quad element
  integer, parameter :: MAXquadHH = (MAXPP+1)**2
  integer, parameter :: MAXquadEE = 2*MAXPP*(MAXPP+1)
  integer, parameter :: MAXquadVV = MAXquadEE
  integer, parameter :: MAXquadQQ = MAXPP**2
  !
  !  ...max number of local dof for a 2D triangular element
  integer, parameter :: MAXtriaHH = (MAXPP+1)*(MAXPP+2)/2
  integer, parameter :: MAXtriaEE = MAXPP*(MAXPP+2)
  integer, parameter :: MAXtriaVV = MAXtriaEE
  integer, parameter :: MAXtriaQQ = MAXPP*(MAXPP+1)/2
  !
  !----------------------------------------------------------------------
  !
  !  ...max number of local dof for a 3D brick element
  integer, parameter :: MAXbrickHH = (MAXPP+1)**3
  integer, parameter :: MAXbrickEE = 3*MAXPP*(MAXPP+1)**2
  integer, parameter :: MAXbrickVV = 3*MAXPP**2*(MAXPP+1)
  integer, parameter :: MAXbrickQQ = MAXPP**3
  !
  !  ...max number of local dof for a 3D prism element
  integer, parameter :: MAXprismHH = MAXtriaHH*(MAXPP+1)
  integer, parameter :: MAXprismEE = MAXtriaEE*(MAXPP+1) + MAXtriaHH*MAXPP
  integer, parameter :: MAXprismVV = MAXtriaVV*MAXPP + MAXtriaQQ*(MAXPP+1)
  integer, parameter :: MAXprismQQ = MAXtriaQQ*MAXPP
  !
  !  ...max number of local dof for a 3D tetrahedral  element
  integer, parameter :: MAXtetraHH = (MAXPP+1)*(MAXPP+2)*(MAXPP+3)/6
  integer, parameter :: MAXtetraEE = MAXPP*(MAXPP+2)*(MAXPP+3)/2
  integer, parameter :: MAXtetraVV = MAXPP*(MAXPP+1)*(MAXPP+3)/2
  integer, parameter :: MAXtetraQQ = MAXPP*(MAXPP+1)*(MAXPP+2)/6
  !
  !  ...max number of local dof for a 3D pyramid  element
  integer, parameter :: MAXpyramHH = 5+8*(MAXPP-1)+(MAXPP-1)**2+(MAXPP-2)*(MAXPP-1)*2+(MAXPP-1)**3
  integer, parameter :: MAXpyramEE = 8*MAXPP+2*MAXPP*(MAXPP-1)+4*MAXPP*(MAXPP-1)+3*(MAXPP-1)**2*MAXPP
  integer, parameter :: MAXpyramVV = MAXquadQQ+4*MAXtriaQQ+3*(MAXPP-1)*MAXPP**2
  integer, parameter :: MAXpyramQQ = MAXPP**3
  !

  contains

    subroutine set_parametersDPG(Nord_add_read)
    integer :: Nord_add_read
    NORD_ADD = Nord_add_read
    end subroutine set_parametersDPG
end module parametersDPG
