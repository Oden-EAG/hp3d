!--------------------------------------------------------------------------------------------
subroutine proj_prism2face(Eta,Ifig, Xi,dXi_dEta)
!--------------------------------------------------------------------------------------------
! REMARKS:
!   1. Eta's  - block LOCAL coordinates
!      Xi's   - face  LOCAL coordinates
!      Zeta   - edge  LOCAL coordiantes
!   2. routine heavily relies on the way orientations are set for the prism!!!
!--------------------------------------------------------------------------------------------
  use kinds
  use element_data
!--------------------------------------------------------------------------------------------
! DUMMY ARGUMENTS
  real*8, dimension(3),   intent(in)  :: Eta
  integer,                  intent(in)  :: Ifig
  real*8, dimension(2),   intent(out) :: Xi
  real*8, dimension(2,3), intent(out) :: dXi_dEta
!--------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  real*8, dimension(3)   :: vshapt
  real*8, dimension(2,3) :: dvshapt
  integer                  :: iv1,iv2
!--------------------------------------------------------------------------------------------
!
! ..affine coordinates for the master triangle
    vshapt(1) = 1.d0 - Eta(1) - Eta(2)
    dvshapt(1,1) = -1.d0; dvshapt(2,1) = -1.d0
    vshapt(2) = Eta(1)
    dvshapt(1,2) = 1.d0; dvshapt(2,2) = 0.d0
    vshapt(3) = Eta(2)
    dvshapt(1,3) = 0.d0; dvshapt(2,3) = 1.d0
! ..project on vertical edge
    Xi(2) = Eta(3)
    dXi_dEta(2,1:3) = (/0.d0, 0.d0, 1.d0/)
! ..project on horizontal edge
    iv1 =  TRIAN_EDGE_TO_VERT(1,Ifig - 2)
    iv2 =  TRIAN_EDGE_TO_VERT(2,Ifig - 2)
    call proj_t2e(Eta(1:2),iv1,iv2,vshapt,dvshapt, Xi(1),dXi_dEta(1,1:2))
    dXi_dEta(1,3) = 0.d0
!
end subroutine proj_prism2face
!--------------------------------------------------------------------------------------------
