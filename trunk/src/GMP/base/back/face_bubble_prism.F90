!--------------------------------------------------------------------------------------------
subroutine face_bubble_prism(Eta,Ifig,Nr,Orie, Phi,dPhi_dEta)
!--------------------------------------------------------------------------------------------
 use element_data
!--------------------------------------------------------------------------------------------
! DUMMY ARGUMENTS
  real*8, dimension(3),   intent(in)  :: Eta
  integer,                  intent(in)  :: Ifig
  integer,                  intent(in)  :: Nr
  integer,                  intent(in)  :: Orie
  real*8, dimension(3,1), intent(out) :: Phi
  real*8, dimension(3,3), intent(out) :: dPhi_dEta
!--------------------------------------------------------------------------------------------  
! LOCAL VARIABLES
  real*8                              :: blend
  real*8, dimension(1,3)              :: dblend
  real*8, dimension(3)                :: vshapt
  real*8, dimension(2,3)              :: dvshapt
  integer                               :: iv1,iv2
  real*8, dimension(2)                :: xi
  real*8, dimension(2,3)              :: dxi_dEta
  real*8, dimension(3,1)              :: psi
  real*8, dimension(3,2)              :: dpsi_dxi
  real*8, dimension(3,3)              :: dpsi_dEta
!--------------------------------------------------------------------------------------------
!
! ..affine coordinates for the master triangle
    vshapt(1) = 1.d0 - Eta(1) - Eta(2)
    dvshapt(1,1) = -1.d0;  dvshapt(2,1) = -1.d0
    vshapt(2) = Eta(1)
    dvshapt(1,2) = 1.d0;   dvshapt(2,2) = 0.d0
    vshapt(3) = Eta(2)
    dvshapt(1,3) = 0.d0;   dvshapt(2,3) = 1.d0
! ..get the edge vertices specifying the local edge orientation
    iv1 =  TRIAN_EDGE_TO_VERT(1,Ifig - 2)
    iv2 =  TRIAN_EDGE_TO_VERT(2,Ifig - 2)
! ..blending function
    blend = vshapt(iv1)*vshapt(iv2)
    dblend(1,1:2) = dvshapt(1:2,iv1)*vshapt(iv1) + vshapt(iv1)*dvshapt(1:2,iv2)
    dblend(1,3)   = 0.d0
! ..project onto face
    call proj_prism2face(Eta,Ifig, xi,dxi_dEta)
! ..multiply kernel function by blending function
    call face_kernel_prism(xi,Nr,Orie, psi,dpsi_dxi)
    dpsi_dEta = matmul(dpsi_dxi, dxi_dEta)
    Phi = psi*blend
    dPhi_dEta = dpsi_dEta*blend + matmul(psi, dblend)
!
end subroutine face_bubble_prism
!--------------------------------------------------------------------------------------------
