!--------------------------------------------------------------------------------------------
subroutine face_kernel_prism(Xi,Nr,Orie, Psi,dPsi_dXi)
!--------------------------------------------------------------------------------------------
 use kinds
!--------------------------------------------------------------------------------------------
! DUMMY ARGUMENTS
  real*8, dimension(3),   intent(in)  :: Xi
  integer,                  intent(in)  :: Nr
  integer,                  intent(in)  :: Orie
  real*8, dimension(3),   intent(out) :: Psi
  real*8, dimension(3,2), intent(out) :: dPsi_dXi
!--------------------------------------------------------------------------------------------
!
    call quadB(Nr,Xi,Orie, Psi,dPsi_dXi)
    dPsi_dXi(1:3,1) = dPsi_dXi(1:3,1)*(1.d0 - Xi(1))*Xi(1) - Psi(1:3)*(1.d0 - 2.d0*Xi(1))/ &
                      ((1.d0 - Xi(1))*Xi(i))**2
    dPsi_dXi(1:3,2) = dPsi_dXi(1:3,2)/((1.d0 - Xi(1))*Xi(1))
    Psi = Psi/((1.d0 - Xi(1))*Xi(2))
!
end subroutine face_kernel_prism
!--------------------------------------------------------------------------------------------

