!--------------------------------------------------------------------------------
subroutine proj_quad2edge(Eta,Ie, Zeta,dZeta_dEta)
!--------------------------------------------------------------------------------
! REMARK: dirty routine!! Needs to be cleaned up!!
!--------------------------------------------------------------------------------
  use kinds
!!!  use element_data
!--------------------------------------------------------------------------------
! DUMMY ARGUMENTS
  real*8, dimension(2), intent(in)  :: Eta
  integer,                intent(in)  :: Ie
  real*8,               intent(out) :: Zeta
  real*8, dimension(2), intent(out) :: dZeta_dEta
!--------------------------------------------------------------------------------
!
  select case(Ie)
    case(1)
      Zeta = Eta(1)
      dZeta_dEta = (/1.d0, 0.d0/)
    case(2)
      Zeta = Eta(2)
      dZeta_dEta = (/0.d0, 1.d0/)
    case(3)
      Zeta = Eta(1)
      dZeta_dEta = (/1.d0, 0.d0/)
    case(4)
      Zeta = Eta(2)
      dZeta_dEta = (/0.d0, 1.d0/)
    case default
      write(*,*)'proj_quad2edge: number of edges exceeded, Ie = ',Ie
      stop
  end select
!
end subroutine proj_quad2edge
!--------------------------------------------------------------------------------
