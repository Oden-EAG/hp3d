!
!> Purpose   : Correct value of Dp in Poisson, 
!              for local test space to be greater than local trial space
!! @param[in]  Mdle           - element identifier
!! @param[in]  Nord_add_ini   - Default value of Dp (usually NORD_ADD)
!! @param[in]  Nord           - Nominal approximation order in current element
!! @param[in]  Nrf            - number of faces in element Mdle
!! @param[out] Nord_add_final - corrected Dp, at most MAXPP (set in module parametersDPG)
!
subroutine correct_dp(Mdle,Nord_add_ini,Nord,Nrf,Nord_add_final)
use parametersDPG, only : MAXPP
!
implicit none
!
integer, intent(in)  :: Mdle,Nord_add_ini,Nord,Nrf
integer, intent(out) :: Nord_add_final
! 
integer              :: ntrial,ntrial_0,ntrial_hat,ntest,nh,nv,nq,nu,nf, nord_test,k
!
! Compute number of degrees of freedom for space components:
! Element L^2
nq= Nord*(Nord+1)*(Nord+2)/6
! Face H^1
nu= (Nord+1)*(Nord+2)/2
! Face L^2
nf=   Nord  *(Nord+1)/2
! Sum trial functions accounting for number of variables and faces
ntrial_0= 12*nq
ntrial_hat = 3*Nrf*(nu+nf) 
ntrial = ntrial_0 + ntrial_hat

! initialize nord_add increase k, number of test functions and enriched order
k=0
ntest = -1
nord_test=0
!
do while (ntrial.ge.ntest.and.nord_test.lt.MAXPP)
! do while ((ntrial_0.ge.ntest.or.ntrial_hat.ge.ntest).and.nord_test.lt.MAXPP)
  ! update Dp final
  Nord_add_final = Nord_add_ini + k
  ! update enriched order
  nord_test=Nord+Nord_add_final
  ! compute number of test functions
  ! H^1
  nh= (nord_test+1)*(nord_test+2)*(nord_test+3)/6
  ! H(div)
  nv=   nord_test  *(nord_test+1)*(nord_test+3)/2
  ! sum test functions (just one each component)
  ntest = 3*(nh+nv)
  ! update increase
  k=k+1
enddo
!
! check if nord_test = MAXPP; then check whether ntest is greater than ntrial. If not, report it
if (nord_test.eq.MAXPP) then
  write(*,*) '-----------------------------------------------------------------'
  write(*,*) 'correct_dp :  MAXPP reached!'
  write(*,*) 'Element = ',mdle
  if (ntrial.ge.ntest) then
    write(*,*) 'Dim of test space LESS than dim of trial space !'
    write(*,*) 'ntrial, ntest = ',ntrial, ntest
  endif
  write(*,*) ' '
endif
!
! end, return Nord_add_final
end subroutine correct_dp