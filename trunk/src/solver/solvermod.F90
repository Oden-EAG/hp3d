!--------------------------------------------------------
!> Purpose : contains the environmental variables for UHM
module solvermod
  integer, parameter :: INCLUDE_ALL_DOMAINS = -1
  integer, save :: NR_COMPUTE_DOMAINS = -1
  integer, save, dimension(20) :: COMPUTE_DOMAINS

  double precision :: ANALYSIS_TIME, ASSEMBLY_TIME, FACTOR_TIME, SOLVE_TIME
  integer :: FACTOR_MEM

contains
  !--------------------------------------------------------
  subroutine set_compute_domains(Nrdomain)
    implicit none
    integer, intent(in) :: Nrdomain
    if (Nrdomain.eq.INCLUDE_ALL_DOMAINS) then
       NR_COMPUTE_DOMAINS = INCLUDE_ALL_DOMAINS
       COMPUTE_DOMAINS(1:20) = 0
    else
       NR_COMPUTE_DOMAINS = Nrdomain
    end if
  end subroutine set_compute_domains

  logical function is_compute_domain(Ndom)
    implicit none
    integer, intent(in) :: Ndom
    integer :: loc
    if (NR_COMPUTE_DOMAINS.eq.INCLUDE_ALL_DOMAINS) then
       is_compute_domain = .TRUE.
    else
       call locate(Ndom, COMPUTE_DOMAINS, NR_COMPUTE_DOMAINS, loc)
       if (loc.gt.0) then
          is_compute_domain = .TRUE.
       else
          is_compute_domain = .FALSE.
       end if
    end if
  end function is_compute_domain
  
  integer function boolean2integer(In)
    implicit none
    logical, intent(in) :: In
    if (In) then
       boolean2integer = 1
    else
       boolean2integer = 0
    end if
  end function boolean2integer

end module solvermod
