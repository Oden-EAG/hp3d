!----------------------------------------------------------------------------
!> @brief Defines variables for initialization
!> @date Mar 2023
!----------------------------------------------------------------------------
!
module environment
!
      implicit none
!
      save
!
!     command line inputs
      integer, parameter :: MAX_ENV_PARAMETER_ = 128
      character(len=128), dimension(MAX_ENV_PARAMETER_) :: ARGS_
      integer :: ARGC_
!
!     refinement lookup location
      character(len=128) :: FILE_REFINE = '../../files/ref'
!
!     basic environments
      character(len=128) :: FILE_CONTROL,  &
                            FILE_HISTORY,  &
                            FILE_GEOM,     &
                            FILE_PHYS,     &
                            FILE_ERR,      &
                            PREFIX
!
!     operation mode
      logical :: IVERBOSE_  = .FALSE.
      logical :: IDRY_      = .FALSE.
      logical :: QUIET_MODE = .FALSE.
      logical :: L2PROJ     = .FALSE.
      logical :: L2GEOM     = .FALSE.
!
!
   contains
!
!
!----------------------------------------------------------------------------
!> @brief Determines value of a BOOLEAN option variable (.TRUE. if present)
!> @note  If option is not present, value is set to a user-provided default
!!
!! @param[in ] Aopt  - name of BOOLEAN option
!! @param[in ] Atext - text explanation of option
!! @param[in ] Ndef  - default value of BOOLEAN option
!! @param[out] Nval  - actual value of BOOLEAN option
!!
!> @date Mar 2023
!----------------------------------------------------------------------------
  subroutine get_option_bool(Aopt,Atext,Ndef, Nval)
 !
    implicit none
    character(len=*), intent(in ) :: Aopt, Atext
    logical,          intent(in ) :: Ndef
    logical,          intent(out) :: Nval
!
    character(len=20 ) :: str1
    character(len=25 ) :: str2
    character(len=100) :: str3
    integer :: i, ifound
!----------------------------------------------------------------------------
 !
    if (IDRY_) then
       write(str1,'(a20)' ) Aopt
       write(str2,'(l1)'  ) Ndef
       write(str3,'(a100)') Atext
       write(*,1) adjustl(str1), adjustl(str2), adjustl(str3)
 1     format(1x,a20,' Boolean ; default = ',a25,' ; ',a100)
    else

       ! set to default value
       Nval = Ndef

       ! scan input arguments and update value
       ifound=0
       do i=1,ARGC_
          if (LLT(trim(ARGS_(i)),trim(Aopt))) then
             cycle
          else if (LGT(ARGS_(i),Aopt)) then
             cycle
          else
             ! if found, set to .TRUE.
             Nval = .TRUE. ; ifound = i
             exit
          endif
       enddo
       !
       if (IVERBOSE_) then
          if (ifound.gt.0) then
             write(*,*) &
                  Aopt, ' is found at ', &
                  ifound, ' and it is set ', Nval, ' ,',Atext
          end if
       end if
    end if

  end subroutine get_option_bool
!
!
!----------------------------------------------------------------------------
!> @brief Determines value of an INTEGER option variable
!> @note  If option is not present, value is set to a user-provided default
!!
!! @param[in ] Aopt  - name of INTEGER option
!! @param[in ] Atext - text explanation of option
!! @param[in ] Ndef  - default value of INTEGER option
!! @param[out] Nval  - actual value of INTEGER option
!!
!> @date Mar 2023
!----------------------------------------------------------------------------
  subroutine get_option_int(Aopt, Atext, Ndef, Nval)
!
    implicit none
    character(len=*),intent(in ) :: Aopt, Atext
    integer,         intent(in ) :: Ndef
    integer,         intent(out) :: Nval
!
    character(len=20 ) :: str1
    character(len=25 ) :: str2
    character(len=100) :: str3
    integer :: i, ifound
!----------------------------------------------------------------------------
!
    if (IDRY_) then
       write(str1,'(a20)' ) Aopt
       write(str2,'(i8)'  ) Ndef
       write(str3,'(a100)') Atext
       write(*,1) adjustl(str1), adjustl(str2), adjustl(str3)
 1     format(1x,a20,' Integer ; default = ',a25,' ; ',a100)
    else

       ! set to default value
       Nval = Ndef

       ! scan input arguments to determine value
       ifound=0
       do i=1,ARGC_
          if (LLT(trim(ARGS_(i)),trim(Aopt))) then
             cycle
          else if (LGT(ARGS_(i),Aopt)) then
             cycle
          else
             ! value is contained in the following argument
             read(ARGS_(i+1),'(i10)') Nval ; ifound = i
             exit
          endif
       enddo
       !
       if (IVERBOSE_) then
          if (ifound.gt.0) then
             write(*,*) &
                  Aopt, ' is found at ', &
                  ifound, ' and it is set ', Nval, ' ,',Atext
          end if
       end if
    end if

  end subroutine get_option_int
!
!
!----------------------------------------------------------------------------
!> @brief Determines value of a REAL option variable
!> @note  If option is not present, value is set to a user-provided default
!!
!! @param[in ] Aopt  - name of REAL option
!! @param[in ] Atext - text explanation of option
!! @param[in ] Def   - default value of REAL option
!! @param[out] Val   - actual value of REAL option
!!
!> @date Mar 2023
!----------------------------------------------------------------------------
  subroutine get_option_real(Aopt,Atext,Def, Val)
!
    implicit none
    character(len=*),intent(in ) :: Aopt, Atext
    real(8),         intent(in ) :: Def
    real(8),         intent(out) :: Val
!
    character(len=20 ) :: str1
    character(len=25 ) :: str2
    character(len=100) :: str3
    integer :: i, ifound
!----------------------------------------------------------------------------

    if (IDRY_) then
       write(str1,'(a20)'  ) Aopt
       write(str2,'(e12.5)') Def
       write(str3,'(a100)' ) Atext
       write(*,1) adjustl(str1), adjustl(str2), adjustl(str3)
 1     format(1x,a20,' Real    ; default = ',a25,' ; ',a100)
    else

       ! set to default value
       Val = Def

       ! scan input arguments to determine value
       ifound = 0
       do i=1,ARGC_
          if (LLT(trim(ARGS_(i)),trim(Aopt))) then
             cycle
          else if (LGT(ARGS_(i),Aopt)) then
             cycle
          else
             ! value is contained in the following argument
             read(ARGS_(i+1),'(f20.15)') Val ; ifound = i
             exit
          endif
       enddo
       !
       if (IVERBOSE_) then
          if (ifound.gt.0) then
             write(*,*) &
                  Aopt, ' is found at ', &
                  ifound, ' and it is set ', Val, ' ,',Atext
          end if
       end if
    end if

  end subroutine get_option_real
!
!
!----------------------------------------------------------------------------
!> @brief Determines value of a COMPLEX option variable
!> @note  If option is not present, value is set to a user-provided default
!!
!! @param[in ] Aopt  - name of COMPLEX option
!! @param[in ] Atext - text explanation of option
!! @param[in ] Def   - default value of COMPLEX option
!! @param[out] Val   - actual value of COMPLEX option
!!
!> @date Mar 2023
!----------------------------------------------------------------------------
  subroutine get_option_comp(Aopt,Atext,Def, Val)

    implicit none
    character(len=*),intent(in ) :: Aopt, Atext
    complex(8),      intent(in ) :: Def
    complex(8),      intent(out) :: Val

    integer :: i, ifound
!
    if (IDRY_) then
       write(*,*) Aopt, ' Real has default value = ', Def, ' ,', Atext
    else

       ! set to default value
       Val = Def

       ! scan input arguments to determine value
       ifound = 0
       do i=1,ARGC_
          if (LLT(trim(ARGS_(i)),trim(Aopt))) then
             cycle
          else if (LGT(ARGS_(i),Aopt)) then
             cycle
          else
             ! value is contained in the following argument
             read(ARGS_(i+1),'(f20.15,f20.15)') Val ; ifound = i
             exit
          endif
       enddo
       !
       if (IVERBOSE_) then
          if (ifound.gt.0) then
             write(*,*) &
                  Aopt, ' is found at ', &
                  ifound, ' and it is set ', Val, ' ,',Atext
          end if
       end if
    end if

  end subroutine get_option_comp
!
!
!----------------------------------------------------------------------------
!> @brief Determines value of a STRING option variable
!> @note  If option is not present, value is set to a user-provided default
!!
!! @param[in ] Aopt  - name of STRING option
!! @param[in ] Atext - text explanation of option
!! @param[in ] Adef  - defaul value of STRING option
!! @param[out] Aval  - actual value of STRING option
!!
!> @date Mar 2023
!----------------------------------------------------------------------------
  subroutine get_option_string(Aopt, Atext, Adef, Aval)
!
    implicit none
    character(len=*),intent(in ) :: Aopt, Atext, Adef
    character(len=*),intent(out) :: Aval
!
    character(len=20 ) :: str1
    character(len=25 ) :: str2
    character(len=100) :: str3
    integer :: i, ifound
!
    if (IDRY_) then
       write(str1,'(a20)' ) Aopt
       write(str2,'(a25)' ) Adef
       write(str3,'(a100)') Atext
       write(*,1) adjustl(str1), adjustl(str2), adjustl(str3)
 1     format(1x,a20,' String  ; default = ',a25,' ; ',a100)
    else

       ! set to default value
       Aval = Adef

       ! scan input arguments to determine value
       ifound=0
       do i=1,ARGC_
          if (LLT(trim(ARGS_(i)),trim(Aopt))) then
             cycle
          else if (LGT(ARGS_(i),Aopt)) then
             cycle
          else
             ! value is contained in the following argument
             Aval = ARGS_(i+1) ; ifound = i
             exit
          endif
       enddo
       !
       if (IVERBOSE_) then
          if (ifound.gt.0) then
             write(*,*) &
                  Aopt, ' is found at ', &
                  ifound, ' and it is set ', Aval, ' ,',Atext
          end if
       end if
    end if

  end subroutine get_option_string
!
!
!----------------------------------------------------------------------------
!> @brief Reads in options arguments
!> @date Mar 2023
!----------------------------------------------------------------------------
  subroutine begin_environment
!
    implicit none
    integer :: i
!
    ! number of arguments
    ARGC_ = command_argument_count()

    ! read in arguments
    do i=1,ARGC_
       !call getarg(i, ARGS_(i))
       call get_command_argument(i, ARGS_(i))
    end do

    ! set up options : -help , -verbose
    call get_option_bool('-help'   ,'Dry run'       ,.FALSE., IDRY_)
    call get_option_bool('-verbose','Verbose output',.FALSE., IVERBOSE_)

  end subroutine begin_environment
!
!
!----------------------------------------------------------------------------
!> @date Mar 2023
!----------------------------------------------------------------------------
  subroutine end_environment
    implicit none
    if (IDRY_) then
       stop 1
    end if
  end subroutine end_environment
!
!
!----------------------------------------------------------------------------
!> @date Mar 2023
!----------------------------------------------------------------------------
  subroutine disp_env(Iout)
    implicit none
    integer, intent(in) :: Iout
    integer :: i
    !
    write(Iout,*) '-- Input Arguments List --'
    do i=1,ARGC_
       write(Iout,*) 'Idx = ', i, '  Content = ', ARGS_(i)
    end do
    !
  end subroutine disp_env
!
!
end module environment
