!----------------------------------------------------------------------------
!> Purpose : define initialization variables
!!
!! @date Dec 14
!----------------------------------------------------------------------------
!
module environment
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
contains
!
!
!
!----------------------------------------------------------------------------
!> Purpose : determines value of a BOOLEAN option variable (.TRUE. if
!!           present). If option is not present, valued is set to a
!!           user-provided default value.
!!
!! @param[in ] Aopt  - name of BOOLEAN option
!! @param[in ] Atext - text explanation of option
!! @param[in ] Ndef  - defaul value of BOOLEAN option
!! @param[out] Nval  - actual value of BOOLEAN option
!!
!> rev@Dec 2012
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

  endsubroutine get_option_bool



!----------------------------------------------------------------------------
!> Purpose : determines value of an INTEGER option variable. If option is
!!           not present, valued is set to a user-provided default value.
!!
!! @param[in ] Aopt  - name of INTEGER option
!! @param[in ] Atext - text explanation of option
!! @param[in ] Ndef  - defaul value of INTEGER option
!! @param[out] Nval  - actual value of INTEGER option
!!
!> rev@Dec 2012
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

  endsubroutine get_option_int



!----------------------------------------------------------------------------
!> Purpose : determines value of a REAL option variable. If option is
!!           not present, valued is set to a user-provided default value.
!!
!! @param[in ] Aopt  - name of REAL option
!! @param[in ] Atext - text explanation of option
!! @param[in ] Def   - defaul value of REAL option
!! @param[out] Val   - actual value of REAL option
!!
!> rev@Dec 2012
!----------------------------------------------------------------------------
  subroutine get_option_real(Aopt,Atext,Def, Val)
!
    implicit none
    character(len=*),intent(in ) :: Aopt, Atext
    real*8,          intent(in ) :: Def
    real*8,          intent(out) :: Val
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

  endsubroutine get_option_real



!----------------------------------------------------------------------------
!> Purpose : determines value of a COMPLEX option variable. If option is
!!           not present, valued is set to a user-provided default value.
!!
!! @param[in ] Aopt  - name of COMPLEX option
!! @param[in ] Atext - text explanation of option
!! @param[in ] Def   - defaul value of COMPLEX option
!! @param[out] Val   - actual value of COMPLEX option
!!
!> rev@Mar 2013
!----------------------------------------------------------------------------
  subroutine get_option_comp(Aopt,Atext,Def, Val)

    implicit none
    character(len=*),intent(in ) :: Aopt, Atext
    complex*16,      intent(in ) :: Def
    complex*16,      intent(out) :: Val

    integer :: i, ifound
!----------------------------------------------------------------------------

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

  endsubroutine get_option_comp



!----------------------------------------------------------------------------
!> Purpose : determines value of a STRING option variable. If option is not
!!           present, valued is set to a user-provided default value.
!!
!! @param[in ] Aopt  - name of STRING option
!! @param[in ] Atext - text explanation of option
!! @param[in ] Adef  - defaul value of STRING option
!! @param[out] Aval  - actual value of STRING option
!!
!> rev@Dec 2012
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
!----------------------------------------------------------------------------

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

  endsubroutine get_option_string



!----------------------------------------------------------------------------
!> Purpose : reads in options arguments
!----------------------------------------------------------------------------
  subroutine begin_environment

    implicit none
    integer :: i
!----------------------------------------------------------------------------

    ! number of arguments
    ARGC_ = iargc()

    ! read in arguments
    do i=1,ARGC_
       call getarg(i, ARGS_(i))
    end do

    ! set up options : -help , -verbose
    call get_option_bool('-help'   ,'Dry run'       ,.FALSE., IDRY_)
    call get_option_bool('-verbose','Verbose output',.FALSE., IVERBOSE_)

  endsubroutine begin_environment



  subroutine end_environment
    implicit none
    if (IDRY_) then
       call exit(0)
    end if
  end subroutine end_environment

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
  endsubroutine disp_env
  !
endmodule environment
