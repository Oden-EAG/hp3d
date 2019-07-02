!> Purpose : define initialization variables
!! rev@Feb.2012
!----------------------------------------------------------
module environment
  save
  integer, parameter :: MAX_ENV_PARAMETER_ = 128
  character(len=128), dimension(MAX_ENV_PARAMETER_) :: ARGS_
  integer :: ARGC_
  !
  logical :: IVERBOSE_ = .FALSE., IDRY_ = .FALSE.
  !
contains
  !----------------------------------------------------------
  subroutine get_option_bool(Aopt, Atext, Ndef, Nval)
    implicit none
    character(len=*), intent(in) :: Aopt, Atext
    logical, intent(in) :: Ndef
    logical, intent(out) :: Nval
    integer :: i, ifound

    if (IDRY_) then
       write(*,*) Aopt, ' Bool has default value = ', Ndef, ' ,', Atext
    else
       ! ** set default value
       Nval = Ndef
       !
       ! ** scan input arguments
       ifound = 0
       do i=1,ARGC_
          if (LLT(trim(ARGS_(i)),trim(Aopt))) then
             cycle
          else if (LGT(ARGS_(i),Aopt)) then
             cycle
          else
             Nval = .TRUE.
             ifound = i
             exit
          end if
       end do
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

  subroutine get_option_int(Aopt, Atext, Ndef, Nval)
    implicit none
    character(len=*), intent(in) :: Aopt, Atext
    integer, intent(in) :: Ndef
    integer, intent(out) :: Nval
    integer :: i, ifound

    if (IDRY_) then
       write(*,*) Aopt, ' Integer has default value = ', Ndef, ' ,', Atext
    else
       ! ** set default value
       Nval = Ndef
       !
       ! ** scan input arguments
       ifound = 0
       do i=1,ARGC_
          if (LLT(trim(ARGS_(i)),trim(Aopt))) then
             cycle
          else if (LGT(ARGS_(i),Aopt)) then
             cycle
          else
             read(ARGS_(i+1),'(i10)') Nval
             ifound = i
             exit
          end if
       end do
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

  subroutine get_option_real(Aopt, Atext, Def, Val)
    implicit none
    character(len=*), intent(in) :: Aopt, Atext
    real*8, intent(in) :: Def
    real*8, intent(out) :: Val
    integer :: i, ifound

    if (IDRY_) then
       write(*,*) Aopt, ' Real has default value = ', Def, ' ,', Atext
    else
       ! ** set default value
       Val = Def
       !
       ! ** scan input arguments
       ifound = 0
       do i=1,ARGC_
          if (LLT(trim(ARGS_(i)),trim(Aopt))) then
             cycle
          else if (LGT(ARGS_(i),Aopt)) then
             cycle
          else
             read(ARGS_(i+1),'(f20.15)') Val
             ifound = i
             exit
          end if
       end do
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

  subroutine get_option_string(Aopt, Atext, Adef, Aval)
    implicit none
    character(len=*), intent(in) :: Aopt, Atext, Adef
    character(len=*), intent(out) :: Aval
    integer :: i, ifound

    if (IDRY_) then
       write(*,*) Aopt, ' String has default value = ', Adef, ' ,', Atext
    else
       ! ** set default value
       Aval = Adef
       !
       ! ** scan input arguments
       ifound = 0
       do i=1,ARGC_
          if (LLT(trim(ARGS_(i)),trim(Aopt))) then
             cycle
          else if (LGT(ARGS_(i),Aopt)) then
             cycle
          else
             Aval = ARGS_(i+1)
             ifound = i
             exit
          end if
       end do
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

  subroutine begin_environment
    implicit none 
    integer :: i
    ARGC_ = iargc()
    do i=1,ARGC_
       call getarg(i, ARGS_(i))
    end do
    call get_option_bool( & 
         '-dry', 'Dry run if it is .TRUE.', &
         .FALSE.,IDRY_)
    call get_option_bool( & 
         '-verbose', 'Verbose output if it is .TRUE.', &
         .FALSE.,IVERBOSE_)
  end subroutine begin_environment

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
  end subroutine disp_env
  !
end module environment
