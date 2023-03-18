!---------------------------------------------------------------------------------------
!> Purpose : computes norm of the error and the exact solution
!!
!! @param[out]  Err   - norm of the error
!! @param[out]  Rnorm - norm of the exact solution
!---------------------------------------------------------------------------------------
subroutine exact_error
  use control          , only : NEXACT
  use data_structure3D
  use environment      , only : QUIET_MODE,L2PROJ,FILE_ERR
  use physics
  use common_prob_data
  implicit none

  integer, dimension(NR_PHYSA) :: flag
!
  real*8 :: errorH,errorE,errorV,errorQ,errorHEVQ,derrorH,derrorE,derrorV,derrorQ
  real*8 :: rnormH,rnormE,rnormV,rnormQ,rnormHEVQ,drnormH,drnormE,drnormV,drnormQ
  real*8 :: errorH_rel,errorE_rel,errorV_rel,errorQ_rel,errorHEVQ_rel
  real*8 :: rateH,rateE,rateQ,rateHEVQ
!
  integer, parameter :: nin = 13
  integer, parameter :: maxvis =2000
!
! for computing error rate
  integer :: ndofH,ndofE,ndofV,ndofQ
  integer, save :: ivis = 0
  integer, save :: nrdof_tot_save
  real*8 , save :: errorH_save,errorE_save,errorV_save,errorQ_save,errorHEVQ_save
  real*8 , dimension(maxvis,10), save :: rwork
  integer, dimension(maxvis,10), save :: iwork
!
! miscellaneous
  integer :: mdle,i,iattr,nrdof_tot,ic
!
! printing flag
  integer :: iprint
!
!---------------------------------------------------------------------------------------
!    I N I T I A L I Z E
!---------------------------------------------------------------------------------------
!
  iprint=0
!
  select case(IERROR_PROB)
  case(IERROR_L2)
    L2PROJ = .TRUE.
  case(IERROR_NATURAL)
    L2PROJ = .FALSE.
  case default
    write(*,*) 'exact_error : Error calculation type not supported'
  end select
!
! check that exact solution is indeed known
  if (NEXACT == 0) then
    write(*,*) 'exact_error: UNKNOWN exact solution!'
    return
  endif
!
! initialize global quantities
  errorH=0.d0 ; rnormH=0.d0
  errorE=0.d0 ; rnormE=0.d0
  errorV=0.d0 ; rnormV=0.d0
  errorQ=0.d0 ; rnormQ=0.d0
!
  write(*,7000)
  7000 format('Declare the attribute to calculate the error of: ',  &
              '   1)Displacement, 2)Stress, 3)Combined, 4)Lagrange Multiplier')
  read(*,*) IERROR_ATTR
!
!---------------------------------------------------------------------------------------
!    A C C U M U L A T E    E L E M E N T    E R R O R S
!---------------------------------------------------------------------------------------
!
      nrdof_tot=0
!
!     loop over active elements
      mdle=0
      do i=1,NRELES
        call nelcon(mdle,mdle)
!
!       set flags
        select case(NODES(Mdle)%case)
        ! PRIMAL
        case(24)
          select case(IERROR_ATTR)
          case(DISPLACEMENT)
            flag = (/1,0,0,0,0/)
          case(STRESS)
            flag = (/0,0,0,0,0/)
          case(COMBINED)
            flag = (/1,0,0,0,0/)
          case(LAGRANGE)
            flag = (/0,0,0,0,0/)
          end select
        ! ULTRA-WEAK
        case(31)
          select case(IERROR_ATTR)
          case(DISPLACEMENT)
            flag = (/0,0,1,0,0/)
          case(STRESS)
            flag = (/0,0,0,1,0/)
          case(COMBINED)
            flag = (/0,0,1,1,0/)
          case(LAGRANGE)
            flag = (/0,0,0,0,1/)
          end select
        end select
!
!       compute error on element
        call element_error(mdle, flag, derrorH,derrorE,derrorV,derrorQ, &
                                       drnormH,drnormE,drnormV,drnormQ )
!
!       accumulate
        errorH = errorH + derrorH; rnormH = rnormH + drnormH
        errorE = errorE + derrorE; rnormE = rnormE + drnormE
        errorV = errorV + derrorV; rnormV = rnormV + drnormV
        errorQ = errorQ + derrorQ; rnormQ = rnormQ + drnormQ
!
!       return dof for element
        call find_ndof(mdle, ndofH,ndofE,ndofV,ndofQ)
!
!       loop over physical attributes
        do iattr=1,NR_PHYSA
!
!         skip if the error not calculated
          if (flag(iattr).eq.0) cycle
          select case(D_TYPE(iattr))
          !  This is not perfect because it does not remove the fixed DOF from the boundary conditions
          case(CONTIN); nrdof_tot = nrdof_tot + ndofH*NR_COMP(iattr)
          case(TANGEN); nrdof_tot = nrdof_tot + ndofE*NR_COMP(iattr)
          case(NORMAL); nrdof_tot = nrdof_tot + ndofV*NR_COMP(iattr)
          case(DISCON); nrdof_tot = nrdof_tot + ndofQ*NR_COMP(iattr)
          end select
        enddo
      enddo
!
!     compute total error
      errorHEVQ=sqrt(errorH+errorE+errorV+errorQ)
      errorH   =sqrt(errorH)
      errorE   =sqrt(errorE)
      errorV   =sqrt(errorV)
      errorQ   =sqrt(errorQ)
!
!     compute norm
      rnormHEVQ=sqrt(rnormH+rnormE+rnormV+rnormQ)
      rnormH   =sqrt(rnormH)
      rnormE   =sqrt(rnormE)
      rnormV   =sqrt(rnormV)
      rnormQ   =sqrt(rnormQ)
!
!     compute relative error
      errorHEVQ_rel=0.d0 ; errorH_rel=0.d0 ; errorE_rel=0.d0 ; errorV_rel=0.d0 ; errorQ_rel=0.d0
      if (rnormHEVQ > 0.d0) errorHEVQ_rel=errorHEVQ/rnormHEVQ
      if (rnormH    > 0.d0) errorH_rel   =errorH   /rnormH
      if (rnormE    > 0.d0) errorE_rel   =errorE   /rnormE
      if (rnormV    > 0.d0) errorV_rel   =errorV   /rnormV
      if (rnormQ    > 0.d0) errorQ_rel   =errorQ   /rnormQ
!
!     compute rate
      rateHEVQ=0.d0 ; rateH=0.d0 ; rateE=0.d0 ; rateQ=0.d0
      if (ivis /= 0) then
        if (nrdof_tot > nrdof_tot_save) then
          rateHEVQ = (log(errorHEVQ_save/errorHEVQ))/log(float(nrdof_tot_save)/float(nrdof_tot))
      endif ; endif
!
!     save quantities
      errorHEVQ_save=errorHEVQ
      errorH_save=errorH ; errorE_save=errorE ; errorV_save=errorV ; errorQ_save=errorQ
      nrdof_tot_save=nrdof_tot
!
!     raise visitation flag
      ivis=ivis+1
!
IF (.NOT. QUIET_MODE) THEN
!
!     check
      if (ivis > maxvis) then
        write(*,*) 'exact_error: increase maxvis!'
        stop
      endif
!
!     store
      rwork(ivis, 1)=errorH
      rwork(ivis, 2)=errorE
      rwork(ivis, 3)=errorV
      rwork(ivis, 4)=errorQ
      rwork(ivis, 5)=errorHEVQ
      rwork(ivis, 6)=rateHEVQ
      iwork(ivis, 1)=IERROR_ATTR
!
ENDIF
!
!     printing
!
!     -- 1st visit --
      if (ivis == 1) then
!
!       open file
        open(unit=nin,file=trim(FILE_ERR),form='formatted',access='sequential',status='unknown',iostat=ic)
        if (ic /= 0) then
          write(*,*)'exact_error: COULD NOT OPEN FILE! [0]'
          stop
        endif
!
!       print header
IF (.NOT. L2PROJ) THEN
        write(nin,*)'-- Error Report --'
ELSE
        write(nin,*)'-- Error Report (L2 only)--'
ENDIF
        write(nin,9998)
 9998   format('             H1            //', &
                           ' H(curl)       //', &
                           ' H(div)        //', &
                           ' L2            //', &
                           ' Total         //', &
                           ' Rate       //',    &
                           ' Case tag')
!
!     -- subsequent visits --
      else
!
!       append to file
        open(unit=nin,file=trim(FILE_ERR),form='formatted',access='sequential',status='old',position='append',iostat=ic)
        if (ic /= 0) then
          write(*,*)'exact_error: COULD NOT OPEN FILE! [1]'
          stop
        endif
      endif
!
!     print to file
      write(nin,9999)ivis,errorH,errorE,errorV,errorQ,errorHEVQ,rateHEVQ,IERROR_ATTR
 9999 format(1x,i6,' ; ',2x,5(e12.5,' ; ',2x),f9.6,' ; ',3x,i8)
!
!     print to screen
IF (.NOT.QUIET_MODE) THEN ; write(*,*)''
IF (.NOT. L2PROJ   ) THEN ; write(*,*)'-- Error Report --'
ELSE                      ; write(*,*)'-- Error Report (L2 only)--'
ENDIF
                            write(*,9998)
      do i=1,ivis         ; write(*,9999)i,rwork(i,1:6),iwork(i,1)
      enddo
                            write(*,*)''
ENDIF
!
!     close file
      close(unit=nin,iostat=ic)
      if (ic /= 0) then
        write(*,*)'exact_error: COULD NOT CLOSE FILE!'
        stop
      endif
!
!
end subroutine exact_error
