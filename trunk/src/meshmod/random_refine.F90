!> Purpose : randomly refine mesh
!! @param[in] Per   - percentage to refine
!! @param[in] Niter - number of iteration
!!
!! REAMARK : routine does NOT update gdofs, since random refinements
!! are only meant for testing

subroutine random_refine(Per, Nitr)
  use data_structure3D
  use refinements
  use environment , only : QUIET_MODE
  implicit none
  !
  ! ** Arguments
  !-----------------------------------------------------
  real(8), intent(in) :: Per
  integer, intent(in) :: Nitr
  !
  ! ** Locals
  !-----------------------------------------------------
  real(8) :: x
  integer :: iprint, iseed, iel, i, istat
  integer :: nsize_list, kref, mdle, nref, idx, nelts
  integer, allocatable ::  mdle_list(:)
  logical :: mode_save
  !-----------------------------------------------------
  iprint=0

  !  ...loop over iteration
  do i=1,Nitr

     !  ...allocate a local list
     nsize_list=NRELES
     allocate(mdle_list(nsize_list), stat=istat)
     if (istat.ne.SUCCESS) then
        call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
     endif

     !  ...dump elements into list
     mdle=0
     do iel=1, nsize_list
        call nelcon(mdle, mdle)
        mdle_list(iel) = mdle
     enddo

!    number of refinements to perform (refine at least one element)
     nref=max(int(Per*nsize_list),1)

IF (.NOT. QUIET_MODE) THEN
     write(*,7010) i, nsize_list, nref
7010 format(1x,i3,' ; elements = ',i7,' ; ref. to perform = ',i7)
ENDIF

     !  ...for each iteration, randomly refine elements
     do while (nref.gt.0)

        iseed = nref
        call random_seed(size=iseed)
        call random_number(harvest=x)

        nref  = nref - 1
        iseed = abs(int(x*1000000000))
        idx   = mod(iseed, nsize_list) + 1
        mdle  = mdle_list(idx)
!
!       leaf element
        if (is_leaf(mdle)) then
!
!          pick refinement kind
           select case (NODES(mdle)%type)
           case ('mdln','mdld','mdlb'); call get_isoref(mdle, kref)
              ! Testing all cases
              ! case ('mdlp'); kref = kref_kind(mod(idx,2)+2, 'mdlp')

              ! Testing anisotropic cases
           case ('mdlp'); kref = kref_kind(mod(idx,3)+1, 'mdlp')
              ! kyungjoo recover this after I complete aniso hcurl constrained approx
              kref = 11
           case default
              write(*,9999)NODES(mdle)%type
9999          format(' random_refine: Element type not supported! Type = ',a4)
              call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
           end select
           call refine(mdle, kref)
           nref = nref - 1
           if (iprint.eq.1) then
              write(*,7000) mdle, NODES(mdle)%type, kref
7000          format('mdle =',i6,' ', a5, '  kref=', i3)
           endif
        endif

     enddo

     deallocate(mdle_list, stat=istat)
     if (istat.ne.SUCCESS) then
        call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
     endif

     !  ...close mesh (do not update gdof's. Those refinements are just for
     !     testing)
     call close

     !  ...testing "neig_face" routine
     call verify_neig

     !  ...if last iteration, print
     if (i.eq.Nitr) then
IF (.NOT. QUIET_MODE) THEN
        write(*,7020) NRELES
7020    format(' NRELES = ',i7)
        write(*,*)'WARNING - Routine does NOT update gdofs!'
ENDIF
     endif
!
!    if number of active elements is large, skip visualization
     if (NRELES > 50000)  cycle

!    dump to Paraview for visualization
     mode_save = QUIET_MODE
     call paraview_driver
     QUIET_MODE = mode_save
!
  enddo
!
!
end subroutine random_refine
