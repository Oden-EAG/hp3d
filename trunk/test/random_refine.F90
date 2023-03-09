! !> @date Mar 2023
program test_random_refine
!
   real(8) :: per
   integer :: nitr
!
!..percentage of elements to refine per iteration
   per  = 0.1
!..number of refinements (iterations)
   nitr = 5
!
   call random_refine_test(per,nitr)
!
end program test_random_refine
!
!----------------------------------------------------------------------
! !> @brief Refines mesh randomly
! !> @date Mar 2023
!----------------------------------------------------------------------
   subroutine random_refine_test(Per,Nitr)
!
      use data_structure3D
      use element_data
!
      implicit none
!
      real(8), intent(in) :: Per
      integer, intent(in) :: Nitr
!
      integer :: kref_prism(3) = (/11,10,1/)
      integer :: mdle_list(NRELES)
!
      integer :: i,iel,icnt,idx,kref
      integer :: mdle,nelts,nref,npx,npy,npz
      real(8) :: x,error_norm,resid_norm
!
      integer :: iprint = 0
!
      write(*,*) 'test not working at the moment (needs revision)'
      goto 99
!
      do i=1,Nitr
!
!  .....collect all mdle
        mdle = 0
        do iel=1, NRELES
          call nelcon(mdle, mdle)
          mdle_list(iel) = mdle
        enddo
!
!  .....set number of refinements to be done
        icnt = 0;        kref = 0;
        x = 0.0;
        nref = int(NRELES*Per)
        nelts = NRELES
!
        if (iprint.eq.1) then
           write(*,7010) i, nelts, nref 
 7010      format('random_refine_test: i=',i3,' nelts=',i6,' nref=',i6) 
        endif

        do while (nref.gt.0)

!  .......random number generage from 0.0 to 1.0
          call random_seed
          call random_number(x)

!  .......pick one mdle from list
          idx  = int(nelts*x+1)
          mdle = mdle_list(idx)
          if (NODES(mdle)%ref_kind.eq.0) then
            select case (NODES(mdle)%ntype)
            case (MDLN)
              call get_isoref(mdle, kref)
              call refine(mdle, kref)
            case (MDLP)
              kref = kref_prism(mod(int(x*100),3)+1)
              call refine(mdle, kref)
            case default
              write(*,*) 'random_refine_test : MDLE IS NOT MDLE'
              stop
            end select
            nref = nref - 1
            if (iprint.eq.1) then
              write(*,7000) mdle, S_Type(NODES(mdle)%ntype), kref
 7000         format('mdle =',i6,' ', a5, 'kref=', i3)
            endif
          endif

!  .......loop exit condition to prevent infinite loop
          icnt = icnt + 1
          if (icnt.gt.(nelts*1000)) then
            exit
          endif
        enddo

!  .....close mesh and solve the problem
        call close
        call solve1(1)

!        call exact_error(error_norm,resid_norm)
 7001   format('random_refine_test: npx,npy,npz=',3i2, &
               ' error_norm,resid_norm = ',2e12.5)
        write(*,7001) npx,npy,npz,error_norm,resid_norm
        if ((error_norm/resid_norm).gt.1.d-6) then
          write(*,*) 'random_refine_test : TEST FAILED'
          return
        endif
!
      enddo
      
      99 continue
!
   end subroutine random_refine_test
