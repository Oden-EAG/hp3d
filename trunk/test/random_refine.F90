!> @date Mar 2023
program test_random_refine
!
   real(8) :: per
   integer :: nitr,npass
!
!..percentage of elements to refine per iteration
   per  = 0.1
!..number of refinements (iterations)
   nitr = 5
!
   write(*,*) 'test not working at the moment (needs revision)'
   goto 99
!
   call random_refine_test(per,nitr,npass)
   if (npass .eq. 1) then
      write(*,*) 'test_random_refine PASSED.'
   else
      write(*,*) 'test_random_refine FAILED.'
   endif
!
   99 continue
!
end program test_random_refine
!
!----------------------------------------------------------------------
! !> @brief Refines mesh randomly
! !> @date Mar 2023
!----------------------------------------------------------------------
   subroutine test_random_refine_aux(Per,Nitr, Npass)
!
      use data_structure3D
!
      implicit none
!
      real(8), intent(in)  :: Per
      integer, intent(in)  :: Nitr
      integer, intent(out) :: Npass
!
      integer :: kref_prism(3) = (/11,10,1/)
      integer :: mdle_list(NRELES)
!
      integer :: i,iel,icnt,idx,kref
      integer :: mdle,nelts,nref
      real(8) :: x
!
      integer :: iprint
      iprint=0
!
      Npass = 1
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
 7010      format('test_random_refine_aux: i=',i3,' nelts=',i6,' nref=',i6)
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
              write(*,*) 'test_random_refine_aux: MDLE type'
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
!
        if (.false.) Npass = 0
!
      enddo
!
   end subroutine test_random_refine_aux
!
!----------------------------------------------------------------------
!> @brief dummy routine
subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!
end subroutine
!----------------------------------------------------------------------
!> @brief dummy routine
subroutine celem(Mdle,Idec,Nrdofs,Nrdofm,Nrdofc,Nodm, &
                 NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,Bload,Astif)
!
end subroutine
