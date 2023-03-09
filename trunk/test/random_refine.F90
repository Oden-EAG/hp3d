c
      program test_random_refine
c
         real(8) :: per
         integer :: nitr
c
c     ...percentage of elements to refine per iteration
         per  = 0.1
c     ...number of refinements (iterations)
         nitr = 5
c
         call random_refine_test(per,nitr)
c
      end program
c
c----------------------------------------------------------------------
c
c   routine name       - random_refine_test
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine to refine mesh randomly
c
c----------------------------------------------------------------------
c
      subroutine random_refine_test(Per,Nitr)
c
      use data_structure3D
      use element_data
c
      implicit none
c
      real(8), intent(in) :: Per
      integer, intent(in) :: Nitr
c
      integer :: kref_prism(1:3) = (/11,10,1/)
      integer, allocatable ::  mdle_list(:)
c
      integer :: i,iel,icnt,kref,mdle,nelts,nref
      real(8) :: x
c
      integer :: iprint = 0
c
      write(*,*) 'test not working at the moment (needs revision)'
c
      do i=1,Nitr
c
c  .....collect all mdle
        allocate(mdle_list(NRELES))
        mdle = 0
        do iel=1, NRELES
          call nelcon(mdle, mdle)
          mdle_list(iel) = mdle
        enddo
c
c  .....set number of refinements to be done
        icnt = 0;        kref = 0;
        x = 0.0;
        nref = int(NRELES*Per)
        nelts = NRELES
c
        if (iprint.eq.1) then
           write(*,7010) i, nelts, nref 
 7010      format('random_refine_test: i=',i3,' nelts=',i6,' nref=',i6) 
        endif
c
c        do while (nref.gt.0)
c
c  .......random number generage from 0.0 to 1.0
c          call random_seed
c          call random_number(x)
c  
c  .......pick one mdle from list
c          idx  = int(nelts*x+1)
c          mdle = mdle_list(idx)
c          if (NODES(mdle)%ref_kind.eq.0) then
c            select case (NODES(mdle)%ntype)
c            case (MDLN)
c              call get_isoref(mdle, kref)
c              call refine(mdle, kref)
c            case (MDLP)
c              kref = kref_prism(mod(int(x*100),3)+1)
c              call refine(mdle, kref)
c            case default
c              write(*,*) 'random_refine_test : MDLE IS NOT MDLE'
c              stop
c            end select
c            nref = nref - 1
c            if (iprint.eq.1) then
c              write(*,7000) mdle, S_Type(NODES(mdle)%ntype), kref
c 7000         format('mdle =',i6,' ', a5, 'kref=', i3)
c            endif
c          endif
c
c  .......loop exit condition to prevent infinite loop
c          icnt = icnt + 1
c          if (icnt.gt.(nelts*1000)) then
c            exit
c          endif
c        enddo
c
        deallocate(mdle_list)
c
c  .....close mesh and solve the problem
c        call close
c        call solve1(1)
c
c        call exact_error(error,rnorm)
c 7001   format('random_refine_test: NPX,NPY,NPZ=',3i2,
c     .       ' error,rnorm = ',2e12.5)
c        write(*,7001) NPX,NPY,NPZ,error,rnorm
c        if ((error/rnorm).gt.1.d-6) then
c          write(*,*) 'random_refine_test : TEST FAIL !!!!!!!!!!!!!!!!'
c          call pause
c          call result
c          return
c        endif
c        
      enddo
c     
      end subroutine random_refine_test
