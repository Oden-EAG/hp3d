c----------------------------------------------------------------------
c
c   routine name       - random_refine_test
c
c----------------------------------------------------------------------
c
c   latest revision    - June 10
c
c   purpose            - routine to refine mesh randomly
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine random_refine_test(Per, Nitr)
c
      use data_structure3D
      use element_data
      use TETRA
c     
#include "syscom.blk"
      integer::kref_prism(1:3) = (/11,10,1/)
      integer, allocatable ::  mdle_list(:)
c
      iprint = 1
c
      do i=1, Nitr
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
        do while (nref.gt.0)
c
c  .......random number generage from 0.0 to 1.0
          call random_seed
          call random_number(x)
c  
c  .......pick one mdle from list
          idx  = int(nelts*x+1) 
          mdle = mdle_list(idx)
          if (NODES(mdle)%ref_kind.eq.0) then
            select case (NODES(mdle)%type) 
            case ('mdln ')
              call get_isoref(mdle, kref)
              call refine(mdle, kref)
            case ('mdlp ')
              kref = kref_prism(mod(int(x*100),3)+1)
              call refine(mdle, kref)
            case default
              write(*,*) 'random_refine_test : MDLE IS NOT MDLE'
              stop1
            end select
            nref = nref - 1
            if (iprint.eq.1) then
              write(*,7000) mdle, NODES(mdle)%type, kref
 7000         format('mdle =',i6,' ', a5, 'kref=', i3)
            endif
          endif
c
c  .......loop exit condition to prevent infinite loop
          icnt = icnt + 1
          if (icnt.gt.(nelts*1000)) then
            exit
          endif
        enddo
c
        deallocate(mdle_list)
c
c  .....close mesh and solve the problem
        call close
        call solve1(1)
        call exact_error(error,rnorm)
 7001   format('random_refine_test: NPX,NPY,NPZ=',3i2,
     .       ' error,rnorm = ',2e12.5)
        write(*,7001) NPX,NPY,NPZ,error,rnorm
        if (error/rnorm.gt.1.d-6) then
          write(*,*) 'random_refine_test : TEST FAIL !!!!!!!!!!!!!!!!'
          call pause
          call result
          return
        endif
c        
      enddo
c     
      end 


