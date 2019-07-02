c//////////////////////////////////////////////////////////////////////
c  P. Gatto, Dec 14: prehistoric routine, do NOT use unless you revise
c                    it!
c//////////////////////////////////////////////////////////////////////
c
c----------------------------------------------------------------------
c
c   routine name       - pstep
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 09
c
c   purpose            - routine performs a single p-adaptive
c                        refinement
c
c   arguments :
c     in:
c         Istep        - step number
c         Solver       - solver to use ('front','mumps')
c         Nr_RHS       - number of right-hand sides
c     out:
c                      - changes in data structure
c
c----------------------------------------------------------------------
c
      subroutine pstep(Istep,Solver,Nr_RHS)
c
      use control
      use data_structure3D
#include "syscom.blk"
c
      character (len=4) :: Solver
c
c  ...work space for element_norm
      dimension dnorm(5,NR_PHYSA)
c
c  ...global errors
      save gnorm
      dimension gnorm(5,10,20)
c
      allocatable :: element_error(:)
      allocatable :: new_order(:)
c
c----------------------------------------------------------------------
c
      if (NR_PHYSA.gt.10) then
        write(*,*) 'pstep: INCREASE DIMENSION OF gnorm'
        stop 1
      endif
      iprint=1
c
c  ...solve the problem on the current mesh
      write(*,*) 'pstep: COARSE GRID SOLVE...'
      select case(Solver)
      case('front')
        call solve1(Nr_RHS)
      case('mumps')
        call mumps_solve(Nr_RHS)
      end select
c
c  ...dumpout the current data structure
      call dumpout_hp3d('files/dumpc3Dhp')
c
c  ...shift solution dof
      write(*,*) 'pstep: SHIFTING DOF...'
      call shift_dof(Nr_RHS,Nr_RHS)
c
c  ...perform a global p-refinement
      write(*,*) 'pstep: GLOBAL p-REFINEMENT...'
      allocate(new_order(NRELES))
      call setp(new_order)
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call nodmod(mdle, new_order(iel))
      enddo
      call enforce_min_rule
c
c  ...solve on the new mesh
      write(*,*) 'pstep: FINE GRID SOLVE...'
      select case(Solver)
      case('front')
        call solve1(Nr_RHS)
      case('mumps')
        call mumps_solve(Nr_RHS)
      end select
      call dumpout_hp3d('files/dumpf3Dhp')
c
c  ...compute the difference between the old and the new solutions
      write(*,*) 'pstep: COMPUTING ERROR...'
      allocate(element_error(NRELES))
      gnorm(1:5,1:NR_PHYSA,Istep) = 0.d0
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call element_norm(mdle,Nr_RHS,Nr_RHS,1, dnorm)
        gnorm(1:5,1:NR_PHYSA,Istep) = gnorm(1:5,1:NR_PHYSA,Istep)
     .                              + dnorm(1:5,1:NR_PHYSA)
        element_error(iel) = 0.d0
        do iphys=1,NR_PHYSA
          element_error(iel) = element_error(iel)
     .                       + dnorm(5,iphys)
        enddo
      enddo
c
c  ...determine the maximimum element error
      error_max = 0.d0
      do iel=1,NRELES
        error_max = max(error_max,element_error(iel))
      enddo
c
c  ...dumpin the coarse grid mesh and perform a selective p-refinement
      write(*,*) 'pstep: OPTIMAL p-REFINEMENT...'
      error_threshold = .7d0*error_max
      call dumpin_hp3d('files/dumpc3Dhp')
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        if (element_error(iel).gt.error_threshold) then
          call find_domain(mdle, no)
          write(*,7021) mdle,no,new_order(iel)
 7021     format('pstep: REFINING mdle ',i6,' IN SUBDOMAIN ',i2,
     .           ' NEW ORDER = ',i3)
          call nodmod(mdle, new_order(iel))
        endif
      enddo
      call enforce_min_rule
c
      deallocate(new_order,element_error)
      write(*,*) 'pstep: ERROR BY ATTRIBUTE'
      do i=1,Istep
        if (NEXACT.ne.0) then
          write(*,7011) i,(sqrt(gnorm(3,j,i)/gnorm(2,j,i)),j=1,NR_PHYSA)
          write(*,7012)   (sqrt(gnorm(4,j,i)/gnorm(2,j,i)),j=1,NR_PHYSA)
 7011     format('STEP = ',i2,' FINE GRID ERROR IN %   = ',10e12.5)
 7012     format('       ',2x,' COARSE GRID ERROR IN % = ',10e12.5)
        else
          write(*,7013) i,(gnorm(5,j,i),j=1,NR_PHYSA)
          write(*,7014) i,(sqrt(gnorm(5,j,i)/gnorm(1,j,i)),j=1,NR_PHYSA)
 7013     format('STEP = ',i2,' FINE - COARSE GRID ERROR  = ',10e12.5)
 7014     format('STEP = ',i2,' FINE - COARSE GRID ERROR IN % = ',
     .           10e12.5)
        endif
      enddo
c
c
      end
