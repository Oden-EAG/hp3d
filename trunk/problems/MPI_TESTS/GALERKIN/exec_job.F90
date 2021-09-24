!----------------------------------------------------------------------
! exec_job
!----------------------------------------------------------------------
      subroutine exec_job
!
      use test_param    
      use common_prob_data
      use data_structure3D
      use MPI           , only: MPI_COMM_WORLD
      use mpi_param     , only: RANK,ROOT,NUM_PROCS
      use par_mesh      , only: EXCHANGE_DOF,distr_mesh
      use zoltan_wrapper, only: zoltan_w_set_lb,zoltan_w_eval
!
      implicit none
      common /cexact/ NP1,NP2,NP3,ICOMP
      integer :: NP1,NP2,NP3,ICOMP
!
      integer :: i
!
!     work space for exact_error_mod
      integer :: isolflag(1:4), nrcomp,ptot
      real(8) :: err,rnorm
!
!----------------------------------------------------------------------
!
      EXCHANGE_DOF = .false.
!
      if(RANK .eq. ROOT) then
        write(*,*) '=================='
        write(*,*) 'exec_job: starting'
        write(*,*) '=================='
      endif
!
!  ...distribute mesh initially
      call distr_mesh
!
!  ...polynomial order
      write(*,*) 'exec_job: IP = ',IP
!
!  ...loop through energy spaces
      do i=1,4
        select case(i)
!
!  .....H1 projection
        case(1)
          PHYSAm(1:4) = (/.true.,.false.,.false.,.false./)
          isolflag = (/1,0,0,0/)
          nrcomp=1; ptot = IP
          write(*,*) 'H1 PROJECTIONS........................'
!
!  .....H(curl) projection
        case(2)
          PHYSAm(1:4) = (/.false.,.true.,.false.,.false./)
          isolflag = (/0,1,0,0/)
          nrcomp=3; ptot = IP-1
          write(*,*) 'H(curl) PROJECTIONS...................'
!
!  .....H(div) projection
        case(3)
          PHYSAm(1:4) = (/.false.,.false.,.true.,.false./)
          isolflag = (/0,0,1,0/)
          nrcomp=3; ptot = IP-1
          write(*,*) 'H(div) PROJECTIONS....................'
!
!  .....L2 projection
        case(4)
          PHYSAm(1:4) = (/.false.,.false.,.false.,.true./)
          isolflag = (/0,0,0,1/)
          nrcomp=1; ptot = IP-1
          write(*,*) 'L2 PROJECTIONS........................'
        end select
!
!  .....loop through components
        do ICOMP = 1,nrcomp
!
!  .......loop through orders
          do NP1=0,ptot
            do NP2=0,ptot
              do NP3=0,ptot
                if (NP1+NP2+NP3.ne.ptot) cycle
                select case(i)
                case(1)
                  write(*,*) '     H1-PROJECTION: NP1,NP2,NP3       = ',NP1,NP2,NP3
                case(2)
                  write(*,*) 'H(curl)-PROJECTION: NP1,NP2,NP3,ICOMP = ',NP1,NP2,NP3,ICOMP
                case(3)
                  write(*,*) ' H(div)-PROJECTION: NP1,NP2,NP3,ICOMP = ',NP1,NP2,NP3,ICOMP
                case(4)
                  write(*,*) '     L2-PROJECTION: NP1,NP2,NP3       = ',NP1,NP2,NP3
                end select

!  .............project
                call solve1(1) !par_mumps_sc('G')
!
!  .............compute the error
                call exact_error_mod(isolflag, err,rnorm)
                if (sqrt(err/rnorm).gt.1.d-12) then
                  write(*,*) 'exec_case: VERIFICATION OF SHAPE FUNCTIONS HAS FAILED'
                  write(*,*) '           i,ICOMP,NP1,NP2,NP3 = ',i,ICOMP,NP1,NP2,NP3
                  write(*,*) '           sqrt(err),sqrt(rnorm) = ',sqrt(err),sqrt(rnorm)
                  call pause
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
!
      if(RANK .eq. ROOT) then
        write(*,*)
        write(*,*) '=================='
        write(*,*) 'exec_job: finished'
        write(*,*) '=================='
        write(*,*)
      endif
!
      end subroutine exec_job
