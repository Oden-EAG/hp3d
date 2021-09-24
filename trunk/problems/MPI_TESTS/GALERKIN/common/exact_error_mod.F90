!----------------------------------------------------------------------
!
!    routine            - exact_error_mod
!
!----------------------------------------------------------------------
!
!     latest revision:  - Jul 21
!
!     purpose:          - compute FE error, a modified version
!     parameters:
!       in:
!               Iflag   - a flag indicating for which physical attribute
!                         should error be computed
!       out:
!               Err     - total error (squared)
!               Rnorm   - total norm  (squared)
!
!----------------------------------------------------------------------
!
      subroutine exact_error_mod(Iflag, Err,Rnorm)
!
      use data_structure3D
      use common_prob_data
      use environment
      use assembly_sc, only: NRDOF_TOT,NRDOF_CON
      use par_mesh   , only: DISTRIBUTED,HOST_MESH
      use mpi_param  , only: ROOT,RANK
      use MPI        , only: MPI_SUM,MPI_COMM_WORLD,MPI_REAL8
!
      implicit none
!
      integer :: Iflag(NR_PHYSA)
      real(8) :: Err,Rnorm
!
!  ...workspace for element_error routine
      real(8) :: errorH,rnormH, errorE,rnormE, &
                 errorV,rnormV, errorQ,rnormQ
      real(8) :: error_subd,rnorm_subd
!
      integer :: iel,mdle,subd,count,ierr
!
      integer :: iprint=0
!
!----------------------------------------------------------------------
!
!  ...fetch active elements
      if (DISTRIBUTED .and. (.not. HOST_MESH)) then
        if (RANK .eq. ROOT) then
          write(*,*) 'exact_error_mod: mesh is distributed. ', &
                                      'computing error in parallel...'
        endif
      else
        if (RANK .ne. ROOT) goto 90
        write(*,*) 'exact_error_mod: mesh is not distributed (or on host). ', &
                                    'computing error on host...'
        ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
        NRELES_SUBD = NRELES
      endif
!
      error_subd = 0.d0; rnorm_subd = 0.d0
!
!$OMP PARALLEL DO                            &
!$OMP PRIVATE(errorH,errorE,errorV,errorQ,   &
!$OMP         rnormH,rnormE,rnormV,rnormQ)   &
!$OMP REDUCTION(+:error_subd,rnorm_subd)
      do iel=1,NRELES_SUBD
        call element_error(ELEM_SUBD(iel),Iflag,           &
                           errorH,errorE,errorV,errorQ,    &
                           rnormH,rnormE,rnormV,rnormQ)
        error_subd = error_subd + errorH + errorE + errorV + errorQ
        rnorm_subd = rnorm_subd + rnormH + rnormE + rnormV + rnormQ
      enddo
!$OMP END PARALLEL DO
!
      Err = 0.d0; Rnorm = 0.d0
      if (DISTRIBUTED .and. (.not. HOST_MESH)) then
        count = 1
        call MPI_REDUCE(error_subd,Err  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(rnorm_subd,Rnorm,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      else
        Err   = error_subd
        Rnorm = rnorm_subd
      endif
!
      if (iprint.eq.1) then
        if (RANK .eq. ROOT) then
          write(*,7020) NRDOF_TOT,sqrt(Err),sqrt(Rnorm)
 7020     format('exact_error_mod: NRDOF_TOT, L2 ERROR AND NORM = ',i8,3x,2es12.5)
          write(*,7030) NRDOF_TOT,sqrt(Err/Rnorm)
 7030     format('exact_error_mod: NRDOF_TOT, RELATIVE L2 ERROR = ',i8,3x,es12.5)
        endif
      endif
!
   90 continue
!
      end subroutine exact_error_mod
