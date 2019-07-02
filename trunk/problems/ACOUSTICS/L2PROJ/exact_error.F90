!----------------------------------------------------------------------
!                                                                     
!    routine            - exact_error
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 17
!                                                                     
!     purpose:          - compute and print exact error
!                                                                    
!----------------------------------------------------------------------
!   
   subroutine exact_error
!  
   use data_structure3D
   use common_prob_data
   use environment
   use m_assembly, ONLY: NRDOF_TOT, NRDOF_CON
!   
   implicit none
!
   integer  :: iflag
   integer, allocatable :: mdle_list(:)
!..work space for elem_error
   real*8 :: errorH,rnormH,errorE,rnormE,      &
             errorV,rnormV,errorQ,rnormQ
   real*8 :: error, rnorm
   integer :: iel, mdle
!
!----------------------------------------------------------------------
!
   iflag=1;
!
   allocate(mdle_list(NRELES))
   error = 0.d0; rnorm = 0.d0
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      mdle_list(iel) = mdle
   enddo
!$OMP PARALLEL DEFAULT(PRIVATE)       &
!$OMP SHARED(NRELES,iflag,mdle_list) &
!$OMP REDUCTION(+:error,rnorm)
!$OMP DO
   do iel=1,NRELES
      call element_error(mdle_list(iel),iflag,       &
                         errorH,errorE,errorV,errorQ,   &
                         rnormH,rnormE,rnormV,rnormQ)  
      error = error + errorQ
      rnorm = rnorm + rnormQ
   enddo
!$OMP END DO        
!$OMP END PARALLEL

   write(*,7020) NRDOF_TOT,sqrt(error),sqrt(rnorm)
 7020 format('exact_error: TOTAL NUMBER OF DOF, L2 ERROR AND NORM = ',i8,3x,2e12.5)
   write(*,7030) NRDOF_TOT ,sqrt(error/rnorm)
 7030 format('exact_error: TOTAL NUMBER OF DOF, RELATIVE L2 ERROR = ',i8,3x,f8.4)

   deallocate(mdle_list)
!
   end subroutine exact_error
!
!
!
