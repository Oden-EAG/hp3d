!----------------------------------------------------------------------
!                                                                     
!     routine name      - uweak_acoustics_module
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 17
!                                                                     
!     purpose:          - module setting up the workspace 
!                         for the UW formulation
!                                                                    
!
!----------------------------------------------------------------------
! 
   module uweak_acoustics_module

   use parametersDPG
#if C_MODE
#define V_TYPE complex*16
#else     
#define V_TYPE real*8
#endif
!
!..H1 discontinuous shape functions
   real*8 shapHH(MAXbrickHH),gradHH(3,MAXbrickHH)
!$OMP THREADPRIVATE (shapHH,gradHH)
!..H(div) discontinuous test functions
   real*8 :: shapVV(3,MAXbrickVV),divVV(MAXbrickVV)
!$OMP THREADPRIVATE (shapVV,divVV)
!
!..max dimension of enriched test space
   integer, parameter :: MAXtest = MAXbrickHH+MAXbrickVV
!      
!..stiffnes matrix for the local Gram matrix in packed format
   V_TYPE AP(MAXtest*(MAXtest+1)/2)
!$OMP THREADPRIVATE (AP)      
!
!..load vector for the enriched space
   V_TYPE BLOADHV(MAXtest)
!$OMP THREADPRIVATE (BLOADHV)      
!
!..stiffnes matrices for the enriched test space
   V_TYPE STIFFHV_H(MAXtest,MAXbrickH)
!$OMP THREADPRIVATE (STIFFHV_H)  
!    
   V_TYPE STIFFHV_V(MAXtest,MAXbrickV)
!$OMP THREADPRIVATE (STIFFHV_V)      
!
   V_TYPE STIFFHV_Q(MAXtest,4*MAXbrickQ)
!$OMP THREADPRIVATE (STIFFHV_Q)
!
   V_TYPE STIFF_ALL(MAXtest,MAXbrickH+MAXbrickV+4*MAXbrickQ+1)
!$OMP THREADPRIVATE (STIFF_ALL) 
    
!
!
   end module uweak_acoustics_module
    
