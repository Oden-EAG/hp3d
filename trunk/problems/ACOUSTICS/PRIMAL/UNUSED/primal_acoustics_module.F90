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
   module primal_acoustics_module

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
!
!..stiffnes matrix for the local Riesz H1 matrix in Lapack format
   V_TYPE AP(MAXbrickHH*(MAXbrickHH+1)/2)
!$OMP THREADPRIVATE (AP)      
!
!..load vector for the enriched space
   V_TYPE BLOADH(MAXbrickHH),BLOADHc(MAXbrickHH)
!$OMP THREADPRIVATE (BLOADH,BLOADHc)      
!
!..stiffnes matrices for the enriched test space
   V_TYPE STIFFHH(MAXbrickHH,MAXbrickH)
!$OMP THREADPRIVATE (STIFFHH)  
!    
   V_TYPE STIFFHV(MAXbrickHH,MAXbrickV)
!$OMP THREADPRIVATE (STIFFHV)      
!
   V_TYPE STIFF_ALL(MAXbrickHH,MAXbrickH+MAXbrickV+1)
!$OMP THREADPRIVATE (STIFF_ALL)      
!
!
   end module primal_acoustics_module
    
