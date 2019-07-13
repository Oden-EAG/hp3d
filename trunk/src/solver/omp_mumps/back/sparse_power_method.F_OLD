!
! -----------------------------------------------------------------------
!
!    routine name       - sparse_power_iteration
!
! -----------------------------------------------------------------------
!
!    latest revision    - FEB 17
!
!    purpose            - routine computes the maximum eigenvalue using 
!                         power iteration method
!
!-----------------------------------------------------------------------



   subroutine sparse_power_iteration(SA,IA,JA,nz,n,lamda)
!
   implicit none      
! 
#if C_MODE
   complex*16 :: SA(nz), x(n), y(n)
   real*8     :: dznrm2
#else   
   real*8     :: SA(nz), x(n), y(n)
   real*8     :: dnrm2
#endif
   integer    :: IA(nz), JA(nz), i, nz,n, maxit   
   real*8     :: lamda, lamda_old

!..workspace for mkl_?coomv
   character  :: matdescra(4)
! 
!-----------------------------------------------------------------------
! 
   matdescra(1) = 'G' ;  matdescra(2) = 'U' ;
   matdescra(3) = 'N' ;  matdescra(4) = 'F'     
   maxit = 200

!..initialize x

#if C_MODE
   x = (1.d0,0.d0)
   lamda = dznrm2(n,x,1)
#else
   x = 1.d0
   lamda = dnrm2(n,x,1)
#endif      

   x = x/lamda
   lamda_old = lamda

   do i=1,maxit

#if C_MODE
      call mkl_zcoomv('N',n,n,(1.d0,0.d0),matdescra,SA,IA,JA,nz,x,(0.d0,0.d0),y)
#else
      call mkl_dcoomv('N',n,n,1.d0,matdescra,SA,IA,JA,nz,x,0.d0,y)
#endif      

#if C_MODE
      lamda = dznrm2(n,y,1)
#else
      lamda = dnrm2(n,y,1)

#endif      
      if (abs(lamda_old-lamda)/abs(lamda) .le. 1e-10 )  then
         exit
      endif

      lamda_old = lamda
      x = y/lamda
      ! write(*,*) lamda

   enddo

    ! write(*,*) 'power_iteration: interations = ', i


   end subroutine sparse_power_iteration

!
! -----------------------------------------------------------------------
!
!    routine name       - sparse_inverse_iteration
!
! -----------------------------------------------------------------------
!
!    latest revision    - FEB 17
!
!    purpose            - routine computes the minimum eigenvalue using 
!                         inverse iteration method
!
!-----------------------------------------------------------------------



   subroutine sparse_inverse_iteration(n,lamda)

   use mumps

#if C_MODE
   complex*16 :: temp(n) 
   real*8     :: dznrm2
#else   
   real*8     :: temp(n)
   real*8     :: dnrm2
#endif
   integer    :: i, maxit
   real*8     :: lamda, lamda_old

   maxit = 20
!...store RHS
   temp = mumps_par%RHS
#if C_MODE
   mumps_par%RHS = (1.0d0,0.0d0)
   lamda = dznrm2(n,mumps_par%RHS,1)
#else
   mumps_par%RHS = 1.0d0
   lamda = dnrm2(n,mumps_par%RHS,1)
#endif   
! 
   lamda_old = lamda
   mumps_par%RHS = mumps_par%RHS/lamda

   do i=1,maxit
      ! write(*,*) 'i = ', i

      mumps_par%JOB = 3
#if C_MODE
      call zmumps(mumps_par)
      lamda = dznrm2(n,mumps_par%RHS,1)
#else
      call dmumps(mumps_par)
      lamda = dnrm2(n,mumps_par%RHS,1)
#endif      

      if (abs(lamda_old-lamda)/abs(lamda) .le. 1e-12 ) then 
         exit
      endif   
      lamda_old = lamda
      ! write(*,*) 1.0d0/lamda
      mumps_par%RHS = mumps_par%RHS/lamda
   enddo
      
   mumps_par%RHS = temp
   lamda = 1.0d0/lamda

   ! write(*,*) 'inverse_iteration: interations = ', i

   end subroutine sparse_inverse_iteration






   !
! -----------------------------------------------------------------------
!
!    routine name       - sparse_power_iteration_indefinite
!
! -----------------------------------------------------------------------
!
!    latest revision    - FEB 17
!
!    purpose            - routine computes the maximum eigenvalue using 
!                         power iteration method
!
!-----------------------------------------------------------------------



   subroutine sparse_power_iteration_indefinite(SA,IA,JA,nz,n,lamda)
!
   implicit none      
! 
#if C_MODE
   complex*16 :: SA(nz), x(n), y(n)
   real*8     :: dznrm2
#else   
   real*8     :: SA(nz), x(n), y(n)
   real*8     :: dnrm2
#endif
   integer    :: IA(nz), JA(nz), i, nz,n, maxit   
   real*8     :: lamda, lamda_old

!..workspace for mkl_?coomv
   character  :: matdescra(4)
! 
!-----------------------------------------------------------------------
! 
   matdescra(1) = 'G' ;  matdescra(2) = 'U' ;
   matdescra(3) = 'N' ;  matdescra(4) = 'F'     
   maxit = 50

!..initialize x

#if C_MODE
   x = (1.d0,0.d0)
   lamda = dznrm2(n,x,1)
#else
   x = 1.d0
   lamda = dnrm2(n,x,1)
#endif      

   x = x/lamda
   lamda_old = lamda

   do i=1,maxit

#if C_MODE
      call mkl_zcoomv('N',n,n,(1.d0,0.d0),matdescra,SA,IA,JA,nz,x,(0.d0,0.d0),y)
      call mkl_zcoomv('N',n,n,(1.d0,0.d0),matdescra,conjg(SA),JA,IA,nz,y,(0.d0,0.d0),x)
      y=x
#else
      call mkl_dcoomv('N',n,n,1.d0,matdescra,SA,IA,JA,nz,x,0.d0,y)
#endif      

#if C_MODE
      lamda = dznrm2(n,y,1)
#else
      lamda = dnrm2(n,y,1)

#endif      
      if (abs(lamda_old-lamda)/abs(lamda) .le. 1e-10 )  then
         exit
      endif

      lamda_old = lamda
      x = y/lamda
   enddo


   lamda = sqrt(lamda)



   end subroutine sparse_power_iteration_indefinite

!
! -----------------------------------------------------------------------
!
!    routine name       - sparse_inverse_iteration_indefinite
!
! -----------------------------------------------------------------------
!
!    latest revision    - FEB 17
!
!    purpose            - routine computes the minimum eigenvalue using 
!                         inverse iteration method
!
!-----------------------------------------------------------------------



   subroutine sparse_inverse_iteration_indefinite(n,lamda)

   use mumps

#if C_MODE
   complex*16 :: temp(n) , tempA(mumps_par%NZ),A(mumps_par%NZ)  
   real*8     :: dznrm2
#else   
   real*8     :: temp(n), tempA(mumps_par%NZ),A(mumps_par%NZ) 
   real*8     :: dnrm2
#endif
   integer    :: i, maxit, tempIA(mumps_par%NZ), tempJA(mumps_par%NZ)
   integer    :: IA(mumps_par%NZ), JA(mumps_par%NZ)
   real*8     :: lamda, lamda_old

   maxit = 20
   A = mumps_par%A
   IA = mumps_par%IRN
   JA = mumps_par%JCN

   nz = mumps_par%NZ
   nrdof = mumps_par%N
#if C_MODE   
   tempA = conjg(mumps_par%A)
#else
   tempA = mumps_par%A
#endif      
   tempIA = mumps_par%JCN
   tempJA = mumps_par%IRN

!...store RHS
   temp = mumps_par%RHS

#if C_MODE
   mumps_par%RHS = (1.0d0,0.0d0)
   lamda = dznrm2(n,mumps_par%RHS,1)
#else
   mumps_par%RHS = 1.0d0
   lamda = dnrm2(n,mumps_par%RHS,1)
#endif   
! 
   lamda_old = lamda
   mumps_par%RHS = mumps_par%RHS/lamda

   do i=1,maxit
!
#if C_MODE
      mumps_par%A = tempA
      mumps_par%IRN = tempIA
      mumps_par%JCN = tempJA
      mumps_par%JOB = 6
      call zmumps(mumps_par)
      mumps_par%A = A
      mumps_par%IRN = IA
      mumps_par%JCN = JA
      mumps_par%JOB = 6
!
      call zmumps(mumps_par)
      lamda = dznrm2(n,mumps_par%RHS,1)
#else
      call dmumps(mumps_par)
      lamda = dnrm2(n,mumps_par%RHS,1)
#endif      

      if (abs(lamda_old-lamda)/abs(lamda) .le. 1e-10 ) then 
         exit
      endif   
      lamda_old = lamda
      mumps_par%RHS = mumps_par%RHS/lamda
   enddo
      
   mumps_par%RHS = temp
   lamda = 1.0d0/lamda

   lamda = sqrt(lamda)

   ! write(*,*) 'inverse_iteration: interations = ', i

   end subroutine sparse_inverse_iteration_indefinite