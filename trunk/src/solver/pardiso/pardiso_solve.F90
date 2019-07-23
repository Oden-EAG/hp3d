!------------------------------------------------------------------------
!
!    routine name       - pardiso_solve
!
!------------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - routine solves the global linear system 
!                         using pardiso solver. The matrix is given in 
!                         CSR format
!    arguments
!          in           - Ia     : array of row pointers
!                       - Ja     : array of column indices of the matrix
!                       - A      : array of nonzero matrix values
!                       - Nnz    : number of non zero values 
!                       - N      : actual size of the matrix
!                       - Nrhs   : supported???
!          in/out       - B      : Right hand size/solution array 
!
!------------------------------------------------------------------------
#include "implicit_none.h"
subroutine pardiso_solve(ia,ja,a,type,nnz,n,nrhs, b)
!
   use assembly_sc, only: IPRINT_TIME
!
   implicit none
!
!..problem workspace   
   integer , intent(in)    :: Nnz,N,Nrhs
   integer , intent(in)    :: Ia(N+1),Ja(Nnz)
   VTYPE   , intent(in)    :: A(Nnz)
   VTYPE   , intent(inout) :: B(N,Nrhs)
!
   VTYPE :: x(N,Nrhs)
!
!..pardiso workspace
!..internal solver memory pointer for 64-bit architectures
   integer*8       :: pt(64)
   integer         :: iparm(64)
   integer         :: mtype,maxfct,mnum,phase,error,msglvl
   integer         :: perm(N)
   character*1     :: type
   real*8          :: tm
!
!..timer
   real(8) :: MPI_Wtime,start_time,end_time
!
!------------------------------------------------------------------------
!
#if C_MODE    
   select case(type)
      case('S')  
         mtype   = 6  ! complex symmetric
      case('H')  
         mtype   = 4  ! complex hermitian positive definite
      case('I')  
         mtype   = -4 ! complex hermitian indefinite
      case default
         mtype   = 13  ! complex non-symmetric
   end select
#else 
   select case(type)
      case('S')  
         mtype   = -2  ! real symmetric indefinite
      case('H')  
         mtype   = 2  ! real symmetric positive definite
      case default  
         mtype   = 11  ! real non-symmetric
   end select   
#endif
!
!..initialize internal address pointers
   call pardisoinit(pt, mtype, iparm)
!   
   mnum    = 1
   maxfct  = 1
   msglvl  = 0   ! with statistical no information
!
!------------------------------------------------------------------------
!
   start_time = MPI_Wtime()
!
   phase     = 11  ! analysis
   iparm(8)  = 1   ! max numbers of iterative refinement steps
!
   call pardiso(pt,maxfct,mnum,mtype,phase,N,A,Ia,Ja,perm,Nrhs,         &
               iparm,msglvl,B,x,error)
!
   end_time = MPI_Wtime()
   tm = end_time-start_time
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1000) tm
 1000 format(' Analysis       : ',f12.5,'  seconds')
   endif 
!
!------------------------------------------------------------------------
!
   start_time = MPI_Wtime()
!
   phase     = 22  ! factorization
   call pardiso(pt,maxfct,mnum,mtype,phase,N,A,Ia,Ja,perm,Nrhs,         &
               iparm,msglvl,B,x,error)
!
   end_time = MPI_Wtime()
   tm = end_time-start_time
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1001) tm
 1001 format(' Factorization  : ',f12.5,'  seconds')
   endif 
!
!------------------------------------------------------------------------
!
   start_time = MPI_Wtime()
!
   phase     = 33  ! F/B Substitution
   call pardiso(pt,maxfct,mnum,mtype,phase,N,A,Ia,Ja,perm,Nrhs,         &
               iparm,msglvl,B,x,error)
!   
   end_time = MPI_Wtime()
   tm = end_time-start_time
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1002) tm
 1002 format(' B/F Solve      : ',f12.5,'  seconds')
   endif 

   if (IPRINT_TIME .eq. 1) then
      write(*,1003) real(max(iparm(15),iparm(16)+iparm(17)),8)/1e6
 1003 format(' Memory in GB   : ',2x,1f7.2)
   endif

   phase     = -1    ! release internal memory
   call pardiso(pt,maxfct,mnum,mtype,phase,N,A,Ia,Ja,perm,Nrhs,         &
                iparm,msglvl,B,x,error)
! 
   B = x
!  
end subroutine pardiso_solve
!
!
!------------------------------------------------------------------------
!
!    routine name       - pardiso_solve_vect
!
!------------------------------------------------------------------------
!
!    latest revision    - July 17
!
!    purpose            - routine solves the global linear system 
!                         using pardiso solver. The matrix is given in 
!                         CSR format. Multiple RHS are given in 1D array column-wise
!    arguments
!          in           - ia     : array of row pointers
!                         ja     : array of column indices of the matrix
!                         a      : array of nonzero matrix values
!                         nnz    : number of non zero values 
!                         n      : actual size of the matrix
!          in/out         b      : Right hand size/solution array 
!
!------------------------------------------------------------------------
subroutine pardiso_solve_vect(ia,ja,a,type,nnz,b_vect,n,nrhs,B,bb)
!
   implicit none
!
!..problem workspace   
   integer  :: n, nnz, nrhs, i, j
   integer  :: ia(n+1), ja(nnz)
   VTYPE    :: a(nnz), b_vect(n*nrhs), x(n*nrhs), B(n,nrhs-1), bb(n)
!
!..pardiso workspace
!..internal solver memory pointer for 64-bit architectures
   integer*8    :: pt(64)
   integer      :: iparm(64)
   integer      :: mtype, maxfct, mnum, phase, error, msglvl
   integer      :: perm(n)
   character*1  :: type
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
#if C_MODE    
   select case(type)
   case('S')  
      mtype   = 6  ! complex symmetric
   case('H')  
      mtype   = 4  ! complex hermitian positive definite
   case('I')  
      mtype   = -4 ! complex hermitian indefinite
   case default
      mtype   = 13  ! complex non-symmetric
   end select
#else 
   select case(type)
   case('S')  
      mtype   = -2  ! real symmetric indefinite
   case('H')  
      mtype   = 2  ! real symmetric positive definite
   case default  
      mtype   = 11  ! real non-symmetric
   end select   
#endif
!
!..initialize internal address pointers                            
   call pardisoinit(pt, mtype, iparm)
!   
   mnum    = 1
   maxfct  = 1
   msglvl  = 0       ! with statistical no information
!   
   phase     = 13  ! analysis, factorization, solve
   iparm(8)  = 1   ! max numbers of iterative refinement steps
!   
   call pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,perm,nrhs,         &
               iparm,msglvl,b_vect,x,error)
! 
   do i = 1, nrhs-1
      do j = 1, n
         B(j,i) = x(j+(i-1)*n)
      enddo
   enddo
   bb = x(n*(nrhs-1)+1:n*(nrhs))
! 
   phase     = -1    ! release internal memory
   call pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,perm,nrhs,         &
                iparm,msglvl,b_vect,x,error)
! 
end subroutine pardiso_solve_vect
