! -----------------------------------------------------------------------

!      routine name      - pbisolver3

! -----------------------------------------------------------------------

!      latest revision:  - Sep 05
!      Author:      Jason Kurtz
!      purpose:          - Routine extracts a linear system (specified by
!                          an extraction vector) from a larger system,
!                          computes the Cholesky factorization of the
!                          extracted matrix, and solves for (possibly
!                          multiple) real or complex-valued right hand
!                          sides. When used repeatedly for a nested
!                          sequence of linear systems, this variant
!                          uses the previously computed factorization
!                          and solution.

!      arguments
!      in:   M           - current step
!            Mblock(M+1) - nrdof in the M-th linear system
!            Aglob       - global matrix
!            Ldglob      - leading dimension of Aglob, Bglob in calling
!                          routine
!            Ldwork      - leading dimension of Awork, Bwork in calling
!                          routine
!            Nextract    - extraction vector for extracting Awork, Bwork
!            Bglob       - global rhs
! #if C_MODE
!                          in complex mode, actual argument is complex
!                          and dummy argument is intentionally a
!                          double precision array, twice as long, storing
!                          alternating real and imaginary parts
! #endif
!            Nrhs        - number of rhs
!      out:  Awork       - Cholesky factorization of extracted matrix
!            Bwork       - solution vector(s)
!--------------------------------------------------------------------------------------------------------
#include "typedefs.h"

subroutine pbisolver3(M,Mblock,Aglob,Ldglob,Awork,Ldwork, &
                              Nextract,Bglob,Bwork,Nrhs)

use data_structure3D
implicit none
integer,    intent(in)  ::  M
integer,    intent(in)  ::  Mblock(*)
integer,    intent(in)  ::  Ldglob
integer,    intent(in)  ::  Ldwork
integer,    intent(in)  ::  Nrhs
integer,    intent(in)  :: Nextract(Ldwork)
real(8),dimension(Ldglob,Ldglob),   intent(in)  :: Aglob
real(8),dimension(Ldwork,Ldwork),   intent(out)  :: Awork
real(8),dimension(Ldglob),  intent(in)  :: Bglob
real(8),dimension(Ldwork),  intent(out)  :: Bwork
integer :: k1,k2,l1,l2,icomp,ivar,k,l,imag,kk,ll
integer :: n,info,offset
real(8), dimension(Ldwork * Ldwork) :: Awc

#if C_MODE
   icomp = 2
#else 
   icomp = 1
#endif
!-----------------------------------------------------------------------

!     EXTRACT NEWLY EXPOSED ROWS OF PROJECTION MATRIX AND RHS

!     rectangular blocks:

      Awc = [Awork]
      Awork = ZERO
      offset = 0

      do k2 = 1,Mblock(M)
        Awork(k2:Mblock(M),k2) = Awc(offset + 1: offset + Mblock(M) - k2 + 1)
        offset = k2 * Mblock(M) + k2
      enddo


!     loop through columns
      do k2=1,Mblock(M)
        l2 = Nextract(k2)

        ! loop through rows in newly exposed block
        do k1=Mblock(M)+1,Mblock(M+1)
          l1 = Nextract(k1)

          if (l1>=l2) then
            Awork(k1,k2) = Aglob(l1,l2)
          else
            Awork(k1,k2) = Aglob(l2,l1)
          endif

        enddo
      enddo
!
!  (lower) triangular block:
! 
!  loop through columns
      do k2=Mblock(M)+1,Mblock(M+1)
        l2 = Nextract(k2)
 
        !  loop through rows in lower triangle
        do k1=k2,Mblock(M+1)
          l1 = Nextract(k1)
 
          if (l1>=l2) then
            Awork(k1,k2) = Aglob(l1,l2)
            ! write(*,*) l1,l2,Aglob(l1,l2)
          else
            Awork(k1,k2) = Aglob(l2,l1)
          endif
 
        enddo
      enddo


    
! !     loop through equations
      do ivar=1,Nrhs

        !  loop through new rows
        do k=Mblock(M)+1,Mblock(M+1)
          l = Nextract(k)

          !  loop through components (real and possibly imaginary)
          do imag=1,icomp
            kk = icomp*((ivar-1)*Ldwork+(k-1))+imag
            ll = icomp*((ivar-1)*Ldglob+(l-1))+imag
            Bwork(kk) = Bglob(ll)
          enddo
        enddo
      enddo



!       !-----------------------------------------------------------------------
      
!       !     DO PREVIOUSLY POSTPONED UPDATES
      
            do n=1,M-1
              do k=1,n-1
      
                ! A(m,n) = A(m,n) - A(m,k)*A(n,k)^T
                !  dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                !  C <- ALPHA*op(A)*op(B) + BETA*C, C - M x N
      
                call dgemm('N','T',Mblock(M+1)-Mblock(M), &
                                  Mblock(n+1)-Mblock(n),  &
                                  Mblock(k+1)-Mblock(k),  &
                                 -1.d0,Awork(Mblock(M)+1,Mblock(k)+1),Ldwork, &
                                       Awork(Mblock(n)+1,Mblock(k)+1),Ldwork, &
                                  1.d0,Awork(Mblock(M)+1,Mblock(n)+1),Ldwork)
    
              enddo
      
              !  A(m,n) = A(m,n)*A(n,n)^{-T}
              !  dtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
              !  B <- ALPHA*B*op(A^{-1}), B - M x N
      
              call dtrsm('R','L','T','N',Mblock(M+1)-Mblock(M),&
                                        Mblock(n+1)-Mblock(n),&
                                  1.d0,Awork(Mblock(n)+1,Mblock(n)+1),Ldwork,&
                                       Awork(Mblock(M)+1,Mblock(n)+1),Ldwork)
      
            enddo

            do n=1,M-1
            
                    !  A(m,m) = A(m,m) - A(m,n)*A(m,n)^T
                    !  dsyrk(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
                    !  C <- ALPHA*A*A^T + BETA*C, C - N x N
              
                      call dsyrk('L','N',Mblock(M+1)-Mblock(M),Mblock(n+1)-Mblock(n), &
                                         -1.d0,Awork(Mblock(M)+1,Mblock(n)+1),Ldwork, &
                                          1.d0,Awork(Mblock(M)+1,Mblock(M)+1),Ldwork)
              
            enddo

          ! -----------------------------------------------------------------------
      
          !  FACTOR NEWLY EXPOSED TRIANGULAR BLOCK
      
          !     A(m,m) = L(m,m)*L(m,m)^T
          !     dpotrf(UPLO,N,A,LDA,INFO)
          !    A = L*L^T, A - N x N
            

            
            call dpotrf('L',Mblock(M+1)-Mblock(M), &
                        Awork(Mblock(M)+1,Mblock(M)+1),Ldwork,info)
            if (info/=0) then
              write(*,*)'pbisolver3: DPOTRF RETURNED INFO = ',info
              stop
            endif




        ! -----------------------------------------------------------------------
        
        !  BACKSOLVE
        
          if (Mblock(M)>0) then
    
            !  1. b_1 <- L_11^T*b_1
            !  DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
            !  X <- A*X 
            !  loop through equations (real and imaginary parts)
            do ivar=1,Nrhs
              do imag=1,icomp
                call dtrmv('L','T','N',Mblock(M),Awork(1,1),Ldwork, &
                            Bwork(icomp*((ivar-1)*Ldwork)+imag),icomp)
              enddo
            enddo
    
            !  2. b_2 <- b_2 - L_21*b_1
            !  DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
            !  Y <- ALPHA*A*X + BETA*Y

            do ivar=1,Nrhs
            do imag=1,icomp
              call dgemv('N',Mblock(M+1)-Mblock(M),Mblock(M), &
                  -1.d0,Awork(Mblock(M)+1,1),Ldwork, &
                        Bwork(icomp*((ivar-1)*Ldwork)+imag),icomp, &
                    1.d0,Bwork(icomp*((ivar-1)*Ldwork+Mblock(M))+imag),icomp)
            enddo
            enddo
    
          endif

              !  3. b_2 <- L_22^{-1}b_2
              !  DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
              !  X <- A^{-1}*X
        
                do ivar=1,Nrhs
                do imag=1,icomp
                  call dtrsv('L','N','N',Mblock(M+1)-Mblock(M),     &
                            Awork(Mblock(M)+1,Mblock(M)+1),Ldwork,  &
                            Bwork(icomp*((ivar-1)*Ldwork+Mblock(M))+imag),icomp)
                enddo
                enddo

              !  4. b_2 <- L_22^{-T}b_2
              !  DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
              !  X <- A^{-1}*X
              
                do ivar=1,Nrhs
                do imag=1,icomp
                  call dtrsv('L','T','N',Mblock(M+1)-Mblock(M),     &
                            Awork(Mblock(M)+1,Mblock(M)+1),Ldwork,  &
                            Bwork(icomp*((ivar-1)*Ldwork+Mblock(M))+imag),icomp)
                enddo
                enddo

                if (Mblock(M)>0) then
                  
                        !  5. b_1 <- b_1 - L_21^T*b_2
                        !  DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
                        !  Y <- ALPHA*A*X + BETA*Y
                
                          do ivar=1,Nrhs
                          do imag=1,icomp
                            call dgemv('T',Mblock(M+1)-Mblock(M),Mblock(M),           &
                               -1.d0,Awork(Mblock(M)+1,1),Ldwork,                     &
                                Bwork(icomp*((ivar-1)*Ldwork+Mblock(M))+imag),icomp,  &
                                1.d0,Bwork(icomp*((ivar-1)*Ldwork)+imag),icomp)
                          enddo
                          enddo
                  
                        !  6. b_1 <- L_11^{-T}b_1
                        !  DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
                        !  X <- A^{-1}*X
                  
                          do ivar=1,Nrhs
                          do imag=1,icomp
                            call dtrsv('L','T','N',Mblock(M),   &
                                      Awork(1,1),Ldwork,        &
                                      Bwork(icomp*((ivar-1)*Ldwork)+imag),icomp)
                          enddo
                          enddo
                  
                endif
            

end subroutine pbisolver3