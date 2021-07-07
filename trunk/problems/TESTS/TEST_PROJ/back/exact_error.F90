!----------------------------------------------------------------------
!                                                                     
!     routine name      - exact_error
!                                                                     
!----------------------------------------------------------------------
!                                                                     
!     latest revision:  - Nov 10
!                                                                     
!     purpose:          - routine computes the approximation error
!                         for a problem with known exact solution
!     arguments:                                                     
!                                                                     
!     out:     
!            Err      - the error
!            Rnorm      - norm of the exact solution
!
!-----------------------------------------------------------------------
!
subroutine exact_error(ErrorH,RnormH,ErrorE,RnormE,       &
                       ErrorV,RnormV,ErrorQ,RnormQ,Rate)
!
      use data_structure3D
      use parameters
!
      implicit none
      real*8,dimension(MAXEQNH),intent(out) :: ErrorH,RnormH
      real*8,dimension(MAXEQNE),intent(out) :: ErrorE,RnormE
      real*8,dimension(MAXEQNV),intent(out) :: ErrorV,RnormV
      real*8,dimension(MAXEQNQ),intent(out) :: ErrorQ,RnormQ
      real*8,                   intent(out) :: Rate
!
      real*8,dimension(MAXEQNH) :: derrH
      real*8,dimension(MAXEQNE) :: derrE
      real*8,dimension(MAXEQNV) :: derrV
      real*8,dimension(MAXEQNQ) :: derrQ
      real*8,dimension(MAXEQNH) :: dnorH
      real*8,dimension(MAXEQNE) :: dnorE
      real*8,dimension(MAXEQNV) :: dnorV
      real*8,dimension(MAXEQNQ) :: dnorQ
!      
      integer :: iprint,mdle,iel
      integer, save :: ivis = 0
      integer, save :: nrdofsh_save
      real*8 , save ::   error_save
!-----------------------------------------------------------------------
!
      iprint=0
!
!  ...intialize
      ErrorH(1:MAXEQNH)=0.d0 ; RnormH(1:MAXEQNH)=0.d0
      ErrorE(1:MAXEQNE)=0.d0 ; RnormE(1:MAXEQNE)=0.d0
      ErrorV(1:MAXEQNV)=0.d0 ; RnormV(1:MAXEQNV)=0.d0
      ErrorQ(1:MAXEQNQ)=0.d0 ; RnormQ(1:MAXEQNQ)=0.d0
!
!  ...loop over active elements      
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
!
!  .....compute error and accumulate
        call exact_error_h1(   mdle, derrH(1:MAXEQNH),dnorH(1:MAXEQNH))
        call exact_error_hcurl(mdle, derrE(1:MAXEQNE),dnorE(1:MAXEQNE))
        call exact_error_hdiv( mdle, derrV(1:MAXEQNV),dnorV(1:MAXEQNV))
        call exact_error_l2(   mdle, derrQ(1:MAXEQNQ),dnorQ(1:MAXEQNQ))
        ErrorH=ErrorH+derrH ; RnormH=RnormH+dnorH
        ErrorE=ErrorE+derrE ; RnormE=RnormE+dnorE
        ErrorV=ErrorV+derrV ; RnormV=RnormV+dnorV
        ErrorQ=ErrorQ+derrQ ; RnormQ=RnormQ+dnorQ
!
!!!!  .....printing        
!!!        if (iprint.eq.1) then
!!!          rel_err=0.d0
!!!          do n=1,neq       
!!!!
!!!!  .........compute relative error if significant          
!!!            if (dnorm(n).gt.0.d0) rel_err=derr(n)/dnorm(n)
!!!!            
!!!            write(*,7001) n,mdle,derr(n),dnorm(n),rel_err
!!! 7001       format('exact_error: n,mdle,derr,dnorm,derr/dnorm = ',i4,i6,3e12.5)
!!!          enddo
!!!        endif
!        
!  ...end of loop over active elements        
      enddo
!      
!  ...the much neglected square root!
      ErrorH(1:MAXEQNH)=sqrt(ErrorH(1:MAXEQNH)) ; RnormH(1:MAXEQNH)=sqrt(RnormH(1:MAXEQNH))
      ErrorE(1:MAXEQNE)=sqrt(ErrorE(1:MAXEQNE)) ; RnormE(1:MAXEQNE)=sqrt(RnormE(1:MAXEQNE))
      ErrorV(1:MAXEQNV)=sqrt(ErrorV(1:MAXEQNV)) ; RnormV(1:MAXEQNV)=sqrt(RnormV(1:MAXEQNV))
      ErrorQ(1:MAXEQNQ)=sqrt(ErrorQ(1:MAXEQNQ)) ; RnormQ(1:MAXEQNQ)=sqrt(RnormQ(1:MAXEQNQ))
!      
!!!!  ...if not 1st visit, compute rate
!!!      write(*,*)'NRDOFSH old,new = ',nrdofsh_save,NRDOFSH
!!!      Rate=0.d0
!!!      if (ivis.ne.0) then
!!!        if (Err(1).ne.0.d0) then
!!!          if (NRDOFSH.gt.nrdofsh_save) then
!!!        Rate = (log(error_save/Err(1)))/log(float(nrdofsh_save)/  &
!!!                                            float(NRDOFSH     ) )
!!!      endif ; endif ; endif
!!!!
!!!!  ...save quantities      
!!!      error_save=Err(1) ; nrdofsh_save=NRDOFSH
!
!  ...raise visitation flag      
      ivis=1
      Rate=0.d0
!
!      
endsubroutine exact_error
