!----------------------------------------------------------------------------
!> Purpose : - evaluate legendre polynomial base 
!!
!! @param[in]  T        - local edge coordinate [0,1]
!! @param[in]  LegOrder - order of approximation of polynomial base

!! @param[out] Npoly    - number of legendre polynomials (dimension base)
!! @param[out] Legendre - values of legendre polynomials at coordinate T
!
!----------------------------------------------------------------------
!   latest revision    - Nov2012
!
!----------------------------------------------------------------------
subroutine shap1Q_legendre(T,LegOrder,Npoly,Legendre)

  use parameters, only : MAXP

  implicit none
  !-------------------------------------------------------------------
  ! subroutine  arguments
  double precision, intent(in) :: T
  integer,intent(in) :: LegOrder 
  integer, intent(out) :: Npoly 
  double precision,dimension(MAXP),intent(inout) :: Legendre
  !-------------------------------------------------------------------
  ! local variables
  logical, save :: initialized = .FALSE.
  double precision :: x
  double precision, dimension(2,3:MAXP), save :: d
  integer :: k, iprint
  logical :: normalization = .TRUE.
  double precision, dimension(MAXP), save :: norm_constant 
 
  !-------------------------------------------------------------------

  iprint=0
  
  if (LegOrder < 0) then 
     Npoly=0; Legendre=0.d0
     return
  endif
  if (LegOrder > MAXP) then
     write(*,7001) LegOrder,MAXP
7001 format('shape1Q_legendre: LegOrder,MAXP = ',2i2)
     stop
  endif


!  Legendre=0.d0
  Npoly = LegOrder+1

  if (.not.initialized) then
     ! initialization of coefficients for 3-term recurrence relation
     initialized = .TRUE.
     do k=3,MAXP
        d(1,k) = real(2*k-3,kind=kind(d))/real(k-1,kind=kind(d))
        d(2,k) = real(k-2,kind=kind(d))/real(k-1,kind=kind(d))
     enddo
     ! normalization constant to make an orthonormal base (T in [0,1])
     if (normalization) then
        norm_constant = (/( 2*(k-1)+1 , k=1,MAXP) /)
        norm_constant = sqrt(norm_constant)
     endif


     if (iprint == 1) then
        write(*,*) 'shape1Q_legendre: d(1:2,3:MAXP) ='
        do k=3,MAXP
           write(*,*) 'd(1:2,',k,')=',d(1:2,k)
        enddo
     endif
  endif
  

! evaluation of legendre polynomials

  x = 2.d0*T-1.d0 ! change from T in [0,1] to x in [-1,1]

  if (LegOrder == 0) then 
     Legendre(1)=1.d0
  endif
  if (LegOrder == 1) then 
     Legendre(1)=1.d0
     Legendre(2)=x
  endif
  
  if (LegOrder > 1) then 
     Legendre(1)=1.d0
     Legendre(2)=x

     do k=3,Npoly
        Legendre(k) = d(1,k)*x*Legendre(k-1) - d(2,k)*Legendre(k-2)
     enddo
  endif

  if (iprint == 1) then
     write(*,*) 'shape1Q_legendre: T=',T
     write(*,*) 'shape1Q_legendre: Legendre(1,',Npoly,') ='
     do k=1,Npoly
        write(*,*) 'Legendre(',k,')=',Legendre(k)
     enddo
     write(*,*)
  endif

  ! normalization of the base 
  if (normalization) then
     Legendre(1:Npoly)=Legendre(1:Npoly)*norm_constant(1:Npoly)
  endif

end subroutine shap1Q_legendre



