#include "typedefs.h"
! 
! MODULE FOR HYPERELASTICITY PROCEDURES AND PARAMETERS
! JUNE 2020
! 
module hyperelasticity

real*8, save :: LOAD_FACTOR_PREV =0.d0
real*8, save :: LOAD_FACTOR =1.d0
logical, parameter :: SAVE_LOAD_FACTOR = .false.
character(len=128), save :: FILE_MATERIALS
! 
!  Kronecker's delta
real*8, save, dimension(3,3) :: DEL
real*8, save, dimension(3,3,3,3) :: DELDEL,DELCROSS,DDSYMM
!  Levi-Civita symbol
real*8, save, dimension(3,3,3) :: LEVI
real*8, save, dimension(3,3,3,3,3,3) :: LEVI2

integer,save :: NR_MATERIALS

TYPE hypermaterial
!  material label
   character(len=6) :: LABEL
!  density
   real*8  :: DENSITY
!  type of stored energy function W (i.e., what are their arguments)
   integer :: WTYPE
!  constitutive law identifier
   integer :: CONSTIT
!  number of parameters for this material's constitutive law
   integer :: NR_PARAMS
!  parameters for this material (must be allocated with length NR_PARAMS)
   real*8,allocatable :: PARAMS(:)
!  set this flag true if material not isotropic (i.e., anisotropic)
   logical :: FLAG_ANISO
!  set this flag true if material not homogeneous (i.e., heterogeneous)
   logical :: FLAG_HETER
!  set this flag true if material not compressible (i.e., incompressible)
   logical :: FLAG_INCOM
ENDTYPE hypermaterial
!  maximum number of parameters (limit for NR_PARAMS)
integer, parameter :: MAX_NR_P = 13
! 
! POSSIBLE VALUES FOR hypermaterial%WTYPE:
!  if W depends on F directly
integer, parameter :: W_F = 0
! if W depends on C=F^T*F
integer, parameter :: W_C = 1
! if W depends on the invariants of C
integer, parameter :: W_INVAR = 2
! if W depends on the principal stretches
integer, parameter :: W_STRETCH = 3
! if W depends on the MODIFIED invariants of C
integer, parameter :: W_INVAR_DEV = 4
! if W depends on the MODIFIED principal stretches
integer, parameter :: W_STRETCH_DEV = 5
! if W is the polyconvex function depending on { F , cof(F) , det(F) }
integer, parameter :: W_POLYCVX = 6
!
! POSSIBLE VALUES FOR hypermaterial%CONSTIT:
integer, parameter :: LINEAR = 0
integer, parameter :: MOONEY = 1
integer, parameter :: GENT = 2
integer, parameter :: OGDEN = 3
! 
type(hypermaterial), allocatable, save :: MATERIALS(:)
! 
integer, allocatable, save :: GMP_MAT(:)
! 
contains
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! functions to evaluate energy for different types of W and different constitutive laws
! 







subroutine eval_W_F(Imat,X,F,W,dWdF,d2WdF)
implicit none
integer,intent(in)::Imat
real*8,intent(in) ::X(3),F(3,3)
real*8,intent(out)::W,dWdF(3,3),d2WdF(3,3,3,3)

integer:: np,i,j,k,l
real*8 :: p(MAX_NR_P),gradu(3,3)

np = MATERIALS(Imat)%NR_PARAMS

p = 0.d0
p(1:np) = MATERIALS(Imat)%PARAMS

select case(MATERIALS(Imat)%CONSTIT)

case(LINEAR)

   if (MATERIALS(Imat)%FLAG_ANISO) then
      write(*,*) 'eval_W_F: no anisotropic constitutive laws have been implemented'
      stop
   endif
   gradu = F - DEL

!  p(1) is mu, p(2) is lambda
!  define the elastictiy tensor
   do l=1,3; do k=1,3; do j=1,3; do i=1,3
      d2WdF(i,j,k,l) = p(2) * DELDEL(i,j,k,l) + 2.d0 * p(1) * DDSYMM(i,j,k,l)
   enddo; enddo; enddo; enddo

   dWdF = contraction_4_2_majorsymm(d2WdF,gradu)

   W = 0.5d0*double_dot_prod(dWdF,gradu)

case default
   write(*,*) 'eval_W_F: this form of W does not support constitutive law no.', &
              MATERIALS(Imat)%CONSTIT
end select


end subroutine







subroutine eval_W_C(Imat,X,F,W,dWdC,d2WdC)
implicit none
integer,intent(in)::Imat
real*8,intent(in) ::X(3),F(3,3)
real*8,intent(out)::W,dWdC(3,3),d2WdC(3,3,3,3)

select case(MATERIALS(Imat)%CONSTIT)

case default
   write(*,*) 'eval_W_C: this form of W does not support constitutive law no.', &
              MATERIALS(Imat)%CONSTIT
end select

end subroutine








subroutine eval_W_invar(Imat,X,FI1,FI2,FI3,W,dW,d2W)
implicit none
integer,intent(in)::Imat
real*8,intent(in) ::X(3),FI1,FI2,FI3
real*8,intent(out)::W,dW(3),d2W(3,3)

integer :: np
real*8  :: p(MAX_NR_P),fJ,wJ,dwdJ,d2wdJ,djdi3,d2jdi3


np = MATERIALS(Imat)%NR_PARAMS

p = 0.d0
p(1:np) = MATERIALS(Imat)%PARAMS

! compute the volumetric part of the energy using J=sqrt(FI3) as argument
fJ = sqrt(FI3)
djdi3 = 0.5d0/fJ
d2jdi3 = -0.25d0/fJ**3
call eval_W_vol_term(fJ,wJ,dwdJ,d2wdJ)


select case (MATERIALS(Imat)%CONSTIT)

case(MOONEY)
   W = p(1)*(FI1-3.d0)+p(2)*(FI2-3.d0) + p(3)*wJ
   dW(1) = p(1)
   dW(2) = p(2)
   dW(3) = p(3)*dwdJ*djdi3
   d2W = 0.d0
   d2W(3,3) = p(3)*( dWdJ*d2jdi3 + d2WdJ*djdi3**2 )

case(GENT)
   W = -p(1)/6.d0 * (p(2)-3.d0) * log( 1.d0 - (FI1-3.d0)/(p(2)-3.d0) ) + p(3)*wJ
   dW(1) = p(1)/6.d0 * (p(2)-3.d0) / ( 1.d0 - (FI1-3.d0)/(p(2)-3.d0) ) / (p(2)-3.d0)
   dW(2) = 0.d0
   dW(3) = p(3)*dwdJ*djdi3
   d2W = 0.d0
   d2W(1,1) = p(1)/6.d0 * (p(2)-3.d0) / ( 1.d0 - (FI1-3.d0)/(p(2)-3.d0) )**2 / (p(2)-3.d0)**2
   d2W(3,3) = p(3)*( dWdJ*d2jdi3 + d2WdJ*djdi3**2 )

case default

   write(*,*) 'eval_W_invar: this form of W does not support constitutive law no.', &
              MATERIALS(Imat)%CONSTIT
   
   stop

end select

end subroutine






subroutine eval_W_invar_dev(Imat,X,FI1,FI2,FI3,W,dW,d2W)
implicit none
integer,intent(in)::Imat
real*8,intent(in) ::X(3),FI1,FI2,FI3
real*8,intent(out)::W,dW(3),d2W(3,3)

integer :: np , i, j, k , l
real*8  :: p(MAX_NR_P),fJ,wJ,dwdJ,d2wdJ,djdi3,d2jdi3,                        &
           fI1dev,fI2dev,dWdev(3),d2Wdev(3,3),dIdev(3,3),d2I1dev(3,3),d2I2dev(3,3)


np = MATERIALS(Imat)%NR_PARAMS

p = 0.d0
p(1:np) = MATERIALS(Imat)%PARAMS

! compute the volumetric part of the energy using J=sqrt(FI3) as argument
fJ = sqrt(FI3)
djdi3 = 0.5d0/fJ
d2jdi3 = -0.25d0/fJ**3
call eval_W_vol_term(fJ,wJ,dwdJ,d2wdJ)

! definition of I_i^dev; i=1,2. case i=3 equals original I_3
fI1dev = FI1*FI3**(-0.3333333333333333d0)
fI2dev = FI2*FI3**(-0.6666666666666667d0)

! dIdev is a matrix with (i,j)-th entry = \partial(I_i^dev) / \partial(I_j)
dIdev = 0.d0
! the nonzero entries are
dIdev(1,1) = FI3**(-0.3333333333333333d0)
dIdev(2,2) = FI3**(-0.6666666666666667d0)
dIdev(1,3) = -0.3333333333333333d0*FI3**(-1.333333333333333d0)*FI1
dIdev(2,3) = -0.6666666666666667d0*FI3**(-1.666666666666667d0)*FI2
dIdev(3,3) = 1.d0

! second derivatives of I_1^dev and I_2^dev
d2I1dev = 0.d0
d2I1dev(3,3) =  0.4444444444444444d0*FI3**(-2.333333333333333d0)*FI1
d2I1dev(1,3) = -0.3333333333333333d0*FI3**(-1.333333333333333d0)
d2I1dev(3,1) = d2I1dev(1,3)

d2I2dev = 0.d0
d2I2dev(3,3) =  1.111111111111111d0 *FI3**(-2.666666666666667d0)*FI2
d2I2dev(2,3) = -0.6666666666666667d0*FI3**(-1.666666666666667d0)
d2I2dev(3,2) = d2I2dev(2,3)

select case (MATERIALS(Imat)%CONSTIT)

case(MOONEY)
   W = p(1)*(FI1dev-3.d0)+p(2)*(FI2dev-3.d0) + p(3)*wJ
   dWdev(1) = p(1)
   dWdev(2) = p(2)
   dWdev(3) = p(3)*dwdJ*djdi3
   d2Wdev = 0.d0
   d2Wdev(3,3) = p(3)*( dWdJ*d2jdi3 + d2WdJ*djdi3**2 )

case(GENT)
   W = -p(1)/6.d0 * (p(2)-3.d0) * log( 1.d0 - (FI1dev-3.d0)/(p(2)-3.d0) ) + p(3)*wJ
   dWdev(1) = p(1)/6.d0 * (p(2)-3.d0) / ( 1.d0 - (FI1dev-3.d0)/(p(2)-3.d0) ) / (p(2)-3.d0)
   dWdev(2) = 0.d0
   dWdev(3) = p(3)*dwdJ*djdi3
   d2Wdev = 0.d0
   d2Wdev(1,1) = p(1)/6.d0 * (p(2)-3.d0) / ( 1.d0 - (FI1dev-3.d0)/(p(2)-3.d0) )**2 / (p(2)-3.d0)**2
   d2Wdev(3,3) = p(3)*( dWdJ*d2jdi3 + d2WdJ*djdi3**2 )

case default

   write(*,*) 'eval_W_invar_dev: this form of W does not support constitutive law no.', &
              MATERIALS(Imat)%CONSTIT
   
   stop

end select

! Now we multiply by the interior derivative
dW(1) = dot_product(dWdev,dIdev(:,1))
dW(2) = dot_product(dWdev,dIdev(:,2))
dW(3) = dot_product(dWdev,dIdev(:,3))

! Similarly for the second derivative
d2W = 0.d0
! term 1=  dIdev' * d2Wdev * dIdev 
do l=1,3; do j=1,3; do k=1,3; do i=1,3
   d2W(k,l) = d2W(k,l) + d2Wdev(i,j)*dIdev(i,k)*dIdev(j,l) 
enddo; enddo; enddo; enddo
! term 2= \sum_a dWdev(a)*d2I(a)dev. (be aware d2I3dev = 0)
do l=1,3; do k=1,3
   d2W(k,l) = d2W(k,l) + dWdev(1)*d2I1dev(k,l) + dWdev(2)*d2I2dev(k,l)
enddo; enddo

end subroutine





subroutine eval_W_stretch(Imat,X,S1,S2,S3,W,dW,d2W)
implicit none
integer,intent(in)::Imat
real*8,intent(in) ::X(3),S1,S2,S3
real*8,intent(out)::W,dW(3),d2W(3,3)

integer :: np,n
real*8  :: p(MAX_NR_P),fJ,wJ,dwds(3),d2wds(3,3),djds(3),d2jds(3,3), &
           dwdj,d2wdj,fi1,fi2,fi3,a,b,ca,cb


np = MATERIALS(Imat)%NR_PARAMS

p = 0.d0
p(1:np) = MATERIALS(Imat)%PARAMS

! compute the volumetric part of the energy using J=sqrt(FI3) as argument
fJ = S1*S2*S3
djds(1) = S2*S3
djds(2) = S3*S1
djds(3) = S1*S2
d2jds = 0.d0
d2jds(2,1) = S3
d2jds(3,1) = S2
d2jds(1,2) = S3
d2jds(3,2) = S1
d2jds(1,3) = S2
d2jds(2,3) = S1

call eval_W_vol_term(fJ,wJ,dwdJ,d2wdJ)
dwds(:) = dwdJ*djds(:)
d2wds = dwdJ * d2jds + d2WdJ * outer_product(djds,djds)

fi1 = S1**2 + S2**2 + S3**2
fi2 = S1**2*S2**2 + S1**2*S3**2 + S2**2*S3**2


select case (MATERIALS(Imat)%CONSTIT)

case(MOONEY)
   W = p(1)*(FI1-3.d0)+p(2)*(FI2-3.d0) + p(3)*wJ
   dW(1) = 2.d0*p(1)*S1 + 2.d0*p(2)*S1*(S2**2+S3**2) + p(3)*dwds(1)
   dW(2) = 2.d0*p(1)*S2 + 2.d0*p(2)*S2*(S1**2+S3**2) + p(3)*dwds(2)
   dW(3) = 2.d0*p(1)*S3 + 2.d0*p(2)*S3*(S1**2+S2**2) + p(3)*dwds(3)
   d2W = 0.d0
   d2W(1,1) = 2.d0*p(1) + 2.d0*p(2)*(S2**2+S3**2) + p(3)*d2wds(1,1)
   d2W(2,1) = 4.d0*p(2)*S1*S2 + p(3)*d2wds(2,1)
   d2W(3,1) = 4.d0*p(2)*S1*S2 + p(3)*d2wds(3,1)
   d2W(2,2) = 2.d0*p(1) + 2.d0*p(2)*(S1**2+S3**2) + p(3)*d2wds(2,2)
   d2W(3,2) = 4.d0*p(2)*S3*S2 + p(3)*d2wds(3,2)
   d2W(3,3) = 2.d0*p(1) + 2.d0*p(2)*(S1**2+S2**2) + p(3)*d2wds(3,3)
   d2W(2,3) = d2W(3,2)
   d2W(1,3) = d2W(3,1)
   d2W(1,2) = d2W(2,1)

case(OGDEN)
   W = p(np)*wJ
   dW(:) = p(np)*dwds(:)
   d2W(:,:) = p(np)*d2wds(:,:)
   do n=0,np-5,4
      a = p(n+1)
      ca= p(n+2)
      b = p(n+3)
      cb= p(n+4)
      W = W + (      S1**a +      S2**a +      S3**a - 3.d0 )*ca/a   &
            + ( (S1*S2)**b + (S2*S3)**b + (S1*S3)**b - 3.d0 )*cb/b
      dW(1)   = dW(1)   + S1**(a-1.d0)*ca                            &
                        + S1**(b-1.d0)*( S2**b + S3**b )*cb
      dW(2)   = dW(2)   + S2**(a-1.d0)*ca                            &
                        + S2**(b-1.d0)*( S1**b + S3**b )*cb
      dW(3)   = dW(3)   + S3**(a-1.d0)*ca                            &
                        + S3**(b-1.d0)*( S2**b + S1**b )*cb
      d2W(1,1)= d2W(1,1)+ S1**(a-2.d0)*ca*(a-1.d0)                   &
                        + S1**(b-2.d0)*( S2**b + S3**b )*cb*(b-1.d0)
      d2W(2,1)= d2W(2,1)+ S1**(b-1.d0)*S2**(b-1.d0)*cb*b
      d2W(3,1)= d2W(3,1)+ S1**(b-1.d0)*S3**(b-1.d0)*cb*b
      d2W(2,2)= d2W(2,2)+ S2**(a-2.d0)*ca*(a-1.d0)                   &
                        + S2**(b-2.d0)*( S1**b + S3**b )*cb*(b-1.d0)
      d2W(3,2)= d2W(3,2)+ S2**(b-1.d0)*S3**(b-1.d0)*cb*b
      d2W(3,3)= d2W(3,3)+ S3**(a-2.d0)*ca*(a-1.d0)                   &
                        + S3**(b-2.d0)*( S1**b + S2**b )*cb*(b-1.d0)
   enddo
   d2W(2,3) = d2W(3,2)
   d2W(1,3) = d2W(3,1)
   d2W(1,2) = d2W(2,1)      


case default

   write(*,*) 'eval_W_stretch: this form of W does not support constitutive law no.', &
              MATERIALS(Imat)%CONSTIT
   
   stop

end select


end subroutine










subroutine eval_W_stretch_sq(Imat,X,S1,S2,S3,W,dW,d2W)
implicit none
integer,intent(in)::Imat
real*8,intent(in) ::X(3),S1,S2,S3
real*8,intent(out)::W,dW(3),d2W(3,3)

integer :: np,n
real*8  :: p(MAX_NR_P),fJ,wJ,dwds(3),d2wds(3,3),djds(3),d2jds(3,3), &
           dwdj,d2wdj,fi1,fi2,fi3,a,b,ca,cb


np = MATERIALS(Imat)%NR_PARAMS

p = 0.d0
p(1:np) = MATERIALS(Imat)%PARAMS

! compute the volumetric part of the energy using J=sqrt(FI3) as argument
fJ = sqrt(S1*S2*S3)
djds(1) = 0.5d0*fJ/S1
djds(2) = 0.5d0*fJ/S2
djds(3) = 0.5d0*fJ/S3
d2jds = 0.d0
d2jds(1,1) = -0.25d0*fJ/S1**2.d0
d2jds(2,1) = 0.25d0*fJ/(S1*S2)
d2jds(3,1) = 0.25d0*fJ/(S1*S3)
d2jds(1,2) = d2jds(2,1)
d2jds(2,2) = -0.25d0*fJ/S2**2.d0
d2jds(3,2) = 0.25d0*fJ/(S3*S2)
d2jds(1,3) = d2jds(3,1)
d2jds(2,3) = d2jds(3,2)
d2jds(3,3) = -0.25d0*fJ/S3**2.d0

call eval_W_vol_term(fJ,wJ,dwdJ,d2wdJ)
dwds(:) = dwdJ*djds(:)
d2wds = dwdJ * d2jds + d2WdJ * outer_product(djds,djds)

fi1 = S1 + S2 + S3
fi2 = S1*S2 + S1*S3 + S2*S3


select case (MATERIALS(Imat)%CONSTIT)

case(MOONEY)
   W = p(1)*(FI1-3.d0)+p(2)*(FI2-3.d0) + p(3)*wJ
   dW(1) = p(1) + p(2)*(S2+S3) + p(3)*dwds(1)
   dW(2) = p(1) + p(2)*(S1+S3) + p(3)*dwds(2)
   dW(3) = p(1) + p(2)*(S1+S2) + p(3)*dwds(3)
   d2W = 0.d0
   d2W(1,1) = p(3)*d2wds(1,1)
   d2W(2,1) = p(2) + p(3)*d2wds(2,1)
   d2W(3,1) = p(2) + p(3)*d2wds(3,1)
   d2W(2,2) = p(3)*d2wds(2,2)
   d2W(3,2) = p(2) + p(3)*d2wds(3,2)
   d2W(3,3) = p(3)*d2wds(3,3)
   d2W(2,3) = d2W(3,2)
   d2W(1,3) = d2W(3,1)
   d2W(1,2) = d2W(2,1)

case(OGDEN)
   W = p(np)*wJ
   dW(:) = p(np)*dwds(:)
   d2W(:,:) = p(np)*d2wds(:,:)
   do n=0,np-5,4
      a = p(n+1)*0.5d0
      ca= p(n+2)*0.5d0
      b = p(n+3)*0.5d0
      cb= p(n+4)*0.5d0
      W = W + (      S1**a +      S2**a +      S3**a - 3.d0 )*ca/a   &
            + ( (S1*S2)**b + (S2*S3)**b + (S1*S3)**b - 3.d0 )*cb/b
      dW(1)   = dW(1)   + S1**(a-1.d0)*ca                            &
                        + S1**(b-1.d0)*( S2**b + S3**b )*cb
      dW(2)   = dW(2)   + S2**(a-1.d0)*ca                            &
                        + S2**(b-1.d0)*( S1**b + S3**b )*cb
      dW(3)   = dW(3)   + S3**(a-1.d0)*ca                            &
                        + S3**(b-1.d0)*( S2**b + S1**b )*cb
      d2W(1,1)= d2W(1,1)+ S1**(a-2.d0)*ca*(a-1.d0)                   &
                        + S1**(b-2.d0)*( S2**b + S3**b )*cb*(b-1.d0)
      d2W(2,1)= d2W(2,1)+ S1**(b-1.d0)*S2**(b-1.d0)*cb*b
      d2W(3,1)= d2W(3,1)+ S1**(b-1.d0)*S3**(b-1.d0)*cb*b
      d2W(2,2)= d2W(2,2)+ S2**(a-2.d0)*ca*(a-1.d0)                   &
                        + S2**(b-2.d0)*( S1**b + S3**b )*cb*(b-1.d0)
      d2W(3,2)= d2W(3,2)+ S2**(b-1.d0)*S3**(b-1.d0)*cb*b
      d2W(3,3)= d2W(3,3)+ S3**(a-2.d0)*ca*(a-1.d0)                   &
                        + S3**(b-2.d0)*( S1**b + S2**b )*cb*(b-1.d0)
   enddo
   d2W(2,3) = d2W(3,2)
   d2W(1,3) = d2W(3,1)
   d2W(1,2) = d2W(2,1)      


case default

   write(*,*) 'eval_W_stretch_sq: this form of W does not support constitutive law no.', &
              MATERIALS(Imat)%CONSTIT
   
   stop

end select


end subroutine










subroutine eval_W_stretch_dev(Imat,X,S1,S2,S3,W,dW,d2W)
implicit none
integer,intent(in)::Imat
real*8,intent(in) ::X(3),S1,S2,S3
real*8,intent(out)::W,dW(3),d2W(3,3)

integer :: np,n, i,j,k,l
real*8  :: p(MAX_NR_P),fJ,wJ,dwds(3),d2wds(3,3),djds(3),d2jds(3,3),            &
           dwdj,d2wdj,fi1dev,fi2dev,a,b,ca,cb,dWdev(3),d2Wdev(3,3),            &
           s1dev,s2dev,s3dev,dsdev(3,3),d2s1dev(3,3),d2s2dev(3,3),d2s3dev(3,3)

np = MATERIALS(Imat)%NR_PARAMS

p = 0.d0
p(1:np) = MATERIALS(Imat)%PARAMS

fJ = S1*S2*S3
djds(1) = S2*S3
djds(2) = S3*S1
djds(3) = S1*S2
d2jds = 0.d0
d2jds(2,1) = S3
d2jds(3,1) = S2
d2jds(1,2) = S3
d2jds(3,2) = S1
d2jds(1,3) = S2
d2jds(2,3) = S1

s1dev = S1*fJ**(-0.3333333333333333d0)
s2dev = S2*fJ**(-0.3333333333333333d0)
s3dev = S3*fJ**(-0.3333333333333333d0)

dsdev(1,1) = -S1*0.3333333333333333d0*fJ**(-1.3333333333333333d0)*djds(1) + fJ**(-0.3333333333333333d0)
dsdev(2,1) = -S2*0.3333333333333333d0*fJ**(-1.3333333333333333d0)*djds(1)
dsdev(3,1) = -S3*0.3333333333333333d0*fJ**(-1.3333333333333333d0)*djds(1)
dsdev(1,2) = -S1*0.3333333333333333d0*fJ**(-1.3333333333333333d0)*djds(2)
dsdev(2,2) = -S2*0.3333333333333333d0*fJ**(-1.3333333333333333d0)*djds(2) + fJ**(-0.3333333333333333d0)
dsdev(3,2) = -S3*0.3333333333333333d0*fJ**(-1.3333333333333333d0)*djds(2)
dsdev(1,3) = -S1*0.3333333333333333d0*fJ**(-1.3333333333333333d0)*djds(3)
dsdev(2,3) = -S2*0.3333333333333333d0*fJ**(-1.3333333333333333d0)*djds(3)
dsdev(3,3) = -S3*0.3333333333333333d0*fJ**(-1.3333333333333333d0)*djds(3) + fJ**(-0.3333333333333333d0)

d2s1dev = 0.d0
do k=1,3; do j=1,3
   d2s1dev(j,k) = d2s1dev(j,k) - 0.3333333333333333d0*fJ**(-1.3333333333333333d0)      &
                               * ( DEL(j,1)*djds(k)+DEL(k,1)*djds(j) + S1*d2jds(j,k)   &
                                  -1.3333333333333333d0/fJ*S1*djds(j)*djds(k)       )
enddo;enddo

d2s2dev = 0.d0
do k=1,3; do j=1,3
   d2s2dev(j,k) = d2s2dev(j,k) - 0.3333333333333333d0*fJ**(-1.3333333333333333d0)      &
                               * ( DEL(j,2)*djds(k)+DEL(k,2)*djds(j) + S2*d2jds(j,k)   &
                                  -1.3333333333333333d0/fJ*S2*djds(j)*djds(k)       )
enddo;enddo

d2s3dev = 0.d0
do k=1,3; do j=1,3
   d2s3dev(j,k) = d2s3dev(j,k) - 0.3333333333333333d0*fJ**(-1.3333333333333333d0)      &
                               * ( DEL(j,3)*djds(k)+DEL(k,3)*djds(j) + S3*d2jds(j,k)   &
                                  -1.3333333333333333d0/fJ*S3*djds(j)*djds(k)       )
enddo;enddo

fi1dev = S1dev**2 + S2dev**2 + S3dev**2
fi2dev = S1dev**2*S2dev**2 + S1dev**2*S3dev**2 + S2dev**2*S3dev**2


select case (MATERIALS(Imat)%CONSTIT)

case(MOONEY)
   W = p(1)*(FI1dev-3.d0)+p(2)*(FI2dev-3.d0)
   dWdev(1) = 2.d0*p(1)*S1dev + 2.d0*p(2)*S1dev*(S2dev**2+S3dev**2)
   dWdev(2) = 2.d0*p(1)*S2dev + 2.d0*p(2)*S2dev*(S1dev**2+S3dev**2)
   dWdev(3) = 2.d0*p(1)*S3dev + 2.d0*p(2)*S3dev*(S1dev**2+S2dev**2)
   d2Wdev = 0.d0
   d2Wdev(1,1) = 2.d0*p(1) + 2.d0*p(2)*(S2dev**2+S3dev**2)
   d2Wdev(2,1) = 4.d0*p(2)*S1dev*S2dev
   d2Wdev(3,1) = 4.d0*p(2)*S1dev*S2dev
   d2Wdev(2,2) = 2.d0*p(1) + 2.d0*p(2)*(S1dev**2+S3dev**2)
   d2Wdev(3,2) = 4.d0*p(2)*S3dev*S2dev
   d2Wdev(3,3) = 2.d0*p(1) + 2.d0*p(2)*(S1dev**2+S2dev**2)
   d2Wdev(2,3) = d2Wdev(3,2)
   d2Wdev(1,3) = d2Wdev(3,1)
   d2Wdev(1,2) = d2Wdev(2,1)

case(OGDEN)
   W = 0.d0
   dWdev(:) = 0.d0
   d2Wdev(:,:) = 0.d0
   do n=0,np-5,4
      a = p(n+1)
      ca= p(n+2)
      b = p(n+3)
      cb= p(n+4)
      W = W + (      S1dev**a    +      S2dev**a    +      S3dev**a - 3.d0    )*ca/a  &
            + ( (S1dev*S2dev)**b + (S2dev*S3dev)**b + (S1dev*S3dev)**b - 3.d0 )*cb/b
      dWdev(1)   = dWdev(1)   + S1dev**(a-1.d0)*ca                                    &
                        + S1dev**(b-1.d0)*( S2dev**b + S3dev**b )*cb
      dWdev(2)   = dWdev(2)   + S2dev**(a-1.d0)*ca                                    &
                        + S2dev**(b-1.d0)*( S1dev**b + S3dev**b )*cb
      dWdev(3)   = dWdev(3)   + S3dev**(a-1.d0)*ca                                    &
                        + S3dev**(b-1.d0)*( S2dev**b + S1dev**b )*cb
      d2Wdev(1,1)= d2Wdev(1,1)+ S1dev**(a-2.d0)*ca*(a-1.d0)                           &
                        + S1dev**(b-2.d0)*( S2dev**b + S3dev**b )*cb*(b-1.d0)
      d2Wdev(2,1)= d2Wdev(2,1)+ S1dev**(b-1.d0)*S2dev**(b-1.d0)*cb*b
      d2Wdev(3,1)= d2Wdev(3,1)+ S1dev**(b-1.d0)*S3dev**(b-1.d0)*cb*b
      d2Wdev(2,2)= d2Wdev(2,2)+ S2dev**(a-2.d0)*ca*(a-1.d0)                           &
                        + S2dev**(b-2.d0)*( S1dev**b + S3dev**b )*cb*(b-1.d0)
      d2Wdev(3,2)= d2Wdev(3,2)+ S2dev**(b-1.d0)*S3dev**(b-1.d0)*cb*b
      d2Wdev(3,3)= d2Wdev(3,3)+ S3dev**(a-2.d0)*ca*(a-1.d0)                           &
                        + S3dev**(b-2.d0)*( S1dev**b + S2dev**b )*cb*(b-1.d0)
   enddo
   d2Wdev(2,3) = d2Wdev(3,2)
   d2Wdev(1,3) = d2Wdev(3,1)
   d2Wdev(1,2) = d2Wdev(2,1)      

case default

   write(*,*) 'eval_W_stretch_dev: this form of W does not support constitutive law no.', &
              MATERIALS(Imat)%CONSTIT
   
   stop

end select

! Now we multiply by the interior derivative
dW(1) = dot_product(dWdev,dsdev(:,1))
dW(2) = dot_product(dWdev,dsdev(:,2))
dW(3) = dot_product(dWdev,dsdev(:,3))

! Similarly for the second derivative
d2W = 0.d0
! term 1=  dsdev' * d2Wdev * dsdev 
do l=1,3; do j=1,3; do k=1,3; do i=1,3
   d2W(k,l) = d2W(k,l) + d2Wdev(i,j)*dsdev(i,k)*dsdev(j,l) 
enddo; enddo; enddo; enddo
! term 2= \sum_a dWdev(a)*d2s(a)dev
do l=1,3; do k=1,3
   d2W(k,l) = d2W(k,l) &
            + dWdev(1)*d2s1dev(k,l) + dWdev(2)*d2s2dev(k,l) + dWdev(3)*d2s3dev(k,l)
enddo; enddo

! compute the volumetric part of the energy using J as argument

call eval_W_vol_term(fJ,wJ,dwdJ,d2wdJ)
dwds(:) = dwdJ*djds(:)
d2wds = dwdJ * d2jds + d2WdJ * outer_product(djds,djds)

W = W + p(np)*wJ
dW = dW + p(np)*dwds
d2W = d2W + p(np)*d2wds

end subroutine





! evaluating with respect to the squared stretches (eigenvalues of F^T*F)

subroutine eval_W_stretch_sq_dev(Imat,X,S1,S2,S3,W,dW,d2W)
implicit none
integer,intent(in)::Imat
real*8,intent(in) ::X(3),S1,S2,S3
real*8,intent(out)::W,dW(3),d2W(3,3)

integer :: np,n, i,j,k,l
real*8  :: p(MAX_NR_P),fI3,dI3ds(3),d2I3ds(3,3),fJ,djds(3),d2jds(3,3),            &
           wJ,dwds(3),d2wds(3,3),dwdj,d2wdj,a,b,ca,cb,dWdev(3),d2Wdev(3,3),            &
           fi1dev,fi2dev,s1dev,s2dev,s3dev,dsdev(3,3),d2s1dev(3,3),d2s2dev(3,3),d2s3dev(3,3)

np = MATERIALS(Imat)%NR_PARAMS

p = 0.d0
p(1:np) = MATERIALS(Imat)%PARAMS

fI3 = S1*S2*S3
dI3ds(1) = S2*S3
dI3ds(2) = S3*S1
dI3ds(3) = S1*S2
d2I3ds = 0.d0
d2I3ds(2,1) = S3
d2I3ds(3,1) = S2
d2I3ds(1,2) = S3
d2I3ds(3,2) = S1
d2I3ds(1,3) = S2
d2I3ds(2,3) = S1

fJ = sqrt(fI3)
djds(1) = 0.5d0*fJ/S1
djds(2) = 0.5d0*fJ/S2
djds(3) = 0.5d0*fJ/S3
d2jds = 0.d0
d2jds(1,1) = -0.25d0*fJ/S1**2
d2jds(2,1) = 0.25d0*fJ/(S1*S2)
d2jds(3,1) = 0.25d0*fJ/(S1*S3)
d2jds(1,2) = d2jds(2,1)
d2jds(2,2) = -0.25d0*fJ/S2**2
d2jds(3,2) = 0.25d0*fJ/(S3*S2)
d2jds(1,3) = d2jds(3,1)
d2jds(2,3) = d2jds(3,2)
d2jds(3,3) = -0.25d0*fJ/S3**2

s1dev = S1*fI3**(-0.3333333333333333d0)
s2dev = S2*fI3**(-0.3333333333333333d0)
s3dev = S3*fI3**(-0.3333333333333333d0)

dsdev(1,1) = -S1*0.3333333333333333d0*fI3**(-1.3333333333333333d0)*dI3ds(1) + fI3**(-0.3333333333333333d0)
dsdev(2,1) = -S2*0.3333333333333333d0*fI3**(-1.3333333333333333d0)*dI3ds(1)
dsdev(3,1) = -S3*0.3333333333333333d0*fI3**(-1.3333333333333333d0)*dI3ds(1)
dsdev(1,2) = -S1*0.3333333333333333d0*fI3**(-1.3333333333333333d0)*dI3ds(2)
dsdev(2,2) = -S2*0.3333333333333333d0*fI3**(-1.3333333333333333d0)*dI3ds(2) + fI3**(-0.3333333333333333d0)
dsdev(3,2) = -S3*0.3333333333333333d0*fI3**(-1.3333333333333333d0)*dI3ds(2)
dsdev(1,3) = -S1*0.3333333333333333d0*fI3**(-1.3333333333333333d0)*dI3ds(3)
dsdev(2,3) = -S2*0.3333333333333333d0*fI3**(-1.3333333333333333d0)*dI3ds(3)
dsdev(3,3) = -S3*0.3333333333333333d0*fI3**(-1.3333333333333333d0)*dI3ds(3) + fI3**(-0.3333333333333333d0)

d2s1dev = 0.d0
do k=1,3; do j=1,3
   d2s1dev(j,k) = d2s1dev(j,k) - 0.3333333333333333d0*fI3**(-1.3333333333333333d0)      &
                               * ( DEL(j,1)*dI3ds(k)+DEL(k,1)*dI3ds(j) + S1*d2I3ds(j,k)   &
                                  -1.3333333333333333d0/fI3*S1*dI3ds(j)*dI3ds(k)       )
enddo;enddo

d2s2dev = 0.d0
do k=1,3; do j=1,3
   d2s2dev(j,k) = d2s2dev(j,k) - 0.3333333333333333d0*fI3**(-1.3333333333333333d0)      &
                               * ( DEL(j,2)*dI3ds(k)+DEL(k,2)*dI3ds(j) + S2*d2I3ds(j,k)   &
                                  -1.3333333333333333d0/fI3*S2*dI3ds(j)*dI3ds(k)       )
enddo;enddo

d2s3dev = 0.d0
do k=1,3; do j=1,3
   d2s3dev(j,k) = d2s3dev(j,k) - 0.3333333333333333d0*fI3**(-1.3333333333333333d0)      &
                               * ( DEL(j,3)*dI3ds(k)+DEL(k,3)*dI3ds(j) + S3*d2I3ds(j,k)   &
                                  -1.3333333333333333d0/fI3*S3*dI3ds(j)*dI3ds(k)       )
enddo;enddo

fi1dev = S1dev + S2dev + S3dev
fi2dev = S1dev*S2dev + S1dev*S3dev + S2dev*S3dev


select case (MATERIALS(Imat)%CONSTIT)

case(MOONEY)
   W = p(1)*(FI1dev-3.d0)+p(2)*(FI2dev-3.d0)
   dWdev(1) = 2.d0*p(1)*S1dev + 2.d0*p(2)*S1dev*(S2dev**2+S3dev**2)
   dWdev(2) = 2.d0*p(1)*S2dev + 2.d0*p(2)*S2dev*(S1dev**2+S3dev**2)
   dWdev(3) = 2.d0*p(1)*S3dev + 2.d0*p(2)*S3dev*(S1dev**2+S2dev**2)
   d2Wdev = 0.d0
   d2Wdev(1,1) = 2.d0*p(1) + 2.d0*p(2)*(S2dev**2+S3dev**2)
   d2Wdev(2,1) = 4.d0*p(2)*S1dev*S2dev
   d2Wdev(3,1) = 4.d0*p(2)*S1dev*S2dev
   d2Wdev(2,2) = 2.d0*p(1) + 2.d0*p(2)*(S1dev**2+S3dev**2)
   d2Wdev(3,2) = 4.d0*p(2)*S3dev*S2dev
   d2Wdev(3,3) = 2.d0*p(1) + 2.d0*p(2)*(S1dev**2+S2dev**2)
   d2Wdev(2,3) = d2Wdev(3,2)
   d2Wdev(1,3) = d2Wdev(3,1)
   d2Wdev(1,2) = d2Wdev(2,1)

   W = p(1)*(FI1dev-3.d0)+p(2)*(FI2dev-3.d0)
   dWdev(1) = p(1) + p(2)*(S2dev+S3dev)
   dWdev(2) = p(1) + p(2)*(S1dev+S3dev)
   dWdev(3) = p(1) + p(2)*(S1dev+S2dev)
   d2Wdev = 0.d0
   d2Wdev(2,1) = p(2) 
   d2Wdev(3,1) = p(2) 
   d2Wdev(3,2) = p(2)
   d2Wdev(2,3) = d2Wdev(3,2)
   d2Wdev(1,3) = d2Wdev(3,1)
   d2Wdev(1,2) = d2Wdev(2,1)


case(OGDEN)
   W = 0.d0
   dWdev(:) = 0.d0
   d2Wdev(:,:) = 0.d0
   do n=0,np-5,4
      a = p(n+1)*0.5d0
      ca= p(n+2)*0.5d0
      b = p(n+3)*0.5d0
      cb= p(n+4)*0.5d0
      W = W + (      S1dev**a    +      S2dev**a    +      S3dev**a - 3.d0    )*ca/a  &
            + ( (S1dev*S2dev)**b + (S2dev*S3dev)**b + (S1dev*S3dev)**b - 3.d0 )*cb/b
      dWdev(1)   = dWdev(1)   + S1dev**(a-1.d0)*ca                                    &
                        + S1dev**(b-1.d0)*( S2dev**b + S3dev**b )*cb
      dWdev(2)   = dWdev(2)   + S2dev**(a-1.d0)*ca                                    &
                        + S2dev**(b-1.d0)*( S1dev**b + S3dev**b )*cb
      dWdev(3)   = dWdev(3)   + S3dev**(a-1.d0)*ca                                    &
                        + S3dev**(b-1.d0)*( S2dev**b + S1dev**b )*cb
      d2Wdev(1,1)= d2Wdev(1,1)+ S1dev**(a-2.d0)*ca*(a-1.d0)                           &
                        + S1dev**(b-2.d0)*( S2dev**b + S3dev**b )*cb*(b-1.d0)
      d2Wdev(2,1)= d2Wdev(2,1)+ S1dev**(b-1.d0)*S2dev**(b-1.d0)*cb*b
      d2Wdev(3,1)= d2Wdev(3,1)+ S1dev**(b-1.d0)*S3dev**(b-1.d0)*cb*b
      d2Wdev(2,2)= d2Wdev(2,2)+ S2dev**(a-2.d0)*ca*(a-1.d0)                           &
                        + S2dev**(b-2.d0)*( S1dev**b + S3dev**b )*cb*(b-1.d0)
      d2Wdev(3,2)= d2Wdev(3,2)+ S2dev**(b-1.d0)*S3dev**(b-1.d0)*cb*b
      d2Wdev(3,3)= d2Wdev(3,3)+ S3dev**(a-2.d0)*ca*(a-1.d0)                           &
                        + S3dev**(b-2.d0)*( S1dev**b + S2dev**b )*cb*(b-1.d0)
   enddo
   d2Wdev(2,3) = d2Wdev(3,2)
   d2Wdev(1,3) = d2Wdev(3,1)
   d2Wdev(1,2) = d2Wdev(2,1)      

case default

   write(*,*) 'eval_W_stretch_sq_dev: this form of W does not support constitutive law no.', &
              MATERIALS(Imat)%CONSTIT
   
   stop

end select

! Now we multiply by the interior derivative
dW(1) = dot_product(dWdev,dsdev(:,1))
dW(2) = dot_product(dWdev,dsdev(:,2))
dW(3) = dot_product(dWdev,dsdev(:,3))

! Similarly for the second derivative
d2W = 0.d0
! term 1=  dsdev' * d2Wdev * dsdev 
do l=1,3; do j=1,3; do k=1,3; do i=1,3
   d2W(k,l) = d2W(k,l) + d2Wdev(i,j)*dsdev(i,k)*dsdev(j,l) 
enddo; enddo; enddo; enddo
! term 2= \sum_a dWdev(a)*d2s(a)dev
do l=1,3; do k=1,3
   d2W(k,l) = d2W(k,l) &
            + dWdev(1)*d2s1dev(k,l) + dWdev(2)*d2s2dev(k,l) + dWdev(3)*d2s3dev(k,l)
enddo; enddo

! compute the volumetric part of the energy using J as argument

call eval_W_vol_term(fJ,wJ,dwdJ,d2wdJ)
dwds(:) = dwdJ*djds(:)
d2wds = dwdJ * d2jds + d2WdJ * outer_product(djds,djds)

W = W + p(np)*wJ
dW = dW + p(np)*dwds
d2W = d2W + p(np)*d2wds

end subroutine




! subroutine eval_W_polycvx(Imat,X,F,H,FJ,W,dWdFp,d2WdFp,dWdH,d2WdH,dWdJ,d2WdJ)
! implicit none
! integer,intent(in)::Imat
! real*8,intent(in) ::X(3),S1,S2,S3
! real*8,intent(out)::W,dW(3),d2W(3,3)
! 
! select case(MATERIALS(Imat)%CONSTIT)
! 
! case default
!    write(*,*) 'eval_W_polycvx: this form of W does not support constitutive law no.', &
!               MATERIALS(Imat)%CONSTIT
! end select
! end subroutine







subroutine eval_W_vol_term(FJ,WJ,dWdJ,d2WdJ)
implicit none
real*8,intent(in ) :: FJ
real*8,intent(out) :: WJ,dWdJ,d2WdJ

! WJ = 0.5d0 * ( FJ**2 - 1.d0 ) + ( log( FJ ) )**2
! dWdJ = FJ + 2.d0 * log( FJ ) / FJ
! d2WdJ = 1.d0 + 2.d0 * ( 1.d0-log(FJ) )/ FJ**2

! Simo & Hughes, pp 307
WJ = 0.5d0*( 0.5d0 * ( FJ**2 - 1.d0 ) - log(FJ) )
dWdJ = 0.5d0*( FJ - 1.d0 / FJ )
d2WdJ = 0.5d0*( 1.d0 + 1.d0 / (FJ**2) )

end subroutine







subroutine get_deviat_F(F,FJ,Fdev)
implicit none
real*8,intent(in ) :: F(3,3),FJ
real*8,intent(out) :: Fdev(3,3)

Fdev = F*FJ**(-1.d0/3.d0)

end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! auxiliary procedures for initialization
! 


subroutine init_tensors
implicit none
integer :: i,j,k,l,m,n

DEL = 0.d0
do i=1,3
   DEL(i,i) = 1.d0
enddo
DELDEL   = 0.d0
DELCROSS = 0.d0
! do l=1,3; do k=1,3; do j=1,3; do i=1,3
!    DELDEL  (i,j,k,l) =  DEL(i,j)*DEL(k,l)
!    DELCROSS(i,j,k,l) =  DEL(i,k)*DEL(j,l)
!    DDSYMM  (i,j,k,l) =( DEL(i,k)*DEL(j,l) + DEL(i,l)*DEL(j,k) )*0.5d0
! enddo;enddo;enddo;enddo

DELDEL = tensor_prod_ij_kl(DEL,DEL)
DELCROSS = tensor_prod_ik_jl(DEL,DEL)

do l=1,3; do k=1,3; do j=1,3; do i=1,3
   DDSYMM  (i,j,k,l) =( DELCROSS(i,j,k,l) + DELCROSS(i,j,l,k) )*0.5d0
enddo;enddo;enddo;enddo


LEVI = 0.d0
LEVI(1,2,3)=  1.d0
LEVI(2,3,1)=  1.d0
LEVI(3,1,2)=  1.d0
LEVI(1,3,2)= -1.d0
LEVI(2,1,3)= -1.d0
LEVI(3,2,1)= -1.d0

LEVI2 = 0.d0
do n=1,3; do m=1,3; do l=1,3; do k=1,3; do j=1,3; do i=1,3
   LEVI2(i,j,k,l,m,n)=LEVI(i,k,m)*LEVI(j,l,n)
enddo;enddo;enddo;enddo;enddo;enddo

end subroutine




subroutine read_materials(Fp)
implicit none
character(len=*), intent(in) :: Fp
integer :: imat
integer, parameter :: nin = 22
!
! file open
open(unit=nin,file=Fp, &
    form='formatted',access='sequential',status='old',action='read')
! 
! read number of materials
read(nin,*) NR_MATERIALS

allocate(MATERIALS(NR_MATERIALS))

do imat=1,NR_MATERIALS
   ! write(*,*) 'Material no.',imat

   read(nin,*) MATERIALS(imat)%LABEL
   ! write(*,*) 'LABEL     :',MATERIALS(imat)%LABEL

   read(nin,*) MATERIALS(imat)%DENSITY

   read(nin,*) MATERIALS(imat)%WTYPE
   ! write(*,*) 'WTYPE     :',MATERIALS(imat)%WTYPE

   read(nin,*) MATERIALS(imat)%CONSTIT
   ! write(*,*) 'CONSTIT   :',MATERIALS(imat)%CONSTIT

   read(nin,*) MATERIALS(imat)%NR_PARAMS
   ! write(*,*) 'NR_PARAMS :',MATERIALS(imat)%NR_PARAMS

   if (MATERIALS(imat)%NR_PARAMS.gt.MAX_NR_P) then
      write(*,*) 'read_materials: number of material parameters exceeded. Increase MAX_NR_P'
      write(*,*) 'MAX_NR_P=',MAX_NR_P
      write(*,*) 'imat,MATERIALS(imat)%NR_PARAMS=',imat,MATERIALS(imat)%NR_PARAMS
      stop
   endif

   allocate(MATERIALS(imat)%PARAMS(MATERIALS(imat)%NR_PARAMS))
   MATERIALS(imat)%PARAMS = 0.d0
   read(nin,*) MATERIALS(imat)%PARAMS(1:MATERIALS(imat)%NR_PARAMS)
   ! write(*,*) 'PARAMS    :',MATERIALS(imat)%PARAMS(:)
   
   read(nin,*) MATERIALS(imat)%FLAG_ANISO
   ! write(*,*) 'FLAG_ANISO:',MATERIALS(imat)%FLAG_ANISO
   
   read(nin,*) MATERIALS(imat)%FLAG_HETER
   ! write(*,*) 'FLAG_HETER:',MATERIALS(imat)%FLAG_HETER

   read(nin,*) MATERIALS(imat)%FLAG_INCOM
   ! write(*,*) 'FLAG_INCOM:',MATERIALS(imat)%FLAG_INCOM
enddo

close(nin)

write(*,*) 'read_materials: materials file has been read.'

end subroutine








subroutine check_material_consistency
implicit none
integer :: imat,np,mflag(NR_MATERIALS)

mflag = 0
write(*,*) ''  
write(*,*) 'check_material_consistency: starting...'
do imat=1,NR_MATERIALS
   write(*,*) 'Material no.',imat
   write(*,*) 'LABEL     :',MATERIALS(imat)%LABEL
   write(*,*) 'DENSITY   :',MATERIALS(imat)%DENSITY
   write(*,*) 'WTYPE     :',MATERIALS(imat)%WTYPE
   write(*,*) 'CONSTIT   :',MATERIALS(imat)%CONSTIT
   write(*,*) 'NR_PARAMS :',MATERIALS(imat)%NR_PARAMS
   write(*,*) 'PARAMS    :',MATERIALS(imat)%PARAMS(:)
   write(*,*) 'FLAG_ANISO:',MATERIALS(imat)%FLAG_ANISO
   write(*,*) 'FLAG_HETER:',MATERIALS(imat)%FLAG_HETER
   write(*,*) 'FLAG_INCOM:',MATERIALS(imat)%FLAG_INCOM

   if (MATERIALS(imat)%DENSITY.le.1.d-8) then
      mflag(imat) = 1
      write(*,*) 'Density must be greater than 1E-8'
   endif

   np = MATERIALS(Imat)%NR_PARAMS
   select case(MATERIALS(Imat)%CONSTIT)
   case(LINEAR)
      if (np.ne.2 .and. .not.MATERIALS(Imat)%FLAG_ANISO) then
         mflag(imat) = 1
         write(*,*) 'An isotropic linear material needs exactly 2 parameters (E and nu)'
      endif
      if (np.eq.2 .and. MATERIALS(Imat)%FLAG_ANISO) then
         mflag(imat) = 1
         write(*,*) 'An anisotropic linear material needs more than 2 parameters'
      endif

   case(MOONEY)
      if (np.ne.3 .and. .not.MATERIALS(Imat)%FLAG_INCOM) then
         mflag(imat) = 1
         write(*,*) 'A compressible Mooney material needs exactly 3 parameters (C_10,C_01,K)'
      endif
      if (np.eq.2 .and. MATERIALS(Imat)%FLAG_INCOM) then
         mflag(imat) = 1
         write(*,*) 'An incompressible Mooney material needs just 2 parameters (C_10,C_01)'
      endif         

   case(GENT)
      if (np.ne.3 .and. .not.MATERIALS(Imat)%FLAG_INCOM) then
         mflag(imat) = 1
         write(*,*) 'A compressible Gent material needs exactly 3 parameters (mu,Jm,K)'
      endif
      if (np.eq.2 .and. MATERIALS(Imat)%FLAG_INCOM) then
         mflag(imat) = 1
         write(*,*) 'An incompressible Gent material needs just 2 parameters (mu,Jm)'
      endif

   case(OGDEN)

      if (mod(np,4).ne.1 .and. .not.MATERIALS(Imat)%FLAG_INCOM) then
         mflag(imat) = 1
         write(*,*) 'A compressible Ogden material needs 4n+1 parameters ({a_i,ca_i,b_i,cb_i},i=1,...,n ; K)'
      endif
      if (mod(np,4).eq.0 .and. MATERIALS(Imat)%FLAG_INCOM) then
         mflag(imat) = 1
         write(*,*) 'An incompressible Ogden material needs 4n parameters ({a_i,ca_i,b_i,cb_i},i=1,...,n)'
      endif

   case default
      mflag(imat) = 1
      write(*,*) 'Unknown constitutive law'


   end select



   if (mflag(imat).eq.1) then
      write(*,*) 'Data of material ',MATERIALS(Imat)%LABEL,' is NOT consistent!'
   else
      write(*,*) 'Data of material ',MATERIALS(Imat)%LABEL,' is consistent'
   endif

enddo 

if (sum(mflag).gt.0) then
   write(*,*) ''
   write(*,*) 'Material inconsistencies have been found. Please fix and try again'
   stop
endif
write(*,*) ''
write(*,*) 'check_material_consistency: DONE!'


end subroutine
! 
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! subroutine to find material as a function of element number
! 
subroutine find_material(Mdle,Imat)
use data_structure3D, only : NODES,ELEMS
use GMP, only : Ndomain_gmp
implicit none
integer, intent(in)  :: Mdle
integer, intent(out) :: Imat

integer :: idom,nfath,nson

! write(*,*) 'find_material: Mdle=',Mdle
! nson = Mdle
nfath = Mdle
do while(nfath.gt.0)
   nson = nfath
   nfath = NODES(nson)%father
enddo
! write(*,*) 'find_material: nfath=',nfath
! write(*,*) 'find_material: nson =',nson
Idom = Ndomain_gmp(ELEMS(nson)%GMPblock)

Imat = GMP_MAT(idom)

end subroutine
! 
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! function to transform derivatives of W with respect to other variables to dWdF, d2WdF
! 
subroutine eval_strain_energy_W_F(Imat,X,F,W,dWdF,d2WdF)
implicit none
integer,intent(in)::Imat
real*8,intent(in) ::X(3),F(3,3)
real*8,intent(out)::W,dWdF(3,3),d2WdF(3,3,3,3)
! 
integer :: i,j,ii,jj,kk,ll,alph, beta, iprint
real*8 :: FT(3,3), C(3,3), fi1,fi2,fi3,di1(3,3),di2(3,3),di3(3,3),  &
          d2i1(3,3,3,3),d2i2(3,3,3,3),d2i3(3,3,3,3),dW(3),d2W(3,3), &
          dWdC(3,3),d2WdC(3,3,3,3),QL(3,3),S(3),QR(3,3),            &
          dLdC(3,3,3),dqdC(3,3,3,3),d2LdC(3,3,3,3,3),               &
          CC(6,6),DDI1(6,6),DDI2(6,6),DDI3(6,6),                    &
          Lsq(3),crosstensor(3,3),factor !,                                   &
          ! cof(3,3),det,dWdFp,dWdH,dWdJ,d2WdFp,d2WdH,d2WdJ
! 
iprint = 0
W = 0.d0; dWdF = 0.d0; d2WdF = 0.d0

FT = transpose(F)
C = matmul(FT,F)

select case(MATERIALS(Imat)%WTYPE)

case(W_F)

   call eval_W_F(Imat,X,F,W,dWdF,d2WdF)



case(W_C)
   
   call eval_W_C(Imat,X,C,W,dWdC,d2WdC)

   dWdF = 2.d0*matmul(F,dWdC)

   do jj=1,3; do j=1,3; do ll=1,3; do ii=1,3; do kk=1,3; do i=1,3

      d2WdF(i,ii,j,jj) = d2WdF(i,ii,j,jj)                         &
                       + 4.d0*F(i,kk)*d2WdC(kk,ii,ll,jj)*FT(ll,j)

   enddo;enddo;enddo;enddo;enddo;enddo

   do jj=1,3; do j=1,3; do ii=1,3; do i=1,3
      if (i.eq.j)  d2WdF(i,ii,j,jj) = d2WdF(i,ii,j,jj) + 2.d0*dWdC(ii,jj)
   enddo;enddo;enddo;enddo



case(W_INVAR , W_INVAR_DEV)

   call get_invar_symm(C,fi1,fi2,fi3)
   call get_deriv1_invar_symm(C,di1,di2,di3)
   call get_deriv2_invar_symm(C,d2i1,d2i2,d2i3)


   if (MATERIALS(Imat)%WTYPE.eq.W_INVAR) then
      call eval_W_invar(Imat,X,fi1,fi2,fi3,W,dW,d2W)
   else
      call eval_W_invar_dev(Imat,X,fi1,fi2,fi3,W,dW,d2W)
   endif
   dWdC  = dW(1)*di1 + dW(2)*di2 + dW(3)*di3
   d2WdC = dW(1)*d2i1+ dW(2)*d2i2+ dW(3)*d2i3  &
         + d2W(1,1)*tensor_prod_ij_kl(di1,di1) &
         + d2W(2,1)*tensor_prod_ij_kl(di1,di2) &
         + d2W(3,1)*tensor_prod_ij_kl(di1,di3) &
         + d2W(1,2)*tensor_prod_ij_kl(di2,di1) &
         + d2W(2,2)*tensor_prod_ij_kl(di2,di2) &
         + d2W(3,2)*tensor_prod_ij_kl(di2,di3) &
         + d2W(1,3)*tensor_prod_ij_kl(di3,di1) &
         + d2W(2,3)*tensor_prod_ij_kl(di3,di2) &
         + d2W(3,3)*tensor_prod_ij_kl(di3,di3)

   dWdF = 2.d0*matmul(F,dWdC)

   do jj=1,3; do j=1,3; do ll=1,3; do ii=1,3; do kk=1,3; do i=1,3

      d2WdF(i,ii,j,jj) = d2WdF(i,ii,j,jj)                         &
                       + 4.d0*F(i,kk)*d2WdC(kk,ii,ll,jj)*FT(ll,j)

   enddo;enddo;enddo;enddo;enddo;enddo

   do jj=1,3; do j=1,3; do ii=1,3; do i=1,3
      if (i.eq.j)  d2WdF(i,ii,j,jj) = d2WdF(i,ii,j,jj) + 2.d0*dWdC(ii,jj)
   enddo;enddo;enddo;enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEBUG PRINT STATEMENTS
if (iprint.ne.0) then
   CC = 0.d0
   j = 0
   do jj=1,3;do ll=1,jj
      j = j+1
      i = 0
      do ii=1,3; do kk=1,ii
         i = i+1
         CC(i,j) = 4.d0*d2WdC(kk,ii,ll,jj)
         DDI1(i,j) = d2i1(kk,ii,ll,jj)
         DDI2(i,j) = d2i2(kk,ii,ll,jj)
         DDI3(i,j) = d2i3(kk,ii,ll,jj)
      enddo; enddo
   enddo;enddo
   
   write(*,*) 'F'
   write(*,7760) F
   write(*,*) 'F^T*F'
   write(*,7760) C
   write(*,*) 'QR'
   write(*,*) ''
   write(*,*) 'di1'
   write(*,7760) di1
   write(*,*) 'di2'
   write(*,7760) di2
   write(*,*) 'di3'
   write(*,7760) di3
   write(*,*) 'dW'
   write(*,7760) dW
   write(*,*) 'dWdC'  
   write(*,7760) dWdC
   write(*,*) 'dWdF'
   write(*,7760) dWdF

   write(*,*) 'd2W'
   write(*,7760) d2W(1,:)
   write(*,7760) d2W(2,:)
   write(*,7760) d2W(3,:)

 7760 format(9(es11.3,','))
   write(*,*) 'd2WdE in matrix form:'
   do i=1,6
      write(*,7770) CC(i,:)
 7770 format(6(es11.3,','))
   enddo
   write(*,*) 'd2I1 in matrix form:'
   do i=1,6
      write(*,7770) DDI1(i,:)
   enddo
   write(*,*) 'd2I2 in matrix form:'
   do i=1,6
      write(*,7770) DDI2(i,:)
   enddo
   write(*,*) 'd2I3 in matrix form:'
   do i=1,6
      write(*,7770) DDI3(i,:)
   enddo
endif
   ! call pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



case(W_STRETCH , W_STRETCH_DEV)


!    write(*,*) 'WITH EIGENVALUES OF C'
   call get_eig_sym(C,Lsq,QR)

   if (MATERIALS(Imat)%WTYPE.eq.W_STRETCH) then
      call eval_W_stretch_sq(Imat,X,Lsq(1),Lsq(2),Lsq(3),W,dW,d2W)
   else
      call eval_W_stretch_sq_dev(Imat,X,Lsq(1),Lsq(2),Lsq(3),W,dW,d2W)
   endif

   call get_eig_tensors(S,QR,dLdC)

!  based on equation 6.1.44 of Ogden's book
   dWdC = dW(1)*dLdC(:,:,1) &
        + dW(2)*dLdC(:,:,2) &
        + dW(3)*dLdC(:,:,3)

   d2WdC = 0.d0


   do beta=1,3; do alph=1,3
      d2WdC = d2WdC + d2W(alph,beta) *                                 &
                      tensor_prod_ij_kl(dLdC(:,:,alph),dLdC(:,:,beta))

      if (alph.eq.beta) cycle

      if (abs(Lsq(alph)-Lsq(beta)).gt.1.d-12) then
         factor = (dW(alph)-dW(beta)) / (Lsq(alph)-Lsq(beta))
      else
         factor = d2W(alph,alph) - d2W(beta,alph)
      endif

      crosstensor = outer_product(QR(:,alph),QR(:,beta))

      d2WdC = d2WdC + 0.5d0 * factor *                                &
              tensor_prod_ij_kl(crosstensor,                          &
                                crosstensor + transpose(crosstensor))

   enddo; enddo



   ! write(*,*) 'WITH SINGULAR VALUES OF F'
   ! call get_svd_mat(F,QL,S,QR)

   ! if (MATERIALS(Imat)%WTYPE.eq.W_STRETCH) then
   !    call eval_W_stretch(Imat,X,S(1),S(2),S(3),W,dW,d2W)
   ! else
   !    call eval_W_stretch_dev(Imat,X,S(1),S(2),S(3),W,dW,d2W)
   ! endif

   ! call get_deriv_eig(S,QR,dLdC,d2LdC)

   ! dWdC = dW(1)*dLdC(:,:,1) &
   !      + dW(2)*dLdC(:,:,2) &
   !      + dW(3)*dLdC(:,:,3)

   ! d2WdC = dW(1)*d2LdC(:,:,:,:,1)       &
   !       + dW(2)*d2LdC(:,:,:,:,2)       &
   !       + dW(3)*d2LdC(:,:,:,:,3)       &
   !       + d2W(1,1)*tensor_prod_ij_kl(dLdC(:,:,1),dLdC(:,:,1)) &
   !       + d2W(2,1)*tensor_prod_ij_kl(dLdC(:,:,1),dLdC(:,:,2)) &
   !       + d2W(3,1)*tensor_prod_ij_kl(dLdC(:,:,1),dLdC(:,:,3)) &
   !       + d2W(1,2)*tensor_prod_ij_kl(dLdC(:,:,2),dLdC(:,:,1)) &
   !       + d2W(2,2)*tensor_prod_ij_kl(dLdC(:,:,2),dLdC(:,:,2)) &
   !       + d2W(3,2)*tensor_prod_ij_kl(dLdC(:,:,2),dLdC(:,:,3)) &
   !       + d2W(1,3)*tensor_prod_ij_kl(dLdC(:,:,3),dLdC(:,:,1)) &
   !       + d2W(2,3)*tensor_prod_ij_kl(dLdC(:,:,3),dLdC(:,:,2)) &
   !       + d2W(3,3)*tensor_prod_ij_kl(dLdC(:,:,3),dLdC(:,:,3))



   dWdF = 2.d0*matmul(F,dWdC)

   do jj=1,3; do j=1,3; do ll=1,3; do ii=1,3; do kk=1,3; do i=1,3

      d2WdF(i,ii,j,jj) = d2WdF(i,ii,j,jj)                         &
                       + 4.d0*F(i,kk)*d2WdC(kk,ii,ll,jj)*FT(ll,j)

   enddo;enddo;enddo;enddo;enddo;enddo

   do jj=1,3; do j=1,3; do ii=1,3; do i=1,3
      if (i.eq.j)  d2WdF(i,ii,j,jj) = d2WdF(i,ii,j,jj) + 2.d0*dWdC(ii,jj)
   enddo;enddo;enddo;enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEBUG PRINT STATEMENTS
if (iprint.ne.0) then
   CC = 0.d0
   j = 0
   do jj=1,3;do ll=1,jj
      j = j+1
      i = 0
      do ii=1,3; do kk=1,ii
         i = i+1
         CC(i,j) = 4.d0*d2WdC(kk,ii,ll,jj)
         DDI1(i,j) = d2LdC(kk,ii,ll,jj,1)
         DDI2(i,j) = d2LdC(kk,ii,ll,jj,2)
         DDI3(i,j) = d2LdC(kk,ii,ll,jj,3)
      enddo; enddo
   enddo;enddo
   
   write(*,*) 'F'
   write(*,7760) F
   write(*,*) 'F^T*F'
   write(*,7760) C
   write(*,*) 'QR'
   write(*,7760) QR
   write(*,*) 'dL1'
   write(*,7760) dLdC(:,:,1)
   write(*,*) 'dL2'
   write(*,7760) dLdC(:,:,2)
   write(*,*) 'dL3'
   write(*,7760) dLdC(:,:,3)
   write(*,*) 'dW'
   write(*,7760) dW
   write(*,*) 'dWdC'

   write(*,7760) dWdC

   write(*,*) 'dWdF'
   write(*,7760) dWdF

   write(*,*) 'd2W'
   write(*,7760) d2W(1,:)
   write(*,7760) d2W(2,:)
   write(*,7760) d2W(3,:)

 ! 7760 format(9(es11.3,','))
   write(*,*) 'd2WdE in matrix form:'
   do i=1,6
      write(*,7770) CC(i,:)
 ! 7770 format(6(es11.3,','))
   enddo
   write(*,*) 'd2L1 in matrix form:'
   do i=1,6
      write(*,7770) !DDI1(i,:)
   enddo
   write(*,*) 'd2L2 in matrix form:'
   do i=1,6
      write(*,7770) !DDI2(i,:)
   enddo
   write(*,*) 'd2L3 in matrix form:'
   do i=1,6
      write(*,7770) !DDI3(i,:)
   enddo
endif
   ! call pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



case default

   write(*,*) 'eval_strain_energy_W_F: we still do not support the energy function of type .', &
              MATERIALS(Imat)%WTYPE

! case(W_POLYCVX)
   ! call get_cof_mat(F,cof)

   ! call get_det_mat(F,det)

   ! call eval_W_polycvx(Imat,X,F,cof,det,W,dWdFp,d2WdFp,dWdH,d2WdH,dWdJ,d2WdJ)



end select

end subroutine











subroutine get_eig_tensors(S,QR,dLdC)
implicit none
real*8, intent(in) :: S(3),QR(3,3)
real*8, intent(out):: dLdC(3,3,3)

real*8 :: dqdC(3,3,3,3),Yps(3,3,3),dLsq_dC(3,3,3),denom1,denom2
integer:: i,j,ii,jj,kk,ll,alph,beta

!  the derivative of lambda^2 and lambda w.r.t C
   do alph=1,3
      dLdC(:,:,alph) = outer_product(QR(:,alph),QR(:,alph))
   enddo


end subroutine














! get derivatives of principal stretches w.r.t symmetric tensor C
! based on Magnus(1985), Ellen Kuhl's lecture notes and own calculations

subroutine get_deriv_eig(S,QR,dLdC,d2LdC)
implicit none
real*8, intent(in) :: S(3),QR(3,3)
real*8, intent(out):: dLdC(3,3,3),d2LdC(3,3,3,3,3)

real*8 :: dqdC(3,3,3,3),Yps(3,3,3),dLsq_dC(3,3,3),denom1,denom2
integer:: i,j,ii,jj,kk,ll,alph,beta

!  the derivative of lambda^2 and lambda w.r.t C
   do alph=1,3
      dLsq_dC(:,:,alph) = outer_product(QR(:,alph),QR(:,alph))
      dLdC(:,:,alph) = 0.5d0/S(alph) * dLsq_dC(:,:,alph)
   enddo

!  now, the derivative of the eigenvectors

!  first, the auxiliary matrix Yps (pseudoinverse of Y=lambda^2*I-C)
   Yps = 0.d0
   do alph=1,3; do beta=1,3

      if (beta.eq.alph) cycle

      denom1= S(alph)+S(beta)
      denom2= S(alph)-S(beta)
! we must divide by the denominator S(alph)^2 - S(beta)^2, but we split it
! in two factors denom1*denom2 to check size of denom2=S(alph)-S(beta)
      if (abs(denom2).ge.1.d-8) then
         Yps(:,:,alph) = Yps(:,:,alph)                                    &
                       + ( dLsq_dC(:,:,beta) / denom1) / denom2 
      endif
   enddo; enddo
!use Yps in the formula  (dq/dC)_pIK = (Yps)_pm DDSYMM_mnIK q_n  herein expanded
   dqdC = 0.d0
   do alph=1,3; do kk=1,3; do ii=1,3; do i=1,3
      dqdC(i,ii,kk,alph) = dqdC(i,ii,kk,alph)                             &
                         + 0.5d0*( Yps(i,ii,alph)*QR(kk,alph)             &
                                  +Yps(i,kk,alph)*QR(ii,alph) )
   enddo; enddo; enddo; enddo


!  now, the second derivative of lambda^2 w.r.t C
   d2LdC = 0.d0
   do alph=1,3; do ll=1,3; do jj=1,3; do kk=1,3; do ii=1,3
      d2LdC(ii,kk,jj,ll,alph)  = d2LdC(ii,kk,jj,ll,alph)                  &
                               + 0.5d0*( dqdC(ii,jj,ll,alph)*QR(kk,alph)  &
                                        +dqdC(kk,jj,ll,alph)*QR(ii,alph)  &
                                        +dqdC(jj,ii,kk,alph)*QR(ll,alph)  &
                                        +dqdC(ll,ii,kk,alph)*QR(jj,alph) )
   enddo; enddo; enddo; enddo; enddo

!  now, the second derivative of lambda w.r.t C (we overwrite the same array)
   do alph=1,3; do ll=1,3; do jj=1,3; do kk=1,3; do ii=1,3
      d2LdC(ii,kk,jj,ll,alph) = 1.d0/S(alph)*( 0.5d0*d2LdC(ii,kk,jj,ll,alph)      &
                                              -dLdC(ii,kk,alph)*dLdC(jj,ll,alph) )
   enddo; enddo; enddo; enddo; enddo

end subroutine








subroutine test_stiffness_tensor
! use hyperelasticity
implicit none
integer :: npoints, i,j,k,mm,nn,m,n,imat,ipiv(10),info, kk,ll
real*8 :: dl,fl1,fl2,fl3,fJ,Ftensor(3,3),X(3),W,dW(3,3),d2W(3,3,3,3), &
          C(10,10),A(10,10),work(100),FnormC,FnormA,dcofF(3,3,3,3),pr

real*8 :: CofF(3,3),cf(9),QR(3,3),QL(3,3)

npoints = 10
dl = 0.02d0/real(npoints,8)

do imat=1,NR_MATERIALS

   if (MATERIALS(imat)%CONSTIT.gt.0) then
      write(*,*) 'Testing elasticity tensor for material ',imat

      do k=0,0!npoints
         fl3 = 1.d0!0.96d0+dl*real(k,8)
         do j=0,0!npoints
            fl2 = 1.d0!0.96d0+dl*real(j,8)
            do i=0,npoints
               fl1=0.99d0+dl*real(i,8)

               fJ = fl1*fl2*fl3
               write(*,*) 'i,j,k=',i,j,k
               write(*,3333) fl1,fl2,fl3
               write(*,3334) fJ
 3333 format('Principal stretches='3(es10.3))
 3334 format('det(F)=',es13.6)

               X = 0.d0

               QL(:,1) = (/ cos(0.5d0),sin(0.5d0),0.d0/)
               QL(:,2) = (/-sin(0.5d0),cos(0.5d0),0.d0/)
               QL(:,3) = (/       0.d0,      0.d0,1.d0/)

               QR(:,1) = (/ cos(0.7d0),0.d0,-sin(0.7d0)/)
               QR(:,2) = (/       0.d0,1.d0,       0.d0/)
               QR(:,3) = (/ sin(0.7d0),0.d0, cos(0.7d0)/)

               Ftensor = fl1*outer_product(QL(:,1),QR(:,1)) &
                       + fl2*outer_product(QL(:,2),QR(:,2)) &
                       + fl3*outer_product(QL(:,3),QR(:,3))
               
               Ftensor = 0.d0
               Ftensor(1,1) = fl1; Ftensor(2,2) = fl2; Ftensor(3,3) = fl3;
               ! Ftensor(1,2) = 0.1d0; Ftensor(1,3) = -0.1d0; Ftensor(2,3) = 0.1d0;
               call eval_strain_energy_W_F(Imat,X,Ftensor,W,dW,d2W)
               dcofF = 0.d0
               call get_deriv_cof_mat(Ftensor,dcofF)
               pr=1.d-5
               ! call resize_tensor_4_2(d2W,C)
               C = 0.d0
               do nn=1,3; do n=1,3
                  ll = (nn-1)*3 + n
                  do mm=1,3; do m=1,3
                     kk = (mm-1)*3 + m
                     C(kk,ll) = d2W (m,mm,n,nn) - pr*dcofF(m,mm,n,nn)
                     ! write(*,*) d2W(m,mm,n,nn)!-d2W(n,nn,m,mm)
                  enddo;enddo
               enddo;enddo

               write(*,*) 'Tensor D2WDF:'
               do kk=1,9
                  write(*,3335) C(kk,1:9)
               enddo
               dW = fJ*matmul(dW,transpose(Ftensor))
               write(*,*) 'Cauchy stress'
               write(*,3335) dW
               write(*,*) ''
               

               FnormC = 0.d0
               do nn=1,10; do mm=1,10
                  FnormC = FnormC + C(mm,nn)**2
               enddo;enddo
               FnormC = sqrt(FnormC)
               write(*,*) 'Frobenius norm of D2WDF =',FnormC

               A = 0.d0
               do mm=1,10
                  A(mm,mm) = 1.d0
               enddo

               if (MATERIALS(imat)%FLAG_INCOM) then
                  call get_cof_mat(Ftensor,CofF)
                  cf = -1.d0*reshape(CofF,(/9/))
                  write(*,*) 'Vectorized cofactor'
                  write(*,3335) cf
                  C(1:9,10) = cf(1:9)
                  C(10,1:9) = C(1:9,10)
                  C(10,10)  = 0.d0
                  call dsysv('U',10,10,C,10,ipiv,A,10,work,100,info)
                  if (info.ne.0) then
                  write(*,*) 'Inversion failed. info=',info
                  endif
               else
                  C(10,10) = 1.d0
                  call dsysv('U',10,10,C,10,ipiv,A,10,work,100,info)
                  if (info.ne.0) then
                  write(*,*) 'Inversion failed. info=',info
                  endif
               endif
               ! write(*,*) 'Optimal LWORK=',work(1)


               if (info.eq.0) then
                  write(*,*) 'Inverse of augmented D2WDF:'
                  do kk=1,10
                     write(*,3335) A(:,kk)
                  enddo
 3335 format(10(es11.3,','))
               endif

               FnormA = 0.d0
               do nn=1,10; do mm=1,10
                  FnormA = FnormA + A(mm,nn)**2
               enddo;enddo
               FnormA = sqrt(FnormA)
               write(*,*) 'Frobenius norm of inv(D2WDF) =',FnormA
               call pause
            enddo
         enddo
      enddo

   endif

enddo

end subroutine

! subroutine resize_tensor_4_2(d2W,C)



! do nn=1,3; do n=1,3
!    ll = (nn-1)*3 + n
!    do mm=1,3; do m=1,3
!       kk = (mm-1)*3 + m
!       C(kk,ll) = d2W (m,mm,n,nn)
!    enddo;enddo
! enddo;enddo


! end subroutine

subroutine invert_d2W(d2W,d2Winv)
implicit none
real*8,intent(in) :: d2W(3,3,3,3)
real*8,intent(out):: d2Winv(3,3,3,3)

integer :: m,mm,n,nn,kk,ll,ipiv(9),info
real*8 :: C(9,9),A(9,9),work(81)

   do nn=1,3; do n=1,3
      ll = (nn-1)*3 + n
      do mm=1,3; do m=1,3
         kk = (mm-1)*3 + m
         C(kk,ll) = d2W (m,mm,n,nn)
      enddo;enddo
   enddo;enddo

   A = 0.d0
   do mm=1,9
      A(mm,mm) = 1.d0
   enddo
   call dsysv('U',9,9,C,9,ipiv,A,9,work,81,info)
   if (info.ne.0) then
      write(*,*) 'Inversion failed. info=',info
   endif

   do ll=1,9
      nn = (ll-1)/3 + 1
      n  = ll-3*(nn-1)
      do kk=1,9
         mm = (kk-1)/3 + 1
         m  = kk - 3*(mm-1)
         d2Winv(m,mm,n,nn) = A(kk,ll)
      enddo
   enddo

end subroutine



! 
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! utilities for 3D tensors
! 





subroutine get_eig_sym(Mat,Lsq,Q)

implicit none
real*8,intent(in) :: Mat(3,3)
real*8,intent(out):: Q(3,3),Lsq(3)

integer :: info
real*8  :: work(8),Copy(3,3)


! Copy the matrix to avoid changes by Lapack to the original matrix
Q = Mat
! 
! Lapack double symmetric matrix eigenvalue/vector subroutine argument key
!
! subroutine dsyev  (  character   JOBZ,
! character   UPLO,
! integer  N,
! double precision, dimension( lda, * )  A,
! integer  LDA,
! double precision, dimension( * )    W,
! double precision, dimension( * )    WORK,
! integer  LWORK,
! integer  INFO 
! )  
!
! Call Lapack
call DSYEV('V','U',3,Q,3,Lsq,work,8,info)
! Subroutine DSYEV returns the eigenvectors in matrix Q
!
! check for errors
if (info.ne.0) then
   write(*,*) 'get_eig_sym: error at eigenvalue computation; info =',info
   stop
endif

end subroutine




subroutine get_svd_mat(Mat,U,S,V)

implicit none
real*8,intent(in) :: Mat(3,3)
real*8,intent(out):: U(3,3),S(3),V(3,3)

integer :: info
real*8  :: work(15),Copy(3,3)


! Create a copy to avoid changes by Lapack to the original matrix
Copy = Mat
! 
! Lapack SVD subroutine argument key
!
! dgesvd   (  character   JOBU,
! character   JOBVT,
! integer  M,
! integer  N,
! double precision, dimension( lda, * )  A,
! integer  LDA,
! double precision, dimension( * )    S,
! double precision, dimension( ldu, * )  U,
! integer  LDU,
! double precision, dimension( ldvt, * )    VT,
! integer  LDVT,
! double precision, dimension( * )    WORK,
! integer  LWORK, >=5*max(M,N)
! integer  INFO 
! )
!
! Call Lapack
call DGESVD('A','A',3,3,Copy,3,S,U,3,V,3,work,15,info)
! Subroutine DGESVD returns the transpose of V, but we prefer to have the
! singular vectors as columns
! 
V = transpose (V)
! check for errors
if (info.ne.0) then
   write(*,*) 'get_svd_mat: error at svd computation; info =',info
   stop
endif

end subroutine









subroutine get_cof_mat(F,Cof)
implicit none
real*8,intent(in) :: F(3,3)
real*8,intent(out):: Cof(3,3)

integer, dimension(3) :: next = (/2,3,1/)
integer, dimension(3) :: npre = (/3,1,2/)
integer :: i,j

do j=1,3; do i=1,3
   Cof(i,j) = F(next(i),next(j))*F(npre(i),npre(j)) &
            - F(next(i),npre(j))*F(npre(i),next(j))
enddo; enddo

! Cof(1,1) = F(2,2)*F(3,3)-F(3,2)*F(2,3)
! Cof(2,1) = F(3,2)*F(1,3)-F(1,2)*F(3,3)
! Cof(3,1) = F(1,2)*F(2,3)-F(2,2)*F(1,3)

! Cof(1,2) = F(2,3)*F(3,1)-F(3,3)*F(2,1)
! Cof(2,2) = F(3,3)*F(1,1)-F(1,3)*F(3,1)
! Cof(3,2) = F(1,3)*F(2,1)-F(2,3)*F(1,1)

! Cof(1,3) = F(2,1)*F(3,2)-F(3,1)*F(2,2)
! Cof(2,3) = F(3,1)*F(1,2)-F(1,1)*F(3,2)
! Cof(3,3) = F(1,1)*F(2,2)-F(2,1)*F(1,2)
end subroutine










subroutine get_det_mat(F,Det)
implicit none
real*8, intent(in ) :: F(3,3)
real*8, intent(out) :: Det

Det = F(1,1)*(F(2,2)*F(3,3)-F(3,2)*F(2,3)) &
    + F(2,1)*(F(3,2)*F(1,3)-F(1,2)*F(3,3)) &
    + F(3,1)*(F(1,2)*F(2,3)-F(2,2)*F(1,3))

end subroutine







subroutine get_inv_mat(F,Finv)
implicit none
real*8, intent(in ) :: F(3,3)
real*8, intent(out) :: Finv(3,3)

real*8 :: Cof(3,3),Det

call get_cof_mat(F,Cof)
call get_det_mat(F,Det)

Finv = transpose(Cof)/Det

end subroutine
  






subroutine get_invar_symm(S,FI1,FI2,FI3)
implicit none
real*8, intent(in ) :: S(3,3)
real*8, intent(out) :: FI1,FI2,FI3

call get_trace_mat(S,FI1)

call get_trace_mat(matmul(S,S),FI2)
FI2 = 0.5d0*(FI1**2-FI2)

call get_det_mat(S,FI3)

end subroutine







subroutine get_trace_mat(S,Tr)
implicit none
real*8, intent(in ) :: S(3,3)
real*8, intent(out) :: Tr

Tr = S(1,1) + S(2,2) + S(3,3)

end subroutine







subroutine get_deriv1_invar_symm(S,DI1,DI2,DI3)
implicit none
real*8, intent(in ) :: S(3,3)
real*8, intent(out) :: DI1(3,3),DI2(3,3),DI3(3,3)

real*8 :: trs

DI1 = DEL

call get_trace_mat(S,trs)
DI2 = trs*DEL - S

call get_cof_mat(S,DI3)

end subroutine







subroutine get_deriv2_invar_symm(S,D2I1,D2I2,D2I3)
implicit none
real*8, intent(in ) :: S(3,3)
real*8, intent(out) :: D2I1(3,3,3,3),D2I2(3,3,3,3),D2I3(3,3,3,3)

real*8 :: trs
integer:: i,j,k,l

D2I1 = 0.d0

D2I2 = DELDEL - DDSYMM

call get_trace_mat(S,trs)

do l=1,3; do k=1,3; do j=1,3; do i=1,3
   D2I3(i,j,k,l) = 0.5d0 *( DEL(i,k)*S(j,l)   &
                           +DEL(i,l)*S(j,k)   &
                           +DEL(j,l)*S(i,k)   &
                           +DEL(j,k)*S(i,l) )
enddo;enddo;enddo;enddo

D2I3 = D2I3 - tensor_prod_ij_kl(DEL,S) - tensor_prod_ij_kl(S,DEL) &
     + trs*(DELDEL - DDSYMM)


end subroutine








subroutine get_deriv_cof_mat(F,Dcof)
implicit none
real*8,intent(in) :: F(3,3)
real*8,intent(out):: Dcof(3,3,3,3)

integer :: i,j,k,l,m,n

Dcof = 0.d0

do n=1,3; do m=1,3; do l=1,3; do k=1,3; do j=1,3; do i=1,3
   Dcof(i,j,k,l) = Dcof(i,j,k,l) + LEVI2(i,j,k,l,m,n)*F(m,n)
enddo; enddo; enddo; enddo; enddo; enddo

end subroutine








subroutine get_deriv2_cof_mat(F,D2cof)
implicit none
real*8,intent(in) :: F(3,3)
real*8,intent(out):: D2cof(3,3,3,3,3,3)

D2cof = LEVI2

end subroutine








function outer_product(Vcol,Vrow)result(Prod)
implicit none
real*8,intent(in) :: Vcol(3),Vrow(3)
real*8 :: Prod(3,3)

Prod(:,1) = Vcol(:)*Vrow(1)
Prod(:,2) = Vcol(:)*Vrow(2)
Prod(:,3) = Vcol(:)*Vrow(3)

end function







function double_dot_prod(A,B)result(C)
implicit none
real*8, intent(in ) :: A(3,3),B(3,3)
real*8 :: C

integer :: i

C = 0.d0
do i=1,3
   C = C + dot_product(A(:,i),B(:,i))
enddo

end function







function contraction_4_2_majorsymm(A,B)result(C)
implicit none
real*8, intent(in ) :: A(3,3,3,3),B(3,3)
real*8 :: C(3,3)

integer :: i,j,k

C = 0.d0
do k=1,3; do j=1,3; do i=1,3
   C(j,k) = C(j,k) + dot_product(A(:,i,j,k),B(:,i))
enddo; enddo; enddo

end function






function tensor_prod_ij_kl(A,B)result(C)
implicit none
real*8, intent(in ) :: A(3,3),B(3,3)
real*8 :: C(3,3,3,3)

integer :: i,j,k,l

do l=1,3; do k=1,3; do j=1,3; do i=1,3
   C(i,j,k,l) = A(i,j)*B(k,l)
enddo;enddo;enddo;enddo

end function







function tensor_prod_ik_jl(A,B)result(C)
implicit none
real*8, intent(in ) :: A(3,3),B(3,3)
real*8 :: C(3,3,3,3)

integer :: i,j,k,l

do l=1,3; do k=1,3; do j=1,3; do i=1,3
   C(i,j,k,l) = A(i,k)*B(j,l)
enddo;enddo;enddo;enddo

end function
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! subroutine to modify load factor
subroutine load_increment(N,Step_size)
implicit none
integer, intent(in) :: N
real*8 , intent(in) :: Step_size

LOAD_FACTOR = LOAD_FACTOR_PREV+real(N,8)*Step_size

end subroutine
end module
