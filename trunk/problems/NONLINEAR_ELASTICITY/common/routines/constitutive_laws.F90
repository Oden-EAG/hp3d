







subroutine eval_W_F(Imat,X,F,W,dWdF,d2WdF)
implicit none
integer,intent(in)::Imat
real*8,intent(in) ::X(3),F(3,3)
real*8,intent(out)::W,dWdF(3,3),d2WdF(3,3,3,3)

integer:: np,i,j,k,l
real*8 :: p(MAX_NR_P),gradu(3,3),

np = MATERIALS(Imat)%NR_PARAMS

p = 0.d0
p(1:np) = MATERIALS(Imat)%PARAMS

select case(MATERIALS(Imat)%CONSTIT)

case(LINEAR)

   if (MATERIALS(Imat)%ANISO) then
      write(*,*) 'eval_W_F: no anisotropic constitutive laws have been implemented'
      stop
   endif
   gradu = F - DEL

!  p(1) is mu, p(2) is lambda
!  define the elastictiy tensor
   do l=1,3; do k=1,3; do j=1,3; do i=1,3
      d2WdF(i,j,k,l) = p(2) * DELDEL(i,j,k,l) + 2.d0 * p(1) * DELSYMM(i,j,k,l)
   enddo; enddo; enddo; enddo

   call contraction_4_2_majorsymm(d2WdF,gradu,dWdF)

   call double_dot_prod(dWdF,gradu,W)
   W = 0.5d0*W

   ! do l=1,3; do k=1,3; do j=1,3; do i=1,3
   !    dWdF(i,j) = d2WdF(i,j,k,l)*gradu(k,l) 
   ! enddo; enddo; enddo; enddo

   ! do j=1,3; do i=1,3
   !    W = 0.5d0*dWdF(i,j)*gradu(i,j) 
   ! enddo; enddo

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








subroutine eval_W_stretch(Imat,X,S1,S2,S3,W,dW,d2W)
implicit none
integer,intent(in)::Imat
real*8,intent(in) ::X(3),S1,S2,S3
real*8,intent(out)::W,dW(3),d2W(3,3)

integer :: np
real*8  :: p(MAX_NR_P),fJ,wJ,dwds(3),d2wds(3,3),djds(3),d2jds(3,3),fi1,fi2,fi3


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
! first compute djds*djds' and store temporarily in d2wds
call outer_product(djds,djds,d2wds)
d2wds = dwdJ * d2wds + d2WdJ * d2jds

fi1 = S1**2 + S2**2 + S3**2
fi2 = S1**2*S2**2 + S1**2*S3**2 + S2**2*S3**2


select case (MATERIALS(Imat)%CONSTIT)

case(MOONEY)
   W = p(1)*(FI1-3.d0)+p(2)*(FI2-3.d0) + p(3)*wJ
   dW(1) = 2.d0*p(1)*S1 + 2.d0*p(2)*S1*(S2**2+S3**2) + p(3)*djds(1)
   dW(2) = 2.d0*p(1)*S2 + 2.d0*p(2)*S2*(S1**2+S3**2) + p(3)*djds(2)
   dW(3) = 2.d0*p(1)*S3 + 2.d0*p(2)*S3*(S1**2+S2**2) + p(3)*djds(3)
   d2W = 0.d0
   d2W(1,1) = 2.d0*p(1) + 2.d0*p(2)*(S2**2+S3**2) + p(3)*d2jds(1,1)
   d2W(2,1) = 4.d0*p(2)*S1*S2 + p(3)*d2jds(2,1)
   d2W(3,1) = 4.d0*p(2)*S1*S2 + p(3)*d2jds(2,1)
   d2W(2,2) = 2.d0*p(1) + 2.d0*p(2)*(S1**2+S3**2) + p(3)*d2jds(2,2)
   d2W(3,2) = 4.d0*p(2)*S3*S2 + p(3)*d2jds(3,2)
   d2W(3,3) = 2.d0*p(1) + 2.d0*p(2)*(S1**2+S2**2) + p(3)*d2jds(3,3)
   d2W(2,3) = d2W(3,2)
   d2W(1,3) = d2W(3,1)
   d2W(1,2) = d2W(2,1)

case(OGDEN)
   W = p(np)*wJ
   dW(:) = p(np)*djds(:)
   d2W(:,:) = p(np)*d2jds(:,:)
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
      
   enddo

case default

   write(*,*) 'eval_W_stretch: this form of W does not support constitutive law no.', &
              MATERIALS(Imat)%CONSTIT
   
   stop

end select


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

WJ = 0.5d0 * ( FJ**2 - 1.d0 ) - log(FJ)
dWdJ = FJ - 1.d0 / FJ
d2WdJ = 1.d0 + 1.d0 / (FJ**2)

end subroutine







subroutine get_deviat_F(F,FJ,Fdev)
implicit none
real*8,intent(in ) :: F(3,3),FJ
real*8,intent(out) :: Fdev(3,3)

Fdev = F*FJ**(-1.d0/3.d0)

end subroutine