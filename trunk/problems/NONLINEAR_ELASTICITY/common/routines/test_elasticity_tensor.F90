subroutine test_stiffness_tensor
use hyperelasticity
implicit none
integer :: npoints, i,j,k,mm,nn,imat,ipiv(9),info, kk,ll
real*8 :: dl,fl1,fl2,fl3,fJ,Ftensor(3,3),X(3),W,dW(3,3),d2W(3,3,3,3), &
          C(9,9),A(3,3),work(81),FnormC,FnormA

npoints = 8
dl = 0.2d0/real(npoints,8)

do imat=1,NR_MATERIALS

   if (MATERIALS(imat)%CONSTIT.gt.0) then
      write(*,*) 'Testing elasticity tensor for material ',imat

      do k=0,npoints
         fl3 = 0.9d0+dl*real(k,8)
         do j=0,npoints
            fl2 = 0.9d0+dl*real(j,8)
            do i=0,npoints
               fl1=0.9d0+dl*real(i,8)

               fJ = fl1*fl2*fl3
               write(*,*) 'i,j,k=',i,j,k
               write(*,3333) fl1,fl2,fl3
               write(*,3334) fJ
 3333 format('Principal stretches='3(es10.3))
 3334 format('det(F)=',es13.6)

               X = 0.d0
               Ftensor = 0.d0
               Ftensor(1,1) = fl1; Ftensor(2,2) = fl2; Ftensor(3,3) = fl3;

               call eval_strain_energy_W_F(Imat,X,Ftensor,W,dW,d2W)

               ! call resize_tensor_4_2(d2W,C)
               do nn=1,3; do n=1,3
                  ll = (nn-1)*3 + n
                  do mm=1,3; do m=1,3
                     kk = (mm-1)*3 + m
                     C(kk,ll) = d2W (m,mm,n,nn)
                  enddo;enddo
               enddo;enddo

               FnormC = 0.d0
               do nn=1,9; do mm=1,9
                  FnormC = FnormC + C(mm,nn)**2
               enddo;enddo;enddo;enddo
               FnormC = sqrt(FnormC)
               write(*,*) 'Frobenius norm of C =',FnormC

               A = 0.d0
               do mm=1,9
                  A(mm,mm) = 1.d0
               enddo
               call dsysrv('U',9,9,C,9,ipiv,A,9,work,81,info)
               if (info.ne.0) then
                  write(*,*) 'Inversion failed. info=',info
               enddo

               FnormA = 0.d0
               do nn=1,9; do mm=1,9
                  FnormA = FnormA + A(mm,nn)**2
               enddo;enddo;enddo;enddo
               FnormA = sqrt(FnormA)
               write(*,*) 'Frobenius norm of inv(C) =',FnormA
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

integer :: m,mm,n,nn,kk,ll,ipic(9),info
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
   call dsysrv('U',9,9,C,9,ipiv,A,9,work,81,info)
   if (info.ne.0) then
      write(*,*) 'Inversion failed. info=',info
   enddo

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