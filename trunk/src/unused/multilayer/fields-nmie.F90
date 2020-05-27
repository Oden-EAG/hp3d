
!................................................................*
! NB! In order to treat very large spheres,                      *
!     one needs to enlarge the parameter NTERMS.                 *
!................................................................*
!This is a modified version of the original code((c) 1999 Sobolev*
! Astronomical Institute, St. Petersburg Univ).----available at  *
! http://www.astro.spbu.ru/staff/ilin2/SOFTWARE/nmie0.html..based*
! on the following publication:                                  *
! Recursive algorithms of Wu & Wang (Radio Sci. 26, 1393, 1991)  *
! created by N.V. Voshchinnikov                                  *
!
!*****************************************************************
!
module mie_parameters
  implicit none
  save
  integer,parameter::n_layers=10,nterms=100,dp=kind(0.0d0)
  real(kind=dp),parameter::pi=4.0_dp*atan(1.0_dp),vlite=2.99792458E008_dp
  integer::n_l,exc
  logical::Jdipole_exc,Mdipole_exc,plane_wave_exc
  integer,parameter::one=1
end module mie_parameters

module mie_variables
  use mie_parameters 
  implicit none
  save
  complex(kind=dp)::rrbb(n_layers,nterms),rrcc(n_layers,nterms),rrd1(n_layers,nterms),&
       rrd2(n_layers,nterms),rrd3(n_layers,nterms),srbb(n_layers,nterms),&
       srd1(n_layers,nterms),srd2(n_layers,nterms),&
       rd11(nterms),rd3x(nterms),&
       rcx(nterms),ri_n(0:n_layers)
  complex(kind=dp)::eps(0:n_layers),mu(0:n_layers),arg_p
  real(kind=dp)::d1x(nterms),lambda
  real(kind=dp),allocatable::cyr_p(:),cyi_p(:)
  real(kind=dp),dimension(n_layers):: a
  real(kind=dp)::r_p!!distance of the center of dipole from origin
  real(kind=dp)::phi_p,theta_p
!!angle of the dipole center,corresponding cos arg
  real(kind=dp)::posit_dipole(3)!!position vector of the center of dipole
  real(kind=dp)::orient_dipole(3)!!unit vector along the orientation of dipole
  real(kind=dp)::orient_dipole_cart(3)!!unit vector along the orientation of dipole in cartesian coordinates
  real(kind=dp)::exc_dipole!!dipole length multiplied by its current amplitude
  real(kind=dp)::the_angle,phi_angle,ethetad,ephid
  real(kind=dp)::thetastart,thetaend,phimin,phimax
  integer::ntheta,nphi
end module mie_variables

! program RCS  
!   use mie_variables
!   IMPLICIT none                                          
!   real(kind=dp)::x,al,freq
!   real(kind=dp),dimension(n_layers):: aa,aaa
!   real(kind=dp)::dummyr,dummyi
!   integer::i,j,iv
!   real(kind=dp)::curr_mag,len_dipole,sp(3)
!   complex(kind=dp)::tmp
! !44 FORMAT(F18.15)
! 20 FORMAT(2F18.15)
! 21 FORMAT(3X,'SPHERES: n-layers',1X,'/',1x,' MIE series Solution')
! 223 FORMAT(I4)                                          
! 30 FORMAT(3x,'Number of layers=',i4/&
!         3x,'Epsilon and Mu values:')
! 301 FORMAT(6x,'i=',i2,',',3x,'eps_i=',2F8.4,3x,'mu_i=',2F8.4)
! 302 FORMAT(3x,'Relative thickness:')
! 303 FORMAT(6x,'i=',i2,',',3x,'d(a_i)/a=',F7.4)
! 31 FORMAT(F18.15)                                               
! 41 FORMAT(1X,64('.'))
! 42 FORMAT(1X,64('-'))                          
! 44 FORMAT(F21.10)  
!   !                  INPUT
!   open(unit=05,file='mie.inp',status='old',access='sequential')

!   Jdipole_exc=.false.
!   Mdipole_exc=.false.
!   plane_wave_exc=.false.
  
!   ! Input
!   READ (5,*) exc

!   select case (exc)
!   case(1)
!      plane_wave_exc=.true.
!   case(2)
!      Jdipole_exc=.true.
!   case(3)
!      Mdipole_exc=.true.
!   case default
!      print*,"The input for excitation can take up values from 1 to 3 only!"
!      stop
!   end select

!   if(.not.plane_wave_exc) then
     
!      read(5,*) curr_mag,len_dipole
!      exc_dipole=curr_mag*len_dipole
!      read(5,*) posit_dipole(1),posit_dipole(2),posit_dipole(3)
 
!      read(5,*) orient_dipole(1),orient_dipole(2),orient_dipole(3)
    
!      orient_dipole=orient_dipole/sqrt(dot_product(orient_dipole,orient_dipole))
!      orient_dipole_cart=orient_dipole
!      call cart2sph(posit_dipole,sp)
!      r_p=sp(1)
!      theta_p=sp(2)
!      phi_p=sp(3)
     
!      sp(1)=sin(theta_p)*cos(phi_p)*orient_dipole(1)+&
!           sin(theta_p)*sin(phi_p)*orient_dipole(2)+&
!           cos(theta_p)*orient_dipole(3)
     
!      sp(2)=cos(theta_p)*cos(phi_p)*orient_dipole(1)+&
!           cos(theta_p)*sin(phi_p)*orient_dipole(2)-&
!           sin(theta_p)*orient_dipole(3)
     
!      sp(3)=-sin(phi_p)*orient_dipole(1)+cos(phi_p)*orient_dipole(2)
!      orient_dipole=sp
!   else

!      read(5,*) the_angle,phi_angle,ethetad,ephid
!      the_angle=the_angle*pi/180.0_dp
!      phi_angle=phi_angle*pi/180.0_dp

!   end if

!   read(5,*) thetastart, thetaend
!   read(5,*) phimin,phimax
!   read(5,*) ntheta,nphi
    
!   READ (5,*) n_l
!   print*," Number of layers: ",n_l
!   if(n_l<3) then
!      print *,'n_l < 3', n_l
!      stop
!   end if
!   if(n_l>n_layers) then
!      print *,'n_l > n_layers', n_l, n_layers
!      stop
!   end if
!   !    
!   print*,"Now reading epsilon and mu values" 

!   do i=1,n_l+1
!      read (5,*) dummyr,dummyi
!      eps(i-1)=cmplx(dummyr,dummyi,dp)
!   end do

!   do i=1,n_l+1
!      read (5,*) dummyr,dummyi
!      mu(i-1)=cmplx(dummyr,dummyi,dp) 
!   end do

!   do i=1,n_l
!      read (5,*) a(i)
!   end do

!   aa(1)=a(1)/a(n_l)
!   do i=1,n_l-1
!      aa(i+1)=(a(i+1)-a(i))/a(n_l)
!   end do
!   al = 0.0_dp
!   do i = 1, n_l-1
!      al = al + aa(i)
!   end do
!   if(al>1.0_dp) then
!      print *,'Sum a(i) > 1', al
!      stop
!   end if
!   aa(n_l) = 1.0_dp - al
     
!   aaa(1) = aa(1)
!   do i = 2, n_l
!      aaa(i) = 0.0_dp
!      do j = 1, i
!         aaa(i) = aaa(i) + aa(j)
!      end do
!   end do
  
!  !...
!   read(5,*) freq

!   print*,"Frequency is ",freq
!   print*,"The relative permittivity and permeability of the surrounding &
!        medium is: ",eps(0),mu(0)
!   lambda=vlite/freq
!   x=2.0_dp*4.0_dp*atan(1.0_dp)*a(n_l)/lambda

!   write (*,21)
!   write (*,30) n_l
!   write (*,301) (i, eps(i), mu(i), i = 1, n_l)
!   write (*,302)
!   write (*,303) (i, aa(i), i = 1, n_l)
!   write (*,41)
!   write (*,42)
  
!   do i=1,n_l+1
!      ri_n(i-1)=sqrt(conjg(eps(i-1))*conjg(mu(i-1)))
!   end do

!   call shexqn1(aa,x,freq)
!   close(5,status='keep')
! end program RCS

!--------------------------------------------------------------------
! **********   shexqn1 - Spheres: n-layers
!                        Theory: exact
!                        Results: efficiency factors
! March 1999, AI SPbU
!--------------------------------------------------------------------
SUBROUTINE shexqn1(aa,x,freq)
  use mie_variables
  IMPLICIT none
  
  real(kind=dp),dimension(n_layers),intent(in)::aa
  real(kind=dp),intent(in)::x,freq
  integer::i,j,NUM,NUM2,NUM1,npts,dummy,ios,ierr2,nz
  integer, external::NM
  complex(kind=dp)::RA(nterms),RB(nterms),rbb(nterms),rcc(nterms),&
       rd1(nterms),&
       rd2(nterms),rd3(nterms),RC(nterms,n_layers),&
       RD(nterms,n_layers),RF(nterms,n_layers),RG(nterms,n_layers),&
       RA_pw(nterms),RB_pw(nterms),RC_pw(nterms,n_layers),&
       RD_pw(nterms,n_layers),RF_pw(nterms,n_layers),RG_pw(nterms,n_layers) 
  real(kind=dp):: AX,ari,ari1,xx(n_layers),r_obs,theta_obs,phi_obs,&
       cart(3),sp(3),power_abs,power_sca,power_ext,factor,radius,power
  complex(kind=dp)::Ef(3),Ecart(3),Hf(3),Hcart(3),temp_field(3),Hcart_ref(3)
  complex(kind=dp)::E_the_dipole,E_phi_dipole,Ecart_ref(3),E_r_dipole
  complex(kind=dp)::const_dipole,const2
!!excitation amplitude of plane wave due to dipole at farfield
  real(kind=dp)::wave_no
  real(kind=dp)::new_cart(3),new_sp(3)
  real(kind=dp)::dummy_freq,tmp_real,tmp_cmplx
  complex(kind=dp)::factor_p
  real(kind=dp)::dtheta,dphi
  real(kind=dp)::theta,phi
  integer::count_phi,count_the,exc_mode
  
  if(ntheta/=1) then
     dtheta=(thetaend-thetastart)/(ntheta-1.0_dp);
  else
     dtheta=0
  end if
  
  if (nphi/=1) then
     dphi=(phimax-phimin)/(nphi-1.0_dp);
  else
     dphi=0;
  end if
  
  if(.not.plane_wave_exc) then
     wave_no=2.0_dp*pi/lambda
     arg_p=wave_no*r_p
     if(Jdipole_exc) then
        const_dipole=4.0_dp*wave_no**2*&
             vlite*1.0E-7_dp*exc_dipole
     else
        const_dipole=wave_no**2*&
             vlite*1.0E-7_dp*exc_dipole/(3600.0_dp*pi**2)
     end if
!!!const_dipole must be multiplied to the result finally
  end if
  


100 format(I5,6ES23.15)
101 format(6ES23.15)
  AX=1.0_dp/X
  xx(1) = x * aa(1)
  xx(n_l) = x
  do i = 2, n_l - 1
     xx(i) = 0.0_dp
     do j = 1, i
        xx(i) = xx(i) + aa(j)
     end do
     xx(i) = x * xx(i)
  end do

  ! d1(x), rd3(x), rc(x)
  NUM = NM(X)
  print*,'NUM is',NUM
  !NUM=30
 
  CALL aax(AX,NUM,d1x)
  CALL cd3x(X,NUM,d1x,rd3x,rcx)

  ari = abs(RI_n(1))
  do i = 2, n_l
     ari1 = abs(RI_n(i))
     if(ari1>ari) ari = ari1
  end do
  NUM2=NM(ari*X)

 ! NUM2=30
        
  CALL aa1(ri_n(1)*xx(1),NUM2,rd11)
  do i = 2, n_l
     
     CALL bcd(ri_n(i)*xx(i-1),NUM2,rd1,rd2,rd3,rbb,rcc)
     do j = 1, num2
        rrcc(i,j) = rcc(j)
        rrbb(i,j) = rbb(j)
        rrd1(i,j) = rd1(j)
        rrd2(i,j) = rd2(j)
        rrd3(i,j) = rd3(j)
     end do
     
     CALL bcd(ri_n(i)*xx(i),NUM2,rd1,rd2,rd3,rbb,rcc)
     do j = 1, num2
        srbb(i,j) = rbb(j)
        srd1(i,j) = rd1(j)
        srd2(i,j) = rd2(j)
     end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!Since rrbb(1,:), rrdd1(1,:), rrd2(:,1) are not used, here, 
!!!I used them for specifying the ..
!!!arguments corrsponding to ri_n(0)*xx(n_l) , 
!!!where ri_n(0)= ri_n of surrounding medium

  CALL bcd(ri_n(0)*xx(n_l),NUM2,rd1,rd2,rd3,rbb,rcc)
  do j = 1, num2
     rrcc(1,j) = rcc(j)
     rrbb(1,j) = rbb(j)
     rrd1(1,j) = rd1(j)
     rrd2(1,j) = rd2(j)
     rrd3(1,j) = rd3(j)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  radius=a(n_l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!*********************************************************!!!!  
!*******************Computation of coefficients****************!
!!!*********************************************************!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  CALL ABCDn1(NUM,NUM1,RA_pw,RB_pw,RC_pw,RD_pw,RF_pw,RG_pw)
  
  if(.not.plane_wave_exc) then
     if(r_p>radius) then
        !CALL ABCDn1(NUM,NUM1,RA,RB,RC,RD,RF,RG)
        RA=RA_pw
        RB=RB_pw
        RC=RC_pw
        RD=RD_pw
        RF=RF_pw
        RG=RG_pw
     else
        print*,"The case where dipole is inside sphere works for only one &
             layer, running the code with material properties of &
             innermost layer&
             and radius set to that of the outermost layer"
        CALL ABCD_dipole_inside(NUM,NUM1,RA,RB,RC,RD,RF,RG)
     end if
     
     allocate(cyr_p(0:num),cyi_p(0:num))
     if(r_p>radius) then
        call zbesh(arg_p,0.0_dp,0.5_dp,one,2,num+1,cyr_p,cyi_p,nz,ierr2)
        if (ierr2/=0) then
           print*,"error in the calculation of the hankel bessel with error &
                code:",ierr2
           stop
        end if
     else
        arg_p=arg_p*ri_n(1)
        tmp_real=real(arg_p)
        tmp_cmplx=aimag(arg_p)
        call zbesj(tmp_real,tmp_cmplx,0.5_dp,one,num+1,cyr_p,cyi_p,nz,ierr2)
        if (ierr2/=0) then
           print*,"error in the calculation of the hankel bessel with error &
                code:",ierr2
           stop
        end if
        const_dipole=const_dipole*ri_n(1)**3
     end if
     
     factor_p=sqrt((pi*0.5_dp)/arg_p)
     
     cyr_p=cyr_p*factor_p;
     cyi_p=cyi_p*factor_p;
     
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!*********************************************************!!!!  
!*******************Far field/RCS computation******************!
!!!*********************************************************!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   
  if(.not.plane_wave_exc) then
     print*,"Opening file E_farfields.out"
     open(unit=08,file='E_farfields.out',status='unknown')          
     if(Jdipole_exc) then
        const2=-cmplx(0.0_dp,1.0_dp,dp)*&
             wave_no*vlite*1.0E-7_dp*exc_dipole
     else
        const2=-cmplx(0.0_dp,1.0_dp,dp)*&
             wave_no*vlite*1.0E-7_dp*exc_dipole/14400.0_dp/pi
     end if

     if(ntheta==0.and.nphi==0) write(8,*) ""
     
     j=0
     do count_phi=1,nphi
        phi=((phimin+(count_phi-1)*dphi)*pi)/180.0_dp
        do count_the=1,ntheta
           j=j+1
           theta=((thetastart+(count_the-1)*dtheta)*pi)/180.0_dp

           sp(1)=1.0_dp
           sp(2)=theta
           sp(3)=phi

           theta_obs=sp(2)
           phi_obs=sp(3)

!!$!!!!!!!!!!!!!!
           
           !For the evaluation of E_theta
           new_cart(1)=-posit_dipole(1)*cos(theta_obs)*cos(phi_obs)-&
                posit_dipole(2)*cos(theta_obs)*sin(phi_obs)+posit_dipole(3)&
                *sin(theta_obs)
           new_cart(2)=-posit_dipole(1)*sin(phi_obs)+&
                posit_dipole(2)*cos(phi_obs)
           new_cart(3)=-posit_dipole(1)*sin(theta_obs)*cos(phi_obs)-&
                posit_dipole(2)*sin(theta_obs)*sin(phi_obs)-posit_dipole(3)&
                *cos(theta_obs)
           
!!$!!!!!!!!!!!!!
!!$              
           call cart2sph(new_cart,new_sp)
           r_obs=new_sp(1)
           theta_obs=new_sp(2)
           phi_obs=new_sp(3)
           call calc_field(0,RA_pw,RB_pw,RC_pw,RD_pw,RF_pw,RG_pw,&
                num,freq,r_obs,theta_obs,phi_obs,Ef,Hf)
           call sph2cart(theta_obs,phi_obs,Ef,Ecart)
           call sph2cart(theta_obs,phi_obs,Hf,Hcart)
           
           if(Mdipole_exc) then
              !        temp_field(1:3)=Ecart(1:3)
              Ecart(1:3)=-Hcart(1:3) 
              !        Hcart(1:3)=temp_field(1:3)
           end if
!!$
!!$!!!Mapping back to original coordinate system
!!$              
           theta_obs=sp(2)
           phi_obs=sp(3)
           Ecart_ref(1)=-Ecart(1)*cos(theta_obs)*cos(phi_obs)-&
                Ecart(2)*sin(phi_obs)-Ecart(3)*sin(theta_obs)*cos(phi_obs)
           Ecart_ref(2)=-Ecart(1)*cos(theta_obs)*sin(phi_obs)+&
                Ecart(2)*cos(phi_obs)-Ecart(3)*sin(theta_obs)*sin(phi_obs)
           Ecart_ref(3)=Ecart(1)*sin(theta_obs)-&
                Ecart(3)*cos(theta_obs)
!!!!The component of this in the direction of dipole is taken
           r_obs=sp(1)
           E_the_dipole=-const2*&
                dot_product(conjg(Ecart_ref),orient_dipole_cart)*&
                cmplx(cos(wave_no*r_obs),-sin(wave_no*r_obs),dp)/r_obs  
           
!!$              
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           r_obs=sp(1)
           theta_obs=sp(2)
           phi_obs=sp(3)
           !For the evaluation of E_phi
           
           new_cart(1)=posit_dipole(1)*(sin(phi_obs))-&
                posit_dipole(2)*cos(phi_obs)
           new_cart(2)=-posit_dipole(1)*cos(theta_obs)*cos(phi_obs)-&
                posit_dipole(2)*&
                cos(theta_obs)*sin(phi_obs)+posit_dipole(3)*sin(theta_obs)
           new_cart(3)=-posit_dipole(1)*sin(theta_obs)*cos(phi_obs)-&
                posit_dipole(2)*sin(theta_obs)*sin(phi_obs)-posit_dipole(3)&
                *cos(theta_obs)
!!$              
!!$!!!!!!!!!!!!!
!!$              
           call cart2sph(new_cart,new_sp)
           r_obs=new_sp(1)
           theta_obs=new_sp(2)
           phi_obs=new_sp(3)
           
           call calc_field(0,RA_pw,RB_pw,RC_pw,RD_pw,RF_pw,RG_pw,&
                num,freq,r_obs,theta_obs,phi_obs,Ef,Hf)
           call sph2cart(theta_obs,phi_obs,Ef,Ecart)
           call sph2cart(theta_obs,phi_obs,Hf,Hcart)
           
           if(Mdipole_exc) then
              !         temp_field(1:3)=Ecart(1:3)
              Ecart(1:3)=-Hcart(1:3) 
              !         Hcart(1:3)=temp_field(1:3)
           end if
           
!!!!Mapping back to original coordinate system
           theta_obs=sp(2)
           phi_obs=sp(3)
           Ecart_ref(1)=Ecart(1)*sin(phi_obs)-&
                Ecart(2)*cos(theta_obs)*sin(phi_obs)&
                -Ecart(3)*sin(theta_obs)*cos(phi_obs)
           Ecart_ref(2)=-Ecart(1)*cos(phi_obs)-&
                Ecart(2)*cos(theta_obs)*sin(phi_obs)-&
                Ecart(3)*sin(theta_obs)*sin(phi_obs)
           Ecart_ref(3)=Ecart(2)*sin(theta_obs)-&
                Ecart(3)*cos(theta_obs)
!!!!The component of this in the direction of dipole is taken
           r_obs=sp(1)
           
           E_phi_dipole=-const2*&
                dot_product(conjg(Ecart_ref),orient_dipole_cart)*&
                cmplx(cos(wave_no*r_obs),-sin(wave_no*r_obs),dp)/r_obs 
           E_r_dipole=cmplx(0.0_dp,0.0_dp,dp)  
!!$              
           write(8,101) theta*180.0_dp/pi,&
                phi*180.0_dp/pi,abs(E_the_dipole),abs(E_phi_dipole)
        end do
     end do
     print*,"Closing Files ... E_farfields.out"
     close(8,status='keep')
  else
     call calc_rcs(RA_pw,RB_pw,num,freq)
  end if
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!*********************************************************!!!!  
!*******************Absorbed power computation*****************!
!!!*********************************************************!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  print*,"Opening file abs_power.out"
  open(unit=12,file='abs_power.out',status='unknown') 
  if(plane_wave_exc) then
     do i=1,n_l
        call Mie_power(i,num,RA_pw,RB_pw,RC_pw,RD_pw,RF_pw,RG_pw,power)
     end do
  else
      do i=1,n_l
        call Mie_power(i,num,RA,RB,RC,RD,RF,RG,power)
     end do
  end if
  print*,"Closing Files ... abs_power.out"
  close(12,status='keep')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!*********************************************************!!!!  
!***********************Near field computation*****************!
!!!*********************************************************!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  print*,"Opening files Efields.out and Hfields.out "
  open(unit=6,file='points.inp',status='old',access='sequential')
  open(unit=10,file='Efields.out',status='unknown')
  open(unit=11,file='Hfields.out',status='unknown')

  if(plane_wave_exc) then

     RA=RA_pw; RB=RB_pw; RC=RC_pw; RD=RD_pw; RF=RF_pw; RG=RG_pw
     exc_mode=0
  else
     exc_mode=1
  end if


  npts=0
  read(6,*, iostat=ios) dummy
  do
     if (ios == 0) then
        npts = npts + 1
        read(6, *, iostat=ios) dummy
     else
        exit
     end if
  end do
  npts=npts
  print*,'n points',npts
  rewind(6)     
  
  if(npts==0) then
     write(10,*) ""
     write(11,*) ""
  end if


  if(plane_wave_exc) then

     do j=1,npts
        read(6,*) cart(1),cart(2),cart(3)
      
!!$!!!!!!!!!!!!!!
!For the rotation of axes and get results for the reference case 
!
        new_cart(1)=-ethetad*(cart(1)*cos(the_angle)*cos(phi_angle)+&
             cart(2)*cos(the_angle)*sin(phi_angle)-cart(3)*&
             sin(the_angle))
        
        new_cart(1)=new_cart(1)+ephid*(cart(1)*sin(phi_angle)-&
             cart(2)*cos(phi_angle))
        
        new_cart(2)=ethetad*(-cart(1)*sin(phi_angle)+&
             cart(2)*cos(phi_angle))
        new_cart(2)= new_cart(2)-ephid*&
             (cart(1)*cos(the_angle)*cos(phi_angle)+&
             cart(2)*cos(the_angle)*sin(phi_angle)-cart(3)*sin(the_angle))
        
        new_cart(3)=-cart(1)*sin(the_angle)*cos(phi_angle)-&
             cart(2)*sin(the_angle)*sin(phi_angle)-cart(3)&
             *cos(the_angle)
!!$!!!!!!!!!!!!!
        call cart2sph(new_cart,new_sp)
        
        r_obs=new_sp(1)
        theta_obs=pi-new_sp(2)
        phi_obs=pi-new_sp(3)

        call calc_field(exc_mode,RA,RB,RC,RD,RF,RG,num,freq,r_obs,theta_obs,&
             phi_obs,Ef,Hf)
        Ef(1:3)=Ef(1:3)
        Hf(1:3)=Hf(1:3)
        call sph2cart(theta_obs,phi_obs,Ef,Ecart)
        call sph2cart(theta_obs,phi_obs,Hf,Hcart)
        
        Ecart_ref(1)=-ethetad*(Ecart(1)*cos(the_angle)*cos(phi_angle)+&
             Ecart(2)*sin(phi_angle))-Ecart(3)*sin(the_angle)*cos(phi_angle)
        Ecart_ref(2)=ethetad*(-Ecart(1)*cos(the_angle)*sin(phi_angle)+&
             Ecart(2)*cos(phi_angle))-Ecart(3)*sin(the_angle)*sin(phi_angle)
        Ecart_ref(3)=ethetad*Ecart(1)*sin(the_angle)-&
             Ecart(3)*cos(the_angle)
        
        
        Ecart_ref(1)=Ecart_ref(1)+ephid*(Ecart(1)*sin(phi_angle)-&
             Ecart(2)*cos(the_angle)*sin(phi_angle))
        Ecart_ref(2)= Ecart_ref(2)-ephid*(Ecart(1)*cos(phi_angle)+&
             Ecart(2)*cos(the_angle)*sin(phi_angle))
        Ecart_ref(3)= Ecart_ref(3)+ephid*Ecart(2)*sin(the_angle)


        Hcart_ref(1)=-ethetad*(Hcart(1)*cos(the_angle)*cos(phi_angle)+&
             Hcart(2)*sin(phi_angle))-Hcart(3)*sin(the_angle)*cos(phi_angle)
        Hcart_ref(2)=ethetad*(-Hcart(1)*cos(the_angle)*sin(phi_angle)+&
             Hcart(2)*cos(phi_angle))-Hcart(3)*sin(the_angle)*sin(phi_angle)
        Hcart_ref(3)=ethetad*Hcart(1)*sin(the_angle)-&
             Hcart(3)*cos(the_angle)
        
        
        Hcart_ref(1)=Hcart_ref(1)+ephid*(Hcart(1)*sin(phi_angle)-&
             Hcart(2)*cos(the_angle)*sin(phi_angle))
        Hcart_ref(2)= Hcart_ref(2)-ephid*(Hcart(1)*cos(phi_angle)+&
             Hcart(2)*cos(the_angle)*sin(phi_angle))
        Hcart_ref(3)= Hcart_ref(3)+ephid*Hcart(2)*sin(the_angle)
        
        
        write(10,100) j,-real(Ecart_ref(1)),-aimag(Ecart_ref(1)),&
             real(Ecart_ref(2)),aimag(Ecart_ref(2)),&
             -real(Ecart_ref(3)),-aimag(Ecart_ref(3))
        
        write(11,100) j,-real(Hcart_ref(1)),-aimag(Hcart_ref(1)),&
             real(Hcart_ref(2)),aimag(Hcart_ref(2)),&
             -real(Hcart_ref(3)),-aimag(Hcart_ref(3))
        
     end do


  else
     do j=1,npts
        read(6,*) cart(1),cart(2),cart(3)
        call cart2sph(cart,sp)
        r_obs=sp(1)
        theta_obs=sp(2)
        phi_obs=sp(3)
        
        call calc_field(exc_mode,RA,RB,RC,RD,RF,RG,num,freq,r_obs,theta_obs,&
             phi_obs,Ef,Hf)
        Ef(1:3)=Ef(1:3)*const_dipole
        Hf(1:3)=Hf(1:3)*const_dipole
        call sph2cart(theta_obs,phi_obs,Ef,Ecart)
        call sph2cart(theta_obs,phi_obs,Hf,Hcart)
        
        if(Mdipole_exc) then
           temp_field(1:3)=Ecart(1:3)
           Ecart(1:3)=-Hcart(1:3) 
           Hcart(1:3)=temp_field(1:3)
        end if
        
        
        write(10,100) j,real(Ecart(1)),aimag(Ecart(1)),&
             real(Ecart(2)),aimag(Ecart(2)),&
             real(Ecart(3)),aimag(Ecart(3))
        
        write(11,100) j,real(Hcart(1)),aimag(Hcart(1)),&
             real(Hcart(2)),aimag(Hcart(2)),&
             real(Hcart(3)),aimag(Hcart(3))
        
     end do
  end if
  
  if(.not.plane_wave_exc)   deallocate(cyr_p,cyi_p)
  print*,"Closing Files ... Efields.out"
  print*,"Closing Files ... Hfields.out"
  close(6,status='keep')
  close(10,status='keep')
  close(11,status='keep')
    
  return
END SUBROUTINE shexqn1

!--------------------------------------------------------------------
! NM-auxiliary function for AA1 & BESSEL
!    (number NM is calculated using X)
! see: Trudy Astronom. Observ. LGU V.28,P.14,1971
!    for X>1 value of NM was raised
! August 1989, AO LGU
!--------------------------------------------------------------------
integer FUNCTION NM(X)
  use mie_parameters
  implicit none
  real(kind=dp),intent(in)::x

  IF(X<1.0_dp) then
     NM=7.5_dp*X+9.0_dp 
  else IF(X.GT.100) then 
     NM=1.0625_dp*X+28.5_dp 
  else
     NM=1.25_dp*X+15.5_dp
  end if
END FUNCTION NM
!--------------------------------------------------------------------
 ! AA1-subroutine for calculations of the ratio of the derivative
 !    to the function for Bessel functions of half order with
!   the complex argument: J'(N)/J(N).
!    The calculations are given by the recursive expression
!    ``from top to bottom'' beginning from N=NUM.
!    RU-array of results.
!    A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
!    RI - complex refractive index.
! August 1989, AO LGU
!--------------------------------------------------------------------
SUBROUTINE AA1(Rx,NUM,RU)
  use mie_parameters
  IMPLICIT none
  
  complex(kind=dp),intent(in)::Rx
  integer,intent(in)::num
  complex(kind=dp),DIMENSION(num),intent(out)::RU
  complex(kind=dp)::S,S1
  integer::num1,j,i,i1

  S = 1.0_dp/ rx
  RU(NUM)=(NUM+1.0_dp)*S
  NUM1=NUM-1
  do J=1,NUM1
     I=NUM-J
     I1=I+1
     S1=I1*S
     RU(I)=S1-1.0_dp/(RU(I1)+S1)
  end do
  RETURN
END SUBROUTINE AA1
!--------------------------------------------------------------------
! AAx-subroutine for calculations of the ratio of the derivative
!    to the function for Bessel functions of half order with
!    the real argument: J'(N)/J(N).
!    The calculations are given by the recursive expression
!    ``from top to bottom'' beginning from N=NUM.
!    RU-array of results.
!    A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
! March 1999, AI SPbU
!--------------------------------------------------------------------
SUBROUTINE AAx(A,NUM,RU)
  use mie_parameters
  IMPLICIT none

  real(kind=dp),intent(in)::a
  integer,intent(in)::num
  real(kind=dp),DIMENSION(num),intent(out):: RU
  integer::num1,i,i1,j
  real(kind=dp)::s1

  RU(NUM)=(NUM+1.0_dp)*a
  NUM1=NUM-1
  DO J=1,NUM1
     I=NUM-J
     I1=I+1
     S1=I1*a
     RU(I)=S1-1.0_dp/(RU(I1)+S1)
  END DO
  RETURN
END SUBROUTINE AAx
!--------------------------------------------------------------------
! CD3X-subroutine for calculations of the ratio of the derivative
!    to the function for Riccati-Bessel functions of half order with
!    the real argument: zeta'(N)/zeta(N)
!    and the ratio of functions: psi(N)/zeta(N).
!    The calculations are given by the recursive expression
!    ``from bottom to top'' beginning from N=0.
!    rd3x, rcx-arrays of results.
!    X - size parameter
! March 1999, AI SPbU
!--------------------------------------------------------------------
SUBROUTINE cd3x(x,NUM,d1x,rd3x,rcx)
  use mie_parameters
  IMPLICIT none

  real(kind=dp),intent(in)::x
  integer,intent(in)::num
  real(kind=dp),dimension(num),intent(in)::d1x
  complex(kind=dp),dimension(num),intent(out)::rd3x, rcx
  real(kind=dp)::ax,a1
  integer::i
  complex(kind=dp)::s1,rd30,rxy,rc0

  S1 = (0.0_dp,1.0_dp)
  ax = 1.0_dp / x
      
  rd30 = s1
  rxy = cos(2.0_dp*x) + s1 * sin(2.0_dp*x)
  rc0 = -(1.0_dp - rxy) / (2.0_dp * rxy)
  rd3x(1) = -ax + 1.0_dp / (ax - rd30)
  rcx(1) = rc0 * (ax + rd3x(1)) / (ax + d1x(1))
  
  do i = 2, NUM
     a1 = I * ax
     rd3x(i) = -a1 + 1.0_dp / (a1 - rd3x(i-1))
     rcx(i) = rcx(i-1) * (a1 + rd3x(i)) / (a1 + d1x(i))
  end do
  
  RETURN
END SUBROUTINE cd3x
!--------------------------------------------------------------------
! BCD-subroutine for calculations of the ratios of the derivative
!    to the function for Riccati-Bessel functions of half order with
!    the complex argument: psi'(N)/psi(N) and khi'(N)/khi(N)
!    and the ratios of functions: psi(N)/khi(N).
!    The calculations are given by the recursive expression
!    ``from bottom to top'' beginning from N=0.
!    rd1, rd2, rbb, rcc-arrays of results.
!    rx - (refr. index) * (size parameter)
! March 1999, AI SPbU
!--------------------------------------------------------------------
SUBROUTINE bcd(rx,NUM,rd1,rd2,rd3,rbb,rcc)
  use mie_parameters
  IMPLICIT none

  complex(kind=dp),intent(in)::rx
  integer,intent(in)::num
  complex(kind=dp),DIMENSION(num),intent(out):: rd1, rd2, rd3, rbb, rcc
  !complex(kind=dp),dimension(nterms)::rcc
  real(kind=dp)::x,y
  complex(kind=dp)::s1,rd30,rxy,rc0,rb0,rx1,r1
  integer::i

  S1 = (0.0_dp,1.0_dp)
  x=real(rx)
  y=aimag(rx)
 ! y=0.0_dp
  rx1 = 1.0_dp / rx

  CALL aa1(rx,NUM,rd1)

! n = 0
  rd30 = s1
  rxy = (cos(2.0_dp*x) + s1 * sin(2.0_dp*x))*exp(-2.0_dp*y)
  rc0 = -(1.0_dp - rxy) / (2.0_dp * rxy)
  rb0 = s1 * (1.0_dp - rxy) / (1.0_dp + rxy)
! n = 1
  rd3(1) = -rx1 + 1.0_dp / (rx1 - rd30)
  rcc(1) = rc0 * (rx1 + rd3(1)) / (rx1 + rd1(1))
  rd2(1) = (rcc(1) * rd1(1) - rd3(1)) / (rcc(1) - 1d0)
  rbb(1) = rb0 * (rx1 + rd2(1)) / (rx1 + rd1(1))
  
  DO i = 2, NUM
     r1 = I * rx1
     rd3(i) = -r1 + 1.0_dp / (r1 - rd3(i-1))
     rcc(i) = rcc(i-1) * (r1 + rd3(i)) / (r1 + rd1(i))
     rd2(i) = (rcc(i) * rd1(i) - rd3(i)) / (rcc(i) - 1.0_dp)
     rbb(i) = rbb(i-1) * (r1 + rd2(i)) / (r1 + rd1(i))
  end do
  
  RETURN
END SUBROUTINE bcd
!--------------------------------------------------------------------
! ABn1-subroutine for calculations of the complex coefficients
!    A(N), B(N) for n-layered spheres.
!    n_l - number of layers
!    RI_n(i) - complex refractive indices for innermost layer (1),
!    layer2, ... (i = 1, n_l)
!    The coefficients are calculated up to the number NUM1.LE.NUM,
!    for which |A(N)**2+B(N)**2|.LE.10**(-40)
!    RA-array of coefficients A(N), RB-array of coefficients B(N)
! March 1999, AI SPbU
!--------------------------------------------------------------------
SUBROUTINE ABCD_dipole_inside(NUM,NUM1,RA,RB,RC,RD,RF,RG)
  use  mie_variables
  IMPLICIT none

  integer,intent(in)::num
  integer,intent(out)::num1
  complex(kind=dp),dimension(nterms),intent(out)::RA, RB
  complex(kind=dp),dimension(nterms,n_layers),intent(out)::RC,RD,RF,RG
  real(kind=dp),allocatable::cyr(:),cyi(:),crwrk(:),ciwrk(:)
  integer::i,j,ierr,nz
  complex(kind=dp)::arg1,arg2
  complex(kind=dp),allocatable,dimension(:)::bessel1a,hankel1a,&
       bessel1da,hankel1da,bessel1b,hankel1b,&
       bessel1db,hankel1db!,bessel2db,bessel2b,bessel2da,bessel2a
  complex(kind=dp)::factor

  RF=cmplx(0.0_dp,0.0_dp,dp)
  RG=cmplx(0.0_dp,0.0_dp,dp)
  RC=cmplx(0.0_dp,0.0_dp,dp)
  RD=cmplx(0.0_dp,0.0_dp,dp)

  allocate(bessel1a(1:num),&
       hankel1a(1:num),bessel1da(1:num),&
       hankel1da(1:num),bessel1b(1:num),&
       hankel1b(1:num),bessel1db(1:num),&
       hankel1db(1:num),cyr(0:num),cyi(0:num),crwrk(0:num),ciwrk(0:num))
  
  arg1=ri_n(1)*a(n_l)*2.0_dp*pi/lambda
  arg2=ri_n(0)*a(n_l)*2.0_dp*pi/lambda

  factor=sqrt((pi*0.5_dp)*arg1)
  call zbesh(real(arg1),aimag(arg1),0.5_dp,1,1,num+1,cyr,cyi,nz,ierr)
  
  hankel1a(1:num)=cmplx(cyr(1:num),cyi(1:num),dp)*factor
  hankel1da(1:num)=factor*cmplx(cyr(0:num-1),cyi(0:num-1),dp)-&
       (/(j,j=1,num)/)*(hankel1a(1:num)/arg1)
  call zbesj(real(arg1),aimag(arg1),0.5_dp,1,num+1,cyr,cyi,nz,ierr)
  
  bessel1a(1:num)=cmplx(cyr(1:num),cyi(1:num),dp)*factor
  bessel1da(1:num)=cmplx(cyr(0:num-1),cyi(0:num-1),dp)*factor-&
       (/(j,j=1,num)/)*(bessel1a(1:num)/arg1)
  
!!!!!!!  
  factor=sqrt((pi*0.5_dp)*arg2)
  call zbesh(real(arg2),aimag(arg2),0.5_dp,1,1,num+1,cyr,cyi,nz,ierr)
  
  hankel1b(1:num)=cmplx(cyr(1:num),cyi(1:num),dp)*factor
  hankel1db(1:num)=cmplx(cyr(0:num-1),cyi(0:num-1),dp)*factor-&
       (/(j,j=1,num)/)*(hankel1b(1:num)/arg2)
  call zbesj(real(arg2),aimag(arg2),0.5_dp,1,num+1,cyr,cyi,nz,ierr)
  
  bessel1b(1:num)=cmplx(cyr(1:num),cyi(1:num),dp)*factor
  bessel1db(1:num)=cmplx(cyr(0:num-1),cyi(0:num-1),dp)*factor-&
       (/(j,j=1,num)/)*(bessel1b(1:num)/arg2)
  
  call zbesy(real(arg2),aimag(arg2),0.5_dp,1,num+1,cyr,cyi,nz,&
       crwrk,ciwrk,ierr)
  

  loop:DO i = 1, NUM
   
 !--------
 ! calculations of a(n), b(n) for dipole outside sphere:
     
!!$     RA(i) =ri_n(1)*conjg(eps(0))*bessel1a(i)*bessel1db(i)-&
!!$          ri_n(0)*conjg(eps(1))*bessel1da(i)*bessel1b(i)
!!$     RA(i)=RA(i)/(ri_n(0)*conjg(eps(1))*bessel1da(i)*hankel1b(i)-&
!!$          ri_n(1)*conjg(eps(0))*bessel1a(i)*hankel1db(i))
!!$     RB(i) =ri_n(1)*conjg(eps(0))*bessel1da(i)*bessel1b(i)-&
!!$          ri_n(0)*conjg(eps(1))*bessel1db(i)*bessel1a(i)
!!$     RB(i)=RB(i)/(ri_n(0)*conjg(eps(1))*bessel1a(i)*hankel1db(i)-&
!!$          ri_n(1)*conjg(eps(0))*bessel1da(i)*hankel1b(i))    
!!$
!!$     RA(i)=-RA(i)
!!$     RB(i)=-RB(i)
!!$
!!$     RC(i,1) =ri_n(1)*(conjg(eps(1))*hankel1db(i)*bessel1b(i)-&
!!$          conjg(eps(1))*hankel1b(i)*bessel1db(i))
!!$     RC(i,1)=RC(i,1)/(ri_n(0)*conjg(eps(1))*bessel1a(i)*hankel1db(i)-&
!!$          ri_n(1)*conjg(eps(0))*bessel1da(i)*hankel1b(i))
!!$     RD(i,1)=(conjg(eps(1))*hankel1b(i)*bessel1db(i)-&
!!$          conjg(eps(1))*hankel1db(i)*bessel1b(i))*ri_n(1)
!!$     RD(i,1)=RD(i,1)/(ri_n(0)*conjg(eps(1))*bessel1da(i)*hankel1b(i)-&
!!$          ri_n(1)*conjg(eps(0))*bessel1a(i)*hankel1db(i))



   ! calculations of a(n), b(n) for dipole inside sphere:
!!$       
     RA(i) =ri_n(0)*conjg(eps(0))*bessel1a(i)*hankel1da(i)-&
          ri_n(0)*conjg(eps(0))*bessel1da(i)*hankel1a(i)
     RA(i)=RA(i)/(ri_n(1)*conjg(eps(0))*bessel1a(i)*hankel1db(i)-&
          ri_n(0)*conjg(eps(1))*bessel1da(i)*hankel1b(i))

     RB(i) =ri_n(0)*conjg(eps(0))*hankel1a(i)*bessel1da(i)-&
          ri_n(0)*conjg(eps(0))*hankel1da(i)*bessel1a(i)
     RB(i)=RB(i)/(ri_n(1)*conjg(eps(0))*bessel1da(i)*hankel1b(i)-&
          ri_n(0)*conjg(eps(1))*bessel1a(i)*hankel1db(i))


     RA(i)=-RA(i)
     RB(i)=-RB(i)


     RC(i,1)=ri_n(0)*conjg(eps(1))*hankel1db(i)*hankel1a(i)-&
          ri_n(1)*conjg(eps(0))*hankel1b(i)*hankel1da(i)
     RC(i,1)=RC(i,1)/(ri_n(1)*conjg(eps(0))*bessel1da(i)*hankel1b(i)-&
          ri_n(0)*conjg(eps(1))*bessel1a(i)*hankel1db(i))

     RD(i,1) =ri_n(0)*conjg(eps(1))*hankel1b(i)*hankel1da(i)-&
          ri_n(1)*conjg(eps(0))*hankel1db(i)*hankel1a(i)
     RD(i,1)=RD(i,1)/(ri_n(1)*conjg(eps(0))*bessel1a(i)*hankel1db(i)-&
          ri_n(0)*conjg(eps(1))*bessel1da(i)*hankel1b(i))
   
     do j=2,n_l
        RC(i,j)=RC(i,1)
        RD(i,j)=RD(i,1)
     end do

     if(abs(RA(i))+abs(RB(i))<=1D-40) then
        NUM1=i
        exit loop
      end if
!     NUM1=NUM
   end do loop
   
   RETURN
 END SUBROUTINE ABCD_dipole_inside


SUBROUTINE ABCDn1(NUM,NUM1,RA,RB,RC,RD,RF,RG)
  use  mie_variables
  IMPLICIT none

  integer,intent(in)::num
  integer,intent(out)::num1
  complex(kind=dp),dimension(nterms),intent(out)::RA, RB
  complex(kind=dp),dimension(nterms,n_layers),intent(out)::RC,RD,RF,RG
  real(kind=dp),allocatable::cyr(:),cyi(:),crwrk(:),ciwrk(:)
  integer::i,j,ierr,nz
  complex(kind=dp),allocatable,dimension(:)::arg1,arg2
  complex(kind=dp),allocatable,dimension(:,:)::bessel1a,bessel2a,hankel1a,&
       bessel1da,bessel2da,hankel1da,bessel1b,bessel2b,hankel1b,&
       bessel1db,bessel2db,hankel1db
  complex(kind=dp)::factor
  complex(kind=dp)::sa(n_layers), sha(n_layers),&
       sb(n_layers), shb(n_layers)
  complex(kind=dp)::dual_epsmu(0:n_layers)
  
  
  if(Mdipole_exc) then
     dual_epsmu=eps
  else
     dual_epsmu=mu
  end if

  allocate(arg1(1:n_l),arg2(1:n_l),cyr(0:num),cyi(0:num),&
       bessel1a(1:num,1:n_l),bessel2a(1:num,1:n_l),&
       hankel1a(1:num,1:n_l),bessel1da(1:num,1:n_l),bessel2da(1:num,1:n_l),&
       hankel1da(1:num,1:n_l),bessel1b(1:num,1:n_l),bessel2b(1:num,1:n_l),&
       hankel1b(1:num,1:n_l),bessel1db(1:num,1:n_l),bessel2db(1:num,1:n_l),&
       hankel1db(1:num,1:n_l),crwrk(0:num),ciwrk(0:num))
 
  do i=1,n_l-1
     arg1(i)=ri_n(i)*a(i)
     arg2(i)=ri_n(i+1)*a(i)
  end do
  arg1(n_l)=ri_n(n_l)*a(n_l)
  arg2(n_l)=ri_n(0)*a(n_l)
  
  arg1=arg1*2.0_dp*pi/lambda
  arg2=arg2*2.0_dp*pi/lambda
  
  do i=1,n_l
     
     factor=sqrt((pi*0.5_dp)*arg1(i))
     call zbesh(real(arg1(i)),aimag(arg1(i)),0.5_dp,1,1,num+1,cyr,cyi,nz,ierr)
     
     hankel1a(1:num,i)=cmplx(cyr(1:num),cyi(1:num),dp)*factor
     hankel1da(1:num,i)=factor*cmplx(cyr(0:num-1),cyi(0:num-1),dp)-&
          (/(j,j=1,num)/)*(hankel1a(1:num,i)/arg1(i))
     call zbesj(real(arg1(i)),aimag(arg1(i)),0.5_dp,1,num+1,cyr,cyi,nz,ierr)

     bessel1a(1:num,i)=cmplx(cyr(1:num),cyi(1:num),dp)*factor
     bessel1da(1:num,i)=cmplx(cyr(0:num-1),cyi(0:num-1),dp)*factor-&
          (/(j,j=1,num)/)*(bessel1a(1:num,i)/arg1(i))
     
     call zbesy(real(arg1(i)),aimag(arg1(i)),0.5_dp,1,num+1,cyr,cyi,nz,&
          crwrk,ciwrk,ierr)
     
     bessel2a(1:num,i)=-cmplx(cyr(1:num),cyi(1:num),dp)*factor
     bessel2da(1:num,i)=(-cmplx(cyr(0:num-1),cyi(0:num-1),dp)*factor-&
          (/(j,j=1,num)/)*(bessel2a(1:num,i)/arg1(i)))
     
!!!!!!!  
     factor=sqrt((pi*0.5_dp)*arg2(i))
     call zbesh(real(arg2(i)),aimag(arg2(i)),0.5_dp,1,1,num+1,cyr,cyi,nz,ierr)
     
     hankel1b(1:num,i)=cmplx(cyr(1:num),cyi(1:num),dp)*factor
     hankel1db(1:num,i)=cmplx(cyr(0:num-1),cyi(0:num-1),dp)*factor-&
          (/(j,j=1,num)/)*(hankel1b(1:num,i)/arg2(i))
     call zbesj(real(arg2(i)),aimag(arg2(i)),0.5_dp,1,num+1,cyr,cyi,nz,ierr)
     
     bessel1b(1:num,i)=cmplx(cyr(1:num),cyi(1:num),dp)*factor
     bessel1db(1:num,i)=cmplx(cyr(0:num-1),cyi(0:num-1),dp)*factor-&
          (/(j,j=1,num)/)*(bessel1b(1:num,i)/arg2(i))
     
     call zbesy(real(arg2(i)),aimag(arg2(i)),0.5_dp,1,num+1,cyr,cyi,nz,&
          crwrk,ciwrk,ierr)
     
     bessel2b(1:num,i)=-cmplx(cyr(1:num),cyi(1:num),dp)*factor
     bessel2db(1:num,i)=(-cmplx(cyr(0:num-1),cyi(0:num-1),dp)*factor-&
          (/(j,j=1,num)/)*(bessel2b(1:num,i)/arg2(i)))
  end do

  loop:DO I = 1, num
     sa(1) = (0.0_dp,0.0_dp)
     sha(1) = rd11(i)
     sb(1) = (0.0_dp,0.0_dp)
     shb(1) = rd11(i)
     !--------
     DO j = 2, n_l
        if(abs(ri_n(j)*sha(j-1)-ri_n(j-1)*rrd2(j,i))==0.0_dp) then
           sa(j) = rrbb(j,i) * (conjg(dual_epsmu(j-1))*ri_n(j) * sha(j-1) - &
                conjg(dual_epsmu(j))*ri_n(j-1) * rrd1(j,i))&
                / (conjg(dual_epsmu(j-1))*ri_n(j) * sha(j-1) - &
                conjg(dual_epsmu(j))*ri_n(j-1) * rrd2(j,i) + 1d-30)
        else
           sa(j) = rrbb(j,i) * (conjg(dual_epsmu(j-1))*ri_n(j) * sha(j-1) - &
                conjg(dual_epsmu(j))*ri_n(j-1) * rrd1(j,i))&
                / (conjg(dual_epsmu(j-1))*ri_n(j) * sha(j-1) - &
                conjg(dual_epsmu(j))*ri_n(j-1) * rrd2(j,i))
        end if
        
        if(abs(ri_n(j)*shb(j-1)-ri_n(j-1)*rrd2(j,i))==0.0_dp) then
           sb(j) = rrbb(j,i) * (conjg(dual_epsmu(j))*ri_n(j-1) * shb(j-1) - &
                conjg(dual_epsmu(j-1))*ri_n(j) * rrd1(j,i))&
                / (conjg(dual_epsmu(j))*ri_n(j-1) * shb(j-1) - &
                conjg(dual_epsmu(j-1))*ri_n(j) * rrd2(j,i) + 1d-30)
        else
           sb(j) = rrbb(j,i) * (conjg(dual_epsmu(j))*ri_n(j-1) * shb(j-1) - &
                conjg(dual_epsmu(j-1))*ri_n(j) * rrd1(j,i))&
                / (conjg(dual_epsmu(j))*ri_n(j-1) * shb(j-1) - &
                conjg(dual_epsmu(j-1))*ri_n(j) * rrd2(j,i))
        end if

        if(abs(srbb(j,i) - sa(j))==0.0_dp) then
           sha(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sa(j))&
                - sa(j) * srd2(j,i) / (srbb(j,i) - sa(j) + 1d-30)
        else
           sha(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sa(j))&
                - sa(j) * srd2(j,i) / (srbb(j,i) - sa(j))
        end if
        
        if(abs(srbb(j,i) - sb(j))==0.0_dp) then
           shb(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sb(j))&
                - sb(j) * srd2(j,i) / (srbb(j,i) - sb(j) + 1d-30)
        else
           shb(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sb(j))&
                - sb(j) * srd2(j,i) / (srbb(j,i) - sb(j))
        end if

     end do
 !--------
 ! calculations of a(n), b(n)
     
     RA(I) = rrcc(1,i) * (ri_n(0)*conjg(dual_epsmu(n_l))*sha(n_l) - &
          ri_n(n_l)*conjg(dual_epsmu(0))* rrd1(1,i)) /&
          (ri_n(0)*conjg(dual_epsmu(n_l))*sha(n_l) - &
          ri_n(n_l)*conjg(dual_epsmu(0)) * rrd3(1,i))
     
     RB(I) =  rrcc(1,i) * (ri_n(n_l) *conjg(dual_epsmu(0))* shb(n_l) - &
          ri_n(0)*conjg(dual_epsmu(n_l))* rrd1(1,i)) /&
          (ri_n(n_l) *conjg(dual_epsmu(0))*shb(n_l) - &
          ri_n(0)*conjg(dual_epsmu(n_l))*rrd3(1,i))
     
!!!!!!!!!!!!!Solution for d's,g's
     
     RD(I,n_l)=bessel1b(I,n_l)-RA(I)*hankel1b(I,n_l)
     RD(I,n_l)=RD(I,n_l)/(bessel1a(I,n_l)-sa(n_l)*bessel2a(I,n_l))
     RD(I,n_l)= RD(I,n_l)*(conjg(dual_epsmu(n_l))/conjg(dual_epsmu(0)))
     RG(I,n_l)=sa(n_l)*RD(I,n_l)

     do j=n_l-1,2,-1
        RD(I,j)=(RD(I,j+1)*bessel1b(I,j)-RG(I,j+1)*bessel2b(I,j))
        RD(I,j)=RD(I,j)/((bessel1a(I,j)-sa(j)*bessel2a(I,j)))
        RD(I,j)= RD(I,j)*(conjg(dual_epsmu(j-1))/conjg(dual_epsmu(j)))
        RG(I,j)=sa(j)*RD(I,j)
     end do
     RD(I,1)=ri_n(1)*(RD(I,2)*bessel1db(I,1)-RG(I,2)*bessel2db(I,1))
     RD(I,1)=RD(I,1)/(ri_n(2)*bessel1da(I,1))
     RG(I,1)=cmplx(0.0_dp,0.0_dp,dp)


!!!!!!!!!!!!!Solution for c's,f's
     
     RC(I,n_l)=ri_n(n_l)*(bessel1b(I,n_l)-RB(I)*hankel1b(I,n_l))
     RC(I,n_l)=RC(I,n_l)/(bessel1a(I,n_l)-sb(n_l)*bessel2a(I,n_l))
     RC(I,n_l)=RC(I,n_l)/ri_n(0)
     RF(I,n_l)=sb(n_l)*RC(I,n_l)
   
     do j=n_l-1,2,-1
        RC(I,j)=(RC(I,j+1)*bessel1b(I,j)-RF(I,j+1)*bessel2b(I,j))
        RC(I,j)=RC(I,j)*ri_n(j)/((bessel1a(I,j)-sb(j)*bessel2a(I,j))*ri_n(j+1))
        RF(I,j)=sb(j)*RC(I,j)
     end do
     RC(I,1)=ri_n(1)*(RC(I,2)*bessel1b(I,1)-RF(I,2)*bessel2b(I,1))
     RC(I,1)=RC(I,1)/(ri_n(2)*bessel1a(I,1))
     RF(I,1)=cmplx(0.0_dp,0.0_dp,dp)
     
     
     if(abs(RA(I))+abs(RB(I))<=1D-40) then
        NUM1=I
        exit loop
     end if
!     NUM1=NUM
  end do loop

  deallocate(arg1,arg2,bessel1a,bessel2a,hankel1a,&
       bessel1da,bessel2da,hankel1da,bessel1b,bessel2b,hankel1b,&
       bessel1db,bessel2db,hankel1db)
  RETURN
END SUBROUTINE ABCDn1


subroutine qq(n,m,theta_in,pp,cc,dd)
  use mie_parameters
  implicit none

  integer,intent(in)::n,m
  real(kind=dp),intent(in)::theta_in
  real(kind=dp),intent(out)::pp,cc,dd
  real(kind=dp)::p1,p2,x,y,itheta
  integer::ierr,flag,ipp

  flag=0
  itheta=theta_in

  x=cos(itheta)
  y=sin(itheta)

  pp=0.0_dp
  CC=0.0_dp
  DD=0.0_dp
  if(abs(theta_in)<1.0E-10_dp) then
     select case(m)
     case(0)
        CC=0;
        DD=-n*(n+1)*0.5_dp*theta_in
        pp=1.0_dp
     case(1)
        CC=-n*(n+1)*0.5_dp
        DD=CC
     end select
     return
  elseif(abs(abs(theta_in)-pi)<1.0E-10_dp) then
     select case(m)
     case(0)
        CC=0;
        DD=((-1.0_dp)**n)*n*(n+1)*0.5_dp*(pi-theta_in)
        pp=1.0_dp
        if(mod(n+m,2)/=0) pp=-pp
     case(1)
        CC=((-1.0_dp)**n)*n*(n+1)*0.5_dp
        DD=-CC
     end select
     return
  end if
  if(itheta>(pi/2.0_dp)) then
     itheta=pi-itheta
     flag=1
  end if
   
  call dxlegf(real(n,dp),0,m,m,itheta,3,pp,ipp,ierr)
  if((mod(n+m,2)/=0).and.(flag==1)) then
     pp=-pp
  end if
  
  if(ipp/=0.0_dp) then
     print*,pp,ipp
      print*,"ipp is not zero ..."
     STOP
  else
     p1=pp
  end if

  CC=m*pp/y

  call dxlegf(real(n-1,dp),0,m,m,itheta,3,p2,ipp,ierr)
  if((mod(n+m,2)==0).and.(flag==1)) then
     p2=-p2
  end if
  if(ipp/=0) then
     print*,"ipp is not zero ...."
     STOP
  end if
     
  dd=(n*x*p1-(n+m)*p2)/y
  return
end subroutine qq
     
     
subroutine calc_field(exc_mode,RA,RB,RC,RD,RF,RG,num,freq,r_obs,theta_obs,phi_obs,Ef,Hf)  
  use mie_variables
  implicit none
  integer,intent(in)::exc_mode !0 for plane wave, dipole otherwise
  complex(kind=dp),dimension(nterms),intent(in)::RA,RB
  complex(kind=dp),dimension(nterms,n_layers),intent(inout)::RC,RD,RF,RG
  integer,intent(in)::num
  real(kind=dp),intent(in)::freq,r_obs,theta_obs,phi_obs
  complex(kind=dp),intent(out)::Ef(3),Hf(3)
  integer::n,m,i,j,count_phi,ierr,flag,nz,ierr2,ipp,index
  real(kind=dp)::pp,ltheta
  real(kind=dp)::zr,zi,D_fixed
  complex(kind=dp)::ri,k,factor,c1,der,der2,factor_p,der_p
  real(kind=dp)::CC,DD,theta,phi,y
  real(kind=dp),allocatable,dimension(:)::cyr,cyi,cyr2,cyi2,crwrk,ciwrk
  complex(kind=dp)::const,Po,Pe,Qo,Qe
  complex(kind=dp)::Mo(3),Me(3),No(3),Ne(3),Mo_p(3),Me_p(3),No_p(3),Ne_p(3)
  real(kind=dp)::const_den
  real(kind=dp)::eta_dual
  complex(kind=dp)::eta_rel_dual(0:n_l)

  if(Mdipole_exc) then
     eta_dual=1.0_dp/(120.0_dp*pi)
     eta_rel_dual(0:n_l)=sqrt(eps(0:n_l)/mu(0:n_l))
  else
     eta_dual=120.0_dp*pi
     eta_rel_dual(0:n_l)=sqrt(mu(0:n_l)/eps(0:n_l))
  end if

  c1=cmplx(0.0_dp,1.0_dp,dp)
  allocate(cyr(0:num),cyi(0:num),cyr2(0:num),cyi2(0:num),&
       crwrk(0:num),ciwrk(0:num))  

  Mo_p(1)=cmplx(0.0_dp,0.0_dp,dp)
  Me_p(1)=cmplx(0.0_dp,0.0_dp,dp)
  Mo(1)=cmplx(0.0_dp,0.0_dp,dp)
  Me(1)=cmplx(0.0_dp,0.0_dp,dp)

  index=n_l+1
  outer: do 
     if(r_obs>=a(index-1)) then
        exit outer
     end if
     
     index=index-1
     if(index==1) exit outer
  end do outer
 
  Ef(1:3)=cmplx(0.0_dp,0.0_dp,dp)
  Hf(1:3)=cmplx(0.0_dp,0.0_dp,dp)
 
  if(exc_mode/=0) then
     
     if(index==n_l+1) then
        ri=(ri_n(0))
        k=(2.0_dp*pi/lambda)*ri
        zr=real(k*r_obs);zi=-aimag(k*r_obs)
        
        call zbesh(zr,zi,0.5_dp,one,2,num+1,cyr,cyi,nz,ierr2)
        if (ierr2/=0) then
           print*,"error in the calculation of the hankel bessel with error &
                code:",ierr2
           stop
        end if
    
        factor=sqrt((pi*0.5_dp)/cmplx(zr,zi,dp))
        
        do n=num,1,-1
           
           der_p=(arg_p*cmplx(cyr_p(n-1),cyi_p(n-1),dp)-&
                n*cmplx(cyr_p(n),cyi_p(n),dp))
           D_fixed=0.25_dp*((2.0_dp*n)+1.0_dp)/(n*(n+1.0_dp))
           der=(zr*cmplx(cyr(n-1),cyi(n-1),dp)-&
                n*cmplx(cyr(n),cyi(n),dp))
           m=0
           const=D_fixed
           
           call qq(n,m,theta_p,pp,CC,DD)
           
           Mo_p(2)=cmplx(cyr_p(n),cyi_p(n),dp)*CC*cos(m*phi_p)
           Mo_p(3)=cmplx(0.0_dp,0.0_dp,dp)
           Me_p(2)=cmplx(0.0_dp,0.0_dp,dp)
           Me_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*cos(m*phi_p)
           
           No_p(1)=cmplx(0.0_dp,0.0_dp,dp)
           No_p(2)=cmplx(0.0_dp,0.0_dp,dp)
           No_p(3)=der_p*CC*cos(m*phi_p)/arg_p
           Ne_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*cos(m*phi_p)/arg_p
           Ne_p(2)=der_p*DD*cos(m*phi_p)/arg_p
           Ne_p(3)=cmplx(0.0_dp,0.0_dp,dp)
           
           Po=sum(Mo_p(1:3)*orient_dipole(1:3))
           Pe=sum(Me_p(1:3)*orient_dipole(1:3))
           
           Qo=sum(No_p(1:3)*orient_dipole(1:3))
           Qe=sum(Ne_p(1:3)*orient_dipole(1:3))
           
           call qq(n,m,theta_obs,pp,CC,DD)
           Mo(2)=cmplx(cyr(n),cyi(n),dp)*CC*cos(m*phi_obs)*factor
           Mo(3)=-cmplx(cyr(n),cyi(n),dp)*DD*sin(m*phi_obs)*factor
           Me(2)=-cmplx(cyr(n),cyi(n),dp)*CC*sin(m*phi_obs)*factor
           Me(3)=-cmplx(cyr(n),cyi(n),dp)*DD*cos(m*phi_obs)*factor
           
           No(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*sin(m*phi_obs)*factor/zr
           No(2)=der*DD*sin(m*phi_obs)*factor/zr
           No(3)=der*CC*cos(m*phi_obs)*factor/zr
           Ne(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*cos(m*phi_obs)*factor/zr
           Ne(2)=der*DD*cos(m*phi_obs)*factor/zr
           Ne(3)=-der*CC*sin(m*phi_obs)*factor/zr
           
           Ef(1:3)=Ef(1:3)+const*(conjg(RB(n))*(Po*Mo(1:3)+Pe*Me(1:3))+&
                conjg(RA(n))*(Qo*No(1:3)+Qe*Ne(1:3)))
           Hf(1:3)=Hf(1:3)+const*(conjg(RB(n))*(Po*No(1:3)+Pe*Ne(1:3))+&
                conjg(RA(n))*(Qo*Mo(1:3)+Qe*Me(1:3)))
           
           const_den=1.0_dp
           do m=1,n
              
              const_den=const_den*real((n-m+1),dp)*real((n+m),dp)
              
              const=2.0_dp*D_fixed/const_den
              call qq(n,m,theta_p,pp,CC,DD)
              Mo_p(2)=cmplx(cyr_p(n),cyi_p(n),dp)*CC*cos(m*phi_p)
              Mo_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*sin(m*phi_p)
              Me_p(2)=-cmplx(cyr_p(n),cyi_p(n),dp)*CC*sin(m*phi_p)
              Me_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*cos(m*phi_p)
              
              No_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*sin(m*phi_p)&
                   /arg_p
              No_p(2)=der_p*DD*sin(m*phi_p)/arg_p
              No_p(3)=der_p*CC*cos(m*phi_p)/arg_p
              Ne_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*cos(m*phi_p)&
                   /arg_p
              Ne_p(2)=der_p*DD*cos(m*phi_p)/arg_p
              Ne_p(3)=-der_p*CC*sin(m*phi_p)/arg_p
              
              
              Po=sum(Mo_p(1:3)*orient_dipole(1:3))
              Pe=sum(Me_p(1:3)*orient_dipole(1:3))
              
              Qo=sum(No_p(1:3)*orient_dipole(1:3))
              Qe=sum(Ne_p(1:3)*orient_dipole(1:3))
              
              call qq(n,m,theta_obs,pp,CC,DD)
              Mo(2)=cmplx(cyr(n),cyi(n),dp)*CC*cos(m*phi_obs)*factor
              Mo(3)=-cmplx(cyr(n),cyi(n),dp)*DD*sin(m*phi_obs)*factor
              Me(2)=-cmplx(cyr(n),cyi(n),dp)*CC*sin(m*phi_obs)*factor
              Me(3)=-cmplx(cyr(n),cyi(n),dp)*DD*cos(m*phi_obs)*factor
              
              No(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*sin(m*phi_obs)*factor&
                   /zr
              No(2)=der*DD*sin(m*phi_obs)*factor/zr
              No(3)=der*CC*cos(m*phi_obs)*factor/zr
              Ne(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*cos(m*phi_obs)*factor&
                   /zr
              Ne(2)=der*DD*cos(m*phi_obs)*factor/zr
              Ne(3)=-der*CC*sin(m*phi_obs)*factor/zr
              
              Ef(1:3)=Ef(1:3)+const*(conjg(RB(n))*(Po*Mo(1:3)+Pe*Me(1:3))+&
                   conjg(RA(n))*(Qo*No(1:3)+Qe*Ne(1:3)))
              Hf(1:3)=Hf(1:3)+const*(conjg(RB(n))*(Po*No(1:3)+Pe*Ne(1:3))+&
                   conjg(RA(n))*(Qo*Mo(1:3)+Qe*Me(1:3)))
              
           end do
        end do
        Hf(1:3)=Hf(1:3)*c1/(eta_rel_dual(0)*eta_dual)
        
     elseif(index==1) then
        
        ri=(ri_n(index))
        k=(2.0_dp*pi/lambda)*ri
        zr=real(k*r_obs);zi=-aimag(k*r_obs)
        call zbesj(zr,zi,0.5_dp,one,num+1,cyr,cyi,nz,ierr2)
        if (ierr2/=0) then
           print*,"error in the calculation of the hankel bessel with error &
                code:",ierr2
           stop
        end if
        
        factor=sqrt((pi*0.5_dp)/cmplx(zr,zi,dp))
        
        do n=num,1,-1
           der=(cmplx(zr,zi,dp)*cmplx(cyr(n-1),cyi(n-1),dp)-&
                n*cmplx(cyr(n),cyi(n),dp))
           
           der_p=(arg_p*cmplx(cyr_p(n-1),cyi_p(n-1),dp)-&
                n*cmplx(cyr_p(n),cyi_p(n),dp))
           D_fixed=0.25_dp*((2.0_dp*n)+1.0_dp)/(n*(n+1.0_dp))
           
           m=0
           const=D_fixed
           
           call qq(n,m,theta_p,pp,CC,DD)
           
           Mo_p(2)=cmplx(cyr_p(n),cyi_p(n),dp)*CC*cos(m*phi_p)
           Mo_p(3)=cmplx(0.0_dp,0.0_dp,dp)
           Me_p(2)=cmplx(0.0_dp,0.0_dp,dp)
           Me_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*cos(m*phi_p)
           
           No_p(1)=cmplx(0.0_dp,0.0_dp,dp)
           No_p(2)=cmplx(0.0_dp,0.0_dp,dp)
           No_p(3)=der_p*CC*cos(m*phi_p)/arg_p
           Ne_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*cos(m*phi_p)/arg_p
           Ne_p(2)=der_p*DD*cos(m*phi_p)/arg_p
           Ne_p(3)=cmplx(0.0_dp,0.0_dp,dp)
           
           Po=sum(Mo_p(1:3)*orient_dipole(1:3))
           Pe=sum(Me_p(1:3)*orient_dipole(1:3))
           
           Qo=sum(No_p(1:3)*orient_dipole(1:3))
           Qe=sum(Ne_p(1:3)*orient_dipole(1:3))
           
           call qq(n,m,theta_obs,pp,CC,DD)
           Mo(2)=cmplx(cyr(n),cyi(n),dp)*CC*cos(m*phi_obs)*factor
           Mo(3)=-cmplx(cyr(n),cyi(n),dp)*DD*sin(m*phi_obs)*factor
           Me(2)=-cmplx(cyr(n),cyi(n),dp)*CC*sin(m*phi_obs)*factor
           Me(3)=-cmplx(cyr(n),cyi(n),dp)*DD*cos(m*phi_obs)*factor
           
           No(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*sin(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           No(2)=der*DD*sin(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           No(3)=der*CC*cos(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           Ne(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*cos(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           Ne(2)=der*DD*cos(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           Ne(3)=-der*CC*sin(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           
           Ef(1:3)=Ef(1:3)+const*(-conjg(RC(n,1))*(Po*Mo(1:3)+Pe*Me(1:3))-&
                conjg(RD(n,1))*(Qo*No(1:3)+Qe*Ne(1:3)))
           Hf(1:3)=Hf(1:3)+const*(-conjg(RC(n,1))*(Po*No(1:3)+Pe*Ne(1:3))-&
                conjg(RD(n,1))*(Qo*Mo(1:3)+Qe*Me(1:3)))
           
           const_den=1.0_dp
           do m=1,n
              const_den=const_den*real((n-m+1),dp)*real((n+m),dp)
              const=2.0_dp*D_fixed/const_den
              call qq(n,m,theta_p,pp,CC,DD)
              Mo_p(2)=cmplx(cyr_p(n),cyi_p(n),dp)*CC*cos(m*phi_p)
              Mo_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*sin(m*phi_p)
              Me_p(2)=-cmplx(cyr_p(n),cyi_p(n),dp)*CC*sin(m*phi_p)
              Me_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*cos(m*phi_p)
              
              No_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*sin(m*phi_p)&
                   /arg_p
              No_p(2)=der_p*DD*sin(m*phi_p)/arg_p
              No_p(3)=der_p*CC*cos(m*phi_p)/arg_p
              Ne_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*cos(m*phi_p)&
                   /arg_p
              Ne_p(2)=der_p*DD*cos(m*phi_p)/arg_p
              Ne_p(3)=-der_p*CC*sin(m*phi_p)/arg_p
              
              Po=sum(Mo_p(1:3)*orient_dipole(1:3))
              Pe=sum(Me_p(1:3)*orient_dipole(1:3))
              Qo=sum(No_p(1:3)*orient_dipole(1:3))
              Qe=sum(Ne_p(1:3)*orient_dipole(1:3))
              
              call qq(n,m,theta_obs,pp,CC,DD)
              Mo(2)=cmplx(cyr(n),cyi(n),dp)*CC*cos(m*phi_obs)*factor
              Mo(3)=-cmplx(cyr(n),cyi(n),dp)*DD*sin(m*phi_obs)*factor
              Me(2)=-cmplx(cyr(n),cyi(n),dp)*CC*sin(m*phi_obs)*factor
              Me(3)=-cmplx(cyr(n),cyi(n),dp)*DD*cos(m*phi_obs)*factor
              
              No(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*sin(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              No(2)=der*DD*sin(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              No(3)=der*CC*cos(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              Ne(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*cos(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              Ne(2)=der*DD*cos(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              Ne(3)=-der*CC*sin(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              
              Ef(1:3)=Ef(1:3)+const*(-conjg(RC(n,1))*(Po*Mo(1:3)+Pe*Me(1:3))-&
                   conjg(RD(n,1))*(Qo*No(1:3)+Qe*Ne(1:3)))
              Hf(1:3)=Hf(1:3)+const*(-conjg(RC(n,1))*(Po*No(1:3)+Pe*Ne(1:3))-&
                   conjg(RD(n,1))*(Qo*Mo(1:3)+Qe*Me(1:3)))
           end do
           
        end do
        Hf(1:3)=Hf(1:3)*c1/(eta_rel_dual(index)*eta_dual)
     else
        ri=(ri_n(index))
        k=(2.0_dp*pi/lambda)*ri
        zr=real(k*r_obs);zi=-aimag(k*r_obs)
        call zbesj(zr,zi,0.5_dp,one,num+1,cyr,cyi,nz,ierr2)
        if (ierr2/=0) then
           print*,"error in the calculation of the hankel bessel with error &
             code:",ierr2
           stop
        end if
        call zbesy(zr,zi,0.5_dp,one,num+1,cyr2,cyi2,nz,&
             crwrk,ciwrk,ierr)
        factor=sqrt((pi*0.5_dp)/cmplx(zr,zi,dp))
        
        do n=num,1,-1
           der=(cmplx(zr,zi,dp)*cmplx(cyr(n-1),cyi(n-1),dp)-&
             n*cmplx(cyr(n),cyi(n),dp))
           der2=(cmplx(zr,zi,dp)*cmplx(cyr2(n-1),cyi2(n-1),dp)-&
                n*cmplx(cyr2(n),cyi2(n),dp))
           der_p=(arg_p*cmplx(cyr_p(n-1),cyi_p(n-1),dp)-&
                n*cmplx(cyr_p(n),cyi_p(n),dp))
           D_fixed=0.25_dp*((2.0_dp*n)+1.0_dp)/(n*(n+1.0_dp))
           
           m=0
           const=D_fixed
           
           call qq(n,m,theta_p,pp,CC,DD)
           
           Mo_p(2)=cmplx(cyr_p(n),cyi_p(n),dp)*CC*cos(m*phi_p)
           Mo_p(3)=cmplx(0.0_dp,0.0_dp,dp)
           Me_p(2)=cmplx(0.0_dp,0.0_dp,dp)
           Me_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*cos(m*phi_p)
           
           No_p(1)=cmplx(0.0_dp,0.0_dp,dp)
           No_p(2)=cmplx(0.0_dp,0.0_dp,dp)
           No_p(3)=der_p*CC*cos(m*phi_p)/arg_p
           Ne_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*cos(m*phi_p)/arg_p
           Ne_p(2)=der_p*DD*cos(m*phi_p)/arg_p
           Ne_p(3)=cmplx(0.0_dp,0.0_dp,dp)
           
           Po=sum(Mo_p(1:3)*orient_dipole(1:3))
           Pe=sum(Me_p(1:3)*orient_dipole(1:3))
           
           Qo=sum(No_p(1:3)*orient_dipole(1:3))
           Qe=sum(Ne_p(1:3)*orient_dipole(1:3))
           
           call qq(n,m,theta_obs,pp,CC,DD)
           Mo(2)=cmplx(cyr(n),cyi(n),dp)*CC*cos(m*phi_obs)*factor
           Mo(3)=-cmplx(cyr(n),cyi(n),dp)*DD*sin(m*phi_obs)*factor
           Me(2)=-cmplx(cyr(n),cyi(n),dp)*CC*sin(m*phi_obs)*factor
           Me(3)=-cmplx(cyr(n),cyi(n),dp)*DD*cos(m*phi_obs)*factor
           
           No(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*sin(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           No(2)=der*DD*sin(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           No(3)=der*CC*cos(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           Ne(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*cos(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           Ne(2)=der*DD*cos(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           Ne(3)=-der*CC*sin(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           
           Ef(1:3)=Ef(1:3)+const*(-conjg(RC(n,index))*(Po*Mo(1:3)+Pe*Me(1:3))-&
                conjg(RD(n,index))*(Qo*No(1:3)+Qe*Ne(1:3)))
           Hf(1:3)=Hf(1:3)+const*(-conjg(RC(n,index))*(Po*No(1:3)+Pe*Ne(1:3))-&
                conjg(RD(n,index))*(Qo*Mo(1:3)+Qe*Me(1:3)))
           
           Mo(2)=cmplx(cyr2(n),cyi2(n),dp)*CC*cos(m*phi_obs)*factor
           Mo(3)=-cmplx(cyr2(n),cyi2(n),dp)*DD*sin(m*phi_obs)*factor
           Me(2)=-cmplx(cyr2(n),cyi2(n),dp)*CC*sin(m*phi_obs)*factor
           Me(3)=-cmplx(cyr2(n),cyi2(n),dp)*DD*cos(m*phi_obs)*factor
           No(1)=n*(n+1.0_dp)*cmplx(cyr2(n),cyi2(n),dp)*pp*sin(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           No(2)=der2*DD*sin(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           No(3)=der2*CC*cos(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           Ne(1)=n*(n+1.0_dp)*cmplx(cyr2(n),cyi2(n),dp)*pp*cos(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           Ne(2)=der2*DD*cos(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           Ne(3)=-der2*CC*sin(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
           
           Ef(1:3)=Ef(1:3)+const*(-conjg(RF(n,index))*(Po*Mo(1:3)+Pe*Me(1:3))-&
                conjg(RG(n,index))*(Qo*No(1:3)+Qe*Ne(1:3)))
           Hf(1:3)=Hf(1:3)+const*(-conjg(RF(n,index))*(Po*No(1:3)+Pe*Ne(1:3))-&
                conjg(RG(n,index))*(Qo*Mo(1:3)+Qe*Me(1:3)))
           
           const_den=1.0_dp
           do m=1,n
              const_den=const_den*real((n-m+1),dp)*real((n+m),dp)
              const=2.0_dp*D_fixed/const_den
              call qq(n,m,theta_p,pp,CC,DD)
              Mo_p(2)=cmplx(cyr_p(n),cyi_p(n),dp)*CC*cos(m*phi_p)
              Mo_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*sin(m*phi_p)
              Me_p(2)=-cmplx(cyr_p(n),cyi_p(n),dp)*CC*sin(m*phi_p)
              Me_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*cos(m*phi_p)
              
              No_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*sin(m*phi_p)&
                /arg_p
              No_p(2)=der_p*DD*sin(m*phi_p)/arg_p
              No_p(3)=der_p*CC*cos(m*phi_p)/arg_p
              Ne_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*cos(m*phi_p)&
                   /arg_p
              Ne_p(2)=der_p*DD*cos(m*phi_p)/arg_p
              Ne_p(3)=-der_p*CC*sin(m*phi_p)/arg_p

              Po=sum(Mo_p(1:3)*orient_dipole(1:3))
              Pe=sum(Me_p(1:3)*orient_dipole(1:3))
              Qo=sum(No_p(1:3)*orient_dipole(1:3))
              Qe=sum(Ne_p(1:3)*orient_dipole(1:3))
              
              call qq(n,m,theta_obs,pp,CC,DD)
              Mo(2)=cmplx(cyr(n),cyi(n),dp)*CC*cos(m*phi_obs)*factor
              Mo(3)=-cmplx(cyr(n),cyi(n),dp)*DD*sin(m*phi_obs)*factor
              Me(2)=-cmplx(cyr(n),cyi(n),dp)*CC*sin(m*phi_obs)*factor
              Me(3)=-cmplx(cyr(n),cyi(n),dp)*DD*cos(m*phi_obs)*factor
              
              No(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*sin(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              No(2)=der*DD*sin(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
              No(3)=der*CC*cos(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              Ne(1)=n*(n+1.0_dp)*cmplx(cyr(n),cyi(n),dp)*pp*cos(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              Ne(2)=der*DD*cos(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              Ne(3)=-der*CC*sin(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
           
              Ef(1:3)=Ef(1:3)+const*(-conjg(RC(n,index))*(Po*Mo(1:3)+Pe*Me(1:3))-&
                   conjg(RD(n,index))*(Qo*No(1:3)+Qe*Ne(1:3)))
              Hf(1:3)=Hf(1:3)+const*(-conjg(RC(n,index))*(Po*No(1:3)+Pe*Ne(1:3))-&
                   conjg(RD(n,index))*(Qo*Mo(1:3)+Qe*Me(1:3)))
              
              Mo(2)=cmplx(cyr2(n),cyi2(n),dp)*CC*cos(m*phi_obs)*factor
              Mo(3)=-cmplx(cyr2(n),cyi2(n),dp)*DD*sin(m*phi_obs)*factor
              Me(2)=-cmplx(cyr2(n),cyi2(n),dp)*CC*sin(m*phi_obs)*factor
              Me(3)=-cmplx(cyr2(n),cyi2(n),dp)*DD*cos(m*phi_obs)*factor
              
              No(1)=n*(n+1.0_dp)*cmplx(cyr2(n),cyi2(n),dp)*pp*sin(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              No(2)=der2*DD*sin(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              No(3)=der2*CC*cos(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
              Ne(1)=n*(n+1.0_dp)*cmplx(cyr2(n),cyi2(n),dp)*pp*cos(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              Ne(2)=der2*DD*cos(m*phi_obs)/&
                   cmplx(zr,zi,dp)*factor
              Ne(3)=-der2*CC*sin(m*phi_obs)/&
                cmplx(zr,zi,dp)*factor
              
              Ef(1:3)=Ef(1:3)+const*(-conjg(RF(n,index))*(Po*Mo(1:3)+Pe*Me(1:3))-&
                   conjg(RG(n,index))*(Qo*No(1:3)+Qe*Ne(1:3)))
              Hf(1:3)=Hf(1:3)+const*(-conjg(RF(n,index))*(Po*No(1:3)+Pe*Ne(1:3))-&
                   conjg(RG(n,index))*(Qo*Mo(1:3)+Qe*Me(1:3)))
           end do
           
        end do
        Hf(1:3)=Hf(1:3)*c1/(eta_rel_dual(index)*eta_dual)
     end if
  else
     
     if(index==n_l+1) then
        ri=ri_n(0)
        k=(2.0_dp*pi/lambda)*ri
        zr=real(k*r_obs);zi=aimag(k*r_obs)
        
        call zbesh(zr,zi,0.5_dp,one,2,num+1,cyr,cyi,nz,ierr2)
        if (ierr2/=0) then
           print*,"error in the calculation of the hankel bessel with error &
                code:",ierr2
           stop
        end if
        
        factor=sqrt((pi*0.5_dp)/cmplx(zr,zi,dp))
        
        cyr=cyr*factor;
        cyi=cyi*factor;
        
        do n=num,1,-1
           call qq(n,1,theta_obs,pp,CC,DD)
           const=((2.0_dp*n)+1.0_dp)/(n*(n+1.0_dp))*((-c1)**n)
           
           Ef(1)= Ef(1)-c1*(conjg(RA(n)))*n*(n+1.0_dp)*&
                (cmplx(cyr(n),cyi(n),dp)/zr)*&
                CC*cos(phi_obs)*const*sin(theta_obs)
           der=(zr*cmplx(cyr(n-1),cyi(n-1),dp)-&
                n*cmplx(cyr(n),cyi(n),dp))
           Ef(2)=Ef(2)-c1*conjg(RA(n))*(der/zr)*&
                DD*cos(phi_obs)*const
           Ef(2)=(Ef(2)-(conjg(RB(n)))*cmplx(cyr(n),cyi(n),dp)*&
                CC*cos(phi_obs)*const)        
           Ef(3)=Ef(3)+c1*conjg(RA(n))*(der/zr)*&
                CC*sin(phi_obs)*const
           Ef(3)=(Ef(3)+conjg(RB(n))*cmplx(cyr(n),cyi(n),dp)*&
                DD*sin(phi_obs)*const)
           
           Hf(1)=Hf(1)-(conjg(RB(n)))*n*(n+1.0_dp)*&
                (cmplx(cyr(n),cyi(n),dp)/zr)*&
                CC*sin(phi_obs)*const*sin(theta_obs)
           Hf(2)=Hf(2)+c1*conjg(RA(n))*cmplx(cyr(n),cyi(n),dp)*&
                CC*sin(phi_obs)*const
           Hf(2)=(Hf(2)-(conjg(RB(n)))*(der/zr)*&
                DD*sin(phi_obs)*const)
           Hf(3)=Hf(3)+c1*conjg(RA(n))*cmplx(cyr(n),cyi(n),dp)*&
                DD*cos(phi_obs)*const
           Hf(3)=(Hf(3)-conjg(RB(n))*(der/zr)*&
                CC*cos(phi_obs)*const)
        end do
        Hf(1:3)=Hf(1:3)*c1/(eta_rel_dual(0)*eta_dual)
     
     elseif(index==1) then
        
        ri=(ri_n(index))
        k=(2.0_dp*pi/lambda)*ri
        zr=real(k*r_obs);zi=aimag(k*r_obs)
        call zbesj(zr,zi,0.5_dp,one,num+1,cyr,cyi,nz,ierr2)
        if (ierr2/=0) then
           print*,"error in the calculation of the hankel bessel with error &
                code:",ierr2
           stop
        end if
        
        factor=sqrt((pi*0.5_dp)/cmplx(zr,zi,dp))
        
        do n=num,1,-1
           call qq(n,1,theta_obs,pp,CC,DD)
           const=((-c1)**n)*((2.0_dp*n)+1.0_dp)/(n*(n+1.0_dp))
           Ef(1)=Ef(1)+c1*const*conjg(RD(n,1)*n*(n+1.0_dp)*&
                (cmplx(cyr(n),cyi(n),dp)/cmplx(zr,zi,dp))*&
                CC*cos(phi_obs)*sin(theta_obs)*factor)
           
           der=(cmplx(zr,zi,dp)*cmplx(cyr(n-1),cyi(n-1),dp)-&
                n*cmplx(cyr(n),cyi(n),dp))*factor
           Ef(2)=Ef(2)+c1*const*conjg(RD(n,1)*(der/cmplx(zr,zi,dp))*&
                DD*cos(phi_obs))
           Ef(2)=Ef(2)+const*conjg(RC(n,1)*cmplx(cyr(n),cyi(n),dp)*&
                CC*cos(phi_obs)*factor)
           Ef(3)=Ef(3)-c1*const*conjg(RD(n,1)*(der/cmplx(zr,zi,dp))*&
                CC*sin(phi_obs))
           Ef(3)=Ef(3)-const*conjg(RC(n,1)*cmplx(cyr(n),cyi(n),dp)*&
                DD*sin(phi_obs)*factor)

           Hf(1)=Hf(1)-const*conjg(RC(n,1)*n*(n+1.0_dp)*&
                (cmplx(cyr(n),cyi(n),dp)/cmplx(zr,zi,dp))*&
                CC*sin(phi_obs)*sin(theta_obs)*factor)
           Hf(2)=Hf(2)-const*conjg(RC(n,1)*(der/cmplx(zr,zi,dp))*&
                DD*sin(phi_obs))
           Hf(2)=Hf(2)-c1*const*conjg(RD(n,1)*cmplx(cyr(n),cyi(n),dp)*&
                CC*sin(phi_obs)*factor)
           Hf(3)=Hf(3)-const*conjg(RC(n,1)*(der/cmplx(zr,zi,dp))*&
                CC*cos(phi_obs))
           Hf(3)=Hf(3)-c1*const*conjg(RD(n,1)*cmplx(cyr(n),cyi(n),dp)*&
                DD*cos(phi_obs)*factor)
        end do
        Hf(1:3)=Hf(1:3)*c1/(eta_rel_dual(index)*eta_dual)
      
     else
        ri=(ri_n(index))
        k=(2.0_dp*pi/lambda)*ri
        zr=real(k*r_obs);zi=aimag(k*r_obs)
        call zbesj(zr,zi,0.5_dp,one,num+1,cyr,cyi,nz,ierr2)
        if (ierr2/=0) then
           print*,"error in the calculation of the hankel bessel with error &
                code:",ierr2
           stop
        end if
        call zbesy(zr,zi,0.5_dp,one,num+1,cyr2,cyi2,nz,&
             crwrk,ciwrk,ierr)
        factor=sqrt((pi*0.5_dp)/cmplx(zr,zi,dp))
        
        do n=num,1,-1
           call qq(n,1,theta_obs,pp,CC,DD)
        
           const=((-c1)**n)*((2.0_dp*n)+1.0_dp)/(n*(n+1.0_dp))
           Ef(1)=Ef(1)+c1*const*conjg(RD(n,index)*n*(n+1.0_dp)*&
                (cmplx(cyr(n),cyi(n),dp)/cmplx(zr,zi,dp))*&
                CC*cos(phi_obs)*sin(theta_obs)*factor)
           Ef(1)=Ef(1)+c1*const*conjg(RG(n,index)*n*(n+1.0_dp)*&
                (cmplx(cyr2(n),cyi2(n),dp)/cmplx(zr,zi,dp))*&
                CC*cos(phi_obs)*sin(theta_obs)*factor)
           
           
           der=(cmplx(zr,zi,dp)*cmplx(cyr(n-1),cyi(n-1),dp)-&
                n*cmplx(cyr(n),cyi(n),dp))*factor
           Ef(2)=Ef(2)+c1*const*conjg(RD(n,index)*(der/cmplx(zr,zi,dp))*&
                DD*cos(phi_obs))
           Ef(2)=Ef(2)+const*conjg(RC(n,index)*cmplx(cyr(n),cyi(n),dp)*&
                CC*cos(phi_obs)*factor)
           
           Ef(3)=Ef(3)-c1*const*conjg(RD(n,index)*(der/cmplx(zr,zi,dp))*&
                CC*sin(phi_obs))
           Ef(3)=Ef(3)-const*conjg(RC(n,index)*cmplx(cyr(n),cyi(n),dp)*&
                DD*sin(phi_obs)*factor)

           Hf(1)=Hf(1)-const*conjg(RC(n,1)*n*(n+1.0_dp)*&
                (cmplx(cyr(n),cyi(n),dp)/cmplx(zr,zi,dp))*&
                CC*sin(phi_obs)*sin(theta_obs)*factor)

           Hf(1)=Hf(1)-const*conjg(RF(n,1)*n*(n+1.0_dp)*&
                (cmplx(cyr2(n),cyi2(n),dp)/cmplx(zr,zi,dp))*&
                CC*sin(phi_obs)*sin(theta_obs)*factor)

           Hf(2)=Hf(2)-const*conjg(RC(n,1)*(der/cmplx(zr,zi,dp))*&
                DD*sin(phi_obs))
           Hf(2)=Hf(2)-c1*const*conjg(RD(n,1)*cmplx(cyr(n),cyi(n),dp)*&
                CC*sin(phi_obs)*factor)
           Hf(3)=Hf(3)-const*conjg(RC(n,1)*(der/cmplx(zr,zi,dp))*&
                CC*cos(phi_obs))
           Hf(3)=Hf(3)-c1*const*conjg(RD(n,1)*cmplx(cyr(n),cyi(n),dp)*&
                DD*cos(phi_obs)*factor)

           
           der=(cmplx(zr,zi,dp)*cmplx(cyr2(n-1),cyi2(n-1),dp)-&
                n*cmplx(cyr2(n),cyi2(n),dp))*factor
           Ef(2)=Ef(2)+c1*const*conjg(RG(n,index)*(der/cmplx(zr,zi,dp))*&
                DD*cos(phi_obs))
           Ef(2)=Ef(2)+const*conjg(RF(n,index)*cmplx(cyr2(n),cyi2(n),dp)*&
                CC*cos(phi_obs)*factor)
           
           Ef(3)=Ef(3)-c1*const*conjg(RG(n,index)*(der/cmplx(zr,zi,dp))*&
                CC*sin(phi_obs))
           Ef(3)=Ef(3)-const*conjg(RF(n,index)*cmplx(cyr2(n),cyi2(n),dp)*&
                DD*sin(phi_obs)*factor)

           Hf(2)=Hf(2)-const*conjg(RF(n,index)*(der/cmplx(zr,zi,dp))*&
                DD*sin(phi_obs))
           Hf(2)=Hf(2)-c1*const*conjg(RG(n,index)*cmplx(cyr2(n),cyi2(n),dp)*&
                CC*sin(phi_obs)*factor)
           
           Hf(3)=Hf(3)-const*conjg(RF(n,index)*(der/cmplx(zr,zi,dp))*&
                CC*cos(phi_obs))
           Hf(3)=Hf(3)-c1*const*conjg(RG(n,index)*cmplx(cyr2(n),cyi2(n),dp)*&
                DD*cos(phi_obs)*factor)

        end do
        Hf(1:3)=Hf(1:3)*c1/(eta_rel_dual(index)*eta_dual)
     end if
     Ef=-Ef
     Hf=-Hf
  end if
  deallocate(cyr,cyi,cyr2,cyi2,crwrk,ciwrk)
  return
 end subroutine calc_field



subroutine calc_rcs(RA,RB,num,freq)  
  use mie_variables
  implicit none
      
  complex(kind=dp),dimension(nterms),intent(in)::RA,RB
  integer,intent(in)::num
  real(kind=dp),intent(in)::freq
  integer::n,i,j,count_phi,ierr,flag,nz,ierr2,ipp,index,&
       count_the
  real(kind=dp)::pp,dtheta,dphi,theta_obs,phi_obs
  real(kind=dp)::const,CC,DD,sp(3),cart(3),new_cart(3),new_sp(3)
  complex(kind=dp)::s1,s2
  real(kind=dp),allocatable,dimension(:)::theta,phi
  real(kind=dp)::Nrcsv,Nrcsh,ueincd(3),ued(3),upropd(3),&
       uthe(3),uthe_to_cart(3),uphi(3),uphi_to_cart(3),&
       Ecart_ref(3),Ecart(3),Esp(3)
  
  print*,"Opening file frcs_mie.out"
  open(unit=07,file='frcs_mie.out',status='unknown')
43 FORMAT(F18.8,4F23.15)
45 FORMAT(I4)
46 FORMAT(2F10.6)  
 
  allocate(theta(1:ntheta),phi(1:nphi))
        
  s1=cmplx(0.0_dp,0.0_dp,dp)
  s2=cmplx(0.0_dp,0.0_dp,dp)
  if(ntheta/=1) then
     dtheta=(thetaend-thetastart)/(ntheta-1.0_dp);
  else
     dtheta=0
  end if

  if (nphi/=1) then
     dphi=(phimax-phimin)/(nphi-1.0_dp);
  else
     dphi=0;
  end if
 
  do i=1,nphi
     phi(i)=((phimin+(i-1)*dphi)*pi)/180.0_dp
  end do
  
  do i=1,ntheta
     theta(i)=((thetastart+(i-1)*dtheta)*pi)/180.0_dp
  end do
 
  upropd(1)=sin(the_angle)*cos(phi_angle)                                  
  upropd(2)=sin(the_angle)*sin(phi_angle)                                  
  upropd(3)=cos(the_angle)                                 
  
  ued(1)=cos(the_angle)*cos(phi_angle)*ethetad-sin(phi_angle)*ephid         
  ued(2)=cos(the_angle)*sin(phi_angle)*ethetad+cos(phi_angle)*ephid         
  ued(3)=-sin(the_angle)*ethetad

  ueincd(1:3)=ued(1:3)/sqrt(dot_product(ued,ued))

  print*,'The plane wave is travelling  to',upropd(1),upropd(2),upropd(3)
  print*,'Incident Efield components:',ueincd(1),ueincd(2),ueincd(3)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do count_phi=1,nphi
     do count_the=1,ntheta
        sp(1)=1.0_dp
        sp(2)=theta(count_the)
        sp(3)=phi(count_phi)

        theta_obs=sp(2)
        phi_obs=sp(3)

        uthe(1)=0.0_dp
        uthe(2)=1.0_dp
        uthe(3)=0.0_dp
        
        call sph2cart_real(theta_obs,phi_obs,uthe,uthe_to_cart)

        uphi(1)=0.0_dp
        uphi(2)=0.0_dp
        uphi(3)=1.0_dp
        
        call sph2cart_real(theta_obs,phi_obs,uphi,uphi_to_cart)

        call sph2cart_coord(sp,cart)
!!$!!!!!!!!!!!!!!
!For the rotation of axes and get results for the reference case 
!
        new_cart(1)=-ethetad*(cart(1)*cos(the_angle)*cos(phi_angle)+&
             cart(2)*cos(the_angle)*sin(phi_angle)-cart(3)*&
             sin(the_angle))
       
        new_cart(1)=new_cart(1)+ephid*(cart(1)*sin(phi_angle)-&
             cart(2)*cos(phi_angle))
       
        new_cart(2)=ethetad*(-cart(1)*sin(phi_angle)+&
             cart(2)*cos(phi_angle))
        new_cart(2)= new_cart(2)-ephid*&
             (cart(1)*cos(the_angle)*cos(phi_angle)+&
             cart(2)*cos(the_angle)*sin(phi_angle)-cart(3)*sin(the_angle))
        
        new_cart(3)=-cart(1)*sin(the_angle)*cos(phi_angle)-&
             cart(2)*sin(the_angle)*sin(phi_angle)-cart(3)&
             *cos(the_angle)
!!$!!!!!!!!!!!!!
        call cart2sph(new_cart,new_sp)
        
        theta_obs=new_sp(2)
        phi_obs=new_sp(3)
        s1=cmplx(0.0_dp,0.0_dp,dp)
        s2=cmplx(0.0_dp,0.0_dp,dp)
        do n=num,1,-1
           const=((2.0_dp*n)+1.0_dp)/(n*(n+1.0_dp))
           call qq(n,1,pi-theta_obs,pp,CC,DD)
           s1=s1+const*(RA(n)*DD+RB(n)*CC)
           s2=s2+const*(RA(n)*CC+RB(n)*DD)
        end do
        
        Esp(1)=0.0_dp
        Esp(2)= lambda*abs(s1)*cos(phi_obs) 
        Esp(3)= lambda*abs(s2)*sin(phi_obs)
        call sph2cart_real(theta_obs,phi_obs,Esp,Ecart)


        Ecart_ref(1)=-ethetad*(Ecart(1)*cos(the_angle)*cos(phi_angle)+&
             Ecart(2)*sin(phi_angle))-Ecart(3)*sin(the_angle)*cos(phi_angle)
        Ecart_ref(2)=ethetad*(-Ecart(1)*cos(the_angle)*sin(phi_angle)+&
             Ecart(2)*cos(phi_angle))-Ecart(3)*sin(the_angle)*sin(phi_angle)
        Ecart_ref(3)=ethetad*Ecart(1)*sin(the_angle)-&
             Ecart(3)*cos(the_angle)
        
        
        Ecart_ref(1)=Ecart_ref(1)+ephid*(Ecart(1)*sin(phi_angle)-&
             Ecart(2)*cos(the_angle)*sin(phi_angle))
        Ecart_ref(2)= Ecart_ref(2)-ephid*(Ecart(1)*cos(phi_angle)+&
             Ecart(2)*cos(the_angle)*sin(phi_angle))
        Ecart_ref(3)= Ecart_ref(3)+ephid*Ecart(2)*sin(the_angle)
        
        Nrcsv=dot_product(Ecart_ref,uthe_to_cart)
        Nrcsh=dot_product(Ecart_ref,uphi_to_cart)

        Nrcsv= Nrcsv**2/pi
        Nrcsh= Nrcsh**2/pi
        
        write(7,43) freq,(theta(count_the))*180.0_dp/pi,&
             phi(count_phi)*180.0_dp/pi,Nrcsv,Nrcsh
     end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$  do n=num,1,-1
!!$   
!!$     const=((2.0_dp*n)+1.0_dp)/(n*(n+1.0_dp))
!!$     do i=1,ntheta
!!$        call qq(n,1,theta(i),pp,CC,DD)
!!$        s1(i)=s1(i)+const*(RA(n)*DD+RB(n)*CC)
!!$        s2(i)=s2(i)+const*(RA(n)*CC+RB(n)*DD)
!!$     end do
!!$     
!!$  end do
!!$  
!!$  do count_phi=1,nphi
!!$     do i=1,ntheta
!!$        Nrcsv(i)= lambda**2/pi * ((abs(s1(i))*cos(phi(count_phi)))**2) 
!!$        Nrcsh(i)= lambda**2/pi * ((abs(s2(i))*sin(phi(count_phi)))**2)
!!$     end do
!!$     
!!$     do j=ntheta,1,-1
!!$        write(7,43) freq,(pi-theta(j))*180.0_dp/pi,&
!!$             phi(count_phi)*180.0_dp/pi,Nrcsv(j),Nrcsh(j)
!!$     end do
!!$  end do
  deallocate(theta,phi)
  print*,"Closing Files ... frcs_mie.out"
  close(07,status='keep')
  return
  
end subroutine calc_rcs


subroutine Mie_power(l,num,RA,RB,RC,RD,RF,RG,power)
  use mie_variables
  implicit none
  integer,intent(in)::l,num
  complex(kind=dp),dimension(nterms),intent(in)::RA, RB
  complex(kind=dp),dimension(nterms,n_layers),intent(in)::RC,RD,RF,RG
  real(kind=dp),intent(out)::power
  real(kind=dp),allocatable,dimension(:)::cyr,cyi,cyr2,cyi2,cyr2wrk,&
       cyi2wrk,bz2,az2,bz2y,az2y,&
       mn1,nn1,cn2,dn2,mn2,nn2,fn2,gn2,qp,wght
  complex(kind=dp),allocatable,dimension(:)::bz,b1z,az,bzy,b1zy,azy,cnfnc,&
       dngnc,mn12c,nn12c
  integer,allocatable::n1(:),n2(:)
  real(kind=dp)::dr,r_obs,en,z2
  integer::nqp,j,nz,ierr,ierr2,i
  complex(kind=dp)::z,c1,sz,ri
  integer::n,m,count_phi,flag,ipp,index
  real(kind=dp)::const_den
  complex(kind=dp)::Po,Pe,Qo,Qe
  complex(kind=dp)::Mo(3),Me(3),No(3),Ne(3),Mo_p(3),Me_p(3),No_p(3),Ne_p(3) 
  real(kind=dp)::CC,DD,theta,phi,y
  real(kind=dp),allocatable,dimension(:)::crwrk,ciwrk
  real(kind=dp)::pp,const
  real(kind=dp)::zr,zi,D_fixed
  real(kind=dp)::dx,xj,x
  complex(kind=dp)::k,factor,der,der2,factor_p,der_p,const_dipole
 
 

  nqp=32
  allocate(qp(1:nqp),wght(1:nqp))
  call oneD_quadrature(nqp,qp,wght)
  c1=cmplx(0.0_dp,1.0_dp)

  if(plane_wave_exc) then
     ri=(ri_n(l))
     x=(2.0_dp*pi/lambda)*a(n_l)
     allocate(cyr(1:num),cyi(1:num),cyr2(1:num),cyi2(1:num),&
          cyr2wrk(1:num),cyi2wrk(1:num),bz(1:num),bz2(1:num),&
          az2(1:num),b1z(1:num),az(1:num),bzy(1:num),&
          bz2y(1:num),az2y(1:num),b1zy(1:num),azy(1:num),n1(1:num),&
          n2(1:num),&
          mn1(1:num),nn1(1:num),cn2(1:num),dn2(1:num),&
          mn2(1:num),nn2(1:num),fn2(1:num),gn2(1:num),&
          mn12c(1:num),nn12c(1:num),cnfnc(1:num),dngnc(1:num))
     if(l==1) then!!!dx denotes x_upper-x_lower
        dx=(x/a(n_l))*(a(l))!/real(nqp-1,dp)
     else
        dx=(x/a(n_l))*(a(l)-a(l-1))!/real(nqp-1,dp)
     end if
     n1=(/(i,i=1,num)/)
     n1=n1*(/(i+1,i=1,num)/)
     n2=2*(/(2*i+1,i=1,num)/)
     
     cn2(1:num)=(abs(RC(1:num,l)))**2
     dn2(1:num)=(abs(RD(1:num,l)))**2
     fn2(1:num)=(abs(RF(1:num,l)))**2
     gn2(1:num)=(abs(RG(1:num,l)))**2
     cnfnc(1:num)=RC(1:num,l)*conjg(RF(1:num,l))
     dngnc(1:num)=RD(1:num,l)*conjg(RG(1:num,l))
     
     power=0.0_dp
     
     do j=1,nqp
        if(l==1) then
           xj=dx*qp(j)!real(j-1,dp)
        else
           xj=(x/a(n_l))*a(l-1)+dx*qp(j)!*real(j-1,dp)
        end if
        if(xj/=0.0_dp) then
           z=ri*xj
           sz=sqrt(0.5_dp*pi/z)
           call zbesj(real(z),aimag(z),1.5_dp,1,num,cyr,cyi,nz,ierr)
           bz(1:num)=(cmplx(cyr(1:num),cyi(1:num),dp)*(sz))
           
           bz2(1:num)=(abs(bz(1:num)))**2
           b1z(1)=sin(z)/z
           b1z(2:num)=bz(1:num-1)
           
           az(1:num)=b1z(1:num)-(/(i,i=1,num)/)*bz(1:num)/z
           az2=(abs(az(1:num)))**2
           z2=(abs(z))**2
           mn1=bz2*n2
           nn1=n2*(az2+bz2*n1/z2)
           call zbesy(real(z),aimag(z),1.5_dp,1,num,cyr2,cyi2,nz,&
                cyr2wrk,cyi2wrk,ierr)
           bzy(1:num)=(cmplx(cyr2(1:num),cyi2(1:num),dp)*(sz))
           bz2y(1:num)=(abs(bzy(1:num)))**2
           b1zy(1)=-cos(z)/z
           b1zy(2:num)=bzy(1:num-1)
           azy(1:num)=b1zy(1:num)-(/(i,i=1,num)/)*bzy(1:num)/z
           az2y=(abs(azy(1:num)))**2
           
           mn2=bz2y*n2
           nn2=n2*(az2y+bz2y*n1/z2)
           mn12c=bz*conjg(bzy)
           mn12c=mn12c*n2
           
           nn12c=n2*(az*conjg(azy)+bz*conjg(bzy)*n1/z2)
           en=sum(cn2(1:num)*mn1(1:num))+sum(dn2(1:num)*nn1(1:num))+&
                sum(fn2(1:num)*mn2(1:num))+sum(gn2(1:num)*nn2(1:num))+2.0_dp*&
                real(sum((cnfnc(1:num))*(mn12c(1:num)))+&
                sum((dngnc(1:num))*(nn12c(1:num))))
           power=power+wght(j)*en*xj**2
        end if
     end do
     power=power*dx*aimag(ri*ri)
     write(12,*) "Numerically calculated Absorbed power in layer",l,&
          "=",power/240.0_dp/((2.0_dp*pi/lambda)**2)
    
     deallocate(cyr,cyi,cyr2,cyi2,cyr2wrk,cyi2wrk,bz,&
          bz2,az2,b1z,az,n1,n2,mn1,nn1,cn2,dn2,mn2,nn2,fn2,gn2,&
          mn12c,nn12c,dngnc,&
          cnfnc,wght,qp)

  elseif(Jdipole_exc) then
     allocate(cyr(0:num+1),cyi(0:num+1),cyr2(0:num+1),cyi2(0:num+1),&
          crwrk(0:num+1),ciwrk(0:num+1))  
     if(r_p>a(n_l)) then
        const_dipole=4.0_dp*(2.0_dp*pi/lambda)**2*&
             vlite*1.0E-7_dp*exc_dipole
        Mo_p(1)=cmplx(0.0_dp,0.0_dp,dp)
        Me_p(1)=cmplx(0.0_dp,0.0_dp,dp)
        Mo(1)=cmplx(0.0_dp,0.0_dp,dp)
        Me(1)=cmplx(0.0_dp,0.0_dp,dp)
        
        if(l==1) then!!!dr denotes x_upper-x_lower
           dr=a(l)
        else
           dr=a(l)-a(l-1)
        end if
        
        ri=(ri_n(l))
        k=(2.0_dp*pi/lambda)*ri
        power=0.0_dp
        
        if(l==1) then
           do j=1,nqp
              r_obs=dr*qp(j)
              zr=real(k*r_obs);zi=-aimag(k*r_obs)
              call zbesj(zr,zi,0.5_dp,one,num+2,cyr,cyi,nz,ierr2)
              if (ierr2/=0) then
                 print*,"error in the calculation of the hankel bessel with error &
                      code:",ierr2
                 stop
              end if
                          
              factor=sqrt((pi*0.5_dp)/cmplx(zr,zi,dp))
            
              do n=num,1,-1
                 der=(cmplx(zr,zi,dp)*cmplx(cyr(n-1),cyi(n-1),dp)-&
                      n*cmplx(cyr(n),cyi(n),dp))
                 
                 der_p=(arg_p*cmplx(cyr_p(n-1),cyi_p(n-1),dp)-&
                      n*cmplx(cyr_p(n),cyi_p(n),dp))
                 D_fixed=0.25_dp*((2.0_dp*n)+1.0_dp)/(n*(n+1.0_dp))
                 
                 m=0
                 const=D_fixed
                 
                 call qq(n,m,theta_p,pp,CC,DD)
                 
                 Mo_p(2)=cmplx(cyr_p(n),cyi_p(n),dp)*CC*cos(m*phi_p)
                 Mo_p(3)=cmplx(0.0_dp,0.0_dp,dp)
                 Me_p(2)=cmplx(0.0_dp,0.0_dp,dp)
                 Me_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*cos(m*phi_p)
                 
                 No_p(1)=cmplx(0.0_dp,0.0_dp,dp)
                 No_p(2)=cmplx(0.0_dp,0.0_dp,dp)
                 No_p(3)=der_p*CC*cos(m*phi_p)/arg_p
                 Ne_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*cos(m*phi_p)/arg_p
                 Ne_p(2)=der_p*DD*cos(m*phi_p)/arg_p
                 Ne_p(3)=cmplx(0.0_dp,0.0_dp,dp)
                 
                 Po=sum(Mo_p(1:3)*orient_dipole(1:3))
                 Pe=sum(Me_p(1:3)*orient_dipole(1:3))
                 
                 Qo=sum(No_p(1:3)*orient_dipole(1:3))
                 Qe=sum(Ne_p(1:3)*orient_dipole(1:3))
                 
                 power=power+const*dr*wght(j)*r_obs**2*factor*conjg(factor)*&
                      (conjg(RC(n,1))*RC(n,1)*(abs(Po)**2+abs(Pe)**2)*&
                      cmplx(cyr(n),cyi(n),dp)*cmplx(cyr(n),-cyi(n),dp)+&
                      conjg(RD(n,1))*RD(n,1)*(abs(Qo)**2+abs(Qe)**2)*&
                      (1.0_dp/real(2*n+1,dp))*&
                      (real(n+1,dp)*cmplx(cyr(n-1),cyi(n-1),dp)*&
                      cmplx(cyr(n-1),-cyi(n-1),dp)+real(n,dp)*&
                      cmplx(cyr(n+1),cyi(n+1),dp)*cmplx(cyr(n+1),-cyi(n+1),dp)))
                 const_den=1.0_dp
                 do m=1,n
                    const_den=const_den*real((n-m+1),dp)*real((n+m),dp)
                    const=2.0_dp*D_fixed/const_den
                    call qq(n,m,theta_p,pp,CC,DD)
                    Mo_p(2)=cmplx(cyr_p(n),cyi_p(n),dp)*CC*cos(m*phi_p)
                    Mo_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*sin(m*phi_p)
                    Me_p(2)=-cmplx(cyr_p(n),cyi_p(n),dp)*CC*sin(m*phi_p)
                    Me_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*cos(m*phi_p)
                    
                    No_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*sin(m*phi_p)&
                   /arg_p
                    No_p(2)=der_p*DD*sin(m*phi_p)/arg_p
                    No_p(3)=der_p*CC*cos(m*phi_p)/arg_p
                    Ne_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*cos(m*phi_p)&
                         /arg_p
                    Ne_p(2)=der_p*DD*cos(m*phi_p)/arg_p
                    Ne_p(3)=-der_p*CC*sin(m*phi_p)/arg_p
                    
                    Po=sum(Mo_p(1:3)*orient_dipole(1:3))
                    Pe=sum(Me_p(1:3)*orient_dipole(1:3))
                    Qo=sum(No_p(1:3)*orient_dipole(1:3))
                    Qe=sum(Ne_p(1:3)*orient_dipole(1:3))
                    
                    power=power+const*dr*wght(j)*r_obs**2*factor*conjg(factor)*&
                         (conjg(RC(n,1))*RC(n,1)*(abs(Po)**2+abs(Pe)**2)*&
                         cmplx(cyr(n),cyi(n),dp)*cmplx(cyr(n),-cyi(n),dp)+&
                         conjg(RD(n,1))*RD(n,1)*(abs(Qo)**2+abs(Qe)**2)*&
                         (1.0_dp/real(2*n+1,dp))*&
                         (real(n+1,dp)*cmplx(cyr(n-1),cyi(n-1),dp)*&
                         cmplx(cyr(n-1),-cyi(n-1),dp)+real(n,dp)*&
                         cmplx(cyr(n+1),cyi(n+1),dp)*cmplx(cyr(n+1),-cyi(n+1),dp)))
                 end do
              end do
              
           end do
        else
           do j=1,nqp
            
              r_obs=a(l-1)+dr*qp(j)
              zr=real(k*r_obs);zi=-aimag(k*r_obs)
              call zbesj(zr,zi,0.5_dp,one,num+2,cyr,cyi,nz,ierr2)
              if (ierr2/=0) then
                 print*,"error in the calculation of the hankel bessel with error &
                      code:",ierr2
                 stop
              end if
              call zbesy(zr,zi,0.5_dp,one,num+2,cyr2,cyi2,nz,&
                   crwrk,ciwrk,ierr)
              factor=sqrt((pi*0.5_dp)/cmplx(zr,zi,dp))
              
              do n=num,1,-1
                 der=(cmplx(zr,zi,dp)*cmplx(cyr(n-1),cyi(n-1),dp)-&
                      n*cmplx(cyr(n),cyi(n),dp))
                 der2=(cmplx(zr,zi,dp)*cmplx(cyr2(n-1),cyi2(n-1),dp)-&
                      n*cmplx(cyr2(n),cyi2(n),dp))
                 der_p=(arg_p*cmplx(cyr_p(n-1),cyi_p(n-1),dp)-&
                      n*cmplx(cyr_p(n),cyi_p(n),dp))
                 D_fixed=0.25_dp*((2.0_dp*n)+1.0_dp)/(n*(n+1.0_dp))
                 
                 m=0
                 const=D_fixed
                 
                 call qq(n,m,theta_p,pp,CC,DD)
                 
                 Mo_p(2)=cmplx(cyr_p(n),cyi_p(n),dp)*CC*cos(m*phi_p)
                 Mo_p(3)=cmplx(0.0_dp,0.0_dp,dp)
                 Me_p(2)=cmplx(0.0_dp,0.0_dp,dp)
                 Me_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*cos(m*phi_p)
                 
                 No_p(1)=cmplx(0.0_dp,0.0_dp,dp)
                 No_p(2)=cmplx(0.0_dp,0.0_dp,dp)
                 No_p(3)=der_p*CC*cos(m*phi_p)/arg_p
                 Ne_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*cos(m*phi_p)/arg_p
                 Ne_p(2)=der_p*DD*cos(m*phi_p)/arg_p
                 Ne_p(3)=cmplx(0.0_dp,0.0_dp,dp)
                 
                 Po=sum(Mo_p(1:3)*orient_dipole(1:3))
                 Pe=sum(Me_p(1:3)*orient_dipole(1:3))
                 
                 Qo=sum(No_p(1:3)*orient_dipole(1:3))
                 Qe=sum(Ne_p(1:3)*orient_dipole(1:3))
                 
                 power=power+const*dr*wght(j)*r_obs**2*factor*conjg(factor)*&
                      (conjg(RC(n,l))*RC(n,l)*(abs(Po)**2+abs(Pe)**2)*&
                      cmplx(cyr(n),cyi(n),dp)*cmplx(cyr(n),-cyi(n),dp)+&
                      conjg(RD(n,l))*RD(n,l)*(abs(Qo)**2+abs(Qe)**2)*&
                      (1.0_dp/real(2*n+1,dp))*&
                      (real(n+1,dp)*cmplx(cyr(n-1),cyi(n-1),dp)*&
                      cmplx(cyr(n-1),-cyi(n-1),dp)+real(n,dp)*&
                      cmplx(cyr(n+1),cyi(n+1),dp)*cmplx(cyr(n+1),-cyi(n+1),dp)))
                 
                 power=power+const*dr*wght(j)*r_obs**2*factor*conjg(factor)*&
                      (conjg(RF(n,l))*RF(n,l)*(abs(Po)**2+abs(Pe)**2)*&
                      cmplx(cyr2(n),cyi2(n),dp)*cmplx(cyr2(n),-cyi2(n),dp)+&
                      conjg(RG(n,l))*RG(n,l)*(abs(Qo)**2+abs(Qe)**2)*&
                      (1.0_dp/real(2*n+1,dp))*&
                      (real(n+1,dp)*cmplx(cyr2(n-1),cyi2(n-1),dp)*&
                      cmplx(cyr2(n-1),-cyi2(n-1),dp)+real(n,dp)*&
                      cmplx(cyr2(n+1),cyi2(n+1),dp)*cmplx(cyr2(n+1),-cyi2(n+1),dp)))
                 
                 power=power+2.0_dp*real(const*dr*wght(j)*r_obs**2*&
                      factor*conjg(factor)*&
                      (conjg(RC(n,l))*RF(n,l)*(abs(Po)**2+abs(Pe)**2)*&
                      cmplx(cyr(n),cyi(n),dp)*cmplx(cyr2(n),-cyi2(n),dp)+&
                      conjg(RD(n,l))*RG(n,l)*(abs(Qo)**2+abs(Qe)**2)*&
                      (1.0_dp/real(2*n+1,dp))*&
                      (real(n+1,dp)*cmplx(cyr(n-1),cyi(n-1),dp)*&
                      cmplx(cyr2(n-1),-cyi2(n-1),dp)+real(n,dp)*&
                      cmplx(cyr(n+1),cyi(n+1),dp)*cmplx(cyr2(n+1),-cyi2(n+1),dp))))
                 
                 
                 const_den=1.0_dp
                 do m=1,n
                    const_den=const_den*real((n-m+1),dp)*real((n+m),dp)
                    const=2.0_dp*D_fixed/const_den
                    call qq(n,m,theta_p,pp,CC,DD)
                    Mo_p(2)=cmplx(cyr_p(n),cyi_p(n),dp)*CC*cos(m*phi_p)
                    Mo_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*sin(m*phi_p)
                    Me_p(2)=-cmplx(cyr_p(n),cyi_p(n),dp)*CC*sin(m*phi_p)
                    Me_p(3)=-cmplx(cyr_p(n),cyi_p(n),dp)*DD*cos(m*phi_p)
                    
                    No_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*sin(m*phi_p)&
                         /arg_p
                    No_p(2)=der_p*DD*sin(m*phi_p)/arg_p
                    No_p(3)=der_p*CC*cos(m*phi_p)/arg_p
                    Ne_p(1)=n*(n+1.0_dp)*cmplx(cyr_p(n),cyi_p(n),dp)*pp*cos(m*phi_p)&
                         /arg_p
                    Ne_p(2)=der_p*DD*cos(m*phi_p)/arg_p
                    Ne_p(3)=-der_p*CC*sin(m*phi_p)/arg_p
                    
                    Po=sum(Mo_p(1:3)*orient_dipole(1:3))
                    Pe=sum(Me_p(1:3)*orient_dipole(1:3))
                    Qo=sum(No_p(1:3)*orient_dipole(1:3))
                    Qe=sum(Ne_p(1:3)*orient_dipole(1:3))
                    
                    
                    power=power+const*dr*wght(j)*r_obs**2*factor*conjg(factor)*&
                         (conjg(RC(n,l))*RC(n,l)*(abs(Po)**2+abs(Pe)**2)*&
                         cmplx(cyr(n),cyi(n),dp)*cmplx(cyr(n),-cyi(n),dp)+&
                         conjg(RD(n,l))*RD(n,l)*(abs(Qo)**2+abs(Qe)**2)*&
                         (1.0_dp/real(2*n+1,dp))*&
                         (real(n+1,dp)*cmplx(cyr(n-1),cyi(n-1),dp)*&
                         cmplx(cyr(n-1),-cyi(n-1),dp)+real(n,dp)*&
                         cmplx(cyr(n+1),cyi(n+1),dp)*cmplx(cyr(n+1),-cyi(n+1),dp)))
                    
                    power=power+const*dr*wght(j)*r_obs**2*factor*conjg(factor)*&
                         (conjg(RF(n,l))*RF(n,l)*(abs(Po)**2+abs(Pe)**2)*&
                         cmplx(cyr2(n),cyi2(n),dp)*cmplx(cyr2(n),-cyi2(n),dp)+&
                         conjg(RG(n,l))*RG(n,l)*(abs(Qo)**2+abs(Qe)**2)*&
                         (1.0_dp/real(2*n+1,dp))*&
                         (real(n+1,dp)*cmplx(cyr2(n-1),cyi2(n-1),dp)*&
                         cmplx(cyr2(n-1),-cyi2(n-1),dp)+real(n,dp)*&
                         cmplx(cyr2(n+1),cyi2(n+1),dp)*&
                         cmplx(cyr2(n+1),-cyi2(n+1),dp)))
                    
                    power=power+2.0_dp*real(const*dr*wght(j)*r_obs**2*&
                         factor*conjg(factor)*&
                         (conjg(RC(n,l))*RF(n,l)*(abs(Po)**2+abs(Pe)**2)*&
                         cmplx(cyr(n),cyi(n),dp)*cmplx(cyr2(n),-cyi2(n),dp)+&
                         conjg(RD(n,l))*RG(n,l)*(abs(Qo)**2+abs(Qe)**2)*&
                         (1.0_dp/real(2*n+1,dp))*&
                         (real(n+1,dp)*cmplx(cyr(n-1),cyi(n-1),dp)*&
                         cmplx(cyr2(n-1),-cyi2(n-1),dp)+real(n,dp)*&
                         cmplx(cyr(n+1),cyi(n+1),dp)*&
                         cmplx(cyr2(n+1),-cyi2(n+1),dp))))
                 end do
                 
              end do
           end do
        
        end if
        power=power*pi*0.5_dp*aimag(ri*ri)*((0.5_dp)/(lambda*1.0E-7_dp*vlite))*&
             abs(const_dipole)**2
        write(12,*) "Numerically calculated Absorbed power in layer",l,&
             "=",power
     else
        write(12,*) "not coded for dipole inside yet"
     end if
     deallocate(cyr,cyi,cyr2,cyi2,crwrk,ciwrk,wght,qp)
  else
     write(12,*) "not coded the Power computation for Magnetic dipole yet!"
  end if

  return
end subroutine Mie_power



 subroutine sph2cart(theta,phi,Esp,Ecart)
   use mie_parameters
   implicit none
   real(kind=dp),intent(in)::theta,phi
   complex(kind=dp),intent(in)::Esp(3)
   complex(kind=dp),intent(out)::Ecart(3)


   Ecart(1)=(Esp(1)*sin(theta)+Esp(2)*cos(theta))*cos(phi)-Esp(3)*sin(phi)
   Ecart(2)=(Esp(1)*sin(theta)+Esp(2)*cos(theta))*sin(phi)+Esp(3)*cos(phi)
   Ecart(3)=(Esp(1)*cos(theta)-Esp(2)*sin(theta))
   
   return
 end subroutine sph2cart
 
 subroutine sph2cart_real(theta,phi,usp,ucart)
   use mie_parameters
   implicit none
   real(kind=dp),intent(in)::theta,phi
   real(kind=dp),intent(in)::usp(3)
   real(kind=dp),intent(out)::ucart(3)


   ucart(1)=(usp(1)*sin(theta)+usp(2)*cos(theta))*cos(phi)-usp(3)*sin(phi)
   ucart(2)=(usp(1)*sin(theta)+usp(2)*cos(theta))*sin(phi)+usp(3)*cos(phi)
   ucart(3)=(usp(1)*cos(theta)-usp(2)*sin(theta))
   
   return
 end subroutine sph2cart_real

 subroutine sph2cart_coord(rthephi,xyz)
   use mie_parameters
   implicit none
   real(kind=dp),intent(in)::rthephi(3)
   real(kind=dp),intent(out)::xyz(3)

   
   xyz(1)=rthephi(1)*sin(rthephi(2))*cos(rthephi(3))
   xyz(2)=rthephi(1)*sin(rthephi(2))*sin(rthephi(3))
   xyz(3)=rthephi(1)*cos(rthephi(2))
   
   return
 end subroutine sph2cart_coord


 subroutine cart2sph(cart,sp)
   use mie_parameters
   implicit none
   real(kind=dp),intent(in)::cart(3)
   real(kind=dp),intent(out)::sp(3)
   real(kind=dp)::tmp
   sp(1)=sqrt(dot_product(cart,cart))
   sp(2)=acos(cart(3)/sp(1))
   tmp=sqrt(cart(1)**2+cart(2)**2)
   if (tmp>10.0E-12_dp) then
      sp(3)=asin(cart(2)/tmp)
   else
      sp(3)=0.0_dp
   end if
  
   if(cart(1)<0.and.abs(cart(1))>10.0E-12_dp) sp(3)=pi-sp(3)

   return
 end subroutine cart2sph


 subroutine oneD_quadrature(qrule,qp,wght)   
   use mie_variables
   implicit none
   integer,intent(in)::qrule
   real(kind=dp),intent(out),dimension(1:qrule)::qp,wght
   select case (qrule)
   case(1)
      qp(1)=0.5_dp
      wght(1)=1.0_dp
   case(2)
      qp(1:2)=(/0.211324865405187_dp,0.788675134594813_dp/)
      wght(1:2)=(/0.5_dp,0.5_dp/)
   case(3)
      qp(1:3)=(/0.5_dp,0.112701665379258_dp,0.887298334620742_dp/)
      wght(1:3) = (/0.444444444444444_dp,0.277777777777776_dp,&
           0.277777777777776_dp/)
   case(4)
      qp(1:4)= (/0.330009478207572_dp,6.943184420297372E-002_dp,&
           0.669990521792428_dp,0.930568155797026_dp/)
      wght(1:4)=(/0.326072577431273_dp,0.173927422568727_dp,&
           0.326072577431273_dp,0.173927422568727_dp/)
   case(5)
      qp(1:5)=(/0.5_dp,4.691007703066802E-002_dp,0.230765344947158_dp,&
           0.769234655052841_dp,0.953089922969332_dp/)
      wght(1:5)=(/0.284444444444444_dp,0.118463442528095_dp,&
           0.239314335249683_dp,0.239314335249683_dp,&
           0.118463442528095_dp/)
   case(6)
      qp(1:6)=(/0.619309593041598_dp,3.376524289842397E-002_dp,&
           0.169395306766868_dp,0.380690406958402_dp,&
           0.830604693233132_dp,0.966234757101576_dp/)
      wght(1:6)=(/0.233956967286346_dp,8.566224618958525E-002_dp,&
           0.180380786524069_dp,0.233956967286346_dp,&
           0.180380786524069_dp,8.566224618958525E-002_dp/)
   case(7)
      qp(1:7)=(/0.5_dp,2.544604382862076E-002_dp,0.129234407200303_dp,&
           0.297077424311301_dp,0.702922575688699_dp,&
           0.870765592799697_dp,0.974553956171379_dp/)
      wght(1:7)=(/0.208979591836735_dp,6.474248308443483E-002_dp,&
           0.139852695744638_dp,0.190915025252560_dp,&
           0.190915025252560_dp,0.139852695744638_dp,&
           6.474248308443483E-002_dp/)
   case(8)
      qp(1:8)=(/0.408282678752175_dp,1.985507175123186E-002_dp,&
           0.101666761293187_dp,0.237233795041836_dp,&
           0.591717321247825_dp,0.762766204958164_dp,&
           0.898333238706813_dp,0.980144928248768_dp/)
      wght(1:8)=(/0.181341891689181_dp,5.061426814518839E-002_dp,&
           0.111190517226687_dp,0.156853322938944_dp,&
           0.181341891689181_dp,0.156853322938944_dp,&
           0.111190517226687_dp,5.061426814518839E-002_dp/)
   case(9)
      qp(1:9)=(/0.5_dp,1.591988024618696E-002_dp,8.198444633668206E-002_dp,&
           0.193314283649705_dp,0.337873288298096_dp,&
           0.662126711701904_dp,0.806685716350295_dp,&
           0.918015553663318_dp,0.984080119753813_dp/)
      wght(1:9)=(/0.165119677500630_dp,4.063719418078731E-002_dp,&
           9.032408034742879E-002_dp,0.130305348201468_dp,&
           0.156173538520001_dp,0.156173538520001_dp,&
           0.130305348201468_dp,9.032408034742879E-002_dp,&
           4.063719418078731E-002_dp/)
   case(10)
      qp(1:10)=(/0.574437169490816_dp,1.304673574141413E-002_dp,&
           6.746831665550773E-002_dp,0.160295215850488_dp,&
           0.283302302935376_dp,0.425562830509184_dp,&
           0.716697697064624_dp,0.839704784149512_dp,&
           0.932531683344492_dp,0.986953264258586_dp/)
      wght(1:10)=(/0.147762112357376_dp,3.333567215434186E-002_dp,&
           7.472567457529027E-002_dp,0.109543181257991_dp,&
           0.134633359654998_dp,0.147762112357376_dp,&
           0.134633359654998_dp,0.109543181257991_dp,&
           7.472567457529027E-002_dp,3.333567215434186E-002_dp/)
   case(11)
      qp(1:11)=(/1.088567092697151E-002_dp,5.646870011595234E-002_dp,&
           0.134923997212975_dp,0.240451935396594_dp,0.365228422023828_dp,&
           0.500000000000000_dp,0.634771577976172_dp,0.759548064603406_dp,&
           0.865076002787025_dp,0.943531299884048_dp,0.989114329073028_dp/)
      wght(1:11)=(/ 2.783428355808480E-002_dp,6.279018473245236E-002_dp,&
           9.314510546386713E-002_dp,0.116596882295989_dp,0.131402272255123_dp,&
           0.136462543388950_dp,0.131402272255123_dp,0.116596882295989_dp,&
           9.314510546386713E-002_dp,6.279018473245236E-002_dp,&
           2.783428355808480E-002_dp/)
   case(12)
      qp(1:12)=(/9.219682876640378E-003_dp,4.794137181476260E-002_dp,&
           0.115048662902848_dp,0.206341022856691_dp,0.316084250500910_dp,&
           0.437383295744266_dp,0.562616704255734_dp,0.683915749499090_dp,&
           0.793658977143309_dp,0.884951337097152_dp,0.952058628185237_dp,&
           0.990780317123360_dp/)
      wght(1:12)=(/2.358766819325391E-002_dp,5.346966299765909E-002_dp,&
           8.003916427167317E-002_dp,0.101583713361533_dp,0.116746268269177_dp,&
           0.124573522906701_dp,0.124573522906701_dp,0.116746268269177_dp,&
           0.101583713361533_dp,8.003916427167317E-002_dp,&
           5.346966299765909E-002_dp,2.358766819325391E-002_dp/)
   case(13)
      qp(1:13)=(/7.908472640705932E-003_dp,4.120080038851098E-002_dp,&
           9.921095463334506E-002_dp,0.178825330279830_dp,0.275753624481777_dp,&
           0.384770842022433_dp,0.500000000000000_dp,0.615229157977567_dp,&
           0.724246375518223_dp,0.821174669720170_dp,0.900789045366655_dp,&
           0.958799199611489_dp,0.992091527359294_dp/)
      wght(1:13)=(/2.024200238265615E-002_dp,4.606074991886419E-002_dp,&
           6.943675510989360E-002_dp,8.907299038097284E-002_dp,&
           0.103908023768444_dp,0.113141590131449_dp,&
           0.116275776615437_dp,0.113141590131449_dp,&
           0.103908023768444_dp,8.907299038097284E-002_dp,&
           6.943675510989360E-002_dp,4.606074991886419E-002_dp,&
           2.024200238265615E-002_dp/)
   case(14)
      qp(1:14)=(/6.858095651593843E-003_dp,3.578255816821319E-002_dp,&
           8.639934246511749E-002_dp,0.156353547594157_dp,&
           0.242375681820923_dp,&
           0.340443815536055_dp,0.445972525646328_dp,0.554027474353672_dp,&
           0.659556184463945_dp,0.757624318179077_dp,0.843646452405843_dp,&
           0.913600657534883_dp,0.964217441831787_dp,0.993141904348406_dp/)
      wght(1:14)=(/1.755973016587453E-002_dp,4.007904357988019E-002_dp,&
           6.075928534395156E-002_dp,7.860158357909680E-002_dp,&
           9.276919873896822E-002_dp,0.102599231860648_dp,&     
           0.107631926731579_dp,0.107631926731579_dp,&
           0.102599231860648_dp,9.276919873896822E-002_dp,&
           7.860158357909680E-002_dp,6.075928534395156E-002_dp,&
           4.007904357988019E-002_dp,1.755973016587453E-002_dp/)
   case(15)
      qp(1:15)=(/6.003740989757256E-003_dp,3.136330379964702E-002_dp,&
           7.589670829478640E-002_dp,0.137791134319915_dp,&
           0.214513913695731_dp,&
           0.302924326461218_dp,0.399402953001283_dp,0.500000000000000_dp,&
           0.600597046998717_dp,0.697075673538782_dp,0.785486086304269_dp,&
           0.862208865680085_dp,0.924103291705214_dp,0.968636696200353_dp,&
           0.993996259010243_dp/)
      wght(1:15)=(/1.537662099805739E-002_dp,3.518302374405414E-002_dp,&
           5.357961023358602E-002_dp,6.978533896307713E-002_dp,&
           8.313460290849568E-002_dp,9.308050000778112E-002_dp,&
           9.921574266355579E-002_dp,0.101289120962781_dp,&
           9.921574266355579E-002_dp,9.308050000778112E-002_dp,&
           8.313460290849568E-002_dp,6.978533896307713E-002_dp,&
           5.357961023358602E-002_dp,3.518302374405414E-002_dp,&
           1.537662099805739E-002_dp/)
      
   case(16)
      qp(1:16)= (/5.299532504175031E-003_dp,2.771248846338370E-002_dp,&
           6.718439880608412E-002_dp,0.122297795822499_dp,&
           0.191061877798678_dp,0.270991611171386_dp,0.359198224610371_dp,&
           0.452493745081181_dp,0.547506254918819_dp,0.640801775389629_dp,&
           0.729008388828614_dp,0.808938122201322_dp,0.877702204177502_dp,&
           0.932815601193916_dp,0.972287511536616_dp,0.994700467495825_dp/)
      wght(1:16)=(/1.357622970587593E-002_dp,3.112676196932389E-002_dp,&
           4.757925584124645E-002_dp,6.231448562776697E-002_dp,&
           7.479799440828665E-002_dp,8.457825969750121E-002_dp,&
           9.130170752246181E-002_dp,9.472530522753424E-002_dp,&
           9.472530522753424E-002_dp,9.130170752246181E-002_dp,&
           8.457825969750121E-002_dp,7.479799440828665E-002_dp,&
           6.231448562776697E-002_dp,4.757925584124645E-002_dp,&
           3.112676196932389E-002_dp,1.357622970587593E-002_dp/)
      
   case(17)
      qp(1:17)=(/4.712262342791318E-003_dp,2.466223911561610E-002_dp,&
           5.988042313650704E-002_dp,0.109242998051599_dp,&
           0.171164420391655_dp,0.243654731456762_dp,0.324384118273062_dp,&
           0.410757909252076_dp,0.500000000000000_dp,0.589242090747924_dp,&
           0.675615881726938_dp,0.756345268543238_dp,0.828835579608345_dp,&
           0.890757001948401_dp,0.940119576863493_dp,0.975337760884384_dp,&
           0.995287737657209_dp/)
      wght(1:17)=(/1.207415143427312E-002_dp,2.772976468699357E-002_dp,&
           4.251807415858957E-002_dp,5.594192359670196E-002_dp,&
           6.756818423426052E-002_dp,7.702288053840499E-002_dp,&
           8.400205107822495E-002_dp,8.828135268349632E-002_dp,&
           8.972323517810327E-002_dp,8.828135268349632E-002_dp,&
           8.400205107822495E-002_dp,7.702288053840499E-002_dp,&
           6.756818423426052E-002_dp,5.594192359670196E-002_dp,&
           4.251807415858957E-002_dp,2.772976468699357E-002_dp,&
           1.207415143427312E-002_dp/)
      
   case(18)
      qp(1:18)=(/4.217415789534495E-003_dp,2.208802521430114E-002_dp,&
           5.369876675122215E-002_dp,9.814752051373843E-002_dp,&
           0.154156478469823_dp,0.220114584463026_dp,0.294124419268579_dp,&
           0.374056887154247_dp,0.457612493479132_dp,0.542387506520868_dp,&
           0.625943112845753_dp,0.705875580731421_dp,0.779885415536974_dp,&
           0.845843521530177_dp,0.901852479486262_dp,0.946301233248778_dp,&
           0.977911974785699_dp,0.995782584210466_dp/)
      wght(1:18)=(/1.080800676324070E-002_dp,2.485727444748492E-002_dp,&
           3.821286512744459E-002_dp,5.047102205314361E-002_dp,&
           6.127760335573650E-002_dp,7.032145733532516E-002_dp,&
           7.734233756313259E-002_dp,8.213824187291634E-002_dp,&
           8.457119148157177E-002_dp,8.457119148157177E-002_dp,&
           8.213824187291634E-002_dp,7.734233756313259E-002_dp,&
           7.032145733532516E-002_dp,6.127760335573650E-002_dp,&
           5.047102205314361E-002_dp,3.821286512744459E-002_dp,&
           2.485727444748492E-002_dp,1.080800676324070E-002_dp/)

   case(19)
      qp(1:19)=(/3.796578078207824E-003_dp,1.989592393258499E-002_dp,&
           4.842204819259105E-002_dp,8.864267173142859E-002_dp,&
           0.139516911332385_dp,0.199727347669160_dp,0.267714629312020_dp,&
           0.341717950018185_dp,0.419820677179887_dp,0.500000000000000_dp,&
           0.580179322820113_dp,0.658282049981815_dp,0.732285370687980_dp,&
           0.800272652330841_dp,0.860483088667615_dp,0.911357328268571_dp,&
           0.951577951807409_dp,0.980104076067415_dp,0.996203421921792_dp/)
      wght(1:19)=(/9.730894114862422E-003_dp,2.240711338284979E-002_dp,&
           3.452227136882063E-002_dp,4.574501081122507E-002_dp,&
           5.578332277366375E-002_dp,6.437698126966783E-002_dp,&
           7.130335108680329E-002_dp,7.638302103292986E-002_dp,&
           7.948442169697721E-002_dp,8.052722492439185E-002_dp,&
           7.948442169697721E-002_dp,7.638302103292986E-002_dp,&
           7.130335108680329E-002_dp,6.437698126966783E-002_dp,&
           5.578332277366375E-002_dp,4.574501081122507E-002_dp,&
           3.452227136882063E-002_dp,2.240711338284979E-002_dp,&
           9.730894114862422E-003_dp/)

   case(20)
      qp(1:20)=(/3.435700407452558E-003_dp,1.801403636104310E-002_dp,&
           4.388278587433703E-002_dp,8.044151408889055E-002_dp,&
           0.126834046769925_dp,0.181973159636742_dp,0.244566499024586_dp,&
           0.313146955642290_dp,0.386107074429177_dp,0.461736739433251_dp,&
           0.538263260566749_dp,0.613892925570823_dp,0.686853044357710_dp,&
           0.755433500975414_dp,0.818026840363258_dp,0.873165953230075_dp,&
           0.919558485911109_dp,0.956117214125663_dp,&     
           0.981985963638957_dp,0.996564299592547_dp/)
      wght(1:20)=(/8.807003569575289E-003_dp,2.030071490019352E-002_dp,&
           3.133602416705452E-002_dp,4.163837078835237E-002_dp,&
           5.096505990861661E-002_dp,5.909726598075887E-002_dp,&
           6.584431922458830E-002_dp,7.104805465919108E-002_dp,&
           7.458649323630188E-002_dp,7.637669356536299E-002_dp,&
           7.637669356536299E-002_dp,7.458649323630188E-002_dp,&
           7.104805465919108E-002_dp,6.584431922458830E-002_dp,&
           5.909726598075887E-002_dp,5.096505990861661E-002_dp,&
           4.163837078835237E-002_dp,3.133602416705452E-002_dp,&
           2.030071490019352E-002_dp,8.807003569575289E-003_dp/)
      
   case(24)
      qp(1:24)=(/2.406390001489289E-003_dp,1.263572201434526E-002_dp,&
           3.086272399863360E-002_dp,5.679223649779946E-002_dp,&
           8.999900701304853E-002_dp,0.129937904210723_dp,&
           0.175953174031512_dp,0.227289264305580_dp,0.283103246186977_dp,&     
           0.342478660151918_dp,0.404440566263192_dp,0.467971553568697_dp,&
           0.532028446431303_dp,0.595559433736808_dp,0.657521339848082_dp,&     
           0.716896753813023_dp,0.772710735694420_dp,0.824046825968488_dp,&     
           0.870062095789277_dp,0.910000992986951_dp,0.943207763502200_dp,&     
           0.969137276001366_dp,0.987364277985655_dp,0.997593609998511_dp/)
      wght(1:24)=(/6.170614899993272E-003_dp,1.426569431446691E-002_dp,&
           2.213871940870984E-002_dp,2.964929245771833E-002_dp,&
           3.667324070553566E-002_dp,4.309508076597601E-002_dp,&
           4.880932605205684E-002_dp,5.372213505798281E-002_dp,&
           5.775283402686281E-002_dp,6.083523646390167E-002_dp,&
           6.291872817341419E-002_dp,6.396909767337612E-002_dp,&
           6.396909767337612E-002_dp,6.291872817341419E-002_dp,&
           6.083523646390167E-002_dp,5.775283402686281E-002_dp,&
           5.372213505798281E-002_dp,4.880932605205684E-002_dp,&
           4.309508076597601E-002_dp,3.667324070553566E-002_dp,&
           2.964929245771833E-002_dp,2.213871940870984E-002_dp,&
           1.426569431446691E-002_dp,6.170614899993272E-003_dp/)
      
   case(32)
      qp(1:32)=(/1.368069075259215E-003_dp,7.194244227365809E-003_dp,&
           1.761887220624681E-002_dp,3.254696203113017E-002_dp,&
           5.183942211697395E-002_dp,7.531619313371501E-002_dp,&
           0.102758102016029_dp,0.133908940629855_dp,0.168477866534892_dp,&     
           0.206142121379619_dp,0.246550045533885_dp,0.289324361934682_dp,&     
           0.334065698858936_dp,0.380356318873931_dp,0.427764019208602_dp,&     
           0.475846167156131_dp,0.524153832843869_dp,0.572235980791398_dp,&     
           0.619643681126068_dp,0.665934301141064_dp,0.710675638065318_dp,&     
           0.753449954466115_dp,0.793857878620381_dp,0.831522133465108_dp,&     
           0.866091059370145_dp,0.897241897983971_dp,0.924683806866285_dp,&     
           0.948160577883026_dp,0.967453037968870_dp,0.982381127793753_dp,&     
           0.992805755772634_dp,0.998631930924741_dp/)
      wght(1:32)=(/3.509305004734761E-003_dp,8.137197365452854E-003_dp,&
           1.269603265463107E-002_dp,1.713693145651071E-002_dp,&
           2.141794901110915E-002_dp,2.549902963118725E-002_dp,&
           2.934204673926756E-002_dp,3.291111138818090E-002_dp,&
           3.617289705442431E-002_dp,3.909694789353522E-002_dp,&
           4.165596211347336E-002_dp,4.382604650220187E-002_dp,&
           4.558693934788195E-002_dp,4.692219954040221E-002_dp,&
           4.781936003963742E-002_dp,4.827004425736393E-002_dp,&
           4.827004425736393E-002_dp,4.781936003963742E-002_dp,&
           4.692219954040221E-002_dp,4.558693934788195E-002_dp,&
           4.382604650220187E-002_dp,4.165596211347336E-002_dp,&
           3.909694789353522E-002_dp,3.617289705442431E-002_dp,&
           3.291111138818090E-002_dp,2.934204673926756E-002_dp,&
           2.549902963118725E-002_dp,2.141794901110915E-002_dp,&
           1.713693145651071E-002_dp,1.269603265463107E-002_dp,&
           8.137197365452854E-003_dp,3.509305004734761E-003_dp/)

   end select
   return
 end subroutine oneD_quadrature
!!$ 
