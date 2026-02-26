!
!----------------------------------------------------------------------
!
!     routine name      - mfd_solutions
!
!----------------------------------------------------------------------
!
!     latest revision:  - June 18
!
!> @brief         - compute all relevant quantities of the exact
!                         solutions for the Maxwell problem
!
!     arguments:
!
!     in:
!           Mdle        - an element (middle node) number
!             Xp        - a point in physical space
!     out:
!              p        - value of the solution pressure
!          Gradp        - corresponding first derivatives
!         Grad2p        - corresponding second derivatives
!
!----------------------------------------------------------------------
!
   subroutine mfd_solutions(Mdle,Xp, p,Gradp,Grad2p)
!
      use data_structure3D
      use commonParam
      use bessel_evaluation
!
      implicit none
!
      integer,    intent(in)  :: Mdle
      real(8),    intent(in)  :: Xp(3)
      complex(8), intent(out) :: p
      complex(8), intent(out) :: Gradp(3)
      complex(8), intent(out) :: Grad2p(3,3)
!
!  ...intermediate variables
      real(8) :: w0, p0, rk, pi_mod
      real(8) :: theta_x, theta_y
      real(8) :: rad_x, rad_y
      real(8) :: sinx, siny, cosx, cosy
      real(8) :: x, y, z
      real(8) :: Rx(3,3), Ry(3,3), Rxy(3,3)
!
!  ...separable solution - factor functions and their derivatives
      complex(8) :: u, du_r,d2u_r ,znu
      complex(8) :: w, dw_th, d2w_th
      real(8) :: v,dv_x1,d2v_x1
!
      complex(8) :: zf_x, zf_y, zf_z, cn
      complex(8) :: dzf_x,dzf_y, dzf_z
      complex(8) :: ddzf_x, ddzf_y, ddzf_z
!
      complex(8) :: p_x, p_xx, p_xy, p_xz
      complex(8) :: p_y, p_yx, p_yy, p_yz
      complex(8) :: p_z, p_zx, p_zy, p_zz
!
      real(8) :: r_x, r_xx, r_xy, r_xz
      real(8) :: r_y, r_yx, r_yy, r_yz
      real(8) :: r_z, r_zx, r_zy, r_zz
!
      real(8) :: dxdxs, dxdys, dxdzs
      real(8) :: dydxs, dydys, dydzs
      real(8) :: dzdxs, dzdys, dzdzs
!
      real(8) :: x1, x2, x3, xshift, yshift, zshift, alpha, a, b, cf
      real(8) :: r,dr_x2,dr_x3,d2r_x2,d2r_x3,d2r_x2x3
      real(8) :: th,dth_x2,dth_x3,d2th_x2,d2th_x3,d2th_x2x3

      real(8) :: kappa_clad,kappa_core

      real(8) :: wavenum0,wavenum1
      complex(8) :: zcurv_st(3),zdcurv_st(3),zd2curv_st(3)
      complex(8):: zJ(3,3),zJinv(3,3),zJdet
      real(8) :: rQ(3,3), curv(3)
      integer :: activePML
!
!---------------------------------------------------------------------------------------
!  ...initialize outputs
      p = ZERO; Gradp = ZERO; Grad2p = ZERO
!  ...extract coordinates
      x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)
!
      select case (ISOL)
!     Constant function
      case(0)
!
         p = ZONE + ZI
         Gradp = ZERO
         Grad2p = ZERO
!
!  ...affine function depending on x (normal to the bending plane yz)
      case(1)
         cn = cmplx(0.5d0,0.d0 , 8 )
         a = 1.d0;    b =1.d0
         p = cn*(a*x1+b)
!     ...1st order derivatives
         Gradp(1) = cn*a
!     ...second order derivatives are all zero
!
!  ...quadratic bubble depending on x (normal to the bending plane yz)
      case(2)
         cn = 1.d0*ZONE
         p = cn*(x1-x1**2)
!     ...1st order derivatives
         Gradp(1) = cn*(1.d0-2.d0*x1)
!     ...second order derivatives
         Grad2p(1,1) = cn*(-2.d0)
!
!  ...polynomial depending on r^2 = y^2+z^2
      case(3)
         cn = 1.d0*ZONE
         p = cn*(x2**2 + x3**2 - RBEND**2)
!     ...1st order derivatives
         Gradp(2) = cn*2.d0*x2
         Gradp(3) = cn*2.d0*x3
!     ...second order derivatives
         Grad2p(2,2) = cn*2.d0
         Grad2p(3,3) = cn*2.d0
!
!  ...linear w.r.t   r = sqrt(y^2+z^2)
      case(4)
         r = dsqrt(x2**2+x3**2)
         dr_x2 = x2/r
         dr_x3 = x3/r
         d2r_x2 = 1.d0/r-1.d0*x2**2/r**3
         d2r_x3 = 1.d0/r-1.d0*x3**2/r**3
         d2r_x2x3 = -1.d0*x2*x3/r**3
!
         cn = 1.d0*ZONE
         p = cn*r
!     ...1st order derivatives
         Gradp(2) = cn*dr_x2
         Gradp(3) = cn*dr_x3
!     ...second order derivatives
         Grad2p(2,2) = cn*d2r_x2
         Grad2p(3,2) = cn*d2r_x2x3
         Grad2p(2,3) = Grad2p(3,2)
         Grad2p(3,3) = cn*d2r_x3
!  ...rational function depending on r = sqrt(y^2+z^2)
      case(5)
            r = dsqrt(x2**2+x3**2)
            dr_x2 = 2.d0*x2/r
            dr_x3 = 2.d0*x3/r
            d2r_x2 = 2.d0/r-4.d0*x2**2/r**3
            d2r_x3 = 2.d0/r-4.d0*x3**2/r**3
            d2r_x2x3 = -4.d0*x2*x3/r**3
!
            cn = RBEND*ZONE
            u = r**(-1)
            du_r = -r**(-2)
            d2u_r = 2.d0*r**(-3)
            p = cn*u
!     ...1st order derivatives
            Gradp(2) = cn*du_r*dr_x2
            Gradp(3) = cn*du_r*dr_x3
!     ...second order derivatives
            Grad2p(2,2) = cn*(d2u_r*dr_x2**2+du_r*d2r_x2)
            Grad2p(3,2) = cn*(d2u_r*dr_x2*dr_x3+du_r*d2r_x2x3)
            Grad2p(2,3) = Grad2p(3,2)
            Grad2p(3,3) = cn*(d2u_r*dr_x3**2+du_r*d2r_x3)      
!
!  ...Function for partly bent waveguide
      case(6)
         if (x3.ge.0.d0) then
            r = dsqrt(x2**2+x3**2)
            dr_x2 = x2/r
            dr_x3 = x3/r
            d2r_x2 = 1.d0/r-1.d0*x2**2/r**3
            d2r_x3 = 1.d0/r-1.d0*x3**2/r**3
            d2r_x2x3 = -1.d0*x2*x3/r**3
         else
            r = x2
            dr_x2 = 1.d0
            dr_x3 = 0.d0
            d2r_x2 = 0.d0
            d2r_x3 = 0.d0
            d2r_x2x3 = 0.d0
         endif
!
         cn = ZONE
         cf = 0.5d0*PI/RCORE
         u = COS((r-RBEND) * cf)
         du_r = -SIN((r-RBEND)*cf) * cf
         d2u_r = -u * cf**2
         v = 1.d0
         dv_x1 = 0.d0
         d2v_x1 = 0.d0
! 
         p = cn * u * v
!     ...1st order derivatives
         Gradp(1) = cn * u * dv_x1
         Gradp(2) = cn * du_r*dr_x2 * v
         Gradp(3) = cn * du_r*dr_x3 * v
!     ...second order derivatives
         Grad2p(1,1) = cn * u * d2v_x1
         Grad2p(1,2) = cn * du_r*dr_x2 * dv_x1
         Grad2p(1,3) = cn * du_r*dr_x3 * dv_x1
         Grad2p(2,1) = Grad2p(1,2)
         Grad2p(2,2) = cn * (d2u_r*dr_x2**2+du_r*d2r_x2) * v
         Grad2p(2,3) = cn * (d2u_r*dr_x2*dr_x3+du_r*d2r_x2x3) * v
         Grad2p(3,1) = Grad2p(1,3)
         Grad2p(3,2) = Grad2p(2,3)
         Grad2p(3,3) = cn * (d2u_r*dr_x3**2+du_r*d2r_x3) * v
!
!  ...Separable function u(r)*v(x)*w(theta), which is a mode of a bent vacuum waveguide with PEC walls: 11-14
!     available for RBEND = 1300 ONLY
      case(11,12,13,14)
         r = dsqrt(x2**2+x3**2)
         dr_x2 = x2/r
         dr_x3 = x3/r
         d2r_x2 = 1.d0/r-1.d0*x2**2/r**3
         d2r_x3 = 1.d0/r-1.d0*x3**2/r**3
         d2r_x2x3 = -1.d0*x2*x3/r**3
!            
!         ...factor u(r) evaluated with the AAA interpolation of preset Bessel functions
         ! call real_eval_aaa(r - RBEND,u,du_r)
!
!        ...factor u(r) evaluated with direct computation of preset Bessel function
         wavenum0 = OMEGA*sqrt(EPSILON*MU)
!        ...u(r) is an eigenfunction of Bessel's equation with appropriate b.c.s
         call bessel_preset(wavenum0, RBEND, r,    u, du_r, d2u_r)
         ! 
         ! write (*,*) 'mfd_solutions: local wavenum0=',wavenum0
         ! call pause
         ! write(*,*) 'mfd_solutions: u,du_r,d2u_r=',u,du_r,d2u_r

         ! u = ZONE; du_r = ZERO; d2u_r = ZERO

         cf = 0.5d0*PI/RCORE
         v =     1.d0 ! COS(x1*cf)
         dv_x1 = 0.d0 !-SIN(x*cf)*cf
         d2v_x1 = 0.d0


!        ...check if we're within the PML      ARE WE ASSUMING THAT THIS IS EVALUATED AT THE BOTTOM FACE?
         call is_pml(Mdle,Xp,activePML)
         ! if so, modify input coordinate r
         if (activePML.ne.0) then 
            ! get rotation matrix
            call get_local_rotation(Mdle,Xp,activePML,rQ)
            ! get pml stretched coordinates and derivatives
            call get_stretched_coords(Mdle,Xp,ActivePML,Zcurv_st,Zdcurv_st,Zd2curv_st)
            ! get PML Jacobian
            call get_stretch_J(zdcurv_st,rQ,zJ,zJinv,zJdet)
         else
            call cartesian2curvilinear_real(Xp,curv)
            zcurv_st = cmplx(curv,0.d0 , 8)
            zJ = cmplx(IDENTITY,0.d0 , 8)
            zJinv = zJ
            ZJdet = ZONE
         endif

         ! FACTOR w(theta)
         dth_x2    = -x3/r**2
         dth_x3    =  x2/r**2
         d2th_x2   =  x3*2.d0*r*dr_x2/r**4
         d2th_x2x3 = -(r**2-x3*2.d0*r*dr_x3)/r**4
         d2th_x3   = -x2*2.d0*r*dr_x3/r**4

         znu = sqrt(ZLAMBDA_MODE_DP)-ENVELOPEK*RBEND
!        complex-valued theta comes from zcurv_st(3) (stretched curvilinear coordinates)
         w = exp(-ZI*znu* zcurv_st(3) )
         dw_th = -ZI*znu*w
         d2w_th = -znu**2*w
         ! current derivatives were computed w.r.t complex stretched coordinate \tilde{\theta}
         ! pass to derivatives w.r.t physical coordinate theta
         if (activePML.ne.0) then
            d2w_th = d2w_th* zdcurv_st(3)**2 + dw_th * zd2curv_st(3)
            dw_th = dw_th * zdcurv_st(3)
         endif

!     ...mfd solution for the polarized component of E
         cn = 1.d0*ZONE
         p = cn * u * v * w
!     ...1st order derivatives
         Gradp(1) = cn * u * dv_x1 * w
         Gradp(2) = cn * (du_r*dr_x2 * v * w  +  u * v * dw_th*dth_x2 )
         Gradp(3) = cn * (du_r*dr_x3 * v * w  +  u * v * dw_th*dth_x3 )
!     ...second order derivatives
         Grad2p(1,1) = cn * (u * d2v_x1 * w )
         Grad2p(1,2) = cn * (du_r*dr_x2 * dv_x1 * w  +  u * dv_x1 * dw_th*dth_x2)
         Grad2p(1,3) = cn * (du_r*dr_x3 * dv_x1 * w  +  u * dv_x1 * dw_th*dth_x3)
         Grad2p(2,1) = Grad2p(1,2)
         Grad2p(2,2) = cn * ( (d2u_r*dr_x2**2+du_r*d2r_x2) * v * w            &
                              + 2.d0*(du_r*dr_x2 * v * dw_th*dth_x2)           &
                              + u * v * (d2w_th*dth_x2**2+dw_th*d2th_x2) )
         Grad2p(2,3) = cn * ( (d2u_r*dr_x2*dr_x3+du_r*d2r_x2x3) * v * w       &
                              + du_r*dr_x2 * v * dw_th*dth_x3                  &
                              + du_r*dr_x3 * v * dw_th*dth_x2                  &
                              + u * v * (dw_th*dth_x2*dth_x3+d2w_th*d2th_x2x3))
         Grad2p(3,1) = Grad2p(1,3)
         Grad2p(3,2) = Grad2p(2,3)
         Grad2p(3,3) = cn * ( (d2u_r*dr_x3**2+du_r*d2r_x3) * v * w            &
                              + 2.d0*(du_r*dr_x3 * v * dw_th*dth_x3)           &
                              + u * v * (d2w_th*dth_x3**2+dw_th*d2th_x3) )
!
!     Propagating modes of straight step-index fiber: 101,111,1110,121,1210,102.   TO BE EVALUATED ONLY AT THE BOTTOM FACE (x3=0)
      case(101) ! LP01 mode
            kappa_core=3.89635821866549d0
            kappa_clad=7.94630920866814d0
            call get_LP01_transversal(Xp,E_AMPL,kappa_core,kappa_clad, p,Gradp)
! 
      case(111) ! LP11a mode
            kappa_core=6.1398251268689d0
            kappa_clad=6.37400853618783d0
            call get_LP11a_transversal(Xp,E_AMPL,kappa_core,kappa_clad, p,Gradp)! 
      case(1110) ! LP11b mode
            kappa_core=6.1398251268689d0
            kappa_clad=6.37400853618783d0
            call get_LP11b_transversal(Xp,E_AMPL,kappa_core,kappa_clad, p,Gradp)
! 
      case(121) ! LP21a mode
            kappa_core=8.08596454606084d0
            kappa_clad=3.59758457409596d0
            call get_LP21a_transversal(Xp,E_AMPL,kappa_core,kappa_clad, p,Gradp)
      case(1210) ! LP21b mode
            kappa_core=8.08596454606084d0
            kappa_clad=3.59758457409596d0
            call get_LP21b_transversal(Xp,E_AMPL,kappa_core,kappa_clad, p,Gradp)
! 
      case(102) ! LP02 mode
            kappa_core=8.47166634119017d0
            kappa_clad=2.56052861953769d0
            call get_LP02_transversal(Xp,E_AMPL,kappa_core,kappa_clad, p,Gradp)
!
!     Straight Step-index slab waveguide modes: 200-202
      case(200) ! first even mode
            kappa_core=2.55570005662253d0
            kappa_clad=8.47312425428247d0
            call get_step_slab_even(Xp,1.d0,kappa_core,kappa_clad, p,Gradp)
!
      case(201) ! first odd mode
            kappa_core=5.06465178791433d0
            kappa_clad=7.25773653938378d0
            call get_step_slab_odd (Xp,1.d0,kappa_core,kappa_clad, p,Gradp)
!
      case(202) ! second even mode
            kappa_core=7.4313660759331d0
            kappa_clad=4.80627045154566d0
            call get_step_slab_even(Xp,1.d0,kappa_core,kappa_clad, p,Gradp)
!
!     Bent Step-index slab waveguide modes with radiation condition: 300-302. TO BE EVALUATED ONLY AT THE BOTTOM FACE (x3=0)
      case(300,301,302) ! can be called with any NEXACT.
            !
            ! EIGENFUNCTIONS OF THE BENT STEP-INDEX SLAB WAVEGUIDE 
            ! WITH NEUMANN BC ON THE INNER FACE, PML ON OUTER FACE
            !
            ! AVAILABLE FOR RBEND = 1300 AND 2600
            !
            ! MAKE SURE THAT PML STRENGTH COEFFICIENT IS SET TO 400!
            !
            if (SLAB_GUIDE.eq.0) then
               write(*,*) 'mfd_solutions: option only available for the SLAB_GUIDE case! STOP'
               stop
            endif

!        ...check if we're within the PML      ARE WE ASSUMING THAT THIS IS EVALUATED AT THE BOTTOM FACE?
            call is_pml(Mdle,Xp,activePML)
            ! if so, modify input coordinate r
            if (activePML.ne.0) then 
               ! get rotation matrix
               call get_local_rotation(Mdle,Xp,activePML,rQ)
               ! get pml stretched coordinates and derivatives
               call get_stretched_coords(Mdle,Xp,ActivePML,Zcurv_st,Zdcurv_st,Zd2curv_st)
               ! get PML Jacobian
               call get_stretch_J(zdcurv_st,rQ,zJ,zJinv,zJdet)
            else
               call cartesian2curvilinear_real(Xp,curv)
               zcurv_st = cmplx(curv,0.d0 , 8)
               zJ = cmplx(IDENTITY,0.d0 , 8)
               zJinv = zJ
               ZJdet = ZONE
            endif

!           ...Separable function u(r)*v(x)*w(theta)

            ! FACTOR u(r)
            ! u(r) is an eigenfunction of Bessel's equation. Its parameters 
            ! (bessel order, linear combination coefficients) were precomputed,
            ! and are chosen from the options in the bessel_evaluation module.
            r = dsqrt(x2**2+x3**2)
            dr_x2 = x2/r
            dr_x3 = x3/r
            d2r_x2 = 1.d0/r-1.d0*x2**2/r**3
            d2r_x3 = 1.d0/r-1.d0*x3**2/r**3
            d2r_x2x3 = -1.d0*x2*x3/r**3
!
            wavenum0 = OMEGA*sqrt(EPSILON*MU)*REFRCORE
            wavenum1 = OMEGA*sqrt(EPSILON*MU)*REFRCLAD

            ! write(*,*) 'mfd_solutions: xp , r  = ',xp
            ! write(*,*) 'mfd_solutions: activePML = ',activePML
            ! write(*,*) 'mfd_solutions: rQ(:,3) = ',rQ(:,3)
            ! write(*,*) 'mfd_solutions: zJdet   = ',zJdet
            ! write(*,*) 'mfd_solutions: Zcurv_st   = ',zcurv_st
            ! write(*,*) 'mfd_solutions: zcurv_st  = ',zcurv_st
            ! write(*,*) 'mfd_solutions: BEFORE bessel_stepindex_preset'
            ! write(*,*) 'mfd_solutions: Mdle, Xp=',Mdle,Xp
            ! call pause   

!           complex-valued r comes from zcurv_st(2) (stretched rotated coordinates)
            call bessel_stepindex_preset(wavenum0, wavenum1, RBEND, RCORE, RCLAD, zcurv_st(2),    u, du_r, d2u_r)

            ! write(*,*) 'mfd_solutions: zcurv_st(2), u = ',zcurv_st(2), u

            ! write(*,*) 'mfd_solutions: AFTER bessel_stepindex_preset'
            ! call pause  

            ! current derivatives were computed w.r.t complex stretched coordinate \tilde{r}
            ! pass to derivatives w.r.t physical coordinate r
            if (activePML.ne.0) then
               d2u_r = d2u_r* zdcurv_st(2)**2 + du_r * zd2curv_st(2)
               du_r = du_r * zdcurv_st(2)
            endif


            ! FACTOR v(x)
            v =     1.d0 !
            dv_x1 = 0.d0 !
            d2v_x1 = 0.d0


            ! FACTOR w(theta)
            dth_x2    = -x3/r**2
            dth_x3    =  x2/r**2
            d2th_x2   =  x3*2.d0*r*dr_x2/r**4
            d2th_x2x3 = -(r**2-x3*2.d0*r*dr_x3)/r**4
            d2th_x3   = -x2*2.d0*r*dr_x3/r**4

            znu = sqrt(ZLAMBDA_MODE)-ENVELOPEK*RBEND
!           complex-valued theta comes from zcurv_st(3) (stretched curvilinear coordinates)
            w = exp(-ZI*znu* zcurv_st(3) )
            dw_th = -ZI*znu*w
            d2w_th = -znu**2*w
            ! current derivatives were computed w.r.t complex stretched coordinate \tilde{\theta}
            ! pass to derivatives w.r.t physical coordinate theta
            if (activePML.ne.0) then
               d2w_th = d2w_th* zdcurv_st(3)**2 + dw_th * zd2curv_st(3)
               dw_th = dw_th * zdcurv_st(3)
            endif

!     ...mfd solution for the polarized component of E
            cn = 1.d0*ZONE
            p = cn * u * v * w
!     ...1st order derivatives
            Gradp(1) = cn * u * dv_x1 * w
            Gradp(2) = cn * (du_r*dr_x2 * v * w  +  u * v * dw_th*dth_x2 )
            Gradp(3) = cn * (du_r*dr_x3 * v * w  +  u * v * dw_th*dth_x3 )
!     ...second order derivatives
            Grad2p(1,1) = cn * (u * d2v_x1 * w )
            Grad2p(1,2) = cn * (du_r*dr_x2 * dv_x1 * w  +  u * dv_x1 * dw_th*dth_x2)
            Grad2p(1,3) = cn * (du_r*dr_x3 * dv_x1 * w  +  u * dv_x1 * dw_th*dth_x3)
            Grad2p(2,1) = Grad2p(1,2)
            Grad2p(2,2) = cn * ( (d2u_r*dr_x2**2+du_r*d2r_x2) * v * w            &
                                + 2.d0*(du_r*dr_x2 * v * dw_th*dth_x2)           &
                                + u * v * (d2w_th*dth_x2**2+dw_th*d2th_x2) )
            Grad2p(2,3) = cn * ( (d2u_r*dr_x2*dr_x3+du_r*d2r_x2x3) * v * w       &
                                + du_r*dr_x2 * v * dw_th*dth_x3                  &
                                + du_r*dr_x3 * v * dw_th*dth_x2                  &
                                + u * v * (dw_th*dth_x2*dth_x3+d2w_th*d2th_x2x3))
            Grad2p(3,1) = Grad2p(1,3)
            Grad2p(3,2) = Grad2p(2,3)
            Grad2p(3,3) = cn * ( (d2u_r*dr_x3**2+du_r*d2r_x3) * v * w            &
                                + 2.d0*(du_r*dr_x3 * v * dw_th*dth_x3)           &
                                + u * v * (d2w_th*dth_x3**2+dw_th*d2th_x3) )

            ! Gradp = matmul(zJ)

   end select
!
!
   end subroutine mfd_solutions
!
!------------------------------------------------------
! subroutines to evaluate Bessel functions
!------------------------------------------------------
function BESSEL_dJ1(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   fval = BESSEL_J0(x) - BESSEL_J1(x)/x
end function

function BESSEL_J2(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   fval = BESSEL_JN(2, x)
end function

function BESSEL_dJ2(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessJY(x,2.d0, a,b,fval,c)
end function

function BESSEL_K0(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,0.d0, a,fval,b,c)
end function

function BESSEL_dK0(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,0.d0, a,b,c,fval)
end function

function BESSEL_K1(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,1.d0, a,fval,b,c)
end function

function BESSEL_dK1(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,1.d0, a,b,c,fval)
end function

function BESSEL_K2(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,2.d0, a,fval,b,c)
end function

function BESSEL_dK2(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,2.d0, a,b,c,fval)
end function


!------------------------------------------------------
! subroutine get_LP01_transversal
!------------------------------------------------------
subroutine get_LP01_transversal(Xp,Ampl,Kappa_core,Kappa_clad, E,dE)
!
   use commonParam
   use control
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: Ampl, Kappa_core, Kappa_clad
   complex(8), intent(out) :: E, dE(3)
!  
   real(8) :: BESSEL_K0, BESSEL_K1
   real(8) :: x1, x2, x3, r, r_x, r_y, ca, cb
!
!------------------------------------------------------
!
   if (NEXACT.ne.0) then
      write(*,*) 'get_LP01_transversal: Error. Transversal LP01 mode to be used if NEXACT=0 only. Stop'
   endif

   E = ZERO; dE = ZERO
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2) - RBEND; x3 = Xp(3)

!..radial coordinate
   r = sqrt(x1*x1+x2*x2)
!  evaluate mode only on bottom face
   if (x3.lt.GEOM_TOL) then
!
      if (abs(r).lt.GEOM_TOL) then
         r_x = 1.d0
         r_y = 1.d0
      else
         r_x = x1/r
         r_y = x2/r
      endif
!
      if (r .le. RCORE) then
         ca = Ampl/BESSEL_J0(Kappa_core*RCORE)
         E = ca*BESSEL_J0(Kappa_core*r)
         cb = -ca*Kappa_core*BESSEL_J1(Kappa_core*r)
      else
         ca = Ampl/BESSEL_K0(Kappa_clad*RCORE)
         E = ca*BESSEL_K0(Kappa_clad*r)
         cb = -ca*Kappa_clad*BESSEL_K1(Kappa_clad*r)
      endif
      dE(1) = cb*r_x
      dE(2) = cb*r_y
   endif
!
end subroutine get_LP01_transversal
!
!------------------------------------------------------
! subroutine get_LP11a_transversal
!------------------------------------------------------
subroutine get_LP11a_transversal(Xp,Ampl,Kappa_core,Kappa_clad, E,dE)
!
   use commonParam 
   use control, only : GEOM_TOL
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: Ampl, Kappa_core, Kappa_clad
   complex(8)  , intent(out) :: E, dE(3)
!
   real(8) :: x1, x2, x3, r, ca, cb, cc
   real(8) :: BESSEL_dJ1, BESSEL_K1, BESSEL_dK1
!
!------------------------------------------------------
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2) - RBEND; x3 = Xp(3)
!
!..shift source away from zero
   if (abs(x1) .lt. GEOM_TOL) then
      x1 = x1+GEOM_TOL
   endif
   if (abs(x2) .lt. GEOM_TOL) then
      x2 = x2+GEOM_TOL
   endif
   r = sqrt(x1*x1+x2*x2)
!
   if (r .le. RCORE) then
      ca = Ampl/BESSEL_J1(Kappa_core*RCORE)
      E = ca*(x1/r)*BESSEL_J1(Kappa_core*r)
      cb = ca*(((x1/r)**(2.d0))*Kappa_core*BESSEL_dJ1(Kappa_core*r)+ &
               ((x2/r)**(2.d0))*BESSEL_J1(Kappa_core*r)/r)
      cc = ca*(x2/r)*(x1/r)*(Kappa_core*BESSEL_dJ1(Kappa_core*r)-BESSEL_J1(Kappa_core*r)/r)
   else
      ca = Ampl/BESSEL_K1(Kappa_clad*RCORE)
      E = ca*(x1/r)*BESSEL_K1(Kappa_clad*r)
      cb = ca*(((x1/r)**(2.d0))*Kappa_clad*BESSEL_dK1(Kappa_clad*r)+ &
               ((x2/r)**(2.d0))*BESSEL_K1(Kappa_clad*r)/r)
      cc = ca*(x2/r)*(x1/r)*(Kappa_clad*BESSEL_dK1(Kappa_clad*r)-BESSEL_K1(Kappa_clad*r)/r)
   endif
!
   dE(1) = cb
   dE(2) = cc
!
   if (RCLAD - r < GEOM_TOL) then
      E = 0.d0
      dE(1:3) = 0.d0
   endif
!
end subroutine get_LP11a_transversal
!
!------------------------------------------------------
! subroutine get_LP11b_transversal (rotated LP11 mode)
!------------------------------------------------------
subroutine get_LP11b_transversal(Xp,Ampl,Kappa_core,Kappa_clad, E,dE)
!
   use commonParam 
   use control, only : GEOM_TOL
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: Ampl, Kappa_core, Kappa_clad
   complex(8) , intent(out) :: E, dE(3)
!
   real(8) :: x1, x2, x3, r, ca, cb, cc
   real(8) :: BESSEL_dJ1, BESSEL_K1, BESSEL_dK1
!
!------------------------------------------------------
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2) - RBEND; x3 = Xp(3)
!
!..shift source away from zero
   if (abs(x1) .lt. GEOM_TOL) then
      x1 = x1+GEOM_TOL
   endif
   if (abs(x2) .lt. GEOM_TOL) then
      x2 = x2+GEOM_TOL
   endif
   r = sqrt(x1*x1+x2*x2)
!
   if (r .le. RCORE) then
      ca = Ampl/BESSEL_J1(Kappa_core*RCORE)
      E = ca*(x2/r)*BESSEL_J1(Kappa_core*r)
      cb = ca*(x2/r)*(x1/r)*(Kappa_core*BESSEL_dJ1(Kappa_core*r)-BESSEL_J1(Kappa_core*r)/r)
      cc = ca*(((x2/r)**(2.d0))*Kappa_core*BESSEL_dJ1(Kappa_core*r) + &
               ((x1/r)**(2.d0))*BESSEL_J1(Kappa_core*r)/r)
   else
      ca = Ampl/BESSEL_K1(Kappa_clad*RCORE)
      E = ca*(x2/r)*BESSEL_K1(Kappa_clad*r)
      cb = ca*(x2/r)*(x1/r)*(Kappa_clad*BESSEL_dK1(Kappa_clad*r)-BESSEL_K1(Kappa_clad*r)/r)
      cc = ca*(((x2/r)**(2.d0))*Kappa_clad*BESSEL_dK1(Kappa_clad*r) + &
               ((x1/r)**(2.d0))*BESSEL_K1(Kappa_clad*r)/r)
   endif
!
   dE(1) = cb
   dE(2) = cc
!
   if (RCLAD - r < GEOM_TOL) then
      E = 0.d0
      dE(1:3) = 0.d0
   endif
!
end subroutine get_LP11b_transversal
!
!------------------------------------------------------
! subroutine get_LP02_transversal
!------------------------------------------------------
subroutine get_LP02_transversal(Xp,Ampl,Kappa_core,Kappa_clad, E,dE)
!
   use commonParam 
   use control, only : GEOM_TOL
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: Ampl, Kappa_core, Kappa_clad
   complex(8)  , intent(out) :: E, dE(3)
!
   real(8) :: x1, x2, x3, r, r_x, r_y, ca, cb
   real(8) :: BESSEL_K0, BESSEL_K1
!
!------------------------------------------------------
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2) - RBEND; x3 = Xp(3)
!
!..radial coordinate
   r = sqrt(x1*x1+x2*x2)
!
   if (abs(r).lt.GEOM_TOL) then
      r_x = 1.d0
      r_y = 1.d0
   else
      r_x = x1/r
      r_y = x2/r
   endif
!
   if (r .le. RCORE) then
      ca = Ampl/BESSEL_J0(Kappa_core*RCORE)
      E = ca*BESSEL_J0(Kappa_core*r)
      cb = -ca*Kappa_core*BESSEL_J1(Kappa_core*r)
   else
      ca = Ampl/BESSEL_K0(Kappa_clad*RCORE)
      E = ca*BESSEL_K0(Kappa_clad*r)
      cb = -ca*Kappa_clad*BESSEL_K1(Kappa_clad*r)
   endif
   dE(1) = cb*r_x
   dE(2) = cb*r_y
!
   if (RCLAD - r < GEOM_TOL) then
      E = 0.d0
      dE(1:3) = 0.d0
   endif
!
end subroutine get_LP02_transversal
!
!------------------------------------------------------
! subroutine get_LP21a_transversal
!------------------------------------------------------
subroutine get_LP21a_transversal(Xp,Ampl,Kappa_core,Kappa_clad, E,dE)
!
   use commonParam 
   use control, only : GEOM_TOL
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: Ampl, Kappa_core, Kappa_clad
   complex(8)  , intent(out) :: E, dE(3)
!
   real(8) :: x1, x2, x3, r, ca, cb, cc
   real(8) :: BESSEL_J2, BESSEL_dJ2, BESSEL_K2, BESSEL_dK2
!
   real(8) :: cos_t,cos_2t
   real(8) :: sin_t,sin_2t
!
!------------------------------------------------------
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2) - RBEND; x3 = Xp(3)
!
!..shift source away from zero
   if (abs(x1) .lt. GEOM_TOL) then
      x1 = x1+GEOM_TOL
   endif
   if (abs(x2) .lt. GEOM_TOL) then
      x2 = x2+GEOM_TOL
   endif
   r = sqrt(x1*x1+x2*x2)
!
   cos_t  = x1/r
   sin_t  = x2/r
   cos_2t = cos_t**(2.d0) - sin_t**(2.d0)
   sin_2t = 2 * sin_t * cos_t
!
   if (r .le. RCORE) then
      ca = Ampl/BESSEL_J2(Kappa_core*RCORE)
      E = ca*cos_2t*BESSEL_J2(Kappa_core*r)
      cb = ca*(cos_t*cos_2t*Kappa_core*BESSEL_dJ2(Kappa_core*r) + &
               2.d0*sin_t*sin_2t*BESSEL_J2(Kappa_core*r)/r)
      cc = ca*(sin_t*cos_2t*Kappa_core*BESSEL_dJ2(Kappa_core*r) - &
               2.d0*sin_2t*cos_t*BESSEL_J2(Kappa_core*r)/r)
   else
      ca = Ampl/BESSEL_K2(Kappa_clad*RCORE)
      E = ca*cos_2t*BESSEL_K2(Kappa_clad*r)
      cb = ca*(cos_t*cos_2t*Kappa_clad*BESSEL_dK2(Kappa_clad*r) + &
               2.d0*sin_t*sin_2t*BESSEL_K2(Kappa_clad*r)/r)
      cc = ca*(sin_t*cos_2t*Kappa_clad*BESSEL_dK2(Kappa_clad*r) - &
      2.d0*sin_2t*cos_t*BESSEL_K2(Kappa_clad*r)/r)
   endif
!
   dE(1) = cb
   dE(2) = cc
!
   if (RCLAD - r < GEOM_TOL) then
      E = 0.d0
      dE(1:3) = 0.d0
   endif
!
end subroutine get_LP21a_transversal

!------------------------------------------------------
! subroutine get_LP21b_transversal (rotated LP21 mode)
!------------------------------------------------------
subroutine get_LP21b_transversal(Xp,Ampl,Kappa_core,Kappa_clad, E,dE)
!
   use commonParam 
   use control, only : GEOM_TOL
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: Ampl, Kappa_core, Kappa_clad
   complex(8) , intent(out) :: E, dE(3)
!
   real(8) :: x1, x2, x3, r, ca, cb, cc
   real(8) :: BESSEL_J2, BESSEL_dJ2, BESSEL_K2, BESSEL_dK2
!
   real(8) :: cos_t,cos_2t
   real(8) :: sin_t,sin_2t
!
!------------------------------------------------------
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2) - RBEND; x3 = Xp(3)
!
!..shift source away from zero
   if (abs(x1) .lt. GEOM_TOL) then
      x1 = x1+GEOM_TOL
   endif
   if (abs(x2) .lt. GEOM_TOL) then
      x2 = x2+GEOM_TOL
   endif
   r = sqrt(x1*x1+x2*x2)
!
   cos_t  = x1/r
   sin_t  = x2/r
   cos_2t = cos_t**(2.d0) - sin_t**(2.d0)
   sin_2t = 2 * sin_t * cos_t
!
   if (r .le. RCORE) then
      ca = Ampl/BESSEL_J2(Kappa_core*RCORE)
      E = ca*sin_2t*BESSEL_J2(Kappa_core*r)
      cb = ca*(cos_t*sin_2t*Kappa_core*BESSEL_dJ2(Kappa_core*r) - &
               2.d0*sin_t*cos_2t*BESSEL_J2(Kappa_core*r)/r)
      cc = ca*(sin_t*sin_2t*Kappa_core*BESSEL_dJ2(Kappa_core*r) + &
               2.d0*cos_2t*cos_t*BESSEL_J2(Kappa_core*r)/r)
   else
      ca = Ampl/BESSEL_K2(Kappa_clad*RCORE)
      E = ca*sin_2t*BESSEL_K2(Kappa_clad*r)
      cb = ca*(cos_t*cos_2t*Kappa_clad*BESSEL_dK2(Kappa_clad*r) - &
               2.d0*sin_t*cos_2t*BESSEL_K2(Kappa_clad*r)/r)
      cc = ca*(sin_t*sin_2t*Kappa_clad*BESSEL_dK2(Kappa_clad*r) + &
               2.d0*cos_2t*cos_t*BESSEL_K2(Kappa_clad*r)/r)
   endif
!
   dE(1) = cb
   dE(2) = cc
!
   if (RCLAD - r < GEOM_TOL) then
      E = 0.d0
      dE(1:3) = 0.d0
   endif
!
end subroutine get_LP21b_transversal


!------------------------------------------------------
!
!------------------------------------------------------
!
subroutine get_step_slab_even(Xp,Ampl,Kappa_core,Kappa_clad, E,dE)
!
   use commonParam, only: ZERO,RBEND,RCORE
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: Ampl, Kappa_core, Kappa_clad
   complex(8)  , intent(out) :: E, dE(3)
!
   real(8) :: X2
!
!..shift the slab center in y axis
   x2 = Xp(2) - RBEND
!..initialize gradient
   dE = ZERO

   if (abs(x2).le.RCORE) then
      E    =  Ampl*cos(Kappa_core*x2)
      dE(2)= -Ampl*Kappa_core*sin(Kappa_core*x2)
   endif
   if (x2.lt.-RCORE) then
      E    =  Ampl*cos(Kappa_core*RCORE)*exp(Kappa_clad*(x2+RCORE))
      dE(2)=  Ampl*cos(Kappa_core*RCORE)*Kappa_clad*exp( Kappa_clad*(x2+RCORE))
   endif
   if (x2.gt. RCORE) then
      E    =  Ampl*cos(Kappa_core*RCORE)*exp(-Kappa_clad*(x2-RCORE))
      dE(2)= -Ampl*cos(Kappa_core*RCORE)*Kappa_clad*exp(-Kappa_clad*(x2-RCORE))
   endif
end subroutine



!------------------------------------------------------
!
!------------------------------------------------------
!
subroutine get_step_slab_odd(Xp,Ampl,Kappa_core,Kappa_clad, E,dE)
!
   use commonParam, only: ZERO,RBEND,RCORE
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: Ampl, Kappa_core, Kappa_clad
   complex(8)  , intent(out) :: E, dE(3)
!
   real(8) :: x2
!
!------------------------------------------------------
!
!..shift the slab center in y axis
   x2 = Xp(2) - RBEND
!..initialize gradient
   dE = ZERO

   if (abs(x2).le.RCORE) then
      E    = Ampl*sin(Kappa_core*x2)
      dE(2)= Ampl*Kappa_core*cos(Kappa_core*x2)
   endif
   if (x2.lt.-RCORE) then
      E    = -Ampl*sin(Kappa_core*RCORE)*exp(Kappa_clad*(x2+RCORE))
      dE(2)= -Ampl*sin(Kappa_core*RCORE)*Kappa_clad*exp( Kappa_clad*(x2+RCORE))
   endif
   if (x2.gt. RCORE) then
      E    = Ampl*sin(Kappa_core*RCORE)*exp(-Kappa_clad*(x2-RCORE))
      dE(2)= -Ampl*sin(Kappa_core*RCORE)*Kappa_clad*exp(-Kappa_clad*(x2-RCORE))
   endif
end subroutine