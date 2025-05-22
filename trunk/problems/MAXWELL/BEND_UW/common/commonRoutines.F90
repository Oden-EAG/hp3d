!
#include "typedefs.h"
!
!------------------------------------------------------------------------------
!> @brief      Propagate flag from father to son nodes; used to correctly
!!             inherit impedance BCs
!!
!> @param[in]  Icomp  - Physics component on which to propagate BC flag
!> @param[in]  Nflag  - flag to propagate
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine propagate_flag(Icomp,Nflag)
!
      use data_structure3D
      use commonParam, only: IBCFLAG
!
      implicit none
!
      integer, intent(in) :: Icomp,Nflag
!
      integer :: ntype
      integer :: iel,mdle,ifc,nrfn,i,j,nod
!
!  ...element nodes and orientations, face nodes
      integer :: nodesl(27),norientl(27),nface_nodes(9)
!
!  ...element face BC flags, decoded BC flag for a node
      integer :: ibc(6,NRINDEX),nodflag(NRINDEX_HEV)
!
!----------------------------------------------------------------------------
!
      if (IBCFLAG .ne. 3) then
         write(*,*) 'propagate_flag called for IBCFLAG.ne.3, returning...'
         return
      endif
!
      if ((Icomp.lt.1) .or. (Icomp.gt.NRINDEX_HEV)) then
         write(*,*) 'propagate_flag: invalid Icomp = ', Icomp
         return
      endif
!
!..loop through active elements
!$OMP PARALLEL                                     &
!$OMP PRIVATE(ntype,mdle,ifc,nrfn,i,j,nod,nodesl,  &
!$OMP         norientl,nface_nodes,ibc,nodflag)
!
!$OMP DO
      do iel=1,NRELES
         mdle = ELEM_ORDER(iel)
         ntype = NODES(mdle)%ntype
!
!     ...determine element nodes
         call elem_nodes(mdle, nodesl,norientl)
!
!     ...get the element boundary conditions flags
         call find_bc(mdle, ibc)
!
!     ...loop through element faces
         do ifc=1,nface(ntype)
!
!        ...if face has a Dirichlet BC flag on this component,
!           then neither propagate Nflag from this face to its edges/vertices,
!           nor prohibit another face from passing Nflag to the edges/vertices.
            if (ibc(ifc,Icomp).eq.1) cycle
!
!        ...determine face node numbers
            call face_nodes(ntype,ifc, nface_nodes,nrfn)
!
!        ...loop through the face nodes
!$OMP CRITICAL
          do i=1,nrfn !-1
             j = nface_nodes(i)
             nod = nodesl(j)
!
!         ...if node belongs to a face that has impedance BC (Nflag),
!            then propagate the flag unless prohibited by another adjacent face
             if (ibc(ifc,Icomp).eq.Nflag) then
                if (NODES(nod)%visit.ne.-Nflag) then
                   NODES(nod)%visit = Nflag
                endif
!         ...prohibit the flag to be passed to the node
!            (if node belongs to a face that has no impedance or Dirichlet BC)
             else
                NODES(nod)%visit = -Nflag
             endif
          enddo
!$OMP END CRITICAL
         enddo
      enddo
!$OMP END DO
!
!  ...change -Nflag to zero
!$OMP DO
      do nod=1,NRNODS
         if (NODES(nod)%visit.eq.0) cycle
         call decod(NODES(nod)%bcond,2,NRINDEX_HEV, nodflag)
         if (NODES(nod)%visit.eq.-Nflag) then
            nodflag(Icomp) = 0
         elseif (NODES(nod)%visit.eq.Nflag) then
            nodflag(Icomp) = 1
         endif
         call encod(nodflag,2,NRINDEX_HEV, NODES(nod)%bcond)
!
      enddo
!$OMP END DO
!$OMP END PARALLEL
!
      call reset_visit
!
   end subroutine propagate_flag




!------------------------------------------------------------------------------
!> @brief      Determines local permittivity tensor epsilon
!!
!> @param[in]  Mdle     - element middle node number
!> @param[in]  x        - physical point to evaluate permittivity
!!
!> @param[in]  eps      - permittivity tensor at point
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine get_local_epsilon(Mdle,Xp,Zeps)
!
      use commonParam, only: EPSILON, ZERO, ZONE
!
      implicit none
!
      integer,    intent(in)  :: Mdle
      real(8),    intent(in)  :: Xp(3)
      complex(8), intent(out) :: Zeps(3,3)
!
      integer :: i
      complex(8) :: zrefr
!
!------------------------------------------------------------------------------
!
!  ...get refractive index of element's subdomain
      call get_refrac(Mdle,zrefr)
!  ...set permittivity to identity for now.
      Zeps = ZERO
      do i=1,3
         Zeps(i,i) = zrefr**2 * EPSILON
      enddo
!
   end subroutine get_local_epsilon






   subroutine get_local_mu(Mdle,Xp,Zmu)
!
      use commonParam, only: MU, ZERO
!
      implicit none
!
      integer,    intent(in)  :: Mdle
      real(8),    intent(in)  :: Xp(3)
      complex(8), intent(out) :: Zmu(3,3)
!
      integer :: i
!
!------------------------------------------------------------------------------
!
!  ...set permeability to mu_0 for now
      Zmu = ZERO
      do i=1,3
         Zmu(i,i) = cmplx(MU,0.d0,8)
      enddo
!
   end subroutine get_local_mu








   subroutine get_refrac(Mdle,Zrefr)
      use commonParam, only: SLAB_GUIDE, ZONE, OMEGA, EPSILON, MU, REFRCORE, REFRCLAD, REFRCOAT,ATTNCOAT
      implicit none
!
      integer,    intent(in) :: Mdle
      complex(8), intent(out):: Zrefr
!
      integer :: ndom
!
      call find_domain(Mdle, ndom)
      if (SLAB_GUIDE.eq.1) then
         select case(ndom)
         case(1)
            Zrefr = cmplx(REFRCORE,0.0,8)
         case(2,3)
            Zrefr = cmplx(REFRCLAD,0.0,8)
         case(4,5)
            Zrefr = cmplx(REFRCOAT,-ATTNCOAT*sqrt(MU*EPSILON)/(2.d0*OMEGA),8)
         case default
            Zrefr = ZONE
         end select
      else
         select case(ndom)
         case(1,2)
            Zrefr = cmplx(REFRCORE,0.0,8)
         case(3)
            Zrefr = cmplx(REFRCLAD,0.0,8)
         case(4)
            Zrefr = cmplx(REFRCOAT,-ATTNCOAT*sqrt(MU*EPSILON)/(2.d0*OMEGA),8)
         case default
            Zrefr = ZONE
         end select
      endif

   end subroutine








   subroutine get_local_rotation(Mdle,Xp,ActivePML,RQ)
!
      use commonParam, only: TOROIDAL_PML
      implicit none
!
      integer, intent(in)  :: Mdle,ActivePML
      real(8), intent(in)  :: Xp(3)
      real(8), intent(out) :: RQ(3,3)
!
      real(8) :: th,Qtor(3,3)
      integer :: flags(8)
      !
      if (activePML.eq.0) then
         write(*,*) "get_local_rotation: no PML flag! Stop"
         stop
      endif
      !
      !  decode the ActivPML variable
      call decod(ActivePML,2,8,flags)
      !  get rotation to cylindrical basis RQ <-- Qcyl
      call get_cylindrical_rotation(Mdle,Xp,RQ)
      !  check if a PML for toroidal geometry is set, and if Xp lies there
      if (TOROIDAL_PML.eq.1 .and.(flags(7).eq.1 .or. flags(8).eq.1)) then 
         !  get rotation from cylindrical to toroidal basis
         call get_toroidal_rotation(Mdle,Xp,Qtor)
         !  update the full rotation RQ <-- Qtor @ Qcyl
         RQ = matmul(Qtor,RQ)
      endif
!
   end subroutine





   subroutine get_cylindrical_rotation(Mdle,Xp,RQ)
!
      implicit none
!
      integer, intent(in)  :: Mdle
      real(8), intent(in)  :: Xp(3)
      real(8), intent(out) :: RQ(3,3)

      real(8) :: th

      th = atan2(Xp(3),Xp(2))

      RQ(:,:) = 0.d0 

      RQ(1,1) = 1.d0
      RQ(2,2) = cos(th)
      RQ(3,2) =-sin(th)
      RQ(2,3) =-RQ(3,2)
      RQ(3,3) = RQ(2,2)

   end subroutine


   subroutine get_toroidal_rotation(Mdle,Xp,RQ)
!
      use commonParam, only: RBEND
      implicit none
!
      integer, intent(in)  :: Mdle
      real(8), intent(in)  :: Xp(3)
      real(8), intent(out) :: RQ(3,3)

      real(8) :: x,r,phi
      
      r = sqrt(Xp(2)**2+Xp(3)**2)
      x = Xp(1)   
      phi = atan2(x,r-RBEND)

      RQ(:,:) = 0.d0 

      RQ(1,1) = sin(phi)
      RQ(2,1) = cos(phi)
      RQ(1,2) = RQ(2,1)
      RQ(2,2) =-RQ(1,1)
      RQ(3,3) = 1.d0

   end subroutine




   ! subroutine get_pml_theta(Mdle,Xp,Zthst,Zdthst,Zd2thst)

   !    use commonParam, only: THETAEND,PMLPROP,OMEGA
   !    implicit none
   !    integer, intent(in) :: Mdle
   !    real(8), intent(in) :: Xp(3)
   !    complex(8),intent(out) :: Zthst,Zdthst,Zd2thst

   !    real(8) :: c,pn,f,df,d2f,th,th_pml,th_trn,th_dif
   !    integer :: n

   !    th = atan2(Xp(3),Xp(2))

   !    th_pml = THETAEND*PMLPROP
   !    th_trn = THETAEND - th_pml
   !    th_dif = Th - th_trn

   !    if (Th.gt.THETAEND.or.th_dif.lt.0.d0) then
   !       write(*,*) ' get_pml_theta: theta>THETAEND or theta<th_trn'
   !       write(*,*) ' theta,th_trn,THETAEND = ',Th,th_trn,THETAEND
   !       write(*,*) ' stop'
   !       stop
   !    endif

   !    n = 3       !!! NEEDS TO BE AT LEAST 2 !!!
   !    c = 10.d0 /(OMEGA* th_pml**n)
   !    f   = c*th_dif**n
   !    df  = c*th_dif**(n-1) * n
   !    d2f = c*th_dif**(n-2) * (n*(n-1))
   !    if((f.le.0.d0).or.(df.le.0.d0).or.(d2f.le.0)) then
   !       write(*,*) ' get_pml_theta: f, df,d2f are negative. stop.'
   !       stop
   !    endif
   !    zthst   = cmplx( Th , -f  , 8)
   !    zdthst  = cmplx(1.d0, -df , 8)
   !    zd2thst = cmplx(0.d0, -d2f, 8)

   ! end subroutine



!===========================================================================
!  SUROUTINE TO COMPUTE COMPLEX-STRETCHED COORDINATE FOR UNIAXIAL PML
!  INPUTS
!  S       Real, original physical coordinate
!  Stra    Real, value of S at which the transition to PML begins
!  Sbnd    Real, minimum/maximum value of S at the boundary where PML lies
!  OUTPUTS
!  Zsst    Complex, stretched
   subroutine get_pml_stretch(S,Stra,Sbnd,Zsst,Zdsst,Zd2sst)

      use commonParam, only: OMEGA
      implicit none
      real(8), intent(in) :: S,Sbnd,Stra
      complex(8),intent(out) :: Zsst,Zdsst,Zd2sst

      real(8) :: c,pn,f,df,d2f,spml,sdif
      integer :: n

      ! get value relative to Stra
      sdif = S - Stra
      ! get PML length (sign indicates side of domain, don't change!)
      spml = Sbnd - Stra

      ! checking that sdif and spml have the same sign
      If (sdif*spml<0.d0) then 
         write(*,*) "get_pml_stretch: value S does not lie between Stra and Sbnd!"
         write(*,*) "S,Stra,Sbnd=",S,Stra,Sbnd
         stop
      endif

      n = 2       !!! NEEDS TO BE AT LEAST 2 !!!
      c = 30.d0 /(OMEGA* spml**n)
      f   = c*sdif**n
      df  = c*sdif**(n-1) * n
      d2f = c*sdif**(n-2) * (n*(n-1))
      if((f.lt.0.d0).or.(df*spml.lt.0.d0)) then
         write(*,*) ' get_pml_stretch: f,df have the wrong sign. stop.'
         write(*,*) ' get_pml_stretch: f,df=',f,df
         stop
      endif
      zsst   = cmplx(  S , -f  , 8)
      zdsst  = cmplx(1.d0, -df , 8)
      zd2sst = cmplx(0.d0, -d2f, 8)

   end subroutine




   subroutine get_stretch_J(Zdsst,RQ,ZJ,ZJinv,ZJdet)

      implicit none 
      complex(8), intent(in) :: Zdsst(3)
      real(8),    intent(in) :: RQ(3,3)
      complex(8), intent(out):: ZJ(3,3),ZJinv(3,3),ZJdet

      complex(8) :: ZQaux(3,3)

      ! precompute tildeJ * Q
      ZQaux(1,:) = RQ(1,:)*Zdsst(1)
      ZQaux(2,:) = RQ(2,:)*Zdsst(2)
      ZQaux(3,:) = RQ(3,:)*Zdsst(3)
      ! multiply Q^T * J * Q
      ZJ = matmul(transpose(RQ),ZQaux)
      !
      ! precompute tildeJinv * Q
      ZQaux(1,:) = RQ(1,:)/Zdsst(1)
      ZQaux(2,:) = RQ(2,:)/Zdsst(2)
      ZQaux(3,:) = RQ(3,:)/Zdsst(3)
      ! multiply Q^T * Jinv * Q
      ZJinv = matmul(transpose(RQ),ZQaux)
      ! Compute determinant of J = det tildeJ
      ZJdet = Zdsst(1)*Zdsst(2)*Zdsst(3)

   end subroutine





   subroutine is_pml(Mdle,Xp,ActivePML)
      use commonParam
      implicit none 
      integer, intent(in) :: Mdle 
      real(8), intent(in) :: Xp(3)
      integer, intent(out):: ActivePML

      real(8) :: th,x,r,rho

      ActivePML = 0
      !  First, PML in the longitudinal direction
      th = atan2(Xp(3),Xp(2))
      if (th .gt.(1.d0-PMLTHUP )*THUP +PMLTHUP *THLO ) ActivePML = ActivePML + 2**7
      if (th .lt.(1.d0-PMLTHLO )*THLO +PMLTHLO *THUP ) ActivePML = ActivePML + 2**6
      !  Now, PML for the transversal geometry
      if (TOROIDAL_PML.eq.0) then
         r = sqrt(Xp(2)**2+Xp(3)**2)
         x = Xp(1)
         if (r  .gt.(1.d0-PMLRUP  )*RUP  +PMLRUP  *RLO  ) ActivePML = ActivePML + 2**5
         if (r  .lt.(1.d0-PMLRLO  )*RLO  +PMLRLO  *RUP  ) ActivePML = ActivePML + 2**4
         if (x  .gt.(1.d0-PMLXUP  )*XUP  +PMLXUP  *XLO  ) ActivePML = ActivePML + 2**3
         if (x  .lt.(1.d0-PMLXLO  )*XLO  +PMLXLO  *XUP  ) ActivePML = ActivePML + 2**2
      else
         rho = sqrt((r-RBEND)**2+x**2)
         if (rho.gt.(1.d0-PMLRHOUP)*RHOUP+PMLRHOUP*RHOLO) ActivePML = ActivePML + 2**1
         if (rho.lt.(1.d0-PMLRHOLO)*RHOLO+PMLRHOLO*RHOUP) ActivePML = ActivePML + 2**0
      endif
!
   end subroutine


   subroutine get_stretched_coords(Mdle,Xp,ActivePML,Zxpst,Zdxpst,Zd2xpst)

      use commonParam
      implicit none
      integer,   intent(in) :: Mdle, ActivePML
      real(8),   intent(in) :: Xp(3)
      complex(8),intent(out):: Zxpst(3),Zdxpst(3),Zd2xpst(3)

      real(8) :: s, stra, sbnd
      integer :: ic,flags(8)

      Zxpst = cmplx(Xp,0.d0,8)
      Zdxpst = ZONE
      Zd2xpst = ZERO
      !  decode the ActivPML variable
      call decod(ActivePML,2,8,flags)
      
      ! write(*,*) "get_stretched_coords:    ActivePML=",ActivePML
      ! write(*,*) "get_stretched_coords:        flags=",flags
      ! write(*,*) "get_stretched_coords:           Xp=",Xp
      ! write(*,*) "get_stretched_coords:          XUP=",XUP
      ! write(*,*) "get_stretched_coords:          XLO=",XLO
      ! call pause

      !  First, PML in the longitudinal direction
      !  coordinate TH, upper bound
      if (flags(1).eq.1 .and. flags(2).eq.0) then
         ic = 3
         s = atan2(Xp(3),Xp(2))
         stra = (1.d0-PMLTHUP )*THUP +PMLTHUP *THLO
         sbnd = THUP
         call get_pml_stretch(s,stra,sbnd,Zxpst(ic),Zdxpst(ic),Zd2xpst(ic))
      endif
      ! coordinate TH, lower bound
      if (flags(1).eq.0 .and. flags(2).eq.1) then
         ic = 3
         s = atan2(Xp(3),Xp(2))
         stra = (1.d0-PMLTHLO )*THUP +PMLTHUP *THUP
         sbnd = THUP
         call get_pml_stretch(s,stra,sbnd,Zxpst(ic),Zdxpst(ic),Zd2xpst(ic))
      endif
      !  Now, PML for the transversal geometry
      if (TOROIDAL_PML.eq.0) then
         ! coordinate R, upper bound
         if (flags(3).eq.1 .and. flags(4).eq.0) then
            ic = 2
            s = sqrt(Xp(2)**2+Xp(3)**2)
            stra = (1.d0-PMLRUP )*RUP +PMLRUP *RLO
            sbnd = RUP
            call get_pml_stretch(s,stra,sbnd,Zxpst(ic),Zdxpst(ic),Zd2xpst(ic))
         endif
         ! coordinate R, lower bound
         if (flags(3).eq.0 .and. flags(4).eq.1) then
            ic = 2
            s = sqrt(Xp(2)**2+Xp(3)**2)
            stra = (1.d0-PMLRLO )*RLO +PMLRLO *RUP
            sbnd = RLO
            call get_pml_stretch(s,stra,sbnd,Zxpst(ic),Zdxpst(ic),Zd2xpst(ic))
         endif
         ! coordinate X, upper bound
         if (flags(5).eq.1 .and. flags(6).eq.0) then
            ic = 1
            s = Xp(1)
            stra = (1.d0-PMLXUP )*XUP +PMLXUP *XLO
            sbnd = XUP
            call get_pml_stretch(s,stra,sbnd,Zxpst(ic),Zdxpst(ic),Zd2xpst(ic))
         endif
         ! coordinate X, lower bound
         if (flags(5).eq.0 .and. flags(6).eq.1) then
            ic = 1
            s = Xp(1)
            stra = (1.d0-PMLXLO )*XLO +PMLXLO *XUP
            sbnd = XLO
            call get_pml_stretch(s,stra,sbnd,Zxpst(ic),Zdxpst(ic),Zd2xpst(ic))
         endif
      else
         ic = 2
         s = sqrt(( sqrt(Xp(2)**2+Xp(3)**2) -RBEND)**2 + Xp(1)**2)
         ! coordinate RHO, upper bound
         if (flags(7).eq.1 .and. flags(8).eq.0) then
            stra = (1.d0-PMLRHOUP )*XUP +PMLRHOUP *RHOLO
            sbnd = RHOUP
            call get_pml_stretch(s,stra,sbnd,Zxpst(ic),Zdxpst(ic),Zd2xpst(ic))
         endif
         ! coordinate RHO, lower bound
         if (flags(7).eq.0 .and. flags(8).eq.1) then
            stra = (1.d0-PMLRHOLO )*XLO +PMLRHOLO *RHOUP
            sbnd = RHOLO
            call get_pml_stretch(s,stra,sbnd,Zxpst(ic),Zdxpst(ic),Zd2xpst(ic))
         endif

      endif

   end subroutine



!------------------------------------------------------------------------------
!> @brief      Returns transformation matrix K arisen due to envelope ansatz 
!!             exp(-ikR\theta)E
!!
!> @param[in]  Mdle     - element middle node number
!> @param[in]  Xp       - physical point to evaluate result
!!
!> @param[out] RK       - transformation matrix K at point X
!!
!> @date       Oct 2024
!------------------------------------------------------------------------------
   subroutine get_matrixK(mdle,Xp, RK)
!
      use data_structure3D
      use commonParam, only: RBEND, ENVELOPEK
      use parameters, only: ZERO, ZONE
!
      implicit none
!
      integer, intent(in)  :: mdle
      real(8), intent(in)  :: Xp(3)
      real(8), intent(out) :: RK(3,3)
!
      integer :: i
      real(8) :: y,z,rr
!
!------------------------------------------------------------------------------
!!
!  ...initialize matrix and copy values of y and z
      RK = 0.d0
      y = Xp(2); z = Xp(3);

!     ENABLE THIS <<IF, ELSE, ENDIF>> CONTROL STRUCTURE IF PARTLY BENT BEHAVIOR IS WANTED
      if (z.gt.0.d0) then
!     ...Recall matrix K = k * Rbend / (y^2 + z^2) * ( 0   -y   -z )
!                                                    ( y    0    0 )
!                                                    ( z    0    0 )  
!    
!     ...first compute coefficient   k*Rbend / (y^2 + z^2)
         rr = ENVELOPEK * RBEND / (y**2+z**2)
!     ...then fill the non-zero entries
         RK(1,2) = -y*rr
         RK(1,3) = -z*rr
         RK(2,1) =  y*rr
         RK(3,1) =  z*rr
      else
!     ...in this case, K = ( 0 -k  0 )
!                          ( k  0  0 )
!                          ( 0  0  0 )
!     ...then fill the non-zero entries
         RK(1,2) = -ENVELOPEK
         RK(2,1) =  ENVELOPEK
!
      endif
!
   end subroutine get_matrixK
!
!
!------------------------------------------------------------------------------
!> @brief      Returns transformation matrix K.K^T arisen due to envelope ansatz 
!!             exp(-ikR\theta)E
!!
!> @param[in]  Mdle     - element middle node number
!> @param[in]  Xp       - physical point to evaluate result
!!
!> @param[out] RKKT     - transformation matrix K.K^T at point X
!!
!> @date       Oct 2024
!------------------------------------------------------------------------------
   subroutine get_matrixKKT(Mdle,Xp, RKKT)
!
      use data_structure3D
      use commonParam, only: RBEND,ENVELOPEK
      use parameters, only: ZERO, ZONE
!
      implicit none
!
      integer, intent(in)  :: Mdle
      real(8), intent(in)  :: Xp(3)
      real(8), intent(out) :: RKKT(3,3)
!
      integer :: i
      real(8) :: y,z,rr
!
!------------------------------------------------------------------------------
!!
!  ...initialize matrix and copy values of y and z
      RKKT = 0.d0
      y = Xp(2); z = Xp(3);
!  ...Recall matrix K.K^T = k^2 * Rbend^2 / (y^2 + z^2)^2 * ( y^2 + z^2  0     0  )
!                                                           (     0     y^2    0  )
!                                                           (     0      0    z^2 )
! 
!  ...first compute coefficient    k^2*Rbend^2 / (y^2 + z^2)^2
      rr = ENVELOPEK**2 * RBEND**2 / (y**2+z**2)**2
!  ...then fill the non-zero entries
      RKKT(1,1) =  (y**2+z**2) * rr
      RKKT(2,2) =  y**2 * rr
      RKKT(3,3) =  z**2 * rr
!
   end subroutine get_matrixKKT
!
!
!------------------------------------------------------------------------------
!> @brief      Returns complex-valued vector K.E, where K is the transformation  
!!             matrix arisen due to the envelope ansatz   exp(-ikR\theta)*E
!!
!> @param[in]  Mdle     - element middle node number
!> @param[in]  Xp       - physical point to evaluate result
!> @param[in]  ZE       - C^3 vector, input E
!!
!> @param[out] ZKE      - resulting C^3 vector: K.E
!!
!> @date       Oct 2024
!------------------------------------------------------------------------------
   subroutine apply_matrixK(Mdle,Xp,ZE,ZKE)
!
      use data_structure3D
      use commonParam, only: RBEND, ENVELOPEK
      use parameters, only: ZERO, ZONE
!
      implicit none
!
      integer,   intent(in)  :: Mdle
      real(8),   intent(in)  :: Xp(3)
      complex(8),intent(in)  :: ZE(3)
      complex(8),intent(out) :: ZKE(3)
!
      integer :: i
      real(8) :: y,z,rr
!
!
#if HP3D_DEBUG
!..Set iprint = 0/1 (Non-/VERBOSE)
      integer :: iprint
      iprint = 0
      ! if (Mdle.eq.937) iprint = 1
#endif
!------------------------------------------------------------------------------
!!
!  ...initialize matrix and copy values of y and z
      ZKE = ZERO
      y = Xp(2); z = Xp(3);

!     ENABLE THIS <<IF, ELSE, ENDIF>> CONTROL STRUCTURE IF PARTLY BENT BEHAVIOR IS WANTED
      if (z.gt.0.d0) then

!     ...Recall matrix K = k * Rbend / (y^2 + z^2) * ( 0   -y   -z )
!                                                    ( y    0    0 )
!                                                    ( z    0    0 )  
! 
!     ...first compute coefficient   k*Rbend / (y^2 + z^2)
         rr = ENVELOPEK * RBEND / (y**2+z**2)
!     ...then fill the entries of K*E
         ZKE(1) = -y*rr*ZE(2) -z*rr*ZE(3)
         ZKE(2) =  y*rr*ZE(1)
         ZKE(3) =  z*rr*ZE(1)

      else
!     ...in this case, K = ( 0 -k  0 )
!                          ( k  0  0 )
!                          ( 0  0  0 )
!     ...so fill the entries of K*E with
         ZKE(1) = -ENVELOPEK*ZE(2)
         ZKE(2) =  ENVELOPEK*ZE(1)
         ZKE(3) =  ZERO
      endif
#if HP3D_DEBUG
      if (iprint.eq.1) then
         write(*,*) 'apply_matrixK: Xp = ', Xp
         write(*,*) 'apply_matrixK: ZE = ', ZE
         write(*,*) 'apply_matrixK: rr = ', rr 
         write(*,*) 'apply_matrixK: ZKE= ', ZKE
      endif
#endif
!
   end subroutine apply_matrixK
!
!
!
!------------------------------------------------------------------------------
!> @brief      Returns complex-valued vector K*E, where K is the transformation  
!!             matrix arisen due to the envelope ansatz   exp(-ikR\theta)*E
!!
!> @param[in]  Mdle     - element middle node number
!> @param[in]  Xp       - physical point to evaluate result
!> @param[in]  ZE       - C^3 vector field E
!> @param[in]  ZDE      - C^3x3 tensor field, must be the gradient of E
!!
!> @param[out] ZDKE     - resulting C^3x3 tensor, gradient of K.E
!!
!> @date       Oct 2024
!------------------------------------------------------------------------------
   subroutine get_gradKE(mdle,Xp,ZE,ZDE,ZDKE)
!
      use data_structure3D
      use commonParam, only: RBEND, ENVELOPEK
      use parameters, only: ZERO, ZONE
!
      implicit none
!
      integer,   intent(in)  :: mdle
      real(8),   intent(in)  :: Xp(3)
      complex(8),intent(in)  :: ZE(3),ZDE(3,3)
      complex(8),intent(out) :: ZDKE(3,3)
!
      integer :: i
      real(8) :: y,z,rr,drr(3),e_y(3),e_z(3)
!
!------------------------------------------------------------------------------
!!
!  ...initialize matrix, copy values of y and z, and define unit vectors
      ZDKE = ZERO
      y = Xp(2); z = Xp(3);
      e_y = (/0.d0,1.d0,0.d0/); e_z = (/0.d0,0.d0,1.d0/)
!  ...Recall K.E = k * Rbend / (y^2 + z^2) * ( -y.E_y -z.E_z )
!                                            (     yE_x      )
!                                            (     zE_x      )  
! 
!  ...first compute coefficient Rbend / (y^2 + z^2)
      rr = ENVELOPEK * RBEND / (y**2+z**2)
!  ...and its gradient
      drr(1) = 0.d0 
      drr(2) = -2.d0*y * ENVELOPEK * RBEND / (y**2+z**2)**2
      drr(3) = -2.d0*z * ENVELOPEK * RBEND / (y**2+z**2)**2
!  ...then fill the gradient of K*E
      ZDKE(1,:) = drr(:)*(-y*ZE(2)-z*ZE(3))  + rr*(-y*ZDE(2,:)-e_y*ZE(2)   &
                                                   -z*ZDE(3,:)-e_z*ZE(3) )
      ZDKE(2,:) = drr(:)*  y*ZE(1)           + rr*( y*ZDE(1,:)+e_y*ZE(1) )
      ZDKE(3,:) = drr(:)*          z*ZE(1)   + rr*( z*ZDE(1,:)+e_z*ZE(1) )
!
   end subroutine get_gradKE

!------------------------------------------------------------------------------
!> @brief      Returns output of adjoint operator A^* on a pair of
!!             complex-valued H(curl) vectors [F ; G]
!!
!> @param[in]  Mdle     - element middle node number
!> @param[in]  Xp       - physical point to evaluate result
!> @param[in]  F,G      - Vector fields in test space
!> @param[in]  CF,CG    - Curl of F and G
!!
!> @param[out] Astar1    - 1st vector of output A^*[F ; G]
!> @param[out] Astar2    - 2nd vector of output A^*[F ; G]
!!
!> @date       Oct 2024 (last updated May 2025)
!------------------------------------------------------------------------------
   subroutine get_Astar(Mdle,Xp,F,G,CF,CG,Astar1,Astar2)
      use commonParam, only: ZERO, ZONE, ZI, OMEGA, EPSILON, MU
!
      implicit none
!
      integer,    intent(in) :: Mdle
      real(8),    intent(in) :: Xp(3)
      complex(8), intent(in) :: F(3),G(3),CF(3),CG(3)
      complex(8), intent(out):: Astar1(3),Astar2(3)
!
      complex(8):: ziKTF(3),ziKTG(3),zeps(3,3),zmu(3,3)
      complex(8):: zxpst(3),zdxpst(3),zd2xpst(3)
      complex(8):: zJ(3,3),zJinv(3,3),zJdet,zKpmlT(3,3)
      real(8) :: rQ(3,3),rK(3,3)
      integer :: activePML
!
!  ...get local permittivity and permeability, as complex tensors
      call get_local_epsilon(Mdle,Xp,zeps)
      call get_local_mu(Mdle,Xp,zmu)
!      
!  ...get matrix K
      call get_matrixK(Mdle,Xp,rK)
!  ...get i.K^T.F and i.K^T.G; We use skew-symmetry K^T = -K 
      ziKTF = -ZI*matmul(rK,F)
      ziKTG = -ZI*matmul(rK,G)
!
!  ...check if we're within the PML
      call is_pml(Mdle,Xp,activePML)
      ! if so, modify terms with the stretch jacobians
      if (activePML.ne.0) then 
         ! get rotation matrix
         call get_local_rotation(Mdle,Xp,activePML,rQ)
         ! get pml stretched coordinates and derivatives
         call get_stretched_coords(Mdle,Xp,ActivePML,Zxpst,Zdxpst,Zd2xpst)
         ! get PML Jacobian
         call get_stretch_J(zdxpst,rQ,zJ,zJinv,zJdet)
         ! modify tensors eps and mu (these are 3x3 matrices, matmul perhaps suffices)
         !       |J| Jinv @ eps @ Jinv^T
         zeps = zJdet*matmul(zJinv,matmul(zeps,transpose(zJinv)))
         !       |J| Jinv @  mu @ Jinv^T
         zmu  = zJdet*matmul(zJinv,matmul(zmu, transpose(zJinv)))
         ! modify operator K to include effects of stretched coords. Jacobian
         !       Kpml^T = Jinv @ K^T @ Jinv^T ; We use skew-symmetry K^T = -K 
         zKpmlT = matmul(zJinv,matmul(-rK,transpose(zJinv)))
         !
         ! ziKTF = ZERO
         !  get term     i . |J| . Kpml^T @ F
         call ZGEMM('N','N',3,1,3,zJdet,zKpmlT,3,F,3,ZERO,ziKTF,3)
         ! ziKTG = ZERO
         !  get term     i . |J| . Kpml^T @ G
         call ZGEMM('N','N',3,1,3,zJdet,zKpmlT,3,G,3,ZERO,ziKTG,3)
      endif
!  ...compute 1st vector of output A^*
      Astar1 = -conjg(ZI*OMEGA*matmul(transpose(zeps),F)) + CG - conjg(ziKTG)
!  ...compute 2nd vector of output A^*
      Astar2 =  conjg(ZI*OMEGA*matmul(transpose(zmu), G)) + CF - conjg(ziKTF)
!      
   end subroutine get_Astar

!------------------------------------------------------------------------------
!> @brief      Returns output of adjoint operator A^* on multiple pairs of
!!             real-valued H(curl) test functions [F ; 0] and [0 ; F]
!!
!> @param[in]  Mdle    - element middle node number
!> @param[in]  Xp      - physical point to evaluate result
!> @param[in]  NE      - number of test function pairs
!> @param[in]  F       - Hcurl test functions, size 3 by NE
!> @param[in]  CF      - Curl of F (columnwise)
!!
!> @param[out] AstarF1 - 1st block of output A^*[F ; 0]
!> @param[out] AstarF2 - 2nd block of output A^*[F ; 0]
!> @param[out] AstarG1 - 1st block of output A^*[0 ; G]
!> @param[out] AstarG2 - 2nd block of output A^*[0 ; G]
!!
!> @date       Mar 2025 (last updated May 2025)
!------------------------------------------------------------------------------
   subroutine get_Astar_multiple(Mdle,Xp,NE,F,CF,AstarF1,AstarF2,AstarG1,AstarG2)
      use parameters, only: ZERO, ZONE
      use commonParam, only: ZI, OMEGA, EPSILON, MU
!
      implicit none
!
      integer,    intent(in) :: Mdle,NE
      real(8),    intent(in) :: Xp(3),F(3,NE),CF(3,NE)
      complex(8), intent(out):: AstarF1(3,NE),AstarF2(3,NE),AstarG1(3,NE),AstarG2(3,NE)
!
      real(8) :: rKF(3,NE),rK(3,3)
      complex(8):: ziKTF(3,NE),zeps(3,3),zmu(3,3),zKpmlT(3,3)
      complex(8):: zxpst(3),zdxpst(3),zd2xpst(3),zJ(3,3),zJinv(3,3),zJdet
      real(8) :: rQ(3,3)
      integer :: activePML
!
!  ...get local permittivity and permeability, as complex tensors
      call get_local_epsilon(Mdle,Xp,zeps)
      call get_local_mu(Mdle,Xp,zmu)
!  ...get matrix K
      call get_matrixK(Mdle,Xp,rK)
!
!  ...check if we're within the PML
      call is_pml(Mdle,Xp,activePML)
      ! if so, modify eps, mu and K^T with the stretch jacobians
      if (activePML.ne.0) then 
         ! get rotation matrix
         call get_local_rotation(Mdle,Xp,activePML,rQ)
         ! get pml stretched coordinates and derivatives
         call get_stretched_coords(Mdle,Xp,ActivePML,Zxpst,Zdxpst,Zd2xpst)
         ! get PML Jacobian
         call get_stretch_J(zdxpst,rQ,zJ,zJinv,zJdet)
         ! modify tensors eps and mu (these are 3x3 matrices, matmul perhaps suffices)
         !       |J| Jinv @ eps @ Jinv^T
         zeps = zJdet*matmul(zJinv,matmul(zeps,transpose(zJinv)))
         !       |J| Jinv @  mu @ Jinv^T
         zmu  = zJdet*matmul(zJinv,matmul(zmu, transpose(zJinv)))
         ! modify operator K to include effects of stretched coords. Jacobian
         !       Kpml^T = Jinv @ K^T @ Jinv^T ; We use skew-symmetry K^T = -K 
         zKpmlT = matmul(zJinv,matmul(-rK,transpose(zJinv)))
         !
         ziKTF(:,:) = ZERO
         !  get term     i . |J| . Kpml^T @ F
         call ZGEMM('N','N',3,NE,3,ZI*zJdet,zKpmlT,3,cmplx(F,0.d0,8),3,ZERO,ziKTF,3)
      else 
         !      
!     ...apply matrix K to F
         call get_matrixK(Mdle,Xp,rK)
         rKF(:,:) = 0.d0
         call DGEMM('N','N',3,NE,3,1.d0,rK,3,F,3,0.d0,rKF,3)
!     ...get i.K^T @ F ; We use skew-symmetry K^T = -K 
         ziKTF = ZI*(-rKF)
      endif
!
!  ...multiply TRANSPOSED permittivity tensor and F
      call ZGEMM('T','N',3,NE,3,ZONE,zeps,3,cmplx(F,0.d0,8),3,ZERO,AstarF1,3)
!
!  ...multiply TRANSPOSED permeability tensor and G
      call ZGEMM('T','N',3,NE,3,ZONE,zmu, 3,cmplx(F,0.d0,8),3,ZERO,AstarG2,3)
!
!  ...compute 1st block of output A^* (F;0)
      AstarF1 = -conjg(ZI*OMEGA*AstarF1)
!  ...compute 2nd block of output A^* (F;0)
      AstarF2 = CF - conjg(ziKTF)
!  ...compute 1st block of output A^* (0;G)
      AstarG1 = AstarF2                     ! = CG - conjg(ziKTG)
!  ...compute 2nd block of output A^* (0;G)
      AstarG2 = conjg(ZI*OMEGA*AstarG2)
!      
   end subroutine get_Astar_multiple



!------------------------------------------------------------------------------
!> @brief      Returns output of direct operator A on a pair of
!!             complex-valued H(curl) vectors [E ; H]
!!
!> @param[in]  Mdle     - element middle node number
!> @param[in]  Xp       - physical point to evaluate result
!> @param[in]  E,H      - Vector fields in test space
!> @param[in]  CE,CH    - Curl of E and H
!!
!> @param[out] A1    - 1st vector of output A[E ; H]
!> @param[out] A2    - 2nd vector of output A[E ; H]
!!
!> @date       May 2025
!------------------------------------------------------------------------------
   subroutine get_A(Mdle,Xp,E,H,CE,CH,A1,A2)
      use commonParam, only: ZERO, ZONE, ZI, OMEGA, EPSILON, MU
!
      implicit none
!
      integer,    intent(in) :: Mdle
      real(8),    intent(in) :: Xp(3)
      complex(8), intent(in) :: E(3),H(3),CE(3),CH(3)
      complex(8), intent(out):: A1(3),A2(3)
!
      complex(8):: ziKE(3),ziKH(3),zeps(3,3),zmu(3,3)
      complex(8):: zxpst(3),zdxpst(3),zd2xpst(3)
      complex(8):: zJ(3,3),zJinv(3,3),zJdet,zKpml(3,3)
      real(8) :: rQ(3,3),rK(3,3)
      integer :: activePML
!
!  ...get local permittivity and permeability, as complex tensors
      call get_local_epsilon(Mdle,Xp,zeps)
      call get_local_mu(Mdle,Xp,zmu)
!      
!  ...get matrix K
      call get_matrixK(Mdle,Xp,rK)
!  ...get i.K.E and i.K.H
      ziKE = ZI*matmul(rK,E)
      ziKH = ZI*matmul(rK,H)
!
!  ...check if we're within the PML
      call is_pml(Mdle,Xp,activePML)
      ! if so, modify terms with the stretch jacobians
      if (activePML.ne.0) then 
         ! get rotation matrix
         call get_local_rotation(Mdle,Xp,activePML,rQ)
         ! get pml stretched coordinates and derivatives
         call get_stretched_coords(Mdle,Xp,ActivePML,Zxpst,Zdxpst,Zd2xpst)
         ! get PML Jacobian
         call get_stretch_J(zdxpst,rQ,zJ,zJinv,zJdet)
         ! modify tensors eps and mu (these are 3x3 matrices, matmul perhaps suffices)
         !       |J| Jinv @ eps @ Jinv^T
         zeps = zJdet*matmul(zJinv,matmul(zeps,transpose(zJinv)))
         !       |J| Jinv @  mu @ Jinv^T
         zmu  = zJdet*matmul(zJinv,matmul(zmu, transpose(zJinv)))
         ! modify operator K to include effects of stretched coords. Jacobian
         !       Kpml = Jinv @ K @ Jinv^T 
         zKpml = matmul(zJinv,matmul(-rK,transpose(zJinv)))
         !
         ! ziKE(:,:) = ZERO
         !  get term     i . |J| . Kpml @ E
         call ZGEMM('N','N',3,1,3,ZI*zJdet,zKpml,3,E,3,ZERO,ziKE,3)
         ! ziKH(:,:) = ZERO
         !  get term     i . |J| . Kpml @ G
         call ZGEMM('N','N',3,1,3,ZI*zJdet,zKpml,3,H,3,ZERO,ziKH,3)
      endif
!  ...compute 1st vector of output A
      A1 = -ZI*OMEGA*matmul(zeps,E) + CH - ziKH
!  ...compute 2nd vector of output A
      A2 =  ZI*OMEGA*matmul(zmu, H) + CE - ziKE
!      
   end subroutine get_A