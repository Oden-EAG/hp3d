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






   subroutine get_local_rotation(Mdle,Xp,RQ)
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






   subroutine get_pml_beta(Mdle,Xp,Zbeta,Zdbeta,Zd2beta)

      use commonParam, only: THETAEND,PMLPROP,OMEGA
      implicit none
      integer, intent(in) :: Mdle
      real(8), intent(in) :: Xp(3)
      complex(8),intent(out) :: Zbeta,Zdbeta,Zd2beta

      real(8) :: c,pn,f,df,d2f,th,th_pml,th_trn,th_dif
      integer :: n

      th = atan2(Xp(3),Xp(2))

      th_pml = THETAEND*PMLPROP
      th_trn = THETAEND - th_pml
      th_dif = Th - th_trn

      if (Th.gt.THETAEND.or.th_dif.lt.0.d0) then
         write(*,*) ' get_pml_beta: theta>THETAEND or theta<th_trn'
         write(*,*) ' theta,th_trn,THETAEND = ',Th,th_trn,THETAEND
         write(*,*) ' stop'
         stop
      endif

      n = 3       !!! NEEDS TO BE AT LEAST 2 !!!
      c = 50.d0 /(OMEGA* th_pml**n)
      f   = c*th_dif**n
      df  = c*th_dif**(n-1) * n
      d2f = c*th_dif**(n-2) * (n*(n-1))
      if((f.le.0.d0).or.(df.le.0.d0).or.(d2f.le.0)) then
         write(*,*) ' get_pml_beta: f, df,d2f are negative. stop.'
         stop
      endif
      zbeta   = cmplx( Th , -f  , 8)
      zdbeta  = cmplx(1.d0, -df , 8)
      zd2beta = cmplx(0.d0, -d2f, 8)

   end subroutine






   subroutine get_stretch_JQ(Zdbeta,RQ,ZJQ,ZJQinv)

      implicit none 
      complex(8), intent(in) :: Zdbeta
      real(8),    intent(in) :: RQ(3,3)
      complex(8), intent(out):: ZJQ(3,3),ZJQinv(3,3)

      complex(8) :: ZQaux(3,3)

      ! precompute J * Q
      ZQaux(1:2,:) = cmplx(RQ(1:2,:),0.d0,8)
      ZQaux(3,:) = RQ(3,:)*Zdbeta
      ! multiply Q^T * J * Q
      ZJQ = matmul(transpose(RQ),ZQaux)
      !
      ! precompute Jinv * Q
      ! ZQaux(1:2,:) = cmplx(RQ(1:2,:),0.d0,8)
      ZQaux(3,:) = RQ(3,:)/Zdbeta
      ! multiply Q^T * Jinv * Q
      ZJQinv = matmul(transpose(RQ),ZQaux)

   end subroutine





   subroutine is_pml(Mdle,Xp,ActivePML)
      use commonParam, only: PMLPROP,THETAEND
      implicit none 
      integer, intent(in) :: Mdle 
      real(8), intent(in) :: Xp(3)
      logical, intent(out):: ActivePML

      real(8) :: th

      ActivePML = .false.
      th = atan2(Xp(3),Xp(2))
      if (th.ge.(1.d0-PMLPROP)*THETAEND) ActivePML = .true.

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
!  ...Recall matrix K = k * Rbend / (y^2 + z^2) * ( 0   -y   -z )
!                                                 ( y    0    0 )
!                                                 ( z    0    0 )  
! 
!  ...first compute coefficient   k*Rbend / (y^2 + z^2)
      rr = ENVELOPEK * RBEND / (y**2+z**2)
!  ...then fill the non-zero entries
      RK(1,2) = -y*rr
      RK(1,3) = -z*rr
      RK(2,1) =  y*rr
      RK(3,1) =  z*rr
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



   subroutine get_refrac(Mdle,Zrefr)
      use commonParam, only: ZONE, OMEGA, EPSILON, MU, REFRCORE, REFRCLAD, REFRCOAT,ATTNCOAT
      implicit none
!
      integer,    intent(in) :: Mdle
      complex(8), intent(out):: Zrefr
!
      integer :: ndom
!
      call find_domain(Mdle, ndom)
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

   end subroutine
!------------------------------------------------------------------------------
!> @brief      Returns output of adjoint operator A^* on a pair of
!!             complex-valued H(curl) test functions [F ; G]
!!
!> @param[in]  Mdle     - element middle node number
!> @param[in]  Xp       - physical point to evaluate result
!> @param[in]  F,G      - Vector fields in test space
!> @param[in]  CF,CG    - Curl of F and G
!!
!> @param[out] Astar1    - 1st vector of output A^*[F ; G]
!> @param[out] Astar2    - 2nd vector of output A^*[F ; G]
!!
!> @date       Oct 2024
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
      complex(8):: zKF(3),zKG(3),zeps(3,3),zmu(3,3)
      complex(8):: zbeta,zdbeta,zd2beta,zJQ(3,3),zJQinv(3,3)
      real(8) :: rQ(3,3)
      logical :: activePML
!
!  ...get local permittivity and permeability, as complex tensors
      call get_local_epsilon(Mdle,Xp,zeps)
      call get_local_mu(Mdle,Xp,zmu)
!      
!  ...apply matrix K to F and G
      call apply_matrixK(Mdle,Xp,F,zKF)
      call apply_matrixK(Mdle,Xp,G,zKG)
!  ...get i.K^T.F and i.K^T.G; We use skew-symmetry K^T = -K 
      zKF = ZI*(-zKF)
      zKG = ZI*(-zKG)
!
!  ...check if we're within the PML
      call is_pml(Mdle,Xp,activePML)
      ! if so, modify terms with the stretch jacobians
      if (activePML) then 
         call get_pml_beta(Mdle,Xp,zbeta,zdbeta,zd2beta)
         call get_local_rotation(Mdle,Xp,rQ)
         call get_stretch_JQ(zdbeta,rQ,zJQ,zJQinv)
         !
         zeps = zdbeta*matmul(zJQinv,matmul(zeps,transpose(zJQinv)))
         zmu  = zdbeta*matmul(zJQinv,matmul(zmu, transpose(zJQinv)))
         zKF = zKF * zdbeta
         zKG = zKG * zdbeta
      endif
!  ...compute 1st vector of output A^*
      Astar1 = -conjg(ZI*OMEGA*matmul(zeps,F)) + CG - conjg(zKG)
!  ...compute 2nd vector of output A^*
      Astar2 =  conjg(ZI*OMEGA*matmul(zmu, G)) + CF - conjg(zKF)
!      
   end subroutine get_Astar

!------------------------------------------------------------------------------
!> @brief      Returns output of adjoint operator A^* on multiple pairs of
!!             complex-valued H(curl) test functions [F ; 0] and [0 ; F]
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
!> @date       Mar 2025
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
      complex(8):: ziKTF(3,NE),zeps(3,3),zmu(3,3)
      complex(8):: zbeta,zdbeta,zd2beta,zJQ(3,3),zJQinv(3,3)
      real(8) :: rQ(3,3)
      logical :: activePML
!
!  ...get local permittivity and permeability, as complex tensors
      call get_local_epsilon(Mdle,Xp,zeps)
      call get_local_mu(Mdle,Xp,zmu)
!      
!  ...apply matrix K to F
      call get_matrixK(Mdle,Xp,rK)
      rKF(:,:) = 0.d0
      call DGEMM('N','N',3,NE,3,1.d0,rK,3,F,3,0.d0,rKF,3)
!  ...get i.K^T.F ; We use skew-symmetry K^T = -K 
      ziKTF = ZI*(-rKF)
!
!  ...check if we're within the PML
      call is_pml(Mdle,Xp,activePML)
      ! if so, modify zeps, zmu and ziKTF with the stretch jacobians
      if (activePML) then 
         call get_pml_beta(Mdle,Xp,zbeta,zdbeta,zd2beta)
         call get_local_rotation(Mdle,Xp,rQ)
         call get_stretch_JQ(zdbeta,rQ,zJQ,zJQinv)
         !
         zeps = zdbeta*matmul(zJQinv,matmul(zeps,transpose(zJQinv)))
         !
         zmu  = zdbeta*matmul(zJQinv,matmul(zmu, transpose(zJQinv)))
         !
         ziKTF = ziKTF * zdbeta
      endif
!
!  ...multiply permittivity tensor and F
      call ZGEMM('N','N',3,NE,3,ZONE,zeps,3,cmplx(F,0.d0,8),3,ZERO,AstarF1,3)
!
!  ...multiply permeability tensor and G
      call ZGEMM('N','N',3,NE,3,ZONE,zmu, 3,cmplx(F,0.d0,8),3,ZERO,AstarG2,3)
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