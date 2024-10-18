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
!> @brief      Evaluates unconstrained stiffness matrix and load vector
!!
!> @param[in]  Mdle     - element middle node number
!> @param[in]  x        - physical point to evaluate permittivity
!!
!> @param[in]  eps      - permittivity tensor at point
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine get_permittivity(mdle,x, eps)
!
      use data_structure3D
      use commonParam
      use parameters, only: ZERO, ZONE
!
      implicit none
!
      integer,    intent(in)  :: mdle
      real(8),    intent(in)  :: x(3)
      complex(8), intent(out) :: eps(3,3)
!
      integer :: i
!
!------------------------------------------------------------------------------
!
! TODO: Implement your own custom permittivity here
!
!  ...set permittivity to identity for now.
      eps = ZERO
      do i=1,3
         eps(i,i) = ZONE
      enddo
!
   end subroutine get_permittivity




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
      use commonParam, only: RBEND, WAVENUM_K
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
      rr = WAVENUM_K * RBEND / (y**2+z**2)
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
      use commonParam, only: RBEND,WAVENUM_K
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
      rr = WAVENUM_K**2 * RBEND**2 / (y**2+z**2)**2
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
      use commonParam, only: RBEND, WAVENUM_K
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
!------------------------------------------------------------------------------
!!
!  ...initialize matrix and copy values of y and z
      ZKE = ZERO
      y = Xp(2); z = Xp(3);
!  ...Recall matrix K = k * Rbend / (y^2 + z^2) * ( 0   -y   -z )
!                                                 ( y    0    0 )
!                                                 ( z    0    0 )  
! 
!  ...first compute coefficient   k*Rbend / (y^2 + z^2)
      rr = WAVENUM_K * RBEND / (y**2+z**2)
!  ...then fill the entries of K*E
      ZKE(1) = -y*rr*ZE(2) -z*rr*ZE(3)
      ZKE(2) =  y*rr*ZE(1)
      ZKE(3) =  z*rr*ZE(1)
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
      use commonParam, only: RBEND, WAVENUM_K
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
      rr = WAVENUM_K * RBEND / (y**2+z**2)
!  ...and its gradient
      drr(1) = 0.d0 
      drr(2) = -2.d0*y * WAVENUM_K * RBEND / (y**2+z**2)**2
      drr(3) = -2.d0*z * WAVENUM_K * RBEND / (y**2+z**2)**2
!  ...then fill the gradient of K*E
      ZDKE(1,:) = drr(:)*(-y*ZE(2)-z*ZE(3))  + rr*(-y*ZDE(2,:)-e_y*ZE(2)   &
                                                   -z*ZDE(3,:)-e_z*ZE(3) )
      ZDKE(2,:) = drr(:)*  y*ZE(1)           + rr*( y*ZDE(1,:)+e_y*ZE(1) )
      ZDKE(3,:) = drr(:)*          z*ZE(1)   + rr*( z*ZDE(1,:)+e_z*ZE(1) )
!
   end subroutine get_gradKE
!------------------------------------------------------------------------------
!> @brief      Returns output of adjoint operator A^* on a pair of
!!             complex-valued H(curl) test functions F and G
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
      use parameters, only: ZERO, ZONE
      use commonParam, only: ZI, OMEGA, EPSILON, MU
!
      implicit none
!
      integer,    intent(in) :: Mdle
      real(8),    intent(in) :: Xp(3)
      complex(8), intent(in) :: F(3),G(3),CF(3),CG(3)
      complex(8), intent(out):: Astar1(3),Astar2(3)
!
      complex(8):: zKF(3),zKG(3)
!      
!  ...apply matrix K to F and G
      call apply_matrixK(Mdle,Xp,F,zKF)
      call apply_matrixK(Mdle,Xp,G,zKG)
!  ...get i.K^T.F and i.K^T.G; We use skew-symmetry K^T = -K 
      zKF = ZI*(-zKF)
      zKG = ZI*(-zKG)
!  ...compute 1st vector of output A^*
      Astar1 = -conjg(ZI*OMEGA*EPSILON)*F + CG - conjg(zKG)
!  ...compute 2nd vector of output A^*
      Astar2 =  conjg(ZI*OMEGA*MU     )*G + CF - conjg(zKF)
!      
   end subroutine get_Astar

