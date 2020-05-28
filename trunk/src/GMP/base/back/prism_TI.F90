!----------------------------------------------------------------------
!   routine name       - prism_TI
!---------------------------------------------------------------------
!   latest revision    - Jul 09
!
!   purpose            - routine defines the parameterization
!                        for a prism based on a transfinite
!                        interpolation (with linear and singular
!                        blending functions) of the paramaterization
!                        of its faces
!
!   arguments :
!     in:
!               No     - prism number
!               Eta    - reference coordinates of a point
!                        in the reference prism
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        to the reference coordinates
!---------------------------------------------------------------------
      subroutine prism_TI_old(No,Eta, X,Dxdeta)
!
      use control
      use GMP
      use element_data
#include "syscom.blk"
      common /cprism_TI/ iprint
      common /ctrianB/ iprint_trianB
!
      dimension Eta(3), X(3),Dxdeta(3,3)
!
!  ...derivatives of edge coordinate
      dimension dtedeta(3)
!
!  ...affine coordinates for triangular faces
      dimension vshapt(3),dvshapt(2,3)
!
!  ...1D vertex shape functions in the xi_3 direction
      dimension vshap(2),dvshap(2)
!
!  ...blending function
      dimension dblend(2)
!
!  ...edge kernels
      dimension xe(3),dxedt(3)
!
!  ...face kernels
      dimension xf(3),dxfdtf(3,2)
!
      dimension phi(3),dphi_dEta(3,3)
!-----------------------------------------------------------------------
!
      iprint_trianB = iprint
      if (iprint .eq. 1) then
        write(*,7001) No,Eta(1:3)
 7001   format('prism_TI: No,Eta = ',i5,2x,3e12.5)
      endif
!
!  ...affine coordinates for the triangular faces
      vshapt(1) = 1.d0 - Eta(1) - Eta(2)
      dvshapt(1,1) = -1.d0; dvshapt(2,1) = -1.d0
      vshapt(2) = Eta(1)
      dvshapt(1,2) = 1.d0; dvshapt(2,2) = 0.d0
      vshapt(3) = Eta(2)
      dvshapt(1,3) = 0.d0; dvshapt(2,3) = 1.d0
!
!  ...1D shape functions in the xi_3 direction
      vshap(1) = 1.d0 - Eta(3); dvshap(1) = -1.d0
      vshap(2) = Eta(3); dvshap(2) = 1.d0
!
! *********************************  VERTEX INTERPOLANT  ********************************* |
      X(1:3) = 0.d0; Dxdeta(1:3,1:3) = 0.d0
      iv = 0
      do j = 1, 2
        do i = 1, 3
          iv = iv + 1
          np = PRISMS(No)%VertNo(iv)
          X(1:3) = X(1:3) + POINTS(np)%Rdata(1:3)*vshapt(i)*vshap(j)
          do k = 1, 2
            Dxdeta(1:3,k) = Dxdeta(1:3,k) &
                          + POINTS(np)%Rdata(1:3)*dvshapt(k,i)*vshap(j)
          enddo
          Dxdeta(1:3,3) = Dxdeta(1:3,3) &
                        + POINTS(np)%Rdata(1:3)*vshapt(i)*dvshap(j)
        enddo
      enddo
      if (iprint .eq. 1) then
        write(*,*) 'prism_TI: VERTEX INTERPOLANT = '
        do ivar = 1, 3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
 7011     format(e12.5,3x,3e12.5)
        enddo
      endif
!
! *******************************  HORIZONTAL EDGE BUBBLES  ****************************** |
      ie = 0
      do j = 1, 2
        do i = 1, 3
          ie = ie + 1
          nc = PRISMS(No)%EdgeNo(ie);  norient = 0
          if (nc .lt. 0) then
            nc = -nc;  norient = 1
          endif
          if (CURVES(nc)%Type .eq. 'Seglin')  cycle
!  .......get the edge vertices specifying the local edge orientation
          iv1 = TRIAN_EDGE_TO_VERT(1,i);  iv2 = TRIAN_EDGE_TO_VERT(2,i)
!  .......project Eta(1:2) onto the edge
          call proj_t2e(Eta(1:2),iv1,iv2,vshapt,dvshapt, te,dtedeta(1:2))
! ........if edge enpoint cycle
          if ((abs(te) .lt. GEOM_TOL) .or. (abs(1.d0 - te) .lt. GEOM_TOL))  cycle
          if (iprint .eq. 1) then
            write(*,7012) ie,nc,CURVES(nc)%Type
 7012       format('prism_TI: ie,nc,Type = ',i2,i5,2x,a5)
          endif
!  .......evaluate edge kernel function
          call curveK(nc,te,norient, xe,dxedt)
!  .......2D blending function
          blend = vshapt(iv1)*vshapt(iv2)
          dblend(1:2) = dvshapt(1:2,iv1)*vshapt(iv2) + vshapt(iv1)*dvshapt(1:2,iv2)
!  .......add edge contribution
          X(1:3) = X(1:3) + xe(1:3)*blend*vshap(j)
          do k = 1, 2
            Dxdeta(1:3,k) = Dxdeta(1:3,k)                          &
                          + dxedt(1:3)*dtedeta(k)*blend*vshap(j)   &
                          + xe(1:3)*dblend(k)*vshap(j)
          enddo
          Dxdeta(1:3,3) = Dxdeta(1:3,3) + xe(1:3)*blend*dvshap(j)
!  .......printing statement
          if (iprint .eq. 1) then
            write(*,*)'prism_TI: AFTER HORIZONTAL EDGE i,j',i,j
            do ivar = 1, 3
              write(*,7011) X(ivar),Dxdeta(ivar,1:3)
            enddo
          endif
        enddo
      enddo
!
! ********************************  VERTICAL EDGE BUBBLES  ******************************* |
      do i = 1, 3
        ie = ie + 1
        nc = PRISMS(No)%EdgeNo(ie);  norient = 0
        if (nc .lt. 0) then
          nc = -nc;  norient = 1
        endif
        if (CURVES(nc)%Type.eq.'Seglin')  cycle
!  .....evaluate edge bubble function
        call curveB(nc,Eta(3),norient, xe,dxedt)
!  .....2D blending function
        blend = vshapt(i)
        dblend(1:2) = dvshapt(1:2,i)
!  .....add edge contribution
        X(1:3) = X(1:3) + xe(1:3)*blend
        do k = 1, 2
          Dxdeta(1:3,k) = Dxdeta(1:3,k) + xe(1:3)*dblend(k)
        enddo
        Dxdeta(1:3,3) = Dxdeta(1:3,3) + dxedt(1:3)*blend
      enddo
!  ...printing statement
      if (iprint .eq. 1) then
        write(*,*) 'prism_TI: AFTER VERTICAL EDGES = '
        do ivar = 1, 3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
        enddo
      endif
!
! *******************************  HORIZONTAL FACE BUBBLES  ****************************** |
      ifig = 0
      do j = 1, 2
        ifig = ifig + 1
        call decode(PRISMS(No)%FigNo(ifig), nt,norient)
        if ((TRIANGLES(nt)%Type .eq. 'TransTri') .or. (TRIANGLES(nt)%Type .eq. 'PlaneTri'))  cycle
!  .....printing statement
        if (iprint .eq. 1) then
          write(*,7013) ifig,nt,TRIANGLES(nt)%Type
 7013     format('prism_TI: ifig,nt,Type = ',i2,i5,2x,a5)
        endif
!  .....compute the face bubble
        call trianB(nt,Eta(1:2),norient, xf,dxfdtf)
!  .....printing statement
        if (iprint .eq. 1) then
          do ivar = 1, 3
            write(*,7033) ivar,xf(ivar),dxfdtf(ivar,1:2)
 7033       format('prism_TI: ivar,xf,dxdtf = ',i2,2x,e12.5,2x,2e12.5)
          enddo
        endif
!  .....add face contribution
        X(1:3) = X(1:3) + xf(1:3)*vshap(j)
        do k = 1, 2
          Dxdeta(1:3,k) = Dxdeta(1:3,k) + dxfdtf(1:3,k)*vshap(j)
        enddo
        Dxdeta(1:3,3) = Dxdeta(1:3,3) + xf(1:3)*dvshap(j)
      enddo
!  ...printing statement
      if (iprint .eq. 1) then
        write(*,*) 'prism_TI: AFTER HORIZONTAL FACES = '
        do ivar = 1, 3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
        enddo
        call pause
      endif
!
! *******************************  VERTICAL FACE BUBBLES  ****************************** |
      do ii = 1, 3
        ifig = ifig + 1
        call decode(PRISMS(No)%FigNo(ifig), nr,norient)
        if ((RECTANGLES(nr)%Type .eq. 'BilQua') .or. (RECTANGLES(nr)%Type .eq. 'TraQua') .or. &
            (RECTANGLES(nr)%Type .eq. 'PTIRec'))   cycle
!  .....printing statement
        if (iprint .eq. 1) then
          write(*,7014) ifig,nr,RECTANGLES(nr)%Type
 7014     format('prism_TI: ifig,nr,Type = ',i2,i5,2x,a5)
        endif
! ......add face bubble
        call face_bubble_prism(Eta,nr,ifig,norient, phi,dphi_dEta)
        X      = X + phi
        Dxdeta = Dxdeta + dphi_dEta
      enddo
!
end subroutine prism_TI_old
