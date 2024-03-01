!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!                                  EDGES
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
!     routine name      - AncPhiE
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!     purpose:          - compute edge H1 ancillary functions and
!                         their gradients
!
!     arguments:
!
!     in:
!             S         - (s0,s1) affine coordinates associated to edge
!             DS        - gradients of S in R^N
!             Nord      - polynomial order
!             Idec      - Binary flag:
!                         = FALSE  s0+s1 != 1
!                         = TRUE   s0+s1  = 1
!             N         - spatial dimension
!
!     out:
!             PhiE      - values of edge H1 ancillary functions
!             DPhiE     - gradients of edge H1 ancillary functions
!
!----------------------------------------------------------------------
!
   subroutine AncPhiE(S,DS,Nord,Idec,N, PhiE,DPhiE)
!
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:1),DS(1:N,0:1)
      double precision, intent(out) ::   PhiE(2:Nord),DPhiE(1:N,2:Nord)
      integer ::                                              minI,maxI
!
!  ...local parameters
      minI = 2
      maxI = Nord
!
      if (N.lt.1) then
        write(*,7001) N
 7001   format('AncPhiE: N = ',i2)
      endif
!
!  ...these are precisely the homogenized Legendre polynomials
      call HomILegendre(S,DS,Nord,Idec,N, PhiE,DPhiE)
!
   end subroutine AncPhiE
!
!----------------------------------------------------------------------
!
!     routine name      - AncEE
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!     purpose:          - compute edge Hcurl ancillary functions and
!                         their curls
!
!     arguments:
!
!     in:
!             S         - (s0,s1) affine coordinates associated to edge
!             DS        - derivatives of S in R^N
!             Nord      - polynomial order
!             Idec      - Binary flag:
!                         = FALSE  s0+s1 != 1
!                         = TRUE   s0+s1  = 1
!             N         - spatial dimension
!
!     out:
!             EE        - edge Hcurl ancillary functions
!             CurlEE    - curls of edge Hcurl ancillary functions
!
!----------------------------------------------------------------------
!
   subroutine AncEE(S,DS,Nord,Idec,N, EE,CurlEE)
!
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:1),DS(1:N,0:1)
      double precision, intent(out) ::                 EE(1:N,0:Nord-1), &
                                               CurlEE(1:2*N-3,0:Nord-1)
      integer ::                                      minI,maxI,Ncurl,i
      double precision ::    homP(0:Nord-1),whiE(1:N),curlwhiE(1:2*N-3)
!
!  ...local parameters
      minI = 0
      maxI = Nord-1
      Ncurl = 2*N-3
!
      if (N.lt.2) then
        write(*,7001) N
 7001   format('AncEE: N = ',i2)
      endif
!
!  ...extract homogenized Legendre polyomials first
      call HomLegendre(S,maxI, homP)
!
!  ...simplified case
      if (Idec) then
        do i=minI,maxI
          EE(1:N,i) = homP(i)*DS(1:N,1)
        enddo
!    ...no need to compute Whitney function or curl
        CurlEE(1:Ncurl,minI:maxI) = 0.d0
!
!  ...in general
      else
!    ...lowest order Whitney function and its curl
        whiE = S(0)*DS(1:N,1)-S(1)*DS(1:N,0)
        call cross(N,DS(1:N,0),DS(1:N,1), curlwhiE)
!    ...now construct the higher order elements
        do i=minI,maxI
          EE(1:N,i) = homP(i)*whiE
          CurlEE(1:Ncurl,i) = (i+2)*homP(i)*curlwhiE
        enddo
      endif
!
   end subroutine AncEE
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!                         QUADRILATERAL FACES
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
!     routine name      - AncPhiQuad
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!     purpose:          - compute quadrilateral face H1 ancillary
!                         functions and their gradients
!
!     arguments:
!
!     in:
!             ST       - affine coordinates associated to face
!                        2x2 matrix [s0,s1;t0,t1]
!             DST      - gradients of ST
!             Nord     - (NordS,NordT) vector polynomial order
!             Idec     - (IdecS,IdecT) vector binary flag:
!                        IdecS = FALSE if s0+s1 != 1
!                        IdecS = TRUE  if s0+s1  = 1, same with IdecT
!             N        - spatial dimension
!
!     out:
!             PhiQuad  - quad H1 ancillary functions
!             DPhiQuad - gradients of quad H1 ancillary functions
!
!----------------------------------------------------------------------
!
   subroutine AncPhiQuad(ST,DST,Nord,Idec,N, PhiQuad,DPhiQuad)
!
      implicit none
      logical,          intent(in)  ::                          Idec(2)
      integer,          intent(in)  ::                        Nord(2),N
      double precision, intent(in)  ::     ST(0:1,1:2),DST(1:N,0:1,1:2)
      double precision, intent(out) ::     PhiQuad(2:Nord(1),2:Nord(2)), &
                                      DPhiQuad(1:N,2:Nord(1),2:Nord(2))
      integer ::                                minI,maxI,minJ,maxJ,i,j
      double precision ::        phiES(2:Nord(1)),DphiES(1:N,2:Nord(1)), &
                                 phiET(2:Nord(2)),DphiET(1:N,2:Nord(2))
!
!  ...local parameters
      minI = 2; maxI = Nord(1)
      minJ = 2; maxJ = Nord(2)
!
      if (N.lt.2) then
        write(*,7001) N
 7001   format('AncPhiQuad: N = ',i2)
      endif
!
!  ...get PhiE for each coordinate pair
      call AncPhiE(ST(0:1,1),DST(1:N,0:1,1),Nord(1),Idec(1),N, &
                phiES,DphiES)
      call AncPhiE(ST(0:1,2),DST(1:N,0:1,2),Nord(2),Idec(2),N, &
                phiET,DphiET)
!  ...the final result is the product of the two phiE
      do j=minJ,maxJ
        do i=minI,maxI
          PhiQuad(i,j) = phiES(i)*phiET(j)
          DphiQuad(1:N,i,j) = phiES(i)*DphiET(1:N,j) &
                            + phiET(j)*DphiES(1:N,i)
        enddo
      enddo
!
   end subroutine AncPhiQuad
!
!----------------------------------------------------------------------
!
!     routine name      - AncEQuad
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!     purpose:          - compute quadrilateral face Hcurl ancillary
!                         functions and their curls
!
!     arguments:
!
!     in:
!             ST       - affine coordinates associated to face
!                        2x2 matrix [s0,s1;t0,t1]
!             DST      - gradients of ST
!             Nord     - (NordS,NordT) vector polynomial order
!             Idec     - (IdecS,IdecT) vector binary flag:
!                        IdecS = FALSE if s0+s1 != 1
!                        IdecS = TRUE  if s0+s1  = 1, same with IdecT
!             N        - spatial dimension
!
!     out:
!             EQuad     - quad Hcurl ancillary functions
!             CurlEQuad - curls of quad Hcurl ancillary functions
!
!----------------------------------------------------------------------
!
   subroutine AncEQuad(ST,DST,Nord,Idec,N, EQuad,CurlEQuad)
!
      implicit none
      logical,          intent(in)  ::                          Idec(2)
      integer,          intent(in)  ::                        Nord(2),N
      double precision, intent(in)  ::     ST(0:1,1:2),DST(1:N,0:1,1:2)
      double precision, intent(out) :: EQuad(1:N,0:Nord(1)-1,2:Nord(2)), &
                               CurlEQuad(1:2*N-3,0:Nord(1)-1,2:Nord(2))
      integer ::                          minI,maxI,minJ,maxJ,Ncurl,i,j
      double precision :: &
                      EES(1:N,0:Nord(1)-1),curlEES(1:2*N-3,0:Nord(1)-1), &
             phiET(2:Nord(2)),DphiET(1:N,2:Nord(2)),DphiETxEES(1:2*N-3)
!
!  ...local parameters
      minI = 0; maxI = Nord(1)-1
      minJ = 2; maxJ = Nord(2)
      Ncurl = 2*N-3
!
      if (N.lt.2) then
        write(*,7001) N
 7001   format('AncEQuad: N = ',i2)
      endif
!
      call AncEE(ST(0:1,1),DST(1:N,0:1,1),Nord(1),Idec(1),N, &
                 EES,curlEES)
      call AncphiE(ST(0:1,2),DST(1:N,0:1,2),Nord(2),Idec(2),N, &
                 phiET,DphiET)
!
      do j=minJ,maxJ
        do i=minI,maxI
          EQuad(1:N,i,j) = EES(1:N,i)*phiET(j)
!
          call cross(N,DphiET(1:N,j),EES(1:N,i), DphiETxEES)
!
          CurlEQuad(1:Ncurl,i,j) = curlEES(1:Ncurl,i)*phiET(j) &
                                 + DphiETxEES(1:Ncurl)
        enddo
      enddo

   end subroutine AncEQuad
!
!----------------------------------------------------------------------
!
!     routine name      - AncVQuad
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!     purpose:          - compute quadrilateral face Hdiv ancillary
!                         functions and their divergences
!
!     arguments:
!
!     in:
!             ST       - affine coordinates associated to face
!                        2x2 matrix [s0,s1;t0,t1]
!             DST      - gradients of ST
!             Nord     - (NordS,NordT) vector polynomial order
!             Idec     - (IdecS,IdecT) vector binary flag:
!                        IdecS = FALSE if s0+s1 != 1
!                        IdecS = TRUE  if s0+s1  = 1, same with IdecT
!             N        - spatial dimension
!
!     out:
!             VQuad    - quad Hdiv ancillary functions
!             DivVQuad - divs of quad Hdiv ancillary functions
!
!----------------------------------------------------------------------
!
   subroutine AncVQuad(ST,DST,Nord,Idec,N, VQuad,DivVQuad)
!
      implicit none
      logical,          intent(in)  ::                          Idec(2)
      integer,          intent(in)  ::                        Nord(2),N
      double precision, intent(in)  ::     ST(0:1,1:2),DST(1:N,0:1,1:2)
      double precision, intent(out) :: &
                                     VQuad(1:N,0:Nord(1)-1,0:Nord(2)-1), &
                                      DivVQuad(0:Nord(1)-1,0:Nord(2)-1)
      integer ::                          minI,maxI,minJ,maxJ,Ncurl,i,j
      double precision :: &
                EES(1:N,0:Nord(1)-1),curlEES(1:2*N-3,0:Nord(1)-1),prod1, &
                EET(1:N,0:Nord(2)-1),curlEET(1:2*N-3,0:Nord(2)-1),prod2
!
!  ...local parameters
      minI = 0; maxI = Nord(1)-1
      minJ = 0; maxJ = Nord(2)-1
      Ncurl = 2*N-3
!
      if (N.lt.3) then
        write(*,7001) N
 7001   format('AncVQuad: N = ',i2)
      endif
!
      call AncEE(ST(0:1,1),DST(1:N,0:1,1),Nord(1),Idec(1),N, &
                 EES,curlEES)
      call AncEE(ST(0:1,2),DST(1:N,0:1,2),Nord(2),Idec(2),N, &
                 EET,curlEET)
!
!      ...slight speedup when Idec=(.true.,.true.)
      if (Idec(1).and.Idec(2)) then
        do j=minJ,maxJ
          do i=minI,maxI
            call cross(N,EES(1:N,i),EET(1:N,j), VQuad(1:N,i,j))
          enddo
        enddo
        DivVQuad(minI:maxI,minJ:maxJ) = 0.d0
      else
        do j=minJ,maxJ
          do i=minI,maxI
            call cross(N,EES(1:N,i),EET(1:N,j), VQuad(1:N,i,j))
!
            call dot_product(EET(1:N,j),curlEES(1:N,i), prod1)
            call dot_product(EES(1:N,i),curlEET(1:N,j), prod2)
!
            DivVQuad(i,j) = prod1-prod2
          enddo
        enddo
      endif
!
   end subroutine AncVQuad
!
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!                           TRIANGULAR FACES
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
!     routine name      - AncPhiTri
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!     purpose:          - compute triangle face H1 ancillary
!                         functions and their gradients
!
!     arguments:
!
!     in:
!             S        - (s0,s1,s2) affine coordinates associated to
!                        triangle face
!             DS       - derivatives of S0,S1,S2
!             Nord     - polynomial order
!             Idec     - Binary flag:
!                        = FALSE s0+s1+s2 != 1
!                        = TRUE  s0+s1+s2  = 1
!             N        - spatial dimension
!
!     out:
!             PhiTri   - triangle H1 ancillary functions
!             DPhiTri  - grads of triangle H1 ancillary functions
!
!----------------------------------------------------------------------
!
   subroutine AncPhiTri(S,DS,Nord,Idec,N, PhiTri,DPhiTri)
!
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:2),DS(1:N,0:2)
      double precision, intent(out) ::        PhiTri(2:Nord-1,1:Nord-2), &
                                         DPhiTri(1:N,2:Nord-1,1:Nord-2)
      logical ::                                                  IdecE
      integer ::       minI,maxI,minJ,maxJ,minIJ,maxIJ,minalpha,i,j,nij
      double precision ::                          sL(0:1),DsL(1:N,0:1), &
           phiE(2:Nord-1),DphiE(1:N,2:Nord-1),homLal(2:Nord-1,1:Nord-2), &
                                         DhomLal(1:N,2:Nord-1,1:Nord-2)
!
!  ...local parameters
      minI = 2; maxI = Nord-1
      minJ = 1; maxJ = Nord-2
      minIJ = minI+minJ; maxIJ = Nord
      minalpha = 2*minI
      IdecE = .false.
!
      if (N.lt.2) then
        write(*,7001) N
 7001   format('AncPhiTri: N = ',i2)
      endif
!
!  ...get PhiE - this is never a simplified case (IdecE=0)
      call AncPhiE(S(0:1),DS(1:N,0:1),Nord-minJ,IdecE,N, phiE,DphiE)
!
!  ...get homogenized Jacobi polynomials, homLal, and gradients
      sL(0) = S(0)+S(1); sL(1) = S(2)
      DsL(1:N,0) = DS(1:N,0)+DS(1:N,1)
      DsL(1:N,1) = DS(1:N,2)
      call HomIJacobi(sL,DsL,maxJ,minalpha,Idec,N, homLal,DhomLal)
!
!  ...simply complete the required information
      do nij=minIJ,maxIJ
        do i=minI,nij-minJ
          j=nij-i
          PhiTri(i,j) = phiE(i)*homLal(i,j)
          DPhiTri(1:N,i,j) = homLal(i,j)*DphiE(1:N,i) &
                              + phiE(i)*DhomLal(1:N,i,j)
        enddo
      enddo
!
   end subroutine AncPhiTri
!
!----------------------------------------------------------------------
!
!     routine name      - AncETri
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!     purpose:          - compute triangle face Hcurl ancillary
!                         functions and their curls
!
!     arguments:
!
!     in:
!             S        - (s0,s1,s2) affine coordinates associated to
!                        triangle face
!             DS       - derivatives of S0,S1,S2
!             Nord     - polynomial order
!             Idec     - Binary flag:
!                        = FALSE s0+s1+s2 != 1
!                        = TRUE  s0+s1+s2  = 1
!             N        - spatial dimension
!
!     out:
!             ETri     - triangle Hcurl ancillary functions
!             CurlETri - curls of triangle Hcurl ancillary functions
!
!----------------------------------------------------------------------
!
   subroutine AncETri(S,DS,Nord,Idec,N, ETri,CurlETri)
!
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:2),DS(1:N,0:2)
      double precision, intent(out) ::      ETri(1:N,0:Nord-2,1:Nord-1), &
                                    CurlETri(1:2*N-3,0:Nord-2,1:Nord-1)
      logical ::                                                  IdecE
      integer :: minI,maxI,minJ,maxJ,minIJ,maxIJ,minalpha,Ncurl,i,j,nij
      double precision ::     EE(1:N,0:Nord-2),curlEE(1:2*N-3,0:Nord-2), &
                         sL(0:1),DsL(1:N,0:1),homLal(0:Nord-2,1:Nord-1), &
                     DhomLal(1:N,0:Nord-2,1:Nord-1),DhomLalxEE(1:2*N-3)
!
!  ...local parameters
      minI = 0; maxI = Nord-2
      minJ = 1; maxJ = Nord-1
      minIJ = minI+minJ; maxIJ = Nord-1
      minalpha = 2*minI+1
      Ncurl = 2*N-3
      IdecE = .false.
!
      if (N.lt.2) then
        write(*,7001) N
 7001   format('AncETri: N = ',i2)
      endif
!
!  ...get EE - this is never a simplified case (IdecE=0)
      call AncEE(S(0:1),DS(1:N,0:1),Nord-minJ,IdecE,N, EE,curlEE)
!
!  ...get homogenized Integrated Jacobi polynomials, homLal, and gradients
      sL(0) = S(0)+S(1); sL(1) = S(2)
      DsL(1:N,0) = DS(1:N,0)+DS(1:N,1)
      DsL(1:N,1) = DS(1:N,2)
      call HomIJacobi(sL,DsL,maxJ,minalpha,Idec,N, homLal,DhomLal)
!
!  ...simply complete the required information
      do nij=minIJ,maxIJ
        do i=minI,nij-minJ
          j=nij-i
            ETri(1:N,i,j) = EE(1:N,i)*homLal(i,j)
!
            call cross(N,DhomLal(1:N,i,j),EE(1:N,i), DhomLalxEE)
!
            CurlETri(1:Ncurl,i,j) = homLal(i,j)*curlEE(1:Ncurl,i) &
                                  + DhomLalxEE
        enddo
      enddo
!
   end subroutine AncETri
!
!----------------------------------------------------------------------
!
!     routine name      - AncVTri
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!     purpose:          - compute triangle face Hcurl ancillary
!                         functions and their curls (family II)
!
!     arguments:
!
!     in:
!             S        - (s0,s1,s2) affine coordinates associated to
!                        triangle face
!             DS       - derivatives of S0,S1,S2
!             Nord     - polynomial order
!             Idec     - Binary flag:
!                        = FALSE s0+s1+s2 != 1
!                        = TRUE  s0+s1+s2  = 1
!             N        - spatial dimension
!
!     out:
!             VTri     - triangle Hdiv ancillary functions
!             DivVTri  - divs of triangle Hdiv ancillary functions
!
!----------------------------------------------------------------------
!
   subroutine AncVTri(S,DS,Nord,Idec,N, VTri,DivVTri)
!
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:2),DS(1:N,0:2)
      double precision, intent(out) ::      VTri(1:N,0:Nord-1,0:Nord-1), &
                                             DivVTri(0:Nord-1,0:Nord-1)
      integer ::       minI,maxI,minJ,maxJ,minIJ,maxIJ,minalpha,i,j,nij
      double precision ::      homP(0:Nord-1),homPal(0:Nord-1,0:Nord-1), &
              DS0xDS1(N),DS1xDS2(N),DS2xDS0(N),V00(N),tripleprod,psiTri
!
!  ...local parameters
      minI = 0; maxI = Nord-1
      minJ = 0; maxJ = Nord-1
      minIJ = minI+minJ; maxIJ = Nord-1
      minalpha = 2*minI+1
!
      if (N.lt.3) then
        write(*,7001) N
 7001   format('AncVTri: N = ',i2)
      endif
!
!  ...get homogenized Legendre polynomials, homP
      call HomLegendre(S(0:1),Nord-1-minJ, homP)
!
!  ...get homogenized Jacobi polynomials, homPal
      call HomJacobi((/S(0)+S(1),S(2)/),maxJ,minalpha, homPal)
!
!  ...simplified case
      if (Idec) then
!    ...construct V00
        call cross(N,DS(1:N,1),DS(1:N,2), V00)
!    ...loop
        do nij=minIJ,maxIJ
          do i=minI,nij-minJ
            j=nij-i
            VTri(1:N,i,j) = homP(i)*homPal(i,j)*V00(1:N)
          enddo
        enddo
!
        DivVTri = 0.d0
!
!  ...general case
      else
!    ...construct V00
        call cross(N,DS(1:N,0),DS(1:N,1), DS0xDS1)
        call cross(N,DS(1:N,1),DS(1:N,2), DS1xDS2)
        call cross(N,DS(1:N,2),DS(1:N,0), DS2xDS0)
        V00 = S(0)*DS1xDS2+S(1)*DS2xDS0+S(2)*DS0xDS1
!    ...loop
        do nij=minIJ,maxIJ
          do i=minI,nij-minJ
            j=nij-i
            psiTri = homP(i)*homPal(i,j)
!
            VTri(1:N,i,j) = psiTri*V00
            DivVTri(i,j)  = (nij+3)*psiTri
          enddo
        enddo
!
        call dot_product(DS(1:N,0),DS1xDS2, tripleprod)
!
        DivVTri = DivVTri*tripleprod
      endif
!
   end subroutine AncVTri
!
