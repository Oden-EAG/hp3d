c
c-----------------------------------------------------------------------
c
c     routine name      - shape3bEM2
c
c-----------------------------------------------------------------------
c
c     latest revision:  - Mar 06
c
c     purpose:          - evaluate scalar H^1 and vector H(curl) shape
c                         functions for a brick element and derivatives
c                         (gradient for H^1, curl for H(curl)) wrt
c                         master element coordinates
c
c     arguments:
c
c     in:
c             Xi        - master element coordinates
c             Norder    - order of approximation of P-field
c     out:
c             NrdofP    - number of scalar shape functions
c             ShapP     - scalar shape functions
c             GradP     - gradient of scalar shape functions
c             NrdofE    - number of vector shape functions
c             ShapE     - vector shape functions
c             CurlE     - curl of vector shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3bEM2(Xi,Norder, NrdofP,ShapP,GradP,
     .                                 NrdofE,ShapE,CurlE)
c
      use parameters   , only : NPMAX=>MAXP

#if C_MODE
      implicit real*8 (a-h,o-y), complex*16 (z)
#else
      implicit double precision (a-h,o-z)
#endif

      dimension Xi(3),Norder(19)
      dimension ShapP(*),GradP(3,*),ShapE(3,*),CurlE(3,*)
c
      dimension P1(NPMAX+1),D1(NPMAX+1),E1(NPMAX)
      dimension P2(NPMAX+1),D2(NPMAX+1),E2(NPMAX)
      dimension P3(NPMAX+1),D3(NPMAX+1),E3(NPMAX)
c
c-----------------------------------------------------------------------
c
c     evaluate scalar shape functions and gradients
c$$$!!!!!!      call shape3b(Xi,Norder, NrdofP,ShapP,GradP)
c
c     decode middle node order
      call decode(Norder(19), nordh, nord3)
      call decode(nordh, nord1,nord2)
c
c     evaluate 1D shape functions and E shape functions
      call shape1_Jason (Xi(1),nord1, nvoid,P1,D1)
      call shape1E_Jason(Xi(1),nord1, nvoid,E1)
c
      call shape1_Jason (Xi(2),nord2, nvoid,P2,D2)
      call shape1E_Jason(Xi(2),nord2, nvoid,E2)
c
      call shape1_Jason (Xi(3),nord3, nvoid,P3,D3)
      call shape1E_Jason(Xi(3),nord3, nvoid,E3)
c
      n = 0
c
c-----------------------------------------------------------------------
c
c     edge 1 (direction 1)
      do i=1,Norder(1)
        n = n+1
c                       i     1     1                   i     1     1
        ShapE(1,n) = E1(i)*P2(1)*P3(1); CurlE(1,n) = 0.d0
        ShapE(2,n) = 0.d0;              CurlE(2,n) = E1(i)*P2(1)*D3(1)
        ShapE(3,n) = 0.d0;              CurlE(3,n) =-E1(i)*D2(1)*P3(1)
c
c        ShapE(1:3,n) = (/ (1-Xi(2))*(1-Xi(3)), 0.d0, 0.d0 /)
c        CurlE(1:3,n) = (/ 0.d0, Xi(2)-1, 1-Xi(3) /)
      enddo
c
c     edge 2 (direction 2)
      do i=1,Norder(2)
        n = n+1
c                       2     i     1                   2     i     1
        ShapE(1,n) = 0.d0;              CurlE(1,n) =-P1(2)*E2(i)*D3(1)
        ShapE(2,n) = P1(2)*E2(i)*P3(1); CurlE(2,n) = 0.d0
        ShapE(3,n) = 0.d0;              CurlE(3,n) = D1(2)*E2(i)*P3(1)
c
c        ShapE(1:3,n) = (/ 0.d0, Xi(1)*(1-Xi(3)), 0.d0 /)
c        CurlE(1:3,n) = (/ Xi(1), 0.d0, 1-Xi(3) /)
      enddo
c
c     edge 3 (direction 1)
      do i=1,Norder(3)
        n = n+1
c                       i     2     1                   i     2     1
        ShapE(1,n) = E1(i)*P2(2)*P3(1); CurlE(1,n) = 0.d0
        ShapE(2,n) = 0.d0;              CurlE(2,n) = E1(i)*P2(2)*D3(1)
        ShapE(3,n) = 0.d0;              CurlE(3,n) =-E1(i)*D2(2)*P3(1)
c
c        ShapE(1:3,n) = (/ Xi(2)*(1-Xi(3)), 0.d0, 0.d0 /)
c        CurlE(1:3,n) = (/ 0.d0, -Xi(2), Xi(3)-1 /)
      enddo
c
c     edge 4 (direction 2)
      do i=1,Norder(4)
        n = n+1
c                       1     i     1                   1     i     1
        ShapE(1,n) = 0.d0;              CurlE(1,n) =-P1(1)*E2(i)*D3(1)
        ShapE(2,n) = P1(1)*E2(i)*P3(1); CurlE(2,n) = 0.d0
        ShapE(3,n) = 0.d0;              CurlE(3,n) = D1(1)*E2(i)*P3(1)
c
c        ShapE(1:3,n) = (/ 0.d0, (1-Xi(1))*(1-Xi(3)), 0.d0 /)
c        CurlE(1:3,n) = (/ 1-Xi(1), 0.d0, Xi(3)-1 /)
      enddo
c
c-----------------------------------------------------------------------
c
c     edge 5 (direction 1)
      do i=1,Norder(5)
        n = n+1
c                       i     1     2                   i     1     2
        ShapE(1,n) = E1(i)*P2(1)*P3(2); CurlE(1,n) = 0.d0
        ShapE(2,n) = 0.d0;              CurlE(2,n) = E1(i)*P2(1)*D3(2)
        ShapE(3,n) = 0.d0;              CurlE(3,n) =-E1(i)*D2(1)*P3(2)
c
c        ShapE(1:3,n) = (/ (1-Xi(2))*Xi(3), 0.d0, 0.d0 /)
c        CurlE(1:3,n) = (/ 0.d0, 1-Xi(2), Xi(3) /)
      enddo
c
c     edge 6 (direction 2)
      do i=1,Norder(6)
        n = n+1
c                       2     i     2                   2     i     2
        ShapE(1,n) = 0.d0;              CurlE(1,n) =-P1(2)*E2(i)*D3(2)
        ShapE(2,n) = P1(2)*E2(i)*P3(2); CurlE(2,n) = 0.d0
        ShapE(3,n) = 0.d0;              CurlE(3,n) = D1(2)*E2(i)*P3(2)
c
c        ShapE(1:3,n) = (/ 0.d0, Xi(1)*Xi(3), 0.d0 /)
c        CurlE(1:3,n) = (/ -Xi(1), 0.d0, Xi(3) /)
      enddo
c
c     edge 7 (direction 1)
      do i=1,Norder(7)
        n = n+1
c                       i     2     2                   i     2     2
        ShapE(1,n) = E1(i)*P2(2)*P3(2); CurlE(1,n) = 0.d0
        ShapE(2,n) = 0.d0;              CurlE(2,n) = E1(i)*P2(2)*D3(2)
        ShapE(3,n) = 0.d0;              CurlE(3,n) =-E1(i)*D2(2)*P3(2)
c
c        ShapE(1:3,n) = (/ Xi(2)*Xi(3), 0.d0, 0.d0 /)
c        CurlE(1:3,n) = (/ 0.d0, Xi(2), -Xi(3) /)
      enddo
c
c     edge 8 (direction 2)
      do i=1,Norder(8)
        n = n+1
c                       1     i     2                   1     i     2
        ShapE(1,n) = 0.d0;              CurlE(1,n) =-P1(1)*E2(i)*D3(2)
        ShapE(2,n) = P1(1)*E2(i)*P3(2); CurlE(2,n) = 0.d0
        ShapE(3,n) = 0.d0;              CurlE(3,n) = D1(1)*E2(i)*P3(2)
c
c        ShapE(1:3,n) = (/ 0.d0, (1-Xi(1))*Xi(3), 0.d0 /)
c        CurlE(1:3,n) = (/ Xi(1)-1, 0.d0, -Xi(3) /)
      enddo
c
c-----------------------------------------------------------------------
c
c     edge 9 (direction 3)
      do i=1,Norder(9)
        n = n+1
c                       1     1     i                   1     1     i
        ShapE(1,n) = 0.d0;              CurlE(1,n) = P1(1)*D2(1)*E3(i)
        ShapE(2,n) = 0.d0;              CurlE(2,n) =-D1(1)*P2(1)*E3(i)
        ShapE(3,n) = P1(1)*P2(1)*E3(i); CurlE(3,n) = 0.d0
c
c        ShapE(1:3,n) = (/ 0.d0, 0.d0, (1-Xi(1))*(1-Xi(2)) /)
c        CurlE(1:3,n) = (/ Xi(1)-1, 1-Xi(2), 0.d0 /)
      enddo
c
c     edge 10 (direction 3)
      do i=1,Norder(10)
        n = n+1
c                       2     1     i                   2     1     i
        ShapE(1,n) = 0.d0;              CurlE(1,n) = P1(2)*D2(1)*E3(i)
        ShapE(2,n) = 0.d0;              CurlE(2,n) =-D1(2)*P2(1)*E3(i)
        ShapE(3,n) = P1(2)*P2(1)*E3(i); CurlE(3,n) = 0.d0
c
c        ShapE(1:3,n) = (/ 0.d0, 0.d0, Xi(1)*(1-Xi(2)) /)
c        CurlE(1:3,n) = (/ -Xi(1), Xi(2)-1, 0.d0 /)
      enddo
c
c     edge 11 (direction 3)
      do i=1,Norder(11)
        n = n+1
c                       2     2     i                   2     2     i
        ShapE(1,n) = 0.d0;              CurlE(1,n) = P1(2)*D2(2)*E3(i)
        ShapE(2,n) = 0.d0;              CurlE(2,n) =-D1(2)*P2(2)*E3(i)
        ShapE(3,n) = P1(2)*P2(2)*E3(i); CurlE(3,n) = 0.d0
c
c        ShapE(1:3,n) = (/ 0.d0, 0.d0, Xi(1)*Xi(2) /)
c        CurlE(1:3,n) = (/ Xi(1), -Xi(2), 0.d0 /)
      enddo
c
c     edge 12 (direction 3)
      do i=1,Norder(12)
        n = n+1
c                       1     2     i                   1     2     i
        ShapE(1,n) = 0.d0;              CurlE(1,n) = P1(1)*D2(2)*E3(i)
        ShapE(2,n) = 0.d0;              CurlE(2,n) =-D1(1)*P2(2)*E3(i)
        ShapE(3,n) = P1(1)*P2(2)*E3(i); CurlE(3,n) = 0.d0
c
c        ShapE(1:3,n) = (/ 0.d0, 0.d0, (1-Xi(1))*Xi(2) /)
c        CurlE(1:3,n) = (/ 1-Xi(1), Xi(2), 0.d0 /)
      enddo
c
c-----------------------------------------------------------------------
c
c     face 1
      call decode_orderq(Norder(12+1), nord1,nord2)
c     direction 1
      do j=3,nord2+1
      do i=1,nord1
        n = n+1
        ShapE(1,n) = E1(i)*P2(j)*P3(1); CurlE(1,n) = 0.d0
        ShapE(2,n) = 0.d0;              CurlE(2,n) = E1(i)*P2(j)*D3(1)
        ShapE(3,n) = 0.d0;              CurlE(3,n) =-E1(i)*D2(j)*P3(1)
      enddo
      enddo
c     direction 2
      do j=1,nord2
      do i=3,nord1+1
        n = n+1
        ShapE(1,n) = 0.d0;              CurlE(1,n) =-P1(i)*E2(j)*D3(1)
        ShapE(2,n) = P1(i)*E2(j)*P3(1); CurlE(2,n) = 0.d0
        ShapE(3,n) = 0.d0;              CurlE(3,n) = D1(i)*E2(j)*P3(1)
      enddo
      enddo
c
c     face 2
      call decode_orderq(Norder(12+2), nord1,nord2)
c     direction 1
      do j=3,nord2+1
      do i=1,nord1
        n = n+1
        ShapE(1,n) = E1(i)*P2(j)*P3(2); CurlE(1,n) = 0.d0
        ShapE(2,n) = 0.d0;              CurlE(2,n) = E1(i)*P2(j)*D3(2)
        ShapE(3,n) = 0.d0;              CurlE(3,n) =-E1(i)*D2(j)*P3(2)
      enddo
      enddo
c     direction 2
      do j=1,nord2
      do i=3,nord1+1
        n = n+1
        ShapE(1,n) = 0.d0;              CurlE(1,n) =-P1(i)*E2(j)*D3(2)
        ShapE(2,n) = P1(i)*E2(j)*P3(2); CurlE(2,n) = 0.d0
        ShapE(3,n) = 0.d0;              CurlE(3,n) = D1(i)*E2(j)*P3(2)
      enddo
      enddo
c
c     face 3
      call decode_orderq(Norder(12+3), nord1,nord2)
c     direction 1
      do j=3,nord2+1
      do i=1,nord1
        n = n+1
        ShapE(1,n) = E1(i)*P2(1)*P3(j); CurlE(1,n) = 0.d0
        ShapE(2,n) = 0.d0;              CurlE(2,n) = E1(i)*P2(1)*D3(j)
        ShapE(3,n) = 0.d0;              CurlE(3,n) =-E1(i)*D2(1)*P3(j)
      enddo
      enddo
c     direction 3
      do j=1,nord2
      do i=3,nord1+1
        n = n+1
        ShapE(1,n) = 0.d0;              CurlE(1,n) = P1(i)*D2(1)*E3(j)
        ShapE(2,n) = 0.d0;              CurlE(2,n) =-D1(i)*P2(1)*E3(j)
        ShapE(3,n) = P1(i)*P2(1)*E3(j); CurlE(3,n) = 0.d0
      enddo
      enddo
c
c     face 4
      call decode_orderq(Norder(12+4), nord1,nord2)
c     direction 2
      do j=3,nord2+1
      do i=1,nord1
        n = n+1
        ShapE(1,n) = 0.d0;              CurlE(1,n) =-P1(2)*E2(i)*D3(j)
        ShapE(2,n) = P1(2)*E2(i)*P3(j); CurlE(2,n) = 0.d0
        ShapE(3,n) = 0.d0;              CurlE(3,n) = D1(2)*E2(i)*P3(j)
      enddo
      enddo
c     direction 3
      do j=1,nord2
      do i=3,nord1+1
        n = n+1
        ShapE(1,n) = 0.d0;              CurlE(1,n) = P1(2)*D2(i)*E3(j)
        ShapE(2,n) = 0.d0;              CurlE(2,n) =-D1(2)*P2(i)*E3(j)
        ShapE(3,n) = P1(2)*P2(i)*E3(j); CurlE(3,n) = 0.d0
      enddo
      enddo
c
c     face 5
      call decode_orderq(Norder(12+5), nord1,nord2)
c     direction 1
      do j=3,nord2+1
      do i=1,nord1
        n = n+1
        ShapE(1,n) = E1(i)*P2(2)*P3(j); CurlE(1,n) = 0.d0
        ShapE(2,n) = 0.d0;              CurlE(2,n) = E1(i)*P2(2)*D3(j)
        ShapE(3,n) = 0.d0;              CurlE(3,n) =-E1(i)*D2(2)*P3(j)
      enddo
      enddo
c     direction 3
      do j=1,nord2
      do i=3,nord1+1
        n = n+1
        ShapE(1,n) = 0.d0;              CurlE(1,n) = P1(i)*D2(2)*E3(j)
        ShapE(2,n) = 0.d0;              CurlE(2,n) =-D1(i)*P2(2)*E3(j)
        ShapE(3,n) = P1(i)*P2(2)*E3(j); CurlE(3,n) = 0.d0
      enddo
      enddo
c
c     face 6
      call decode_orderq(Norder(12+6), nord1,nord2)
c     direction 2
      do j=3,nord2+1
      do i=1,nord1
        n = n+1
        ShapE(1,n) = 0.d0;              CurlE(1,n) =-P1(1)*E2(i)*D3(j)
        ShapE(2,n) = P1(1)*E2(i)*P3(j); CurlE(2,n) = 0.d0
        ShapE(3,n) = 0.d0;              CurlE(3,n) = D1(1)*E2(i)*P3(j)
      enddo
      enddo
c     direction 3
      do j=1,nord2
      do i=3,nord1+1
        n = n+1
        ShapE(1,n) = 0.d0;              CurlE(1,n) = P1(1)*D2(i)*E3(j)
        ShapE(2,n) = 0.d0;              CurlE(2,n) =-D1(1)*P2(i)*E3(j)
        ShapE(3,n) = P1(1)*P2(i)*E3(j); CurlE(3,n) = 0.d0
      enddo
      enddo
c
c-----------------------------------------------------------------------
c
c     middle node
      call decode(Norder(19), nordh, nord3)
      call decode(nordh, nord1,nord2)

c     direction 1
      do k=3,nord3+1
      do j=3,nord2+1
      do i=1,nord1
        n = n+1
        ShapE(1,n) = E1(i)*P2(j)*P3(k); CurlE(1,n) = 0.d0
        ShapE(2,n) = 0.d0;              CurlE(2,n) = E1(i)*P2(j)*D3(k)
        ShapE(3,n) = 0.d0;              CurlE(3,n) =-E1(i)*D2(j)*P3(k)
      enddo
      enddo
      enddo
c     direction 2
      do k=3,nord3+1
      do j=1,nord2
      do i=3,nord1+1
        n = n+1
        ShapE(1,n) = 0.d0;              CurlE(1,n) =-P1(i)*E2(j)*D3(k)
        ShapE(2,n) = P1(i)*E2(j)*P3(k); CurlE(2,n) = 0.d0
        ShapE(3,n) = 0.d0;              CurlE(3,n) = D1(i)*E2(j)*P3(k)
      enddo
      enddo
      enddo
c     direction 3
      do k=1,nord3
      do j=3,nord2+1
      do i=3,nord1+1
        n = n+1
        ShapE(1,n) = 0.d0;              CurlE(1,n) = P1(i)*D2(j)*E3(k)
        ShapE(2,n) = 0.d0;              CurlE(2,n) =-D1(i)*P2(j)*E3(k)
        ShapE(3,n) = P1(i)*P2(j)*E3(k); CurlE(3,n) = 0.d0
      enddo
      enddo
      enddo
c
      NrdofE = n
c
c

      CONTAINS

c     ----------------------------------------------------
      subroutine decode_orderq(Norderq, Nord1,Nord2)
c
c     RECALL: Norderq = Nord1*10+Nord2
c
      Nord1 = Norderq/10
      Nord2 = Norderq-Nord1*10
c
      if (Nord1<1.or.Nord1>9) then
        write(*,*)'decode_orderq: OUT OF RANGE Nord1 = ',Nord1
        stop
      endif
      if (Nord2<1.or.Nord2>9) then
        write(*,*)'decode_orderq: OUT OF RANGE Nord2 = ',Nord2
        stop
      endif
c
      end subroutine decode_orderq
c     ----------------------------------------------------

      end subroutine shape3bEM2


