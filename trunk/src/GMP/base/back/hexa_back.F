c----------------------------------------------------------------------
c
c   routine name       - hexa
c
c---------------------------------------------------------------------
c
c   latest revision    - Jan 06
c
c   purpose            - routine defines the parameterization
c                        for a hexahedron based on the trivariate
c                        interpolation (with linear blending
c                        functions) of the paramaterization
c                        of its faces and edges
c
c   arguments :
c     in:
c               No     - hexahedron number
c               Eta    - reference coordinates of a point
c                        in the reference hexahedron
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the reference coordinates
c
c---------------------------------------------------------------------
c
      subroutine hexa(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
c
      dimension Eta(3),X(3),Dxdeta(3,3)
c
c
c  ...face figures numbers, the corresponding orientation
      dimension nof(6),norientf(6)
c
c  ...edge curves numbers, the corresponding orientation
      dimension noe(12),noriente(12)
c
c  ...points numbers
      dimension nop(8)
c
c  ...local edge curves  and face rectangles coordinates and their
c     derivatives wrt eta
      dimension detac(3), etas(2),detas(2,3)
c
c  ...values and derivatives of blending functions
      dimension val(26),dval(3,26)
c
c  ...value of point edge or face parametrizations with derivatives
      dimension xyz(3),dedg(3),dsid(3,2)
c
c
c----------------------------------------------------------------------
c
      iprint=0
c
      select case(HEXAS(No)%Type)
      case('Linear','TIhexa')
      case default
        write(*,7000) HEXAS(No)%Type
 7000   format('hexa: TYPE ',a10,' NOT SUPPORTED')
        stop 1
      end select
c
      if (iprint.eq.1) then
        write(*,7001) No,Eta
 7001   format('hexa: No,Eta = ',i4,2x,3f8.3)
      endif
c
c  ...decode the face rectangles numbers and their orientation
      do is=1,6
        call decode(iabs(HEXAS(No)%FigNo(is)), nof(is),norientf(is))
      enddo
c
      if (iprint.eq.1) then
        write(*,*) 'hexa: FACE NUMBERS AND ORIENTATIONS'
        write(*,8001) (nof(i),i=1,6)
 8001   format(12(i8,2x))
        write(*,8001) (norientf(i),i=1,6)
      endif
c
c------------------------------------------------------------------------
c
c  ...get the numbers of edge curves
      do iedg=1,12
        noe(iedg) = iabs(HEXAS(No)%EdgeNo(iedg))
c
c  .....the edge curve number and orientation
        select case(noe(iedg)/HEXAS(No)%EdgeNo(iedg))
        case(-1)
          noriente(iedg) = 1
        case(1)
          noriente(iedg) = 0
        case default
          write(*,*) 'hexa: ERROR 1'
          stop 1
        end select
      enddo
      if (iprint.eq.1) then
        write(*,*) 'hexa: EDGE NUMBERS AND ORIENTATIONS'
        write(*,8001) (noe(i),i=1,12)
        write(*,8001) (noriente(i),i=1,12)
      endif
c
c------------------------------------------------------------------------
c
c  ...get the vertex point numbers
      do iv=1,8
        nop(iv) = HEXAS(No)%VertNo(iv)
      enddo
      if (iprint.eq.1) then
        write(*,*) 'hexa: VERTEX POINTS NUMBER'
        write(*,8001) (nop(i),i=1,8)
      endif
c
c------------------------------------------------------------------------
c
c  ...compute the blending functions with derivatives
      call hexa_blend(Eta, val,dval)
c
c  ...initiate the value of the parametrization
      X(1:3) = 0.d0
      Dxdeta(1:3,1:3) = 0.d0
      k=0
c
c  ...trilinear interpolant of the vertex values
      do i=1,8
        k=k+1
        call pointr(nop(i), xyz)
        X(1:3) = X(1:3) + xyz(1:3)*val(k)
        do j=1,3
          Dxdeta(1:3,j) = Dxdeta(1:3,j) + xyz(1:3)*dval(j,k)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'hexa: AFTER VERTEX INTERPOLATION'
        do ivar=1,3
          write(*,7002) ivar,X(ivar),(Dxdeta(ivar,ivar1),ivar1=1,3)
 7002     format('hexa: ivar,X,Dxdeta = ',i2,2x,f8.3,2x,3f8.3)
        enddo
      endif
      if (HEXAS(No)%Type.eq.'Linear') return
c
c------------------------------------------------------------------------
c
c  ...bilinear interpolant of the edge values
      do i=1,12
        k=k+1
        detac(1:3) = 0.d0
        select case(i)
        case(1,3,5,7)
          etac = Eta(1); detac(1) = 1.d0
        case(2,4,6,8)
          etac = Eta(2); detac(2) = 1.d0
        case(9,10,11,12)
          etac = Eta(3); detac(3) = 1.d0
        end select
        call curve_local(noe(i),noriente(i),etac, xyz,dedg)
        X(1:3) = X(1:3) - xyz(1:3)*val(k)
        do j=1,3
          Dxdeta(1:3,j) = Dxdeta(1:3,j)
     .                  - dedg(1:3)*detac(j)*val(k)
     .                  - xyz(1:3)*dval(j,k)
        enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,*)'hexa: AFTER EDGE INTERPOLATION'
        do ivar=1,3
          write(*,7002) ivar,X(ivar),(Dxdeta(ivar,ivar1),ivar1=1,3)
        enddo
      endif
c
c------------------------------------------------------------------------
c
c  ...linear interpolant of the face values
      do i=1,6
        k=k+1
        detas(1:2,1:3) = 0.d0
        select case(i)
        case(1,2)
          etas(1) = Eta(1); etas(2) = Eta(2)
          detas(1,1) = 1.d0; detas(2,2) = 1.d0
        case(3,5)
          etas(1) = Eta(1); etas(2) = Eta(3)
          detas(1,1) = 1.d0; detas(2,3) = 1.d0
        case(4,6)
          etas(1) = Eta(2); etas(2) = Eta(3)
          detas(1,2) = 1.d0; detas(2,3) = 1.d0
        end select
        call recta_local(nof(i),norientf(i),etas, xyz, dsid)
        X(1:3) = X(1:3) + xyz(1:3)*val(k)
        do j=1,3
          Dxdeta(1:3,j) = Dxdeta(1:3,j)
     .         + (dsid(1:3,1)*detas(1,j)+dsid(1:3,2)*detas(2,j))
     .           *val(k) + xyz(1:3)*dval(j,k)
        enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,*)'hexa: AFTER FACE INTERPOLATION'
        do ivar=1,3
          write(*,7002) ivar,X(ivar),(Dxdeta(ivar,ivar1),ivar1=1,3)
        enddo
        call pause
      endif
c
c
      return
      end
c
c-----------------------------------------------------------------------
c
c   routine name       - hexa_blend
c
c-----------------------------------------------------------------------
c
c
c   computer           - machine independent
c
c   latest revision    - Jan 06
c
c   purpose            - routine evaluates triliner, bilinear and linear
c                        blending functions for the transfinite
c                        interpolation on master hexahedron and their
c                        derivatives
c
c   arguments :
c     in:
c          Eta        - reference coordinates of a point in the
c                       reference hexahedron
c     out:
c          Val        - values of the functions
c          Dval       - derivatives
c
c   required  routines -
c
c-----------------------------------------------------------------------
c
      subroutine hexa_blend(Eta, Val,Dval)
c
#include "syscom.blk"
c
      dimension Eta(3),Val(26),Dval(3,26)
c
c-----------------------------------------------------------------------
c
c  ...vertex trilinear blending functions...
c
c  ...first vertex trilinear funtion
      k=1
      Val(k) = (1.d0-Eta(1))*(1.d0-Eta(2))*(1.d0-Eta(3))
      Dval(1,k) = - (1.d0-Eta(2))*(1.d0-Eta(3))
      Dval(2,k) = - (1.d0-Eta(1))*(1.d0-Eta(3))
      Dval(3,k) = - (1.d0-Eta(1))*(1.d0-Eta(2))
c
c  ...second vertex trilinear funtion
      k=2
      Val(k) = Eta(1)*(1.d0-Eta(2))*(1.d0-Eta(3))
      Dval(1,k) =   (1.d0-Eta(2))*(1.d0-Eta(3))
      Dval(2,k) = - Eta(1)*(1.d0-Eta(3))
      Dval(3,k) = - Eta(1)*(1.d0-Eta(2))
c
c  ...third vertex trilinear funtion
      k=3
      Val(k) = Eta(1)*Eta(2)*(1.d0-Eta(3))
      Dval(1,k) =   Eta(2)*(1.d0-Eta(3))
      Dval(2,k) =   Eta(1)*(1.d0-Eta(3))
      Dval(3,k) = - Eta(1)*Eta(2)
c
c  ...fourth vertex trilinear funtion
      k=4
      Val(k) = (1.d0-Eta(1))* Eta(2)*(1.d0-Eta(3))
      Dval(1,k) = - Eta(2)*(1.d0-Eta(3))
      Dval(2,k) =   (1.d0-Eta(1))*(1.d0-Eta(3))
      Dval(3,k) = - (1.d0-Eta(1))*Eta(2)
c
c  ...fifth vertex trilinear funtion
      k=5
      Val(k) = (1.d0-Eta(1))*(1.d0-Eta(2))*Eta(3)
      Dval(1,k) = - (1.d0-Eta(2))*Eta(3)
      Dval(2,k) = - (1.d0-Eta(1))*Eta(3)
      Dval(3,k) =   (1.d0-Eta(1))*(1.d0-Eta(2))
c
c  ...sixth vertex trilinear funtion
      k=6
      Val(k) = Eta(1)*(1.d0-Eta(2))*Eta(3)
      Dval(1,k) =   (1.d0-Eta(2))*Eta(3)
      Dval(2,k) =  - Eta(1)*Eta(3)
      Dval(3,k) =    Eta(1)*(1.d0-Eta(2))
c
c  ...seventh vertex trilinear funtion
      k=7
      Val(k) = Eta(1)*Eta(2)*Eta(3)
      Dval(1,k) =   Eta(2)*Eta(3)
      Dval(2,k) =   Eta(1)*Eta(3)
      Dval(3,k) =   Eta(1)*Eta(2)
c
c  ...eighth vertex trilinear funtion
      k=8
      Val(k) = (1.d0-Eta(1))*Eta(2)*Eta(3)
      Dval(1,k) = - Eta(2)* Eta(3)
      Dval(2,k) =   (1.d0-Eta(1))*Eta(3)
      Dval(3,k) =   (1.d0-Eta(1))*Eta(2)
c
c-----------------------------------------------------------------------
c
c  ...edge bilinear blending functions..
c
c  ...first edge bilinear funtion
      k=9
      Val(k) = (1.d0-Eta(2))*(1.d0-Eta(3))
      Dval(1,k) =    0.d0
      Dval(2,k) = - (1.d0-Eta(3))
      Dval(3,k) = - (1.d0-Eta(2))
c
c  ...second  edge bilinear funtion
      k=10
      Val(k) = Eta(1)*(1.d0-Eta(3))
      Dval(1,k) =   (1.d0-Eta(3))
      Dval(2,k) =    0.d0
      Dval(3,k) = - Eta(1)
c
c  ...third edge bilinear funtion
      k=11
      Val(k) = Eta(2)*(1.d0-Eta(3))
      Dval(1,k) =   0.d0
      Dval(2,k) =  (1.d0-Eta(3))
      Dval(3,k) = - Eta(2)
c
c  ...fourth edge bilinear funtion
      k=12
      Val(k) = (1.d0-Eta(1))*(1.d0-Eta(3))
      Dval(1,k) = - (1.d0-Eta(3))
      Dval(2,k) =   0.d0
      Dval(3,k) = - (1.d0-Eta(1))
c
c  ...fifth edge bilinear funtion
      k=13
      Val(k) = (1.d0-Eta(2))*Eta(3)
      Dval(1,k) =   0.d0
      Dval(2,k) = - Eta(3)
      Dval(3,k) =   (1.d0-Eta(2))
c
c  ...sixth edge bilinear funtion
      k=14
      Val(k) = Eta(1)*Eta(3)
      Dval(1,k) =   Eta(3)
      Dval(2,k) =   0.d0
      Dval(3,k) =   Eta(1)
c
c  ...seventh edge bilinear funtion
      k=15
      Val(k) = Eta(2)*Eta(3)
      Dval(1,k) =   0.d0
      Dval(2,k) =   Eta(3)
      Dval(3,k) =   Eta(2)
c
c  ...eigth edge bilinear funtion
      k=16
      Val(k) = (1.d0-Eta(1))*Eta(3)
      Dval(1,k) = - Eta(3)
      Dval(2,k) =   0.d0
      Dval(3,k) =   (1.d0-Eta(1))
c
c  ...nineth edge bilinear funtion
      k=17
      Val(k) = (1.d0-Eta(1))*(1.d0-Eta(2))
      Dval(1,k) = - (1.d0-Eta(2))
      Dval(2,k) = - (1.d0-Eta(1))
      Dval(3,k) =   0.d0
c
c  ...tenth  edge bilinear funtion
      k=18
      Val(k) = Eta(1)*(1.d0-Eta(2))
      Dval(1,k) =   (1.d0-Eta(2))
      Dval(2,k) = - Eta(1)
      Dval(3,k) =   0.d0
c
c  ...eleventh edge bilinear funtion
      k=19
      Val(k) = Eta(1)* Eta(2)
      Dval(1,k) =   Eta(2)
      Dval(2,k) =   Eta(1)
      Dval(3,k) =   0.d0
c
c  ...twelve edge bilinear funtion
      k=20
      Val(k) = (1.d0-Eta(1))*Eta(2)
      Dval(1,k) = - Eta(2)
      Dval(2,k) =   (1.d0-Eta(1))
      Dval(3,k) =   0.d0
c
c-----------------------------------------------------------------------
c
c  ...face linear blending functions...
c
c
c  ...first face linear function
      k=21
      Val(k) = 1.d0 - Eta(3)
      Dval(1,k) =   0.d0
      Dval(2,k) =   0.d0
      Dval(3,k) = - 1.d0
c
c  ...second face linear function
      k=22
      Val(k) = Eta(3)
      Dval(1,k) =   0.d0
      Dval(2,k) =   0.d0
      Dval(3,k) =   1.d0
c
c  ...third face linear function
      k=23
      Val(k) = 1.d0 - Eta(2)
      Dval(1,k) =   0.d0
      Dval(2,k) = - 1.d0
      Dval(3,k) =   0.d0
c
c  ...fourth face linear function
      k=24
      Val(k) = Eta(1)
      Dval(1,k) =   1.d0
      Dval(2,k) =   0.d0
      Dval(3,k) =   0.d0
c
c  ...fifth face linear function
      k=25
      Val(k) = Eta(2)
      Dval(1,k) =   0.d0
      Dval(2,k) =   1.d0
      Dval(3,k) =   0.d0
c
c  ...sixth face linear function
      k=26
      Val(k) = 1.d0 - Eta(1)
      Dval(1,k) = - 1.d0
      Dval(2,k) =   0.d0
      Dval(3,k) =   0.d0
c
c
      end





