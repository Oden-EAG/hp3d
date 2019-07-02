c----------------------------------------------------------------------
c
c   routine name       - prism_TI
c
c---------------------------------------------------------------------
c
c
c   computer           - machine independent
c
c   latest revision    - Dec 07
c
c   purpose            - routine defines the parameterization
c                        for a prism based on a transfinite
c                        interpolation (with linear and singular
c                        blending functions) of the paramaterization
c                        of its faces
c
c   arguments :
c     in:
c               No     - prism number
c               Eta    - reference coordinates of a point
c                        in the reference prism
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the reference coordinates
c
c   required  routines -
c
c---------------------------------------------------------------------
c
      subroutine prism_TI_OLD_OLD(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
c
      dimension Eta(3),X(3),Dxdeta(3,3)
c
c
c  ...face figures numbers, the corresponding orientation
      dimension nof(5),norientf(5)
c
c  ...edge curves numbers, the corresponding orientation
      dimension noe(9),noriente(9)
c
c  ...points numbers
      dimension nop(6)
c
c  ...values and derivatives of blending functions
      dimension val(29), dval(3,29)
c
c  ...vertex coordinates
      dimension xvert(3,6)
c
c  ...vertical edges parametrizations
      dimension vedg(3,3),dvedg(3,3)
c
c  ...value of edge or face parametrizations with derivatives
      dimension xyz(3),dedg(3),dsid(3,2)
c
c  ...derivatives of the local edge variable wrt reference coordinates
      dimension detac(3)
c
c  ...local face coordinates with derivatives
      dimension etas(2),detas(2,3)
c
c  ...affine coordinates
      dimension rlambda(3),drlambda(3,2)
c
      data drlambda /-1.d0,-1.d0, 1.d0,0.d0, 0.d0,1.d0/
c
c
      imod(j,mod) = j-(j-1)/mod*mod
c
c----------------------------------------------------------------------
c
      iprint=0
c
c
      if (iprint.eq.1) then
        write(*,7001) No,Eta
 7001   format('prism_TI: No,Eta = ',i4,2x,3f8.3)
      endif
c
c  ...decode the face figures numbers and their orientation
      do is=1,5
        call decode(iabs(PRISMS(No)%FigNo(is)), nof(is),norientf(is))
      enddo
c
      if (iprint.eq.1) then
        write(*,*) 'prism: FACE NUMBERS AND ORIENTATIONS'
        write(*,8001) (nof(i),i=1,6)
 8001   format(12(i8,2x))
        write(*,8001) (norientf(i),i=1,6)
      endif
c
c------------------------------------------------------------------------
c
c  ...get the numbers of edge curves
      do iedg=1,9
        noe(iedg) = iabs(PRISMS(No)%EdgeNo(iedg))
c
c  .....the edge curve number and orientation
        select case(noe(iedg)/PRISMS(No)%EdgeNo(iedg))
        case(-1)
          noriente(iedg) = 1
        case(1)
          noriente(iedg) = 0
        case default
          write(*,*) 'prism: ERROR 1'
          stop 1
        end select
      enddo
      if (iprint.eq.1) then
        write(*,*) 'prism: EDGE NUMBERS AND ORIENTATIONS'
        write(*,8001) (noe(i),i=1,12)
        write(*,8001) (noriente(i),i=1,12)
      endif
c
c------------------------------------------------------------------------
c
c  ...get the vertex point numbers
      do iv=1,6
        nop(iv) = PRISMS(No)%VertNo(iv)
      enddo
      if (iprint.eq.1) then
        write(*,*) 'prism: VERTEX POINTS NUMBER'
        write(*,8001) (nop(i),i=1,8)
      endif
c
c------------------------------------------------------------------------
c
      rlambda(1) = 1-Eta(1)-Eta(2)
      rlambda(2) = Eta(1)
      rlambda(3) = Eta(2)
c
c  ...initiate the value of the parametrization
      X(1:3) = 0.d0
      Dxdeta(1:3,1:3) = 0.d0
c
c  ...compute the blending functions and their derivatives
      call prism_blend(Eta, val,dval)
c
c  ...bilinear interpolant of the vertex values
      do i=1,6
        call pointr(nop(i), xvert(1:3,i))
        X(1:3) = X(1:3) - xvert(1:3,i)*val(i)
        do j=1,3
          Dxdeta(1:3,j) = Dxdeta(1:3,j) - xvert(1:3,i)*dval(j,i)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'prism: AFTER VERTEX INTERPOLATION'
        do ivar=1,3
          write(*,7002) ivar,X(ivar),(Dxdeta(ivar,ivar1),ivar1=1,3)
 7002     format('hexa: ivar,X,Dxdeta = ',i2,2x,f8.3,2x,3f8.3)
        enddo
      endif
c
c------------------------------------------------------------------------
c
c  ...horizontal edges contributions
      k=6
      do iz=1,2
      do i=1,3
        i1 = imod(i+1,3)
        do j=1,2
          k=k+1
c
c  .......find the local edge coordinate
          select case(j)
          case(1)
            etac = rlambda(i1)
            detac(1:2) = drlambda(i1,1:2)
          case(2)
            etac = 1.d0 - rlambda(i)
            detac(1:2) = -drlambda(i,1:2)
          end select
          detac(3) = 0.d0
c
          ie = (iz-1)*2+i
          call curve_local(noe(ie),noriente(ie),etac, xyz,dedg)
c
c  .......compute the egde bubbles and their derivatives
          xyz(1:3) = xyz(1:3)
     .             - (xvert(1:3,i)*(1.d0-etac)+ xvert(1:3,i1)*etac)
          dedg(1:3) = dedg(1:3)
     .              - (xvert(1:3,i1) - xvert(1:3,i))
c
c  .......upgrade the parametrization
          do ivar=1,3
            X(ivar) = X(ivar) - .5d0*xyz(ivar)*val(k)
            Dxdeta(ivar,1:3) = Dxdeta(ivar,1:3)
     .                       - .5d0*( dedg(ivar)*detac(1:3)*val(k)
     .                              + xyz(ivar)*dval(1:3,k))
          enddo
c
c  .....end of loop through the two edge contributions
        enddo
      enddo
      enddo
c
c  ...vertical edges contributions
      do i=1,3
        ie = 6+i
        k=k+1
        etac = Eta(3)
        detac(1:2) = 0.d0
        detac(3) = 1.d0
        call curve_local(noe(ie),noriente(ie),etac,
     .                   vedg(1:3,i),dvedg(1:3,i))
c
c  .....upgrade the parametrization
        do ivar=1,3
          X(ivar) = X(ivar) + vedg(ivar,i)*val(k)
          Dxdeta(ivar,1:3) = Dxdeta(ivar,1:3)
     .                     + dvedg(ivar,i)*detac(1:3)*val(k)
     .                     + vedg(ivar,i)*dval(1:3,k)
        enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,*)'prism: AFTER EDGE INTERPOLATION'
        do ivar=1,3
          write(*,7002) ivar,X(ivar),(Dxdeta(ivar,ivar1),ivar1=1,3)
        enddo
      endif
c
c------------------------------------------------------------------------
c
c  ...bottom and top face
      do if=1,2
        k=k+1
        etas(1:2) = Eta(1:2)
        detas(1:2,1:3) = 0.d0; detas(1,1) = 1.d0; detas(2,2) = 1.d0
        call trian_local(nof(if),norientf(if),etas, xyz,dsid)
c
c  .....upgrade the parametrization
        do ivar=1,3
          X(ivar) = X(ivar) + xyz(ivar)*val(k)
          Dxdeta(ivar,1:3) = Dxdeta(ivar,1:3)
     .    + (dsid(ivar,1)*detas(1,1:3)+dsid(ivar,2)*detas(2,1:3))*val(k)
     .    + xyz(ivar)*dval(1:3,k)
        enddo
      enddo
c
c  ...lateral faces
      do i=1,3
        i1 = imod(i+1,3)
        if=2+i
c
c  .....two contributions per face
        do j=1,2
          k=k+1
c
c  .......find the local face coordinates
          detas(1:2,1:3) = 0.d0
          select case(j)
          case(1)
            etas(1) = rlambda(i1)
            detas(1,1:2) = drlambda(i1,1:2)
          case(2)
            etas(1) = 1.d0 - rlambda(i)
            detas(1,1:2) = -drlambda(i,1:2)
          end select
          detas(2,3) = 1.d0
c
          call recta_local(nof(if),norientf(if),etas, xyz,dsid)
c
c  .......compute the face bubbles and their derivatives
          xyz(1:3) = xyz(1:3)
     .             - (vedg(1:3,i)*(1.d0-etas(1))+ vedg(1:3,i1)*etas(1))
          dsid(1:3,1) = dsid(1:3,1)
     .              - (xvert(1:3,i1) - xvert(1:3,i))
c
c  .......upgrade the parametrization
          do ivar=1,3
            X(ivar) = X(ivar) + xyz(ivar)*val(k)
            Dxdeta(ivar,1:3) = Dxdeta(ivar,1:3)
     .      + ( dsid(ivar,1)*detas(1,1:3)
     .         +dsid(ivar,2)*detas(2,1:3))*val(k)
     .      + xyz(ivar)*dval(1:3,k)
          enddo
        enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,*)'prism: AFTER FACE INTERPOLATION'
        do ivar=1,3
          write(*,7002) ivar,X(ivar),(Dxdeta(ivar,ivar1),ivar1=1,3)
        enddo
        call pause
      endif
c
c
      end
c
c----------------------------------------------------------------------
c
c   routine name       - prism_blend
c
c---------------------------------------------------------------------
c
c
c   computer           - machine independent
c
c   latest revision    - Jan 06
c
c   purpose            - routine evaluates blending functions with
c                        their derivatives for the transfinite
c                        interpolation on master prism; functions
c                        corresponding to horizontal edges and faces
c                        are singular, the remaining functions
c                        are bilinear; there are TWO singular
c                        blending functions for each singular
c                        edge or lateral face
c
c   arguments :
c     in:
c          Eta        - reference coordinates of a point in the
c                       reference hexahedron
c     out:
c          Val        - value of the function
c          Dval       - derivatives
c
c   required  routines -
c
c---------------------------------------------------------------------
c
      subroutine prism_blend(Eta, Val,Dval)
c
#include "syscom.blk"
c
      dimension Eta(3),Val(29),Dval(3,29)
c
c  ...affine coordinates and their derivatives wrt etas(x,y)
      dimension rlambda(3),drlambda(2,3)
c
      data eps/1.d-13/
      data one/1.d0/
      data drlambda /-1.d0,-1.d0, 1.d0,0.d0, 0.d0,1.d0/
c
      imod(j,mod) = j-(j-1)/mod*mod
c
      rz = Eta(3); rzp = one-rz
      rlambda(1) = one-Eta(1)-Eta(2)
      rlambda(2) = Eta(1)
      rlambda(3) = Eta(2)
c
c  ...check input values
      if ((rlambda(1).lt.0.d0).or.(rlambda(2).lt.0.d0).or.
     .   (rlambda(3).lt.0.d0)) then
        write(*,7001) Eta(1:3)
 7001   format('prism_blend: WRONG INPUT, Eta = ',3e)
        stop 1
      endif
c
      k=0
c
c  ...bottom vertex  blending functions
      do i=1,3
        k=k+1
        Val(k) = rlambda(i)*rzp
        Dval(1:2,k) = drlambda(1:2,i)*rzp
        Dval(3,k) = -rlambda(i)
      enddo
c
c  ...top vertex blending functions
      do i=1,3
        k=k+1
        Val(k) = rlambda(i)*rz
        Dval(1:2,k) = drlambda(1:2,i)*rz
        Dval(3,k) =  rlambda(i)
      enddo
c
c  ...lateral faces singular functions
      k=23
      do i=1,3
        i1 = imod(i+1,3)
c
c  .....first singular function
        k=k+1
        if (one - rlambda(i1).lt.eps) then
          Val(k) = 0.d0; Dval(1:3,k) = 0.d0
        else
          Val(k) = rlambda(i)/(one - rlambda(i1))
          Dval(1:2,k) = drlambda(1:2,i)/(one - rlambda(i1))
     .    - rlambda(i)/(one - rlambda(i1))**2*drlambda(1:2,i1)
          Dval(3,k) =  0.d0
        endif
c
c  .....second singular function
        k=k+1
        if (one-rlambda(i).lt.eps) then
          Val(k) = 0.d0; Dval(1:3,k) = 0.d0
        else
          Val(k) = rlambda(i1)/(one - rlambda(i))
          Dval(1:2,k) = drlambda(1:2,i1)/(one - rlambda(i))
     .    - rlambda(i1)/(one - rlambda(i))**2*drlambda(1:2,i)
          Dval(3,k) =  0.d0
        endif
      enddo
c
c  ...bottom horizontal edges singular functions
      k=6
      do j=1,6
        k=k+1
        Val(k) = Val(17+k)*rzp
        Dval(1:2,k) = Dval(1:2,17+k)*rzp
        Dval(3,k) =  -Val(17+k)
      enddo
c
c  ...top horizontal edges singular functions
      k=12
      do j=1,6
        k=k+1
        Val(k) = Val(11+k)*rz
        Dval(1:2,k) = Dval(1:2,11+k)*rz
        Dval(3,k) = Val(11+k)
      enddo
c
c  ...vertical edges shape functions
      k=18
      do i=1,3
        k=k+1
        Val(k) = rlambda(i)
        Dval(1:2,k) = drlambda(1:2,i)
        Dval(3,k) = 0.d0
      enddo
c
c  ...horizontal faces shape functions
      k=k+1
      Val(k) = rzp
      Dval(1:2,k) = 0.d0
      Dval(3,k) = -one
      k=k+1
      Val(k) = rz
      Dval(1:2,k) = 0.d0
      Dval(3,k) = one
c
      end


