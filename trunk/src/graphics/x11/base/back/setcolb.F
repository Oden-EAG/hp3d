c----------------------------------------------------------------------
c
c   routine name       - setcolb
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Dec 05
c
c   purpose            - routine sets the color for subtriangle
c                        of a face of a brick element
c
c   arguments :
c     in:
c          Norder      - order of an element
c          Nside       - face number
c          Nsub1,Nsub2 - subdivision numbers
c          Icase       - number of subtriangles family
c          I,J         - subtriangle local coordinates
c     out:
c          Ncol        - color number
c
c   required  routines - decode
c
c----------------------------------------------------------------------
c
      subroutine setcolb(Norder,Nside,Nsub1,Nsub2,Icase,I,J, Ncol)
c
      use graphmod
      dimension Norder(19)
c
c
      imod(j,mod) = j-(j-1)/mod*mod
c
      iprint=0
      if (iprint.eq.1) then
        write(*,7001) Norder
 7001   format('setcolb: Norder = ',12i2,6i3,i4)
      endif
c
c  ...determine face orders in a local face system of coordinates
c     (consult routine encode_nodes for the definition of the face
c      systems of coordinates)
      select case(Nside)
      case(1)
        nord1 = Norder(4)
        nord2 = Norder(3)
        nord3 = Norder(2)
        nord4 = Norder(1)
        call decode(Norder(13), nordv,nordh)
      case(2)
        nord1 = Norder(5)
        nord2 = Norder(6)
        nord3 = Norder(7)
        nord4 = Norder(8)
        call decode(Norder(14), nordh,nordv)
      case(3)
        nord1 = Norder(1)
        nord2 = Norder(10)
        nord3 = Norder(5)
        nord4 = Norder(9)
        call decode(Norder(15), nordh,nordv)
      case(4)
        nord1 = Norder(2)
        nord2 = Norder(11)
        nord3 = Norder(6)
        nord4 = Norder(10)
        call decode(Norder(16), nordh,nordv)
      case(5)
        nord1 = Norder(3)
        nord2 = Norder(12)
        nord3 = Norder(7)
        nord4 = Norder(11)
        call decode(Norder(17), nordh,nordv)
      case(6)
        nord1 = Norder(4)
        nord2 = Norder(9)
        nord3 = Norder(8)
        nord4 = Norder(12)
        call decode(Norder(18), nordh,nordv)
      end select
c
c  ...determine which face quadrant we are on and set the
c     qudrant colors using a local quadrant system of coordinates
      nhalf1 = Nsub1/2
      if (2*nhalf1.ne.Nsub1) then
        write(*,7002) Nsub1
 7002   format('setcolb: WRONG Nsub1 = ',i3)
        stop 1
      endif
      nhalf2 = Nsub2/2
      if (2*nhalf2.ne.Nsub2) then
        write(*,7003) Nsub2
 7003   format('setcolb: WRONG Nsub2 = ',i3)
        stop 1
      endif
c
c  ...first quadrant
      if ((I.le.nhalf1).and.(J.le.nhalf2)) then
        ii = I; jj = J; ll = Icase
        nc1 = nord1; nc2 = nord4; nc3 = nordh; nc4 = nordv
c
c  ...second quadrant
      elseif ((I.gt.nhalf1).and.(J.le.nhalf2)) then
        ii = J; jj = Nsub1+1-I; ll = imod(Icase+1,2)
        nc1 = nord2; nc2 = nord1; nc3 = nordv; nc4 = nordh
c
c  ...third quadrant
      elseif ((I.gt.nhalf1).and.(J.gt.nhalf2)) then
        ii = Nsub1+1-I; jj = Nsub2+1-J; ll = imod(Icase+1,2)
        nc1 = nord3; nc2 = nord2; nc3 = nordh; nc4 = nordv
c
c  ...fourth quadrant
      elseif ((I.le.nhalf1).and.(J.gt.nhalf2)) then
        ii = Nsub2+1-J; jj = I; ll = Icase
        nc1 = nord4; nc2 = nord3; nc3 = nordv; nc4 = nordh
c
      else
        write(*,7004) I,J,nhalf1,nhalf2
 7004   format('setcolb: I,J,nhalf1,nhalf2 = ',4i4)
        stop 1
      endif
c
c  ...assign color
      if ((jj.eq.1).and.((ii.gt.1).or.(ll.eq.1))) then
        Ncol = nc1
      elseif ((ii.eq.1).and.((jj.gt.1).or.(ll.eq.2))) then
        Ncol = nc2
      elseif ((ii.gt.jj).or.((ii.eq.jj).and.(ll.eq.1))) then
        Ncol = nc3
      elseif ((ii.lt.jj).or.((ii.eq.jj).and.(ll.eq.2))) then
        Ncol = nc4
      else
        write(*,7005) ii,jj,ll,Nsub1,Nsub2
 7005   format('setcolb: ii,jj,ll,Nsub1,Nsub2 = ',5i3)
        stop 1
      endif
c
c
      return
      end
