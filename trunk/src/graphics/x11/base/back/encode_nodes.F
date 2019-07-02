c----------------------------------------------------------------------
c
c   routine name       - encode_nodes
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jul 98
c
c   purpose            - encode node numbers for a side of a 3D element
c
c   usage              - call encode_nodes(Type,Isid,Nodesl, Nsidnod)
c
c   arguments :
c     in:
c              Type    - element type
c              Isid    - side number
c              Nodesl  - element nodes
c     out:
c              Nsidnod - the side nodes numbers
c
c   required  routines - none
c
c----------------------------------------------------------------------
c
      subroutine encode_nodes(Type,Isid,Nodesl, Nsidnod)
c
c
      character*4 Type
      dimension Nodesl(27),Nsidnod(9)
c
      if (Type.eq.'pris') then
        if (Isid.eq.1) then
          Nsidnod(1) = Nodesl(1)
          Nsidnod(2) = Nodesl(2)
          Nsidnod(3) = Nodesl(3)
        elseif (Isid.eq.2) then
          Nsidnod(1) = Nodesl(4)
          Nsidnod(2) = Nodesl(5)
          Nsidnod(3) = Nodesl(6)
        elseif (Isid.eq.3) then
          Nsidnod(1) = Nodesl(1)
          Nsidnod(2) = Nodesl(2)
          Nsidnod(3) = Nodesl(5)
          Nsidnod(4) = Nodesl(4)
        elseif (Isid.eq.4) then
          Nsidnod(1) = Nodesl(2)
          Nsidnod(2) = Nodesl(3)
          Nsidnod(3) = Nodesl(6)
          Nsidnod(4) = Nodesl(5)
        elseif (Isid.eq.5) then
          Nsidnod(1) = Nodesl(3)
          Nsidnod(2) = Nodesl(1)
          Nsidnod(3) = Nodesl(4)
          Nsidnod(4) = Nodesl(6)
        else
          write(*,*) 'encode_nodes: Isid = ',Isid
          stop 1
        endif
      elseif (Type.eq.'bric') then
        if (Isid.eq.1) then
          Nsidnod(1) = Nodesl(1)
          Nsidnod(2) = Nodesl(4)
          Nsidnod(3) = Nodesl(3)
          Nsidnod(4) = Nodesl(2)
          Nsidnod(5) = Nodesl(12)
          Nsidnod(6) = Nodesl(11)
          Nsidnod(7) = Nodesl(10)
          Nsidnod(8) = Nodesl(9)
          Nsidnod(9) = Nodesl(21)
        elseif (Isid.eq.2) then
          Nsidnod(1) = Nodesl(5)
          Nsidnod(2) = Nodesl(6)
          Nsidnod(3) = Nodesl(7)
          Nsidnod(4) = Nodesl(8)
          Nsidnod(5) = Nodesl(13)
          Nsidnod(6) = Nodesl(14)
          Nsidnod(7) = Nodesl(15)
          Nsidnod(8) = Nodesl(16)
          Nsidnod(9) = Nodesl(22)
        elseif (Isid.eq.3) then
          Nsidnod(1) = Nodesl(1)
          Nsidnod(2) = Nodesl(2)
          Nsidnod(3) = Nodesl(6)
          Nsidnod(4) = Nodesl(5)
          Nsidnod(5) = Nodesl(9)
          Nsidnod(6) = Nodesl(18)
          Nsidnod(7) = Nodesl(13)
          Nsidnod(8) = Nodesl(17)
          Nsidnod(9) = Nodesl(23)
        elseif (Isid.eq.4) then
          Nsidnod(1) = Nodesl(2)
          Nsidnod(2) = Nodesl(3)
          Nsidnod(3) = Nodesl(7)
          Nsidnod(4) = Nodesl(6)
          Nsidnod(5) = Nodesl(10)
          Nsidnod(6) = Nodesl(19)
          Nsidnod(7) = Nodesl(14)
          Nsidnod(8) = Nodesl(18)
          Nsidnod(9) = Nodesl(24)
        elseif (Isid.eq.5) then
          Nsidnod(1) = Nodesl(3)
          Nsidnod(2) = Nodesl(4)
          Nsidnod(3) = Nodesl(8)
          Nsidnod(4) = Nodesl(7)
          Nsidnod(5) = Nodesl(11)
          Nsidnod(6) = Nodesl(20)
          Nsidnod(7) = Nodesl(15)
          Nsidnod(8) = Nodesl(19)
          Nsidnod(9) = Nodesl(25)
        elseif (Isid.eq.6) then
          Nsidnod(1) = Nodesl(4)
          Nsidnod(2) = Nodesl(1)
          Nsidnod(3) = Nodesl(5)
          Nsidnod(4) = Nodesl(8)
          Nsidnod(5) = Nodesl(12)
          Nsidnod(6) = Nodesl(17)
          Nsidnod(7) = Nodesl(16)
          Nsidnod(8) = Nodesl(20)
          Nsidnod(9) = Nodesl(26)
        else
          write(*,*) 'encode_nodes: Isid = ',Isid
          stop 1
        endif
      else
        write(*,*) 'encode_nodes: Type = ',Type
        stop 1
      endif
c
c
      return
      end




