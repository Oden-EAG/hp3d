c----------------------------------------------------------------------
c
c   routine name       - fincut
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine constructs small triangles resulting
c                        from cutting a plane triangle with a cutting
c                        plane
c
c   arguments :
c     in :
c              Xtro    - triangle vertex coordinates
c              Usol    - values of the cutting plane equation at
c                        vertices
c     out:
c              Nstrl   - number of small triangles
c              Costr   - small triangles vertex coordinates
c              Isflag  = 0   -no changes,
c                        1   -new triangles
c
c---------------------------------------------------------------------
      subroutine fincut(Xtro,Usol, Nstrl,Costr,Isflag)
c
      implicit none
c
      integer :: Nstrl,Isflag
      real(8) :: Xtro(3,3),Usol(3),Costr(2,3,3)
c
      real(8) :: slev,U1,U2,U3
      integer :: icart
c
      U1=Usol(1)
      U2=Usol(2)
      U3=Usol(3)
      slev=0.d0
c
c  ...perform multistage checking
      Nstrl=0
      if (U1.gt.slev) then
        if (U2.gt.slev) then
          if (U3.gt.slev) then
            Nstrl=0
            Isflag=0
          else
            Nstrl=1
            do icart=1,3
              Costr(1,icart,1) = Xtro(icart,1) +
     .        (Xtro(icart,3)-Xtro(icart,1)) * (slev-U1)/(U3-U1)
              Costr(1,icart,2) = Xtro(icart,2) +
     .        (Xtro(icart,3)-Xtro(icart,2)) * (slev-U2)/(U3-U2)
              Costr(1,icart,3) = Xtro(icart,3)
            enddo
            Isflag=1
          endif
        else
          if (U3.gt.slev) then
            Nstrl=1
            do icart=1,3
              Costr(1,icart,2) = Xtro(icart,1) +
     .        (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
              Costr(1,icart,3) = Xtro(icart,2)
              Costr(1,icart,1) = Xtro(icart,2) +
     .        (Xtro(icart,3)-Xtro(icart,2)) * (slev-U2)/(U3-U2)
            enddo
            Isflag=1
          else
            Nstrl=2
            do icart=1,3
              Costr(1,icart,1) = Xtro(icart,1) +
     .        (Xtro(icart,3)-Xtro(icart,1)) * (slev-U1)/(U3-U1)
              Costr(1,icart,2) = Xtro(icart,1) +
     .        (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
              Costr(1,icart,3) = Xtro(icart,3)
              Costr(2,icart,1) = Xtro(icart,1) +
     .        (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
              Costr(2,icart,2) = Xtro(icart,2)
              Costr(2,icart,3) = Xtro(icart,3)
            enddo
            Isflag=1
          endif
        endif
      else
        if (U2.gt.slev) then
          if (U3.gt.slev) then
            Nstrl=1
            do icart=1,3
              Costr(1,icart,3) = Xtro(icart,1)
              Costr(1,icart,1) = Xtro(icart,1) +
     .        (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
              Costr(1,icart,2) = Xtro(icart,1) +
     .        (Xtro(icart,3)-Xtro(icart,1)) * (slev-U1)/(U3-U1)
            enddo
            Isflag=1
          else
            Nstrl=2
            do icart=1,3
              Costr(2,icart,1) = Xtro(icart,1)
              Costr(2,icart,2) = Xtro(icart,1) +
     .        (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
              Costr(2,icart,3) = Xtro(icart,3)
              Costr(1,icart,1) = Xtro(icart,1) +
     .        (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
              Costr(1,icart,2) = Xtro(icart,2) +
     .        (Xtro(icart,3)-Xtro(icart,2)) * (slev-U2)/(U3-U2)
              Costr(1,icart,3) = Xtro(icart,3)
            enddo
            Isflag=1
          endif
        else
          if (U3.gt.slev) then
            Nstrl=2
            do icart=1,3
              Costr(1,icart,3) = Xtro(icart,1)
              Costr(1,icart,1) = Xtro(icart,2) +
     .        (Xtro(icart,3)-Xtro(icart,2)) * (slev-U2)/(U3-U2)
              Costr(1,icart,2) = Xtro(icart,1) +
     .        (Xtro(icart,3)-Xtro(icart,1)) * (slev-U1)/(U3-U1)
              Costr(2,icart,1) = Xtro(icart,1)
              Costr(2,icart,2) = Xtro(icart,2)
              Costr(2,icart,3) = Xtro(icart,2) +
     .        (Xtro(icart,3)-Xtro(icart,2)) * (slev-U2)/(U3-U2)
            enddo
            Isflag=1
          else
            Nstrl=1
            Isflag=0
          endif
        endif
      endif
c
c
      end subroutine fincut
