#if HP3D_USE_X11

c----------------------------------------------------------------------
c
c   routine name       - finstr
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine constructs small triangles
c                        where the solution is greater then Slev
c                        within the bigger triangle with vertices
c                        coordinates Xtro and solution values Usol
c
c   arguments :
c     in:
c              Xtro      - bigger triangle vertices coordinates
c              Usol      - values of the solution at vertices
c              Slev      - limiting level of solution
c              Slevup    - next value on the scale
c     out:
c              Nstrl     - number of small triangles
c              Costr     - small triangles vertices coordinates
c              Nsid_flag - a flag indicating which edge should be
c                          be drawn (B/W version)
c
c----------------------------------------------------------------------
c
      subroutine finstr(Xtro,Usol,Slev,Slevup, Nstrl,Costr,Nsid_flag)
c
      implicit none
c
      integer :: Nstrl
      integer :: Nsid_flag(2)
      real(8) :: Xtro(3,3),Usol(3),Costr(2,3,3)
      real(8) :: Slev,Slevup
c
      integer :: icart
c
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
      if (iprint.eq.1) then
         write(*,*) 'in finstr'
         write(*,*) 'Usol',Usol
         write(*,*) 'Slev,Slevup',Slev,Slevup
         call pause
      endif
#endif
c
      Nsid_flag=0
c
c  ...perform multistage checking
      Nstrl=0
      if (Usol(1).le.Slev) then
        if (Usol(2).le.Slev) then
          if (Usol(3).le.Slev) then
            Nstrl=0
            return
          else
            Nstrl=1
            Nsid_flag(1)=-1
            do 10 icart=1,3
              Costr(1,icart,1) = Xtro(icart,1) +
     . (Xtro(icart,3)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(3)-Usol(1))
              Costr(1,icart,2) = Xtro(icart,2) +
     . (Xtro(icart,3)-Xtro(icart,2)) * (Slev-Usol(2))/(Usol(3)-Usol(2))
              Costr(1,icart,3) = Xtro(icart,3)
   10       continue
            return
          endif
        else
          if (Usol(3).le.Slev) then
            Nstrl=1
            Nsid_flag(1)=-3
            do 20 icart=1,3
              Costr(1,icart,1) = Xtro(icart,1) +
     . (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(1,icart,2) = Xtro(icart,2)
              Costr(1,icart,3) = Xtro(icart,2) +
     . (Xtro(icart,3)-Xtro(icart,2)) * (Slev-Usol(2))/(Usol(3)-Usol(2))
   20       continue
            return
          else
            Nstrl=2
            Nsid_flag(1)=-1
            do 30 icart=1,3
              Costr(1,icart,1) = Xtro(icart,1) +
     . (Xtro(icart,3)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(3)-Usol(1))
              Costr(1,icart,2) = Xtro(icart,1) +
     . (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(1,icart,3) = Xtro(icart,3)
              Costr(2,icart,1) = Xtro(icart,1) +
     . (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(2,icart,2) = Xtro(icart,2)
              Costr(2,icart,3) = Xtro(icart,3)
   30       continue
            return
          endif
        endif
      else
        if (Usol(2).le.Slev) then
          if (Usol(3).le.Slev) then
            Nstrl=1
            Nsid_flag(1)=-2
            do 40 icart=1,3
              Costr(1,icart,1) = Xtro(icart,1)
              Costr(1,icart,2) = Xtro(icart,1) +
     . (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(1,icart,3) = Xtro(icart,1) +
     . (Xtro(icart,3)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(3)-Usol(1))
   40       continue
            return
          else
            Nstrl=2
            Nsid_flag(2)=-1
            do 50 icart=1,3
              Costr(1,icart,1) = Xtro(icart,1)
              Costr(1,icart,2) = Xtro(icart,1) +
     . (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(1,icart,3) = Xtro(icart,3)
              Costr(2,icart,1) = Xtro(icart,1) +
     . (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(2,icart,2) = Xtro(icart,2) +
     . (Xtro(icart,3)-Xtro(icart,2)) * (Slev-Usol(2))/(Usol(3)-Usol(2))
              Costr(2,icart,3) = Xtro(icart,3)
   50       continue
            return
          endif
        else
          if (Usol(3).le.Slev) then
            Nstrl=2
            Nsid_flag(1)=-2
            do 60 icart=1,3
              Costr(1,icart,1) = Xtro(icart,1)
              Costr(1,icart,2) = Xtro(icart,2) +
     . (Xtro(icart,3)-Xtro(icart,2)) * (Slev-Usol(2))/(Usol(3)-Usol(2))
              Costr(1,icart,3) = Xtro(icart,1) +
     . (Xtro(icart,3)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(3)-Usol(1))
              Costr(2,icart,1) = Xtro(icart,1)
              Costr(2,icart,2) = Xtro(icart,2)
              Costr(2,icart,3) = Xtro(icart,2) +
     . (Xtro(icart,3)-Xtro(icart,2)) * (Slev-Usol(2))/(Usol(3)-Usol(2))
   60       continue
            return
          else
c            iup=1
c            do 70 iver=1,3
c              if (Usol(iver).lt.Slevup) iup=0
c   70       continue
c            if (iup.eq.1) then
c              Nstrl=0
c              return
c            else
              Nstrl=1
              do 80 icart=1,3
                Costr(1,icart,1) = Xtro(icart,1)
                Costr(1,icart,2) = Xtro(icart,2)
                Costr(1,icart,3) = Xtro(icart,3)
   80         continue
              return
c            endif
          endif
        endif
      endif
c
      end subroutine finstr

#endif
