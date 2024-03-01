#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - finstr
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine constructs small triangles
!                        where the solution is greater then Slev
!                        within the bigger triangle with vertices
!                        coordinates Xtro and solution values Usol
!
!   arguments :
!     in:
!              Xtro      - bigger triangle vertices coordinates
!              Usol      - values of the solution at vertices
!              Slev      - limiting level of solution
!              Slevup    - next value on the scale
!     out:
!              Nstrl     - number of small triangles
!              Costr     - small triangles vertices coordinates
!              Nsid_flag - a flag indicating which edge should be
!                          be drawn (B/W version)
!
!----------------------------------------------------------------------
!
   subroutine finstr(Xtro,Usol,Slev,Slevup, Nstrl,Costr,Nsid_flag)
!
      implicit none
!
      integer :: Nstrl
      integer :: Nsid_flag(2)
      real(8) :: Xtro(3,3),Usol(3),Costr(2,3,3)
      real(8) :: Slev,Slevup
!
      integer :: icart
!
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
!
      Nsid_flag=0
!
!  ...perform multistage checking
      Nstrl=0
      if (Usol(1).le.Slev) then
        if (Usol(2).le.Slev) then
          if (Usol(3).le.Slev) then
            Nstrl=0
            return
          else
            Nstrl=1
            Nsid_flag(1)=-1
            do icart=1,3
              Costr(1,icart,1) = Xtro(icart,1) + &
       (Xtro(icart,3)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(3)-Usol(1))
              Costr(1,icart,2) = Xtro(icart,2) + &
       (Xtro(icart,3)-Xtro(icart,2)) * (Slev-Usol(2))/(Usol(3)-Usol(2))
              Costr(1,icart,3) = Xtro(icart,3)
            enddo
            return
          endif
        else
          if (Usol(3).le.Slev) then
            Nstrl=1
            Nsid_flag(1)=-3
            do icart=1,3
              Costr(1,icart,1) = Xtro(icart,1) + &
       (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(1,icart,2) = Xtro(icart,2)
              Costr(1,icart,3) = Xtro(icart,2) + &
       (Xtro(icart,3)-Xtro(icart,2)) * (Slev-Usol(2))/(Usol(3)-Usol(2))
            enddo
            return
          else
            Nstrl=2
            Nsid_flag(1)=-1
            do icart=1,3
              Costr(1,icart,1) = Xtro(icart,1) + &
       (Xtro(icart,3)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(3)-Usol(1))
              Costr(1,icart,2) = Xtro(icart,1) + &
       (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(1,icart,3) = Xtro(icart,3)
              Costr(2,icart,1) = Xtro(icart,1) + &
       (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(2,icart,2) = Xtro(icart,2)
              Costr(2,icart,3) = Xtro(icart,3)
            enddo
            return
          endif
        endif
      else
        if (Usol(2).le.Slev) then
          if (Usol(3).le.Slev) then
            Nstrl=1
            Nsid_flag(1)=-2
            do icart=1,3
              Costr(1,icart,1) = Xtro(icart,1)
              Costr(1,icart,2) = Xtro(icart,1) + &
       (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(1,icart,3) = Xtro(icart,1) + &
       (Xtro(icart,3)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(3)-Usol(1))
            enddo
            return
          else
            Nstrl=2
            Nsid_flag(2)=-1
            do icart=1,3
              Costr(1,icart,1) = Xtro(icart,1)
              Costr(1,icart,2) = Xtro(icart,1) + &
       (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(1,icart,3) = Xtro(icart,3)
              Costr(2,icart,1) = Xtro(icart,1) + &
       (Xtro(icart,2)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(2)-Usol(1))
              Costr(2,icart,2) = Xtro(icart,2) + &
       (Xtro(icart,3)-Xtro(icart,2)) * (Slev-Usol(2))/(Usol(3)-Usol(2))
              Costr(2,icart,3) = Xtro(icart,3)
            enddo
            return
          endif
        else
          if (Usol(3).le.Slev) then
            Nstrl=2
            Nsid_flag(1)=-2
            do icart=1,3
              Costr(1,icart,1) = Xtro(icart,1)
              Costr(1,icart,2) = Xtro(icart,2) + &
       (Xtro(icart,3)-Xtro(icart,2)) * (Slev-Usol(2))/(Usol(3)-Usol(2))
              Costr(1,icart,3) = Xtro(icart,1) + &
       (Xtro(icart,3)-Xtro(icart,1)) * (Slev-Usol(1))/(Usol(3)-Usol(1))
              Costr(2,icart,1) = Xtro(icart,1)
              Costr(2,icart,2) = Xtro(icart,2)
              Costr(2,icart,3) = Xtro(icart,2) + &
       (Xtro(icart,3)-Xtro(icart,2)) * (Slev-Usol(2))/(Usol(3)-Usol(2))
            enddo
            return
          else
!            iup=1
!            do iver=1,3
!              if (Usol(iver).lt.Slevup) iup=0
!            enddo
!            if (iup.eq.1) then
!              Nstrl=0
!              return
!            else
              Nstrl=1
              do icart=1,3
                Costr(1,icart,1) = Xtro(icart,1)
                Costr(1,icart,2) = Xtro(icart,2)
                Costr(1,icart,3) = Xtro(icart,3)
              enddo
              return
!            endif
          endif
        endif
      endif
!
   end subroutine finstr

#endif
