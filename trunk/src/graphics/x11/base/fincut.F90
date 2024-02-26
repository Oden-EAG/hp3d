#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - fincut
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine constructs small triangles resulting
!                        from cutting a plane triangle with a cutting
!                        plane
!
!   arguments :
!     in :
!              Xtro    - triangle vertex coordinates
!              Usol    - values of the cutting plane equation at
!                        vertices
!     out:
!              Nstrl   - number of small triangles
!              Costr   - small triangles vertex coordinates
!              Isflag  = 0   -no changes,
!                        1   -new triangles
!
!---------------------------------------------------------------------
      subroutine fincut(Xtro,Usol, Nstrl,Costr,Isflag)
!
      implicit none
!
      integer :: Nstrl,Isflag
      real(8) :: Xtro(3,3),Usol(3),Costr(2,3,3)
!
      real(8) :: slev,U1,U2,U3
      integer :: icart
!
      U1=Usol(1)
      U2=Usol(2)
      U3=Usol(3)
      slev=0.d0
!
!  ...perform multistage checking
      Nstrl=0
      if (U1.gt.slev) then
        if (U2.gt.slev) then
          if (U3.gt.slev) then
            Nstrl=0
            Isflag=0
          else
            Nstrl=1
            do icart=1,3
              Costr(1,icart,1) = Xtro(icart,1) +   &
              (Xtro(icart,3)-Xtro(icart,1)) * (slev-U1)/(U3-U1)
              Costr(1,icart,2) = Xtro(icart,2) +   &
              (Xtro(icart,3)-Xtro(icart,2)) * (slev-U2)/(U3-U2)
              Costr(1,icart,3) = Xtro(icart,3)
            enddo
            Isflag=1
          endif
        else
          if (U3.gt.slev) then
            Nstrl=1
            do icart=1,3
              Costr(1,icart,2) = Xtro(icart,1) +   &
              (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
              Costr(1,icart,3) = Xtro(icart,2)
              Costr(1,icart,1) = Xtro(icart,2) +   &
              (Xtro(icart,3)-Xtro(icart,2)) * (slev-U2)/(U3-U2)
            enddo
            Isflag=1
          else
            Nstrl=2
            do icart=1,3
              Costr(1,icart,1) = Xtro(icart,1) +   &
              (Xtro(icart,3)-Xtro(icart,1)) * (slev-U1)/(U3-U1)
              Costr(1,icart,2) = Xtro(icart,1) +   &
              (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
              Costr(1,icart,3) = Xtro(icart,3)
              Costr(2,icart,1) = Xtro(icart,1) +   &
              (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
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
              Costr(1,icart,1) = Xtro(icart,1) +   &
              (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
              Costr(1,icart,2) = Xtro(icart,1) +   &
              (Xtro(icart,3)-Xtro(icart,1)) * (slev-U1)/(U3-U1)
            enddo
            Isflag=1
          else
            Nstrl=2
            do icart=1,3
              Costr(2,icart,1) = Xtro(icart,1)
              Costr(2,icart,2) = Xtro(icart,1) +   &
              (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
              Costr(2,icart,3) = Xtro(icart,3)
              Costr(1,icart,1) = Xtro(icart,1) +   &
              (Xtro(icart,2)-Xtro(icart,1)) * (slev-U1)/(U2-U1)
              Costr(1,icart,2) = Xtro(icart,2) +   &
              (Xtro(icart,3)-Xtro(icart,2)) * (slev-U2)/(U3-U2)
              Costr(1,icart,3) = Xtro(icart,3)
            enddo
            Isflag=1
          endif
        else
          if (U3.gt.slev) then
            Nstrl=2
            do icart=1,3
              Costr(1,icart,3) = Xtro(icart,1)
              Costr(1,icart,1) = Xtro(icart,2) +   &
              (Xtro(icart,3)-Xtro(icart,2)) * (slev-U2)/(U3-U2)
              Costr(1,icart,2) = Xtro(icart,1) +   &
              (Xtro(icart,3)-Xtro(icart,1)) * (slev-U1)/(U3-U1)
              Costr(2,icart,1) = Xtro(icart,1)
              Costr(2,icart,2) = Xtro(icart,2)
              Costr(2,icart,3) = Xtro(icart,2) +   &
              (Xtro(icart,3)-Xtro(icart,2)) * (slev-U2)/(U3-U2)
            enddo
            Isflag=1
          else
            Nstrl=1
            Isflag=0
          endif
        endif
      endif
!
!
      end subroutine fincut

#endif
