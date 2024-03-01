#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - int_coord
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine determines integer coordinates of
!                        vertices for a subtriangle of a 2D element,
!                        the corresponding color flag, and flags
!                        indicating adjacency to element sides
!
!   arguments :
!     in:
!             Numlev   - a flag = 0 for mesh display,
!                               > 0 for contour maps
!             Nsub     - (even,.ge.4) number of subdivisions
!             Ntype    - type of element
!             Norder   - order of approximation (4 edges + 1 middle)
!             I,J      - integer indices of the field
!     out:
!             Ixi1     - integer coordinates of the lower triangle
!             Ixi2     - integer coordinates of the upper triangle
!             Np1      - color code for the lower triangle
!             Np2      - color code for the upper triangle
!             Nsid1    - adjacency info for the lower triangle
!             Nsid2    - adjacency info for the upper triangle
!
!----------------------------------------------------------------------
!
   subroutine int_coord(Numlev,Nsub,Ntype,Norder,I,J, &
                           Ixi1,Ixi2,Np1,Np2,Nsid1,Nsid2)
!
      use node_types
!
      implicit none
!
      integer, intent(in)  :: Numlev,Nsub,Ntype
      integer, intent(in)  :: Norder(5)
      integer, intent(in)  :: I,J
      integer, intent(out) :: Ixi1(2,3),Ixi2(2,3)
      integer, intent(out) :: Np1,Np2,Nsid1,Nsid2
!
      integer :: nhalf,nordh,nordv
!
#if HP3D_DEBUG
      integer :: k
      integer :: iprint
      iprint=0
#endif
!
!  ...initialize
      Ixi1 = 0; Ixi2 = 0
!
!  ...check input
      nhalf = Nsub/2
      if ((Nsub.lt.4).or.(Nsub.ne.2*nhalf)) then
        write(*,*) 'int_coord:  WRONG INPUT !! Nsub = ',Nsub
        stop 1
      endif
      select case(Ntype)
      case(TRIA)
        if ((I.le.0).or.(I.gt.Nsub).or.(J.le.0).or.(J.gt.Nsub).or.   &
            (I+J.gt.nsub+1)) then
          write(*,*)  'int_coord:  WRONG INPUT (T) !! I,J = ', I,J
          stop 1
        endif
      case(RECT,QUAD)
        if ((I.le.0).or.(I.gt.Nsub).or.(J.le.0).or.(J.gt.Nsub)) then
          write(*,*)  'int_coord:  WRONG INPUT (Q) !! I,J = ', I,J
          stop 1
        endif
      case default
        write(*,*) 'int_coord:  WRONG INPUT !! Type = ',S_Type(Ntype)
        stop 1
      end select
!
      Nsid1=0
      Nsid2=0
      select case(Ntype)
      case(TRIA,MDLT)
        if (I.eq.1) then
          if (J.eq.1) then
            if (Numlev.eq.0) then
              Ixi1(1,1)=0   ; Ixi1(1,2)=1   ; Ixi1(1,3)=1
              Ixi1(2,1)=0   ; Ixi1(2,2)=0   ; Ixi1(2,3)=1
              Np1=Norder(1)
              Nsid1=1*4+0*2+0
              Ixi2(1,1)=0   ; Ixi2(1,2)=1   ; Ixi2(1,3)=0
              Ixi2(2,1)=0   ; Ixi2(2,2)=1   ; Ixi2(2,3)=1
              Np2=Norder(3)
              Nsid2=0*4+0*2+1
            else
              Ixi1(1,1)=0   ; Ixi1(1,2)=1   ; Ixi1(1,3)=0
              Ixi1(2,1)=0   ; Ixi1(2,2)=0   ; Ixi1(2,3)=1
              Np1=Norder(4)
              Nsid1=1*4+0*2+1
              Ixi2(1,1)=0   ; Ixi2(1,2)=1   ; Ixi2(1,3)=1
              Ixi2(2,1)=1   ; Ixi2(2,2)=0   ; Ixi2(2,3)=1
              Np2=Norder(4)
              Nsid2=0*4+0*2+0
            endif
          elseif (J.lt.Nsub-1) then
            Ixi1(1,1)=0   ; Ixi1(1,2)=1   ; Ixi1(1,3)=0
            Ixi1(2,1)=J-1 ; Ixi1(2,2)=J-1 ; Ixi1(2,3)=J
            Np1=Norder(3)
            Nsid1=0*4+0*2+1
            Ixi2(1,1)=0   ; Ixi2(1,2)=1   ; Ixi2(1,3)=1
            Ixi2(2,1)=J   ; Ixi2(2,2)=J-1 ; Ixi2(2,3)=J
            Np2=Norder(3)
          elseif (J.eq.Nsub-1) then
            if (Numlev.eq.0) then
              Ixi1(1,1)=0   ; Ixi1(1,2)=1   ; Ixi1(1,3)=0
              Ixi1(2,1)=J-1 ; Ixi1(2,2)=J-1 ; Ixi1(2,3)=J
              Np1=Norder(3)
              Nsid1=0*4+0*2+1
              Ixi2=0
              Np2=0
            else
              Ixi1(1,1)=0   ; Ixi1(1,2)=1   ; Ixi1(1,3)=0
              Ixi1(2,1)=J-1 ; Ixi1(2,2)=J-1 ; Ixi1(2,3)=J
              Np1=Norder(3)
              Nsid1=0*4+0*2+1
              Ixi2(1,1)=0   ; Ixi2(1,2)=1   ; Ixi2(1,3)=1
              Ixi2(2,1)=J   ; Ixi2(2,2)=J-1 ; Ixi2(2,3)=J
              Np2=Norder(4)
              Nsid2=0
            endif
          elseif (J.eq.Nsub) then
            if (Numlev.eq.0) then
              Ixi1(1,1)=0   ; Ixi1(1,2)=1   ; Ixi1(1,3)=0
              Ixi1(2,1)=J-1 ; Ixi1(2,2)=J-2 ; Ixi1(2,3)=J
              Np1=Norder(3)
              Nsid1=0*4+0*2+1
              Ixi2(1,1)=0   ; Ixi2(1,2)=1   ; Ixi2(1,3)=1
              Ixi2(2,1)=J   ; Ixi2(2,2)=J-2 ; Ixi2(2,3)=J-1
              Np2=Norder(2)
              Nsid2=0*4+0*2+1
            else
              Ixi1(1,1)=0   ; Ixi1(1,2)=1   ; Ixi1(1,3)=0
              Ixi1(2,1)=J-1 ; Ixi1(2,2)=J-1 ; Ixi1(2,3)=J
              Np1=Norder(4)
              Nsid1=0*4+1*2+1
              Ixi2 = 0
              Np2=0
            endif
          endif
        elseif (I.lt.Nsub-1) then
          if (J.eq.1) then
            Ixi1(1,1)=I-1 ; Ixi1(1,2)=I   ; Ixi1(1,3)=I-1
            Ixi1(2,1)=0   ; Ixi1(2,2)=0   ; Ixi1(2,3)=1
            Np1=Norder(1)
            Nsid1=1*4+0*2+0
            Ixi2(1,1)=I-1 ; Ixi2(1,2)=I   ; Ixi2(1,3)=I
            Ixi2(2,1)=J   ; Ixi2(2,2)=J-1 ; Ixi2(2,3)=J
            Np2=Norder(1)
          elseif (J.lt.Nsub-I) then
            Ixi1(1,1)=I-1 ; Ixi1(1,2)=I   ; Ixi1(1,3)=I-1
            Ixi1(2,1)=J-1 ; Ixi1(2,2)=J-1 ; Ixi1(2,3)=J
            Np1=Norder(4)
            Ixi2(1,1)=I-1 ; Ixi2(1,2)=I   ; Ixi2(1,3)=I
            Ixi2(2,1)=J   ; Ixi2(2,2)=J-1 ; Ixi2(2,3)=J
            Np2=Norder(4)
          elseif (J.eq.Nsub-I) then
            Ixi1(1,1)=I-1 ; Ixi1(1,2)=I   ; Ixi1(1,3)=I-1
            Ixi1(2,1)=J-1 ; Ixi1(2,2)=J-1 ; Ixi1(2,3)=J
            Np1=Norder(4)
            Ixi2(1,1)=I-1 ; Ixi2(1,2)=I   ; Ixi2(1,3)=I
            Ixi2(2,1)=J   ; Ixi2(2,2)=J-1 ; Ixi2(2,3)=J
            Np2=Norder(2)
          elseif (J.eq.Nsub-I+1) then
            Ixi1(1,1)=I-1 ; Ixi1(1,2)=I   ; Ixi1(1,3)=I-1
            Ixi1(2,1)=J-1 ; Ixi1(2,2)=J-1 ; Ixi1(2,3)=J
            Np1=Norder(2)
            Nsid1=0*4+1*2+0
            Ixi2=0
            Np2=0
          endif
        elseif (I.eq.Nsub-1) then
          if (J.eq.1) then
            if (Numlev.eq.0) then
              Ixi1(1,1)=I-1 ; Ixi1(1,2)=I   ; Ixi1(1,3)=I-1
              Ixi1(2,1)=0   ; Ixi1(2,2)=0   ; Ixi1(2,3)=1
              Np1=Norder(1)
              Nsid1=1*4+0*2+0
              Ixi2=0
              Np2=0
            else
              Ixi1(1,1)=I-1 ; Ixi1(1,2)=I   ; Ixi1(1,3)=I-1
              Ixi1(2,1)=0   ; Ixi1(2,2)=0   ; Ixi1(2,3)=1
              Np1=Norder(1)
              Nsid1=1*4+0*2+0
              Ixi2(1,1)=I-1 ; Ixi2(1,2)=I   ; Ixi2(1,3)=I
              Ixi2(2,1)=J   ; Ixi2(2,2)=J-1 ; Ixi2(2,3)=J
              Np2=Norder(4)
            endif
          elseif (J.eq.2) then
            Ixi1(1,1)=I-1 ; Ixi1(1,2)=I   ; Ixi1(1,3)=I-1
            Ixi1(2,1)=1   ; Ixi1(2,2)=1   ; Ixi1(2,3)=2
            Np1=Norder(2)
            Nsid1=0*4+1*2+0
            Ixi2=0
            Np2=0
          endif
        elseif (I.eq.Nsub) then
          if (Numlev.eq.0) then
            Ixi1(1,1)=I-1 ; Ixi1(1,2)=I   ; Ixi1(1,3)=I-2
            Ixi1(2,1)=0   ; Ixi1(2,2)=0   ; Ixi1(2,3)=1
            Np1=Norder(1)
            Nsid1=1*4+0*2+0
            Ixi2(1,1)=I-2 ; Ixi2(1,2)=I   ; Ixi2(1,3)=I-1
            Ixi2(2,1)=1   ; Ixi2(2,2)=0   ; Ixi2(2,3)=1
            Np2=Norder(2)
            Nsid2=0*4+1*2+0
          else
            Ixi1(1,1)=I-1 ; Ixi1(1,2)=I   ; Ixi1(1,3)=I-1
            Ixi1(2,1)=0   ; Ixi1(2,2)=0   ; Ixi1(2,3)=1
            Np1=Norder(4)
            Nsid1=1*4+1*2+0
            Ixi2=0
            Np2=0
          endif
        endif
!
      case(QUAD,RECT,MDLQ)
!
!   ....determine the coordinates
        if ( (I.le.nhalf .and. J.le.nhalf)   &
        .or. (I.gt.nhalf .and. J.gt.nhalf)) then
          Ixi1(1,1)=I-1;   Ixi1(1,2)=I  ;   Ixi1(1,3)=I
          Ixi1(2,1)=J-1;   Ixi1(2,2)=J-1;   Ixi1(2,3)=J
!
          Ixi2(1,1)=I-1;   Ixi2(1,2)=I  ;   Ixi2(1,3)=I-1
          Ixi2(2,1)=J-1;   Ixi2(2,2)=J  ;   Ixi2(2,3)=J
        else
          Ixi1(1,1)=I-1;   Ixi1(1,2)=I  ;   Ixi1(1,3)=I-1
          Ixi1(2,1)=J-1;   Ixi1(2,2)=J-1;   Ixi1(2,3)=J
!
          Ixi2(1,1)=I-1;   Ixi2(1,2)=I  ;   Ixi2(1,3)=I
          Ixi2(2,1)=J  ;   Ixi2(2,2)=J-1;   Ixi2(2,3)=J
        endif
        call decode(Norder(5), nordh,nordv)
!
!  .....determine the color flags
        if (I.eq.1) then
          if (J.eq.1) then
            Np1=Norder(1)
            Np2=Norder(4)
          elseif (J.lt.Nsub) then
            Np1=Norder(4)
            Np2=Norder(4)
          elseif (J.eq.Nsub) then
            Np1=Norder(4)
            Np2=Norder(3)
          endif
        elseif (I.le.nhalf) then
          if (J.eq.1) then
            Np1=Norder(1)
            Np2=Norder(1)
          elseif (J.lt.I) then
            Np1=nordh
            Np2=nordh
          elseif (J.eq.I) then
            Np1=nordh
            Np2=nordv
          elseif (J.lt.Nsub-I+1) then
            Np1=nordv
            Np2=nordv
          elseif (J.eq.Nsub-I+1) then
            Np1=nordv
            Np2=nordh
          elseif (J.lt.Nsub) then
            Np1=nordh
            Np2=nordh
          elseif (J.eq.Nsub) then
            Np1=Norder(3)
            Np2=Norder(3)
          endif
        elseif (I.lt.Nsub) then
          if (J.eq.1) then
            Np1=Norder(1)
            Np2=Norder(1)
          elseif (J.lt.Nsub-I+1) then
            Np1=nordh
            Np2=nordh
          elseif (J.eq.Nsub-I+1) then
            Np1=nordh
            Np2=nordv
          elseif (J.lt.I) then
            Np1=nordv
            Np2=nordv
          elseif (J.eq.I) then
            Np1=nordv
            Np2=nordh
          elseif (J.lt.Nsub) then
            Np1=nordh
            Np2=nordh
          elseif (J.eq.Nsub) then
            Np1=Norder(3)
            Np2=Norder(3)
          endif
        elseif (I.eq.Nsub) then
          if (J.eq.1) then
            Np1=Norder(1)
            Np2=Norder(2)
          elseif (J.lt.Nsub) then
            Np1=Norder(2)
            Np2=Norder(2)
          elseif (J.eq.Nsub) then
            Np1=Norder(2)
            Np2=Norder(3)
          endif
        endif
!
!  .....determine the adjacency flags
        if (J.eq.1) then
          Nsid1 = 1*4+0*2+0
        endif
        if (I.eq.Nsub) then
          if (J.le.nhalf) then
            Nsid2 = 0*4+1*2+0
          else
            Nsid1 = 0*4+1*2+0
          endif
        endif
        if (J.eq.Nsub) then
          if (I.le.nhalf) then
            Nsid2 = 0*4+0*2+1
          else
            Nsid2 = 0*4+1*2+0
          endif
        endif
        if (I.eq.1) then
          if (J.le.nhalf) then
            Nsid2 = 0*4+0*2+1
          else
            Nsid1 = 0*4+0*2+1
          endif
        endif
!
      case default
        write(*,*) 'int_coord: Type = ',S_Type(Ntype)
        stop 1
      end select
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) I,J
 7001   format('int_coord: I,J = ',2i3)
        write(*,7002) (Ixi1(1:2,k),k=1,3),Nsid1
 7002   format('Ixi1,Nsid1 = ',3(2i3,2x),3x,i3)
        write(*,7003) (Ixi2(1:2,k),k=1,3),Nsid2
 7003   format('Ixi2,Nsid2 = ',3(2i3,2x),3x,i3)
        call pause
      endif
#endif
!
   end subroutine int_coord

#endif
