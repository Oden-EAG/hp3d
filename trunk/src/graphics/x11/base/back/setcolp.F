      subroutine setcolp(Norder,Nside,Nsub1,Nsub2,Icase,I,J,Ncol)
c
c----------------------------------------------------------------------
c
c   routine name       - setcolp
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - jan, 92
c
c   purpose            - routine sets the colour for subtriangle
c                        of a prism element
c
c   usage              - call setcolp(Norder,Nside,Nsub1,
c                                    Nsub2,Icase,I,J,Ncol)
c
c   arguments :
c     in:  Norder      - order of an element
c          Nside       - side number
c          Nsub1,Nsub2 - subdivision numbers
c          Icase       - number of subtriangles family
c          I,J         - subtriangle local coordinates
c     out: Ncol        - color number
c
c   required  routines - decode
c
c----------------------------------------------------------------------
c
      use graphmod
      dimension Norder(15),nord(15), nordh(15)
c . ..only nordh(12-14) is meanningfull
c
      iprint = 0
c
      do 5 k=1,15
        call decode(Norder(k), nordh(k),nord(k))
    5 continue
      do 10 k=1,15
        nord(k)  = max0(nord(k),1)
        nordh(k) = max0(nordh(k),1)
   10 continue
      minp1=2
      minp2=2
      do 20 k=1,3
        minp1=min0(minp1,nord(k))
        minp2=min0(minp2,nord(k+3))
   20 continue
      nord(10) = max0(minp1,nord(10))
      nord(11) = max0(minp2,nord(11))
c
c  ...prepering for sides
      if( Nside>2 )then
        nhalf = Nsub1/2
        if( 2*nhalf/=Nsub2 )then
           print*,'SETCOLB: Nsub1 /= Nsub2',Nsub1, '/=', Nsub2
           stop
        endif

c    ...first quarter
        if( I<=nhalf .and. J<=nhalf )then
           jcase = 1; ii = I; jj = J
           norh = nordh(9+Nside);  norv = nord(9+Nside)

c    ...third quarter
        elseif( I>nhalf .and. J>nhalf )then
           jcase = 1; ii = I-nhalf; jj = J-nhalf
           norv = nordh(9+Nside);  norh = nord(9+Nside)

c    ...fourth quarter
        elseif( I<=nhalf .and. J>nhalf )then
           jcase = 2; ii = I; jj = J-nhalf
           norh = nord(9+Nside);  norv = nordh(9+Nside)

c    ...second quarter
        elseif( I>nhalf .and. J<=nhalf )then
           jcase = 2; ii = I-nhalf; jj = J
           norv = nord(9+Nside);  norh = nordh(9+Nside)
        else
           print*,'SETCOLB: I,J,Nsub1 = ',I,J,Nsub1
           stop
        endif
      endif

c
      if (Nside.eq.1) then
        if (Icase.eq.1) then
          if ((J.eq.1).and.(I.gt.1).and.(I.lt.Nsub1)) then
            Ncol = (nord(3))
          elseif ((J.eq.Nsub2).and.(I.gt.1).and.(I.lt.Nsub1)) then
            Ncol = (nord(2))
          elseif ((I.eq.1).and.(J.gt.1).and.(J.lt.Nsub2)) then
            Ncol = (nord(1))
          else
            Ncol = (nord(10))
          endif
        else
          if ((J.eq.1).and.(I.gt.1).and.(I.lt.Nsub1-1)) then
            Ncol = (nord(3))
          elseif ((J.eq.Nsub2-1).and.(I.gt.1).and.(I.lt.Nsub1-1)) then
            Ncol = (nord(2))
          elseif ((I.eq.1).and.(J.gt.1).and.(J.lt.Nsub2-1)) then
            Ncol = (nord(1))
          else
            Ncol = (nord(10))
          endif
        endif
      endif
c
      if (Nside.eq.2) then
        if (Icase.eq.1) then
          if ((J.eq.1).and.(I.gt.1).and.(I.lt.Nsub1)) then
            Ncol = (nord(4))
          elseif ((J.eq.Nsub2).and.(I.gt.1).and.(I.lt.Nsub1)) then
            Ncol = (nord(5))
          elseif ((I.eq.1).and.(J.gt.1).and.(J.lt.Nsub2)) then
            Ncol = (nord(6))
          else
            Ncol = (nord(11))
          endif
        else
          if ((J.eq.1).and.(I.gt.1).and.(I.lt.Nsub1-1)) then
            Ncol = (nord(4))
          elseif ((J.eq.Nsub2-1).and.(I.gt.1).and.(I.lt.Nsub1-1)) then
            Ncol = (nord(5))
          elseif ((I.eq.1).and.(J.gt.1).and.(J.lt.Nsub2-1)) then
            Ncol = (nord(6))
          else
            Ncol = (nord(11))
          endif
        endif
      endif
c
      if (Nside.eq.3) then
        if ((J.eq.1).and.(I.gt.1 .or. Icase==1)
     .              .and.(I.lt.Nsub1 .or. Icase==1)) then
          Ncol = (nord(1))
        elseif ((I.eq.Nsub1).and.(J.gt.1 .or. Icase==2)
     .                      .and.(J.lt.Nsub2 .or. Icase==1)) then
          Ncol = (nord(8))
        elseif ((J.eq.Nsub2).and.(I.gt.1 .or. Icase==2)
     .                      .and.(I.lt.Nsub1 .or. Icase==2)) then
          Ncol = (nord(4))
        elseif ((I.eq.1).and.(J.gt.1 .or. Icase==2)
     .                  .and.(J.lt.Nsub2 .or. Icase==1)) then
          Ncol = (nord(7))
        else
          go to 100
        endif
      endif
c
      if (Nside.eq.4) then
        if ((J.eq.1).and.(I.gt.1 .or. Icase==1)
     .              .and.(I.lt.Nsub1 .or. Icase==1)) then
          Ncol = (nord(2))
        elseif ((I.eq.Nsub1).and.(J.gt.1 .or. Icase==2)
     .                      .and.(J.lt.Nsub2 .or. Icase==1)) then
          Ncol = (nord(9))
        elseif ((J.eq.Nsub2).and.(I.gt.1 .or. Icase==2)
     .                      .and.(I.lt.Nsub1 .or. Icase==2)) then
          Ncol = (nord(5))
        elseif ((I.eq.1).and.(J.gt.1 .or. Icase==2)
     .                  .and.(J.lt.Nsub2 .or. Icase==1)) then
          Ncol = (nord(8))
        else
          go to 100
       endif
      endif
c
      if (Nside.eq.5) then
        if ((J.eq.1).and.(I.gt.1 .or. Icase==1)
     .              .and.(I.lt.Nsub1 .or. Icase==1)) then
          Ncol = (nord(3))
        elseif ((I.eq.Nsub1).and.(J.gt.1 .or. Icase==2)
     .                      .and.(J.lt.Nsub2 .or. Icase==1)) then
          Ncol = (nord(7))
        elseif ((J.eq.Nsub2).and.(I.gt.1 .or. Icase==2)
     .                      .and.(I.lt.Nsub1 .or. Icase==2)) then
          Ncol = (nord(6))
        elseif ((I.eq.1).and.(J.gt.1 .or. Icase==2)
     .                  .and.(J.lt.Nsub2 .or. Icase==1)) then
          Ncol = (nord(9))
        else
          go to 100
        endif
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'SETCOLP: Icase,Nside,Nsub1,Nsub2,I,J,Ncol = '
        write(*,*)  Icase,Nside,Nsub1,Nsub2,I,J,Ncol
        call pause
      endif
c
      return

 100  continue
      if( jcase==1 )then
         if( (ii>jj) .or. (ii==jj .and. Icase==1) )then
            Ncol = norh
         else
            Ncol = norv
         endif

      elseif( jcase==2 )then
       if( (ii+jj)<=nhalf .or. ( (ii+jj-1)==nhalf .and. Icase==1 ) )then
            Ncol = norh
         else
            Ncol = norv
         endif
      else
        write(*,*) 'SETCOLP: jcase = ',jcase
        stop
      endif

      return
      end

