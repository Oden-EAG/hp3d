#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - check_edge
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine identifies sharp curves
!
!   arguments
!     in:
!              Ncu     - a curve
!     out:
!              Idec    = 2 if the curve is on a sharp edge
!                      = 0,1 if not
!
!----------------------------------------------------------------------
!
      subroutine check_edge(Ncu, Idec)
!
      use GMP
      use control
!
      implicit none
!
      integer :: Ncu,Idec
!
!  ...list of attached triangles
      integer :: itrian(2)
      real(8) :: s,aux1(3),aux2(3),aux3(3),aux4(3)
!
      integer :: i,ifound,it,j,lab,nc,np,np1,np2,np3,np4
      integer :: nr_trian,nt1,nt2
!
!  ...list of sharp edges
      integer :: nsharp(1000)
!
      integer :: ivisit,nr_sharp
!
!  ...save visitation flag and list of sharp edges
      save ivisit,nsharp,nr_sharp
!
!  ...visitation flag
      ivisit = 0
!
!-----------------------------------------------------------------------
!
!  ...select routine vesion you want to use (100 - default)
      goto 100

!=======================================================================
!  VERSION 1: use surface numbers to determine sharp edges
!
!
!  STEP 0: if 1st visit determine list of sharp edges
!
 100  continue

      Idec=0
!
      if (ivisit.eq.0) then
        ivisit = 1
!
        nr_sharp = 0
!  .....loop over curves
        do nc = 1,NRCURVE
!  .......skip if not an Hermite curve
          if (CURVES(nc)%Type.ne.'HermCur') cycle
!  .......build up list of attached 'G1RecTri'
          itrian = 0; nr_trian = 0
!  .......loop over attached figures
          do i = 1,CURVES(nc)%NrFig
            call decode(CURVES(nc)%FigNo(i), it,lab)
            it = abs(it)
            if (TRIANGLES(it)%Type.ne.'G1RecTri') cycle
            call locate(it,itrian,nr_trian, j)
!  .........if triangle is not on the list
            if (j.eq.0) then
              nr_trian = nr_trian + 1
              itrian(nr_trian) = it
            endif
!  .......end of loop over attached figures
          enddo
          if (nr_trian.ne.2) then
            write(*,*)'check_edge: have found an Hermite curve', &
                      ' without 2 adjacent G1RecTri'
            write(*,*)'nc,nr_trian = ',nc,nr_trian
            write(*,*)'     itrian = ',itrian
            call print_GMP
            stop
          endif
          nt1 = itrian(1); nt2 = itrian(2)
          if (TRIANGLES(nt1)%Idata(1).ne.TRIANGLES(nt2)%Idata(1)) then
            nr_sharp = nr_sharp + 1
            nsharp(nr_sharp) = nc
          endif
!  .....end of loop over curves
        enddo
!
        if (nr_sharp.gt.0) then
          write(*,*)'-------------------------------------------'
          write(*,*)'SHARP EDGES:'
          do i = 1,nr_sharp
            write(*,*)'i,nsharp(i) = ',i,nsharp(i)
          enddo
          write(*,*)'-------------------------------------------'
          call pause
        endif
!
!  ...end of 1st visit
      endif
!
!
!  STEP 1: determine wheter curve is on the list of sharp edges
!
      Idec = 0
      call locate(Ncu,nsharp,nr_sharp, ifound)
      if (ifound.ne.0) Idec = 2
!
      return
!
!=======================================================================
!  VERSION 2: use dot product b/w normals to detect sharp edges
!
!
 200  continue
!
      if (ivisit.eq.0) then
        ivisit = 1
!
        nr_sharp = 0
!  .....loop over curves
        do nc = 1,NRCURVE
!  .......skip if not an Hermite curve
          if (CURVES(nc)%Type.ne.'HermCur') cycle
!  .......build up list of attached 'G1RecTri'
          itrian = 0; nr_trian = 0
!  .......loop over attached figures
          do i = 1,CURVES(nc)%NrFig
            call decode(CURVES(nc)%FigNo(i), it,lab)
            it = abs(it)
            if (TRIANGLES(it)%Type.ne.'G1RecTri') cycle
            call locate(it,itrian,nr_trian, j)
!  .........if triangle is not on the list
            if (j.eq.0) then
              nr_trian = nr_trian + 1
              itrian(nr_trian) = it
            endif
!  .......end of loop over attached figures
          enddo
          if (nr_trian.ne.2) then
            write(*,*)'check_edge: have found an Hermite curve' &
                      ' without 2 adjacent G1RecTri'
            write(*,*)'nc,nr_trian = ',nc,nr_trian
            write(*,*)'     itrian = ',itrian
            call print_GMP
            stop
          endif
          nt1 = itrian(1); nt2 = itrian(2)
          np1=CURVES(nc)%EndPoNo(1) ; np2=CURVES(nc)%EndPoNo(2)
          do i=1,3
            np=TRIANGLES(nt1)%VertNo(i)
            if ((np.ne.np1).and.(np.ne.np2)) np3=np
          enddo
          do i=1,3
            np=TRIANGLES(nt2)%VertNo(i)
            if ((np.ne.np1).and.(np.ne.np2)) np4=np
          enddo
          aux1(1:3)=POINTS(np3)%Rdata(1:3)-POINTS(np1)%Rdata(1:3)
          aux2(1:3)=POINTS(np2)%Rdata(1:3)-POINTS(np1)%Rdata(1:3)
          call cross_product(aux2,aux1, aux3)
          call normalize(aux3)
          aux1=aux3
          aux3(1:3)=POINTS(np4)%Rdata(1:3)-POINTS(np1)%Rdata(1:3)
          call cross_product(aux3,aux2, aux4)
          call normalize(aux4)
          call scalar_product(aux1,aux4, s)
!
!         SET DOT PRODUCT HERE ------------------------|
          if (s.lt.0.6d0) then !  <--------------------|
!
            nr_sharp = nr_sharp + 1
            nsharp(nr_sharp) = nc
          endif
!  .....end of loop over curves
        enddo
!
        if (nr_sharp.gt.0) then
          write(*,*)'-------------------------------------------'
          write(*,*)'SHARP EDGES:'
          do i = 1,nr_sharp
            write(*,*)'i,nsharp(i) = ',i,nsharp(i)
          enddo
          write(*,*)'-------------------------------------------'
          call pause
        endif
!
!  ...end of 1st visit
      endif
!
      Idec = 0
      call locate(Ncu,nsharp,nr_sharp, ifound)
      if (ifound.ne.0) Idec = 2
!
      return
!
!=======================================================================
!  VERSION 3: explicitly list sharp edges
!
 300  continue
      nr_sharp=194

      nsharp = 0
      nsharp( 1) =      85
      nsharp( 2) =      86
      nsharp( 3) =      87
      nsharp( 4) =     314
      nsharp( 5) =     350
      nsharp( 6) =     351
      nsharp( 7) =     416
      nsharp( 8) =     740
      nsharp( 9) =     836
      nsharp(10) =     933
      nsharp(11) =     955
      nsharp(12) =    1310
      nsharp(13) =    1320
      nsharp(14) =    1321
      nsharp(15) =    1499
      nsharp(16) =    1831
      nsharp(17) =    1996
      nsharp(18) =    2419
      nsharp(19) =    2502
      nsharp(20) =    2614
      nsharp(21) =    2980
      nsharp(22) =    3477
      nsharp(23) =    3528
      nsharp(24) =    3554
      nsharp(25) =    3619
      nsharp(26) =    4331
      nsharp(27) =    4380
      nsharp(28) =    4473
      nsharp(29) =    4496
      nsharp(30) =    4502
      nsharp(31) =    4844
      nsharp(32) =    4906
      nsharp(33) =    4986
      nsharp(34) =    5207
      nsharp(35) =    5235
      nsharp(36) =    5387
      nsharp(37) =    5434
      nsharp(38) =    6659
      nsharp(39) =    6806
      nsharp(40) =    6853
      nsharp(41) =    6855
      nsharp(42) =    6873
      nsharp(43) =    6942
      nsharp(44) =    7077
      nsharp(45) =    7078
      nsharp(46) =    7117
      nsharp(47) =    7134
      nsharp(48) =    7153
      nsharp(49) =    7154
      nsharp(50) =    7682
      nsharp(51) =    7736
      nsharp(52) =    7754
      nsharp(53) =    7781
      nsharp(54) =    7800
      nsharp(55) =    7830
      nsharp(56) =    7962
      nsharp(57) =    8357
      nsharp(58) =    8434
      nsharp(59) =    8445
      nsharp(60) =    8488
      nsharp(61) =    8581
      nsharp(62) =    8784
      nsharp(63) =    8992
      nsharp(64) =    9089
      nsharp(65) =    9171
      nsharp(66) =    9411
      nsharp(67) =    9424
      nsharp(68) =    9437
      nsharp(69) =    9440
      nsharp(70) =    9464
      nsharp(71) =    9565
      nsharp(72) =    9643
      nsharp(73) =    9644
      nsharp(74) =    9781
      nsharp(75) =    9805
      nsharp(76) =   10008
      nsharp(77) =   10116
      nsharp(78) =   10204
      nsharp(79) =   10209
      nsharp(80) =   10451
      nsharp(81) =   10509
      nsharp(82) =   10541
      nsharp(83) =   10913
      nsharp(84) =   10951
      nsharp(85) =   10952
      nsharp(86) =   10979
      nsharp(87) =   11369
      nsharp(88) =   11513
      nsharp(89) =   11559
      nsharp(90) =   11565
      nsharp(91) =   11566
      nsharp(92) =   11684
      nsharp(93) =   11822
      nsharp(94) =   11853
      nsharp(95) =   11925
      nsharp(96) =   12096
      nsharp(97) =   12171
      nsharp(98) =   12253
      nsharp(99) =   12528
      nsharp(100) =   12565
      nsharp(101) =   12895
      nsharp(102) =   13086
      nsharp(103) =   13087
      nsharp(104) =   13184
      nsharp(105) =   13212
      nsharp(106) =   13496
      nsharp(107) =   13507
      nsharp(108) =   13534
      nsharp(109) =   13549
      nsharp(110) =   13628
      nsharp(111) =   14145
      nsharp(112) =   14465
      nsharp(113) =   14790
      nsharp(114) =   15073
      nsharp(115) =   15138
      nsharp(116) =   15159
      nsharp(117) =   15396
      nsharp(118) =   15490
      nsharp(119) =   15777
      nsharp(120) =   15779
      nsharp(121) =   15807
      nsharp(122) =   16008
      nsharp(123) =   16564
      nsharp(124) =   16836
      nsharp(125) =   17468
      nsharp(126) =   17573
      nsharp(127) =   17607
      nsharp(128) =   17752
      nsharp(129) =   17878
      nsharp(130) =   18118
      nsharp(131) =   18486
      nsharp(132) =   18536
      nsharp(133) =   18537
      nsharp(134) =   19223
      nsharp(135) =   19607
      nsharp(136) =   19943
      nsharp(137) =   20002
      nsharp(138) =   20048
      nsharp(139) =   20322
      nsharp(140) =   20476
      nsharp(141) =   20574
      nsharp(142) =   20723
      nsharp(143) =   20830
      nsharp(144) =   20973
      nsharp(145) =   20974
      nsharp(146) =   21102
      nsharp(147) =   21103
      nsharp(148) =   21211
      nsharp(149) =   21568
      nsharp(150) =   21822
      nsharp(151) =   21873
      nsharp(152) =   21943
      nsharp(153) =   22161
      nsharp(154) =   22162
      nsharp(155) =   22223
      nsharp(156) =   22304
      nsharp(157) =   22305
      nsharp(158) =   22769
      nsharp(159) =   22770
      nsharp(160) =   23259
      nsharp(161) =   23263
      nsharp(162) =   23365
      nsharp(163) =   23449
      nsharp(164) =   23450
      nsharp(165) =   23852
      nsharp(166) =   24190
      nsharp(167) =   24450
      nsharp(168) =   24599
      nsharp(169) =   25002
      nsharp(170) =   25151
      nsharp(171) =   25182
      nsharp(172) =   25293
      nsharp(173) =   25912
      nsharp(174) =   26625
      nsharp(175) =   26853
      nsharp(176) =   27092
      nsharp(177) =   27487
      nsharp(178) =   27686
      nsharp(179) =   28198
      nsharp(180) =   28498
      nsharp(181) =   28506
      nsharp(182) =   28652
      nsharp(183) =   30254
      nsharp(184) =   30344
      nsharp(185) =   32543
      nsharp(186) =   35106
      nsharp(187) =   36202
      nsharp(188) =   36242
      nsharp(189) =   37164
      nsharp(190) =   40077
      nsharp(191) =   42058
      nsharp(192) =   43066
      nsharp(193) =   48371
      nsharp(194) =   50711
!
!
      Idec = 0
      call locate(Ncu,nsharp,nr_sharp, ifound)
      if (ifound.ne.0) Idec = 2
!
!
      end subroutine check_edge

#endif
