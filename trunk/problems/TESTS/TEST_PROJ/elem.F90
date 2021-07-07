subroutine elem(Mdle, Itest,Jtrial)
!
      use parameters , only : ZERO
      use physics    , only : NR_PHYSA
      use assembly   , only : ALOC,BLOC,NR_RHS
!
      implicit none
      integer,                    intent(in ) :: Mdle
      integer,dimension(NR_PHYSA),intent(out) :: Itest,Jtrial
!--------------------------------------------------------------------------
!
!  ...initialize to 0
      Itest(1:NR_PHYSA)=0 ; Jtrial(1:NR_PHYSA)=0
!
!==========================================================================
!     H1
!==========================================================================
      Itest(1)=1 ; Jtrial(1)=1
      call elem_h1(   Mdle, BLOC(1  )%array,BLOC(1  )%nrow,NR_RHS, &
                            ALOC(1,1)%array,ALOC(1,1)%ncol         )
!  ...off-diagonal (unused) blocks                             
      ALOC(1,2)%array=ZERO ; ALOC(1,3)%array=ZERO ; ALOC(1,4)%array=ZERO
!                       
!
!==========================================================================
!     H(curl)
!==========================================================================
      Itest(2)=1 ; Jtrial(2)=1
      call elem_hcurl(Mdle, BLOC(2  )%array,BLOC(2  )%nrow,NR_RHS, &
                            ALOC(2,2)%array,ALOC(2,2)%ncol         )
!  ...off-diagonal (unused) blocks                             
      ALOC(2,1)%array=ZERO ; ALOC(2,3)%array=ZERO ; ALOC(2,4)%array=ZERO
!
!
!==========================================================================
!     H(div)
!==========================================================================
      Itest(3)=1 ; Jtrial(3)=1
      call elem_hdiv( Mdle, BLOC(3  )%array,BLOC(3  )%nrow,NR_RHS, &
                            ALOC(3,3)%array,ALOC(3,3)%ncol         )
!  ...off-diagonal (unused) blocks                             
      ALOC(3,1)%array=ZERO ; ALOC(3,2)%array=ZERO ; ALOC(3,4)%array=ZERO
!
!
!==========================================================================
!     L2
!==========================================================================
      Itest(4)=1 ; Jtrial(4)=1
      call elem_l2(   Mdle, BLOC(4  )%array,BLOC(4  )%nrow,NR_RHS, &
                            ALOC(4,4)%array,ALOC(4,4)%ncol         )
!  ...off-diagonal (unused) blocks                             
      ALOC(4,1)%array=ZERO ; ALOC(4,2)%array=ZERO ; ALOC(4,3)%array=ZERO
!
!
endsubroutine elem
