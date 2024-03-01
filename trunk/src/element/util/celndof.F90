!----------------------------------------------------------------------
!
!   routine name       - celndof
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine determines number of single dof
!                        for a 3D element
!
!   arguments :
!     in:
!             Ntype    - middle node type
!             Norder   - order of approximation for the element
!     out:
!             NrdofH   - number of H1 dof
!             NrdofE   - number of H(curl) dof
!             NrdofV   - number of H(div) dof
!             NrdofQ   - number of L2 dof
!
!----------------------------------------------------------------------
!
   subroutine celndof(Ntype,Norder, NrdofH,NrdofE,NrdofV,NrdofQ)
!
      use element_data
!
      implicit none
!
      integer, intent(in)  :: Ntype,Norder(19)
      integer, intent(out) :: NrdofH,NrdofE,NrdofV,NrdofQ
!
      integer :: i,j
      integer :: nredg
      integer :: ndofH,ndofE,ndofV,ndofQ
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
      NrdofH=0; NrdofE=0; NrdofV=0; NrdofQ=0
!
!  ...V E R T I C E S..................................................
      do i=1,nvert(Ntype)
        call ndof_nod(VERT,1, ndofH,ndofE,ndofV,ndofQ)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7001) 'vert',1,ndofH,ndofE,ndofV,ndofQ
 7001     format('celndof: type,nord,ndofH,ndofE,ndofV,ndofQ = ', &
                           a5,i4,4i6)
        endif
#endif
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
      enddo
!
!  ...E D G E S........................................................
      nredg = nedge(Ntype)
      do i=1,nredg
        call ndof_nod(MEDG,Norder(i), ndofH,ndofE,ndofV,ndofQ)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7001) 'medg',Norder(i),ndofH,ndofE,ndofV,ndofQ
        endif
#endif
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
      enddo
      j=nredg
!
!  ...F A C E S   A N D   M I D D L E..................................
      select case(Ntype)
!  ...HEXA
      case(BRIC,MDLB)
        do i=1,6
          j=j+1
          call ndof_nod(MDLQ,Norder(j), ndofH,ndofE,ndofV,ndofQ)
          NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
          NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        enddo
        j=j+1
        call ndof_nod(MDLB,Norder(j), ndofH,ndofE,ndofV,ndofQ)
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
!  ...PRISM
      case(PRIS,MDLP)
        do i=1,2
          j=j+1
          call ndof_nod(MDLT,Norder(j), ndofH,ndofE,ndofV,ndofQ)
#if HP3D_DEBUG
          if (iprint.eq.1) then
            write(*,7001) 'mdlt',Norder(j),ndofH,ndofE,ndofV,ndofQ
          endif
#endif
          NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
          NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        enddo
        do i=1,3
          j=j+1
          call ndof_nod(MDLQ,Norder(j), ndofH,ndofE,ndofV,ndofQ)
#if HP3D_DEBUG
          if (iprint.eq.1) then
            write(*,7001) 'mdlq',Norder(j),ndofH,ndofE,ndofV,ndofQ
          endif
#endif
          NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
          NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        enddo
        j=j+1
        call ndof_nod(MDLP,Norder(j), ndofH,ndofE,ndofV,ndofQ)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7001) 'mdlp',Norder(j),ndofH,ndofE,ndofV,ndofQ
        endif
#endif
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
!  ...TETRA
      case(TETR,MDLN)
        do i=1,4
          j=j+1
          call ndof_nod(MDLT,Norder(j), ndofH,ndofE,ndofV,ndofQ)
          NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
          NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        enddo
        j=j+1
        call ndof_nod(MDLN,Norder(j), ndofH,ndofE,ndofV,ndofQ)
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
!  ...PYRAMID
      case(PYRA,MDLD)
        j=j+1
        call ndof_nod(MDLQ,Norder(j), ndofH,ndofE,ndofV,ndofQ)
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        do i=1,4
          j=j+1
          call ndof_nod(MDLT,Norder(j), ndofH,ndofE,ndofV,ndofQ)
          NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
          NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        enddo
        j=j+1
        call ndof_nod(MDLD,Norder(j), ndofH,ndofE,ndofV,ndofQ)
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
      end select
!
!
   end subroutine celndof
