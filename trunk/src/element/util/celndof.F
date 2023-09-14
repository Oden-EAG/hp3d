c----------------------------------------------------------------------
c
c   routine name       - celndof
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine determines number of single dof
c                        for a 3D element
c
c   arguments :
c     in:
c             Ntype    - middle node type
c             Norder   - order of approximation for the element
c     out:
c             NrdofH   - number of H1 dof
c             NrdofE   - number of H(curl) dof
c             NrdofV   - number of H(div) dof
c             NrdofQ   - number of L2 dof
c
c----------------------------------------------------------------------
c
      subroutine celndof(Ntype,Norder, NrdofH,NrdofE,NrdofV,NrdofQ)
c
      use element_data
c
      implicit none
c
      integer, intent(in)  :: Ntype,Norder(19)
      integer, intent(out) :: NrdofH,NrdofE,NrdofV,NrdofQ
c
      integer :: i,j
      integer :: nredg
      integer :: ndofH,ndofE,ndofV,ndofQ
c
#if DEBUG_MODE
      integer :: iprint
      iprint=0
#endif
      NrdofH=0; NrdofE=0; NrdofV=0; NrdofQ=0
c
c  ...V E R T I C E S..................................................
      do i=1,NVERT(Ntype)
        call ndof_nod(VERT,1, ndofH,ndofE,ndofV,ndofQ)
#if DEBUG_MODE
        if (iprint.eq.1) then
          write(*,7001) 'vert',1,ndofH,ndofE,ndofV,ndofQ
 7001     format('celndof: type,nord,ndofH,ndofE,ndofV,ndofQ = ',
     .                     a5,i4,4i6)
        endif
#endif
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
      enddo
c
c  ...E D G E S........................................................
      nredg = NEDGE(Ntype)
      do i=1,nredg
        call ndof_nod(MEDG,Norder(i), ndofH,ndofE,ndofV,ndofQ)
#if DEBUG_MODE
        if (iprint.eq.1) then
          write(*,7001) 'medg',Norder(i),ndofH,ndofE,ndofV,ndofQ
        endif
#endif
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
      enddo
      j=nredg
c
c  ...F A C E S   A N D   M I D D L E..................................
      select case(Ntype)
c  ...HEXA
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
c  ...PRISM
      case(PRIS,MDLP)
        do i=1,2
          j=j+1
          call ndof_nod(MDLT,Norder(j), ndofH,ndofE,ndofV,ndofQ)
#if DEBUG_MODE
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
#if DEBUG_MODE
          if (iprint.eq.1) then
            write(*,7001) 'mdlq',Norder(j),ndofH,ndofE,ndofV,ndofQ
          endif
#endif
          NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
          NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        enddo
        j=j+1
        call ndof_nod(MDLP,Norder(j), ndofH,ndofE,ndofV,ndofQ)
#if DEBUG_MODE
        if (iprint.eq.1) then
          write(*,7001) 'mdlp',Norder(j),ndofH,ndofE,ndofV,ndofQ
        endif
#endif
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
c  ...TETRA
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
c  ...PYRAMID
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
c
c
      end subroutine celndof
