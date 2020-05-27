c----------------------------------------------------------------------
c
c   routine name       - celndof
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jan 08
c
c   purpose            - routine determines number of single dof
c                        for a 3D element
c
c   arguments :
c     in:
c             Type     - middle node type
c             Norder   - order of approximation for the element
c     out:
c             NrdofH   - number of H1 dof
c             NrdofE   - number of H(curl) dof
c             NrdofV   - number of H(div) dof
c             NrdofQ   - number of L2 dof
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine celndof(Type,Norder, NrdofH,NrdofE,NrdofV,NrdofQ)
c
      use element_data
c
      character (len=4) :: Type
c  ...use brick to dimentionalize
      dimension Norder(19)
c
      iprint=0
c
      NrdofH=0; NrdofE=0; NrdofV=0; NrdofQ=0
c
c  ...V E R T I C E S..................................................
      do i=1,nvert(Type)
        call ndof_nod('vert',1, ndofH,ndofE,ndofV,ndofQ)
        if (iprint.eq.1) then
          write(*,7001) 'vert',1,ndofH,ndofE,ndofV,ndofQ
 7001     format('celndof: type,nord,ndofH,ndofE,ndofV,ndofQ = ',
     .                     a5,i4,4i6)
        endif
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
      enddo
c
c  ...E D G E S........................................................
      nredg = nedge(Type)
      do i=1,nredg
        call ndof_nod('medg',Norder(i), ndofH,ndofE,ndofV,ndofQ)
        if (iprint.eq.1) then
          write(*,7001) 'medg',Norder(i),ndofH,ndofE,ndofV,ndofQ
        endif
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
      enddo
      j=nredg
c
c  ...F A C E S   A N D   M I D D L E..................................
      select case(Type)
c  ...PRISM
      case('pris','mdlp')
        do i=1,2
          j=j+1
          call ndof_nod('mdlt',Norder(j), ndofH,ndofE,ndofV,ndofQ)
          if (iprint.eq.1) then
            write(*,7001) 'mdlt',Norder(j),ndofH,ndofE,ndofV,ndofQ
          endif
          NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
          NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        enddo
        do i=1,3
          j=j+1
          call ndof_nod('mdlq',Norder(j), ndofH,ndofE,ndofV,ndofQ)
          if (iprint.eq.1) then
            write(*,7001) 'mdlq',Norder(j),ndofH,ndofE,ndofV,ndofQ
          endif
          NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
          NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        enddo
        j=j+1
        call ndof_nod('mdlp',Norder(j), ndofH,ndofE,ndofV,ndofQ)
        if (iprint.eq.1) then
          write(*,7001) 'mdlp',Norder(j),ndofH,ndofE,ndofV,ndofQ
        endif
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
c  ...HEXA
      case('bric','mdlb')
        do i=1,6
          j=j+1
          call ndof_nod('mdlq',Norder(j), ndofH,ndofE,ndofV,ndofQ)
          NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
          NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        enddo
        j=j+1
        call ndof_nod('mdlb',Norder(j), ndofH,ndofE,ndofV,ndofQ)
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
c  ...TETRA
      case('tetr','mdln')
        do i=1,4
          j=j+1
          call ndof_nod('mdlt',Norder(j), ndofH,ndofE,ndofV,ndofQ)
          NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
          NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        enddo
        j=j+1
        call ndof_nod('mdln',Norder(j), ndofH,ndofE,ndofV,ndofQ)
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
c  ...PYRAMID
      case('pyra','mdld')
        j=j+1
        call ndof_nod('mdlq',Norder(j), ndofH,ndofE,ndofV,ndofQ)
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        do i=1,4
          j=j+1
          call ndof_nod('mdlt',Norder(j), ndofH,ndofE,ndofV,ndofQ)
          NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
          NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
        enddo
        j=j+1
        call ndof_nod('mdld',Norder(j), ndofH,ndofE,ndofV,ndofQ)
        NrdofH=NrdofH+ndofH; NrdofE=NrdofE+ndofE
        NrdofV=NrdofV+ndofV; NrdofQ=NrdofQ+ndofQ
      end select
c
c
      end
