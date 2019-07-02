c----------------------------------------------------------------------
c
c   routine name       - logic
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 08
c
c   purpose            - routine establishes list of nodes
c                        for a modified element, and the corresponding
c                        constraints coefficients
c
c   arguments :
c     in:
c            Mdle      - an element number, same as the middle node
c                        number
c            Idec      = 1 if only the element nodes are to be
c                          determined
c     out:
c            Nodm      - actual (unconstrained) nodes returned in the
c                        standard order: vertex, mid-edge, mid-face
c                        and middle nodes
c            NdofmH,NdofmE,NdofmV - the corresponding number of
c                        H1, H(curl) and H(div) dof
c            Nrnodm    - number of modified element nodes
c            NrconH    - number of parent dof's for an i-th local H1 dof
c            NacH(j,i),j=1,...,nrcon(i) - list of parent dof's 
c                        for i-th local dof
c            ConstrH(j,i),j=1,...,nrcon(i) - the corresponding
c                        constraints coefficients
c            NrconE,NacE,ConstrE - same for H(curl) dof
c            NrconV,NacV,ConstrV - same for H(div) dof
c
c----------------------------------------------------------------------
c
      subroutine logic(Mdle,Idec,
     .                 Nodm,NdofmH,NdofmE,NdofmV,Nrnodm,
     .                 NrconH,NacH,ConstrH,
     .                 NrconE,NacE,ConstrE,
     .                 NrconV,NacV,ConstrV)
c
      use element_data
      use data_structure3D
#include "syscom.blk"
c      
      dimension Nodm(MAXNODM),NdofmH(MAXNODM),NdofmE(MAXNODM),
     .                        NdofmV(MAXNODM),NdofmQ(MAXNODM)
c
      dimension nrconH(MAXbrickH),nrconE(MAXbrickE),nrconV(MAXbrickV),
     .          nacH(NACDIM,MAXbrickH),constrH(NACDIM,MAXbrickH),
     .          nacE(NACDIM,MAXbrickE),constrE(NACDIM,MAXbrickE),
     .          nacV(NACDIM,MAXbrickV),constrV(NACDIM,MAXbrickV)
c
c  ...element type
      character(len=4) :: type
c
c  ...order for the element
      dimension norder(19), norientl(27)
c
c----------------------------------------------------------------------
c
      iprint=0
c
c  ...this is a temporary version for unconstrained elements only
c
      call elem_nodes(Mdle, Nodm, norientl)
      call find_order(Mdle, norder)
      type = NODES(Mdle)%type

ccc      write(*,*)'type = ',type

      Nrnodm = nvert(type)+nedge(type)+nface(type)+1
      do i=1,nvert(type)
        call ndof_nod('vert',1, 
     .    NdofmH(i),NdofmE(i),NdofmV(i),NdofmQ(i))
      enddo
      do j=1,nedge(type)+nface(type)+1
        i=nvert(type)+j

ccc        write(*,*)'i,Nodm(i) = ',i,Nodm(i)

        nod = Nodm(i)
        call ndof_nod(NODES(nod)%type,norder(j), 
     .    NdofmH(i),NdofmE(i),NdofmV(i),NdofmQ(i))
      enddo
c      
      if (iprint.eq.1) then
        write(*,7001) Mdle,NODES(Mdle)%type,Nrnodm
 7001   format('logic: Mdle, type, Nrnodm = ',i8,2x,a5,2x,i4)
        do i=1,Nrnodm
          write(*,7002)i,Nodm(i),NdofmH(i),NdofmE(i),NdofmV(i),NdofmQ(i)
 7002     format(7x,'i,Nodm(i),NdofmH(i),NdofmE(i),NdofmV(i),NdofmQ(i)',
     .           ' = ',i3,i8,2x,4i4)
        enddo
        call pause
      endif
      if (Idec.eq.1) return
c
      call celndof(type,Norder, nrdofH,nrdofE,nrdofV,nrdofQ)
c      
      do k=1,nrdofH
        nrconH(k)=1; nacH(1,k)=k; constrH(1,k) = 1.d0
      enddo
c
      do k=1,nrdofE
        nrconE(k)=1; nacE(1,k)=k; constrE(1,k) = 1.d0
      enddo
c
      do k=1,nrdofV
        nrconV(k)=1; nacV(1,k)=k; constrV(1,k) = 1.d0
      enddo
c
c
      end
      
