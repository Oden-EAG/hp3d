c---------------------------------------------------------------------
c
c   routine name       - dofin
c
c---------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Sep 07
c
c   purpose            - copy dof for a node from local arrays into
c                        data structure arrays
c
c   arguments
c     in:
c          Nod         - a node number
c          ZvalH,ZvalE,ZvalV,ZvalQ - H1,H(curl),H(div) and L2 dof
c                        for the node in the expanded mode
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine dof_in(Nod,ZvalH,ZvalE,ZvalV,ZvalQ)
c
      use data_structure3D
c
#include "syscom.blk"
c
      dimension ZvalH(MAXEQNH,*), ZvalE(MAXEQNE,*),
     .          ZvalV(MAXEQNE,*), ZvalQ(MAXEQNQ,*)
c
c  ...node physics attributes flags
      dimension ncase(NR_PHYSA)
c
c
c  ...decode the physical attributes of the node
      call decod(NODES(Nod)%case,2,NR_PHYSA, ncase)
c
c  ...determine number of dof for the node
      call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
c
c  ...initiate the node component counters
      ivarH=0; ivarE=0; ivarV=0; ivarQ=0
c
c  ...loop through multiple copies of variables
      do j=1,NRCOMS
c
c  .....loop through physical attributes of the node
        do i=1,NR_PHYSA
          if (ncase(i).eq.1) then
c
c  .........loop through the components of the attribute
            do k=1,NR_COMP(i)
              select case(DTYPE(i))
              case('contin')
                nvar = (j-1)*NRHVAR+ADRES(i)+k
                ivarH = ivarH+1
                NODES(Nod)%zdofH(ivarH,1:ndofH)
     .            = ZvalH(nvar,1:ndofH)
              case('tangen')
                nvar = (j-1)*NREVAR+ADRES(i)+k
                ivarE = ivarE+1
                NODES(Nod)%zdofE(ivarE,1:ndofE)
     .            = ZvalE(nvar,1:ndofE)
              case('normal')
                nvar = (j-1)*NRVVAR+ADRES(i)+k
                ivarV = ivarV+1
                NODES(Nod)%zdofV(ivarV,1:ndofV)
     .            = ZvalV(nvar,1:ndofV)
              case('discon')
                nvar = (j-1)*NRQVAR+ADRES(i)+k
                ivarQ = ivarQ+1
                NODES(Nod)%zdofQ(ivarQ,1:ndofQ)
     .            = ZvalQ(nvar,1:ndofQ)
              case default
                write(*,*) 'dofin:INCONSISTENCY'; stop 1
              end select
c
c  .........end of loop through components
            enddo
c
c  .......if the node owns the physical attribute
          endif
c
c  .....end of loop through physical attributes
        enddo
c
c  ...end of loop through multiple copies
      enddo
c
c
      end




