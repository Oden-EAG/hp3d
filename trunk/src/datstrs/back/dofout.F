c---------------------------------------------------------------------
c
c   routine name       - dofout
c
c---------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - copy dof for a node from data structure
c                        arrays into local arrays
c
c   arguments
c     in:
c          Nod         - a node number
c     in/out:
c          KdofH,KdofE,KdofV,KdofQ - dof counters for the
c                        local arrays
c     out:
c          ZvalH,ZvalE,ZvalV,ZvalQ - H1,H(curl),H(div) and L2 dof
c                        extracted from the data structure in the
c                        expanded mode
c
c----------------------------------------------------------------------
c
      subroutine dof_out(Nod, KdofH,KdofE,KdofV,KdofQ,
     .                   ZvalH,ZvalE,ZvalV,ZvalQ)
c
      use data_structure3D
#include "syscom.blk"
c
      dimension ZvalH(MAXEQNH,*), ZvalE(MAXEQNE,*),
     .          ZvalV(MAXEQNE,*), ZvalQ(MAXEQNQ,*)
c
c  ...node physics attributes flags
      dimension ncase(NR_PHYSA)
c
      iprint=1
c
c-----------------------------------------------------------------------
c
c  ...decode the physical attributes of the node
      call decod(NODES(Nod)%case,2,NR_PHYSA, ncase)
c
c  ...determine number of dof for the node
      call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
      if (iprint.eq.1) then
        write(*,*) 'dof_out: Nod,ndofH,ndofE = ',Nod,ndofH,ndofE
      endif
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
                ivarH = ivarH+1
                nvar = (j-1)*NRHVAR+ADRES(i)+k
                if (ndofH.gt.0) then
                  ZvalH(nvar,KdofH+1:KdofH+ndofH)
     .              = NODES(Nod)%zdofH(ivarH,1:ndofH)
                  if (iprint.eq.1) then
                    write(*,7001) DTYPE(i),
     .                ivarH,nvar,KdofH+1,KdofH+ndofH,
     .                ZvalH(nvar,KdofH+1:KdofH+ndofH)
                  endif
                endif
              case('tangen')
                ivarE = ivarE+1
                nvar = (j-1)*NREVAR+ADRES(i)+k
                if (ndofE.gt.0) then
                  ZvalE(nvar,KdofE+1:KdofE+ndofE)
     .              = NODES(Nod)%zdofE(ivarE,1:ndofE)
                  if (iprint.eq.1) then
                    write(*,*) 'kj ',ZvalE(nvar, KdofE+1:KdofE+ndofE)
                    write(*,7001) DTYPE(i),
     .                ivarE,nvar,KdofE+1,KdofE+ndofE,
     .                ZvalE(nvar,KdofE+1:KdofE+ndofE)
                  endif
                endif
              case('normal')
                ivarV = ivarV+1
                nvar = (j-1)*NRVVAR+ADRES(i)+k
                ZvalV(nvar,KdofV+1:KdofV+ndofV)
     .            = NODES(Nod)%zdofV(ivarV,1:ndofV)
              case('discon')
                ivarQ = ivarQ+1
                nvar = (j-1)*NRQVAR+ADRES(i)+k
                ZvalQ(nvar,KdofQ+1:KdofQ+ndofQ)
     .            = NODES(Nod)%zdofQ(ivarQ,1:ndofQ)
              case default
                write(*,*) 'dofout:INCONSISTENCY'; stop 1
              end select
#if C_MODE
 7001         format('dof_out: ',a6,i3,
     .          ' Zval(',i3,',',i3,':',i3,') = ',
     .          10(/,5(2e12.5,2x)))
#else
 7001         format('dof_out: ',a6,i3,
     .          ' Zval(',i3,',',i3,':',i3,') = ',
     .          10(/,10e12.5))
#endif

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
c  ...update the local dof counters
      KdofH= KdofH+ndofH
      KdofE= KdofE+ndofE
      KdofV= KdofV+ndofV
      KdofQ= KdofQ+ndofQ
c
c
      end




