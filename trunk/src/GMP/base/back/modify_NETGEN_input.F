c//////////////////////////////////////////////////////////////////////
c P. Gatto, Nov 12 : since routine input_NETGEN was modified to
c                    automatically determine surface numbers, this
c                    routine is not really needed anymore.
c//////////////////////////////////////////////////////////////////////
c
c
c----------------------------------------------------------------------
c
c   routine name       - modify_NETGEN_input
c
c----------------------------------------------------------------------
c
c   latest revision    - Jul 08
c
c   purpose            - Routine scans the NETGEN input file for
c                        duplicated surfaces and modifies the data
c                        to avoid it
c
c   arguments            none
c
c---------------------------------------------------------------------
c
      subroutine modify_NETGEN_input
c
      use GMP
#include "syscom.blk"
ccc #include "cinout.blk" (removed common block)
c
c  ...defines "old" to "new" surface number function
      integer, dimension(:), pointer:: no_surf
c
c  ...work space
      dimension rdata(10),nvert(3)
c
c  ...surface type
      character(len=10) :: type
c
c---------------------------------------------------------------------
c
      iprint=0
c
c  ...set up the dimension of the problem
      NDIM=3; MANDIM=3
c
c  ...set up the unused parameters
      NRRECTA=0 ; NRHEXAS=0 ; NRPRISM=0 ; NRPYRAM=0
c
c  ...allocate memory for GMP data structure
      call alloc_GMP
c
c  ...open the modified input file
      kout=30
      open(unit=kout,file='files/input_tets1',form='formatted',
     .     access='sequential',status='unknown')
c
c---------------------------------------------------------------------
c
c  Step 1: Input data on algebraic surfaces
c
c
c  ...read in number of surfaces
      read(KIN,*) NRSURFS
      if (MAXSU.lt.NRSURFS) then
        write(*,7010)
 7010   format(' modify_NETGEN_input: INCREASE MAXSU')
        stop 1
      endif
      allocate(no_surf(NRSURFS))
      no_surf(1:NRSURFS)=0
      if (iprint.eq.1) then
        write(*,7011)
 7011   format(' modify_NETGEN_input: READING SURFACES...')
      endif
c
c  ...initiate the number of surfaces after elimination of duplicated
c     surfaces
      ns=0
c
c  ...loop through surfaces in the input file
      do js=1,NRSURFS
        read(KIN,*) type
        write(*,*) ' js,type = ',js,type
c
        select case(type)
c
c  .......a plane normal to a given vector and passing through
c         a point, orientation specified by the normal
          case('VecPt')
            read(KIN,*) rdata(1:3)
            read(KIN,*) rdata(4:6)
c
c  .........check for duplication
            do is=1,ns
              iflag=0
              if (SURFACES(is)%Type.eq.'VecPt') then
                iflag=1
                do k=1,6
                  if (SURFACES(is)%Rdata(k).ne.rdata(k)) iflag=0
                enddo
                if (iprint.eq.1) then
                  write(*,*)'js,is = ',js,is
                  write(*,6001) rdata(1:6)
                  write(*,6001) SURFACES(is)%Rdata(1:6)
 6001             format(10e12.5)
                  write(*,*) 'iflag = ',iflag
                  call pause
                endif
              endif
              if (iflag.eq.1) then
                no_surf(js)=is
                write(*,*) 'modify_NETGEN_input: ELIMINATE PLANE ',js
                call pause
              endif
            enddo
            if (no_surf(js).eq.0) then
c
c  ...........new surface
              ns=ns+1; no_surf(js)=ns
              SURFACES(ns)%Type = type
              allocate(SURFACES(ns)%Rdata(6))
              SURFACES(ns)%Rdata(1:6) = rdata(1:6)
            endif
c
c  .......a sphere with radius and central point
          case('Sphere')
            read(KIN,*) rdata(1:3)
            read(KIN,*) rdata(4)
c
c  .........check for duplication
            do is=1,ns
              iflag=0
              if (SURFACES(is)%Type.eq.'Sphere') then
                iflag=1
                do k=1,4
                  if (SURFACES(is)%Rdata(k).ne.rdata(k)) iflag=0
                enddo
              endif
              if (iflag.eq.1) then
                no_surf(js)=is
                write(*,*) 'modify_NETGEN_input: ELIMINATE SPHERE ',js
                call pause
              endif
            enddo
            if (no_surf(js).eq.0) then
c
c  ...........new surface
              ns=ns+1; no_surf(js)=ns
              SURFACES(ns)%Type = type
              allocate(SURFACES(ns)%Rdata(4))
              SURFACES(ns)%Rdata(1:4) = rdata(1:4)
            endif
c
c    .....a cylinder defined with an endpoint of its axis, direction
c         vector and radius
          case('Cylinder')
            read(KIN,*) rdata(1:3)
            read(KIN,*) rdata(4:6)
            read(KIN,*) rdata(7)
c
c  .........check for duplication
            do is=1,ns
              iflag=0
              if (SURFACES(is)%Type.eq.'Cylinder') then
                write(*,*) 'is = ',is
                write(*,*) 'SURFACES(is)%Rdata = ',
     .                      SURFACES(is)%Rdata(1:7)
                write(*,*) '             rdata = ',rdata(1:7)
                call pause
                iflag=1
                do k=1,7
                  if (SURFACES(is)%Rdata(k).ne.rdata(k)) iflag=0
                enddo
              endif
              if (iflag.eq.1) then
                no_surf(js)=is
                write(*,*) 'modify_NETGEN_input: ELIMINATE CYLINDER ',js
                call pause
              endif
            enddo
            if (no_surf(js).eq.0) then
c
c  ...........new surface
              ns=ns+1; no_surf(js)=ns
              write(*,*) 'STORING DATA FOR js = ',js
              SURFACES(ns)%Type = type
              allocate(SURFACES(ns)%Rdata(7))
              SURFACES(ns)%Rdata(1:7) = rdata(1:7)
            endif
c
c  .......cone
          case('Cone')
            read(KIN,*) rdata(1:3)
            read(KIN,*) rdata(4:6)
            read(KIN,*) rdata(7)
c
c  .........check for duplication
            do is=1,ns
              iflag=0
              if (SURFACES(is)%Type.eq.'Cone') then
                iflag=1
                do k=1,7
                  if (SURFACES(is)%Rdata(k).ne.rdata(k)) iflag=0
                enddo
              endif
              if (iflag.eq.1) then
                no_surf(js)=is
                write(*,*) 'modify_NETGEN_input: ELIMINATE CONE ',js
                call pause
              endif
            enddo
            if (no_surf(js).eq.0) then
c
c  ...........new surface
              ns=ns+1; no_surf(js)=ns
              SURFACES(ns)%Type = type
              allocate(SURFACES(ns)%Rdata(7))
              SURFACES(ns)%Rdata(1:7) = rdata(1:7)
            endif
c
c  .......reconstructed surface
          case('Rec')
c
c  .......other cases
          case default
            write(*,7012)
 7012       format('modify_NETGEN_input: UNKNOWN SURFACE TYPE')
            write(*,*) type
            stop 1
        end select
      enddo
c
c  ...input number of subdomains
      read(KIN,*) NRDOMAIN
      write(kout,*) NRDOMAIN
      write(*,7500) ns,no_surf(1:NRSURFS)
 7500 format('modify_NETGEN_input: ns = ',i3,' no_surf = ',20i4)
      call pause
c
c---------------------------------------------------------------------
c
c  Step 2: Input points
c
c
c  ...read in number of points
      read(KIN,*) NRPOINT
      write(kout,*) NRPOINT
      if (MAXNP.lt.NRPOINT) then
        write(*,*) 'input_NETGEN:INCREASE MAXNP'
        stop 1
      endif
ccc      allocate(POINTS(MAXNP))
      do np=1,NRPOINT
        POINTS(np)%Type  = 'Regular'
        allocate(POINTS(np)%Rdata(3))
        read(KIN,*) POINTS(np)%Rdata(1:3)
        write(kout,8001) POINTS(np)%Rdata(1:3)
      enddo
c
c---------------------------------------------------------------------
c
c  Step 3: Input tetrahedra
c
c  ...read in number of tets
      read(KIN,*) NRTETRA
      write(kout,*) NRTETRA
      if (NRTETRA.gt.MAXTE) then
        write(*,7013) NRTETRA, MAXTE
 7013   format('input_NETGEN: INCREASE MAXTE, NRTETRA,MAXTE = ',2i8)
        stop 1
      endif
      do ntet=1,NRTETRA
        TETRAS(ntet)%Type = 'Linear'
        TETRAS(ntet)%FigNo(1:4) = 0
        TETRAS(ntet)%EdgeNo(1:6) = 0
        read(KIN,*) TETRAS(ntet)%domain,TETRAS(ntet)%VertNo(1:4)
        write(kout,8002) TETRAS(ntet)%domain,TETRAS(ntet)%VertNo(1:4)
      enddo
c
c----------------------------------------------------------------------
c
c  Step 4: for each surface (both algebraic and reconstructed)
c          read in the surface triangles and change the surface number
c          accordingly
c
c
c  ...read the number of surface triangles
      read(KIN,*) no_surface_triangles
      write(kout,*) no_surface_triangles
c
c  ...loop through surface triangles
      do itsurf=1,no_surface_triangles
        read(KIN,*) nos, nvert(1:3)
        write(kout,8002) no_surf(nos), nvert(1:3)
      enddo
c
c     NOUT was defined in now-removed common block "cinout.blk"
      close(NOUT)
      stop
c
 8001 format(3e15.6)
 8002 format(i4,5i9)
c
      end subroutine modify_NETGEN_input

