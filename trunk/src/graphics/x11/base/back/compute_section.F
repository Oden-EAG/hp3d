c----------------------------------------------------------------------
c
c   routine name       - compute section
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jan 09
c
c   purpose            - routine computes values of a desired quantity
c                        along a straight line
c
c   arguments :
c     in:
c             Xpoint   - a point within the FE domain
c             Vector   - components of a vector specifying the direction
c                        of the line
c             Nrval    - number of quantities to be drwan
c             Na       - dimension of Array below
c     out:
c             Nrsegm   - number of straight line segments to be drawn
c             Array    - contains data to be drawn, for an n-th segment
c                        Array(1:2,n) contain endpoints
c                          coordinates (measured along the specified
c                          line)
c                        Array(3:Nrval-1,n) contain the corresponding
c                          values of the functions to be drawn
c
c----------------------------------------------------------------------
c
      subroutine compute_section(Xpoint,Vector,Nrval,Na, Nrsegm,Array)
c
      use control
      use physics
      use data_structure3D
#include "syscom.blk"
      common /cnodcorb/ iprint_nodcorb
      common /cinvmap/ iprint_invmap
      common /csoleval/ iprint_soleval
c
c  ...real and imaginary part
      double precision, external :: dreal_part,dimag_part
c
      dimension Xpoint(3),Vector(3),Array(2+Nrval,Na)
c
c  ...geometry dof
      dimension xnod(NDIMEN,MAXbrickH)
c
c  ...solution dof
      dimension zdofH(MAXEQNH,MAXbrickH)
      dimension zdofE(MAXEQNE,MAXbrickE)
      dimension zdofV(MAXEQNV,MAXbrickV)
      dimension zdofQ(MAXEQNQ,MAXbrickQ)
c
c  ...master element coordinates of endpoints
      dimension xi(3,2)
c  ...points in element local and global coordinates systems
      dimension xi_moved(3,2),xi_new(3,2),xi_org(3,2) !LCS
      dimension x_moved(3,2),x_moved2(3,2),x_org(3,2) !GCS
c
c  ...physical coordinates and solution values
      dimension x(3,2),dxdxi(3,3),dxidx(3,3)
c
c  ...a point, max of elements in two directions
      dimension xp(3), jmax(2), x_save(3,2)
c
c  ...shape functions
      dimension shapH(MAXbrickH), gradH(3,MAXbrickH)
c
c  ...work space
      parameter (nwork = 1000)
      double precision ,allocatable :: work(:,:,:)
c
      double precision ,allocatable :: work2(:,:)
c
c  ...element order
      dimension norder(19)
c
c  ...node orientations
      dimension nedge_orient(12), nface_orient(6)
c
c  ...middle nodes of two elements
      dimension mdles(2)
c
c  ...approximate solution at a point
      dimension zsolH(MAXEQNH),zgradH(MAXEQNH,3),
     .          zsolE(MAXEQNE,3),zcurlE(MAXEQNE,3),
     .          zsolV(MAXEQNV,3),zdivV(MAXEQNV),
     .          zsolQ(MAXEQNQ)
c
      dimension zpinc(MAXNRHS),zdpinc(3,MAXNRHS)
c
c  ...pressure gradient
      dimension zpressure_grad(3)
c
c-----------------------------------------------------------------------
c
      iprint_nodcorb=0
      iprint=0 !2 for bisection debugging
      if (iprint.eq.1) then
        write(*,7001) Xpoint(1:3), Vector(1:3)
 7001   format('comp_section: Xpoint,Vector = ',3f8.3,2x,3f8.3)
        call pause
      endif
c
      write(*,*)'SELECT QUANTITY TO PLOT'
      write(*,*)'REAL PRESSURE............1'
      write(*,*)'IMAG PRESSURE............2'
      write(*,*)'MODULE OF PRESSURE.......3'
      write(*,*)'MODULE OF PRESSURE [DB]..4'
      read(*,*)ichoose

c
      allocate(work(1:2,1:2,-nwork:nwork))
      allocate(work2(1:10,-nwork:nwork))
c
      call normalize(Vector)
c
      step  = 0.01d0
c
c----------------------------------------------------------------------
c
c  ...find element containing point Xpoint
 6    mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call find_xi(mdle,Xpoint(1:3),Vector,
     .    idec,xi(1:3,1),xi(1:3,2))
        if (idec.eq.1) go to 11
      enddo
      write(*,*) 'comp_section: HAVE NOT FOUND THE ELEMENT FOR ',
     .           'x_moved(1:3,',idir,') = ', x_moved(1:3,idir)
      write(*,*) '              PERTURB THE POINT',Xpoint
      read(*,*) Xpoint(1:3)
      go to 6
 11   continue
      call find_domain(mdle, ndom)
c
c  ...skip PML elements
ccc      if(ndom.eq.6)then
ccc        write(*,*) 'comp_section: POINT INSIDE PML DOMAIN'
ccc        write(*,*) '              RESET THE POINT',Xpoint
ccc        read(*,*) Xpoint(1:3)
ccc        go to 6
ccc      endif
c
c  ...compute value at the point
c
c  ...determine order of approximation
      call find_order(mdle, norder)
c
c  ...determine the node orientation
      call find_orient(mdle, nedge_orient,nface_orient)
c
c  ...determine nodes coordinates
      call nodcor(mdle, xnod)
c
c  ...determine solution dof
      call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
c
c  ...determine physical coordinates at xi_moved(1:3,idir)
      call invmap(NODES(mdle)%type,Xpoint(1:3),
     .    norder,nedge_orient,nface_orient,
     .    xnod, idec,xi(1:3,1))
      if (idec.ne.1) then
        write(*,*) 'compute_section: INCONSISTENCY'
        call pause
      endif
c
c  ...determine physical coordinates and solution dof for two points
      call soleval(NODES(mdle)%type,mdle,xi(1:3,1),
     .             nedge_orient,nface_orient,norder,
     .             xnod,zdofH,
     .             1, x(1:3,idir),dxdxi,zsolH,zgradH)
c
      work2(1,0) = 0.d0
c
      select case(NODES(mdle)%case)
c
c  ...elasticity
      case(2)
c
c  .....dilatation (divergence)
        zdil = zgradH(1,1) + zgradH(2,2) + zgradH(3,3)
c
c  .....get material data
        call getmat_elast(mdle, omeg,zEYoung,rni, rhos,ufactor)
        zmu  = zEYoung/2.d0/(1.d0+rni)
        zla  = zEYoung*rni/(1.d0+rni)/(1.d0-2.d0*rni)
        zpressure = - (2.d0/3.d0*zmu+zla)*zdil*ufactor
        ztraction = (2.d0*zmu*zgradH(1,1)+zla*zdil)*ufactor
        if (x(1,idir).gt.0.d0) then
          zuradial = zsolH(1)*ufactor
        else
          zuradial = -zsolH(1)*ufactor
        endif
c
c  ...acoustics
      case(1)
c
        call getmat_acous(mdle, ifPML,omeg,csound,rhof,ufactor)
        call incident_pressure(x(1:3,idir),zpinc,zdpinc)
        zpressure = zpinc(1) + zsolH(4)
        ztraction = -zpressure
        zpressure_grad(1:3) = zdpinc(1:3,1) + zgradH(4,1:3)
        if (x(1,idir).gt.0.d0) then
          zuradial = zpressure_grad(1)/rhof/omeg**2
        else
          zuradial = -zpressure_grad(1)/rhof/omeg**2
        endif
c
      case default
        write(*,*) 'comp_section: iel, ELEMS(iel)%physics = ',
     .                            iel, ELEMS(iel)%physics
        stop 1
      end select
c
      work2(2,0)= dreal_part(zpressure)
      work2(3,0)= dimag_part(zpressure)
      work2(4,0)= abs(zpressure)
      work2(5,0)= dreal_part(zuradial)
      work2(6,0)= dimag_part(zuradial)
      work2(7,0)= abs(zuradial)
      work2(8,0) = dreal_part(ztraction)
      work2(9,0) = dimag_part(ztraction)
      work2(10,0) = abs(ztraction)
c
c-------------------------------------------------------------------
c
      mdle=0
c
c  ...0-th element
      j=0
c
c  ...two endpoints for the first section
 5    x_moved(1:3,1) = Xpoint(1:3) - step * Vector(1:3)
      x_moved(1:3,2) = Xpoint(1:3) + step * Vector(1:3)
c
      do idir=1,2
c
c  .....look for the element containg x_moved(1:3,idir)
        if(iprint.eq.1)then
          write(*,*)'************************************'
          write(*,*)'compute_section:looking for element with ',
     .      x_moved(1:3,idir)
        endif
        mdle=0
        do iel=1,NRELES
          call nelcon(mdle, mdle)
          call find_xi(mdle,x_moved(1:3,idir),Vector,
     .      idec,xi(1:3,1),xi(1:3,2))
          if (idec.eq.1) go to 10
        enddo
        write(*,*) 'comp_section: HAVE NOT FOUND THE ELEMENT FOR ',
     .             'x_moved(1:3,',idir,') = ', x_moved(1:3,idir)
        write(*,*) '              PERTURB THE POINT',Xpoint
        read(*,*) Xpoint(1:3)
        go to 5
 10     continue
        call find_domain(mdle, ndom)
c
c  .....skip PML elements
ccc        if(ndom.eq.6)then
ccc          write(*,*) 'comp_section: POINT INSIDE PML DOMAIN'
ccc          write(*,*) '              RESET THE POINT',Xpoint
ccc          read(*,*) Xpoint(1:3)
ccc          go to 5
ccc        endif
c
        if (iprint.eq.1) then
          write(*,7002) mdle
 7002     format('comp_section HAVE FOUND mdle = ',i7)
          write(*,*)'************************************'
          call pause
        endif
c
c  .....determine order of approximation
        call find_order(mdle, norder)
c
c  .....determine the node orientation
        call find_orient(mdle, nedge_orient,nface_orient)
c
c  .....determine nodes coordinates
        call nodcor(mdle, xnod)
c
c  .....determine physical coordinates at xi_moved(1:3,idir)
        call invmap(NODES(mdle)%type,x_moved(1:3,idir),
     .    norder,nedge_orient,nface_orient,
     .    xnod, idec,xi_moved(1:3,idir))
        if (idec.ne.1) then
          write(*,*) 'compute_section: INCONSISTENCY'
          call pause
        endif
c
c  .....determine solution dof
        call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
c
c  .....determine physical coordinates and solution dof for two points
        call soleval(NODES(mdle)%type,mdle,xi_moved(1:3,idir),
     .               nedge_orient,nface_orient,norder,
     .               xnod,zdofH,
     .               1, x(1:3,idir),dxdxi,zsolH,zgradH)
c
c  .....compute the line coordinate
        call norm(x(1:3,idir)-Xpoint(1:3), d)
        work(1,idir,j) = d*(-1.d0)**idir
c
        select case(NODES(mdle)%case)
c
c  .....elasticity
        case(2)
c
c  .......dilatation (divergence)
          zdil = zgradH(1,1) + zgradH(2,2) + zgradH(3,3)
c
c  .......get material data
          call getmat_elast(mdle, omeg,zEYoung,rni, rhof,ufactor)
          zmu  = zEYoung/2.d0/(1.d0+rni)
          zla  = zEYoung*rni/(1.d0+rni)/(1.d0-2.d0*rni)
          zpressure = - (2.d0/3.d0*zmu+zla)*zdil*ufactor
          ztraction = (2.d0*zmu*zgradH(1,1)+zla*zdil)*ufactor
          if (x(1,idir).gt.0.d0) then
            zuradial = zsolH(1)*ufactor
          else
            zuradial = -zsolH(1)*ufactor
          endif
          select case(ichoose)
          case(1)
            work(2,idir,j) = real(zpressure)
          case(2)
            work(2,idir,j) = dimag_part(zpressure)
          case(3)
            work(2,idir,j) = abs(zpressure)
          case(4)
            work(2,idir,j) = 20.d0*log10(
     .                      max(abs(zpressure),10.d0**(-5.d0/2.d0)))
          end select
c
c  .....acoustics
        case(1)
c
          call getmat_acous(mdle, ifPML,omeg,csound,rhof,ufactor)
          call incident_pressure(x(1:3,idir),zpinc,zdpinc)
          zpressure = zpinc(1) + zsolH(4)
          ztraction = -zpressure
          zpressure_grad(1:3) = zdpinc(1:3,1) + zgradH(4,1:3)
          if (x(1,idir).gt.0.d0) then
            zuradial = zpressure_grad(1)/rhof/omeg**2
          else
            zuradial = -zpressure_grad(1)/rhof/omeg**2
          endif
          select case(ichoose)
          case(1)
            work(2,idir,j) = real(zpressure)
            work(2,idir,j) = real(zpressure-zpinc(1))
          case(2)
            work(2,idir,j) = dimag_part(zpressure)
          case(3)
            work(2,idir,j) = abs(zpressure)
          case(4)
            work(2,idir,j) = 20.d0*log10(
     .                       max(abs(zpressure),10.d0**(-5.d0/2.d0)))
          end select
c
        case default
          write(*,*) 'comp_section: iel, ELEMS(iel)%physics = ',
     .                              iel, ELEMS(iel)%physics
          stop 1
        end select
c
        select case (idir)
        case(1)
          iplace=-1
          work2(1,iplace)= -d
        case(2)
          iplace=1
          work2(1,iplace)= d
        end select
        work2(2,iplace)= real(zpressure)
        work2(3,iplace)= dimag_part(zpressure)
        work2(4,iplace)= abs(zpressure)
        work2(5,iplace)= real(zuradial)
        work2(6,iplace)= dimag_part(zuradial)
        work2(7,iplace)= abs(zuradial)
        work2(8,iplace) = real(ztraction)
        work2(9,iplace) = dimag_part(ztraction)
        work2(10,iplace) = abs(ztraction)
c
c  ...end of loop through endpoints
      enddo

      if(iprint.eq.1)then
        write(*,*)'compute_section:found:',work(1,1:2,0),work(2,1:2,0)
        call pause
      endif

c  ...AT THIS POINT:
c       Xpoint (GCS) is the point inside the first element mdle
c       xi_moved (LCS) contains two shifted points
c       x_moved (GCS) contains two shifted points
c       x = x_moved
c       mdle contains middle point for the first element
c
c  ...store the original points (first moved points)
      x_org(1:3,2) = x(1:3,2)

      if (iprint.eq.1) then
        write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        write(*,*)'compute_section:LOOP WRT DIRECTIONS, step:',step
        call pause
      endif

c  ...loop wrt the directions
      do idir=1,2
c
        if (iprint.eq.1) then
          write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
          write(*,*)'compute_section:DIRECTION',idir
          call pause
        endif
c
c  .....take the original points
        x(1:3,1:2) = x_org(1:3,1:2)
c
c  .....reset the element number
        j=0
c
        do istep=1,100000
c
          select case(idir)
          case(1)
            j=j-1
            eps = -step
          case(2)
            j=j+1
            eps = step
          end select
c
          x_moved(1:3,idir)=x(1:3,idir) + eps * Vector(1:3)
c
c  .......look for the element containg x_moved(:,idir)
          if (iprint.eq.1) then
            write(*,*)'************************************'
            write(*,*)'compute_section:looking for element with ',
     .        'moved point ',x_moved(1:3,idir)
          endif
          mdle=0
          do iel=1,NRELES
            call nelcon(mdle, mdle)
            call find_xi(mdle,x_moved(1:3,idir),Vector,
     .                   idec,xi_new(1:3,1),xi_new(1:3,2))
            if (idec.eq.1) go to 20
          enddo
          write(*,*) 'comp_section: HAVE NOT FOUND THE ELEMENT FOR ',
     .               'xi_moved(:,',idir,')= ', x_moved(1:3,idir)
          if(iprint.eq.1)call pause
          go to 100 !exit the loop
 20       if (iprint.eq.1) then
            write(*,7002) mdle
            call pause
          endif
          call find_domain(mdle, ndom)
c
c  .......skip PML elements
ccc          if(ndom.eq.6)go to 100 !exit the loop

c   ......AT THIS POINT
c         x(1:3,idir) (GCS) contains previous points
c         x_moved(1:3.idir) (GCS) contains moved points
c         xi_new (LCS) contains crosssections with element faces
c
c  .......determine order of approximation
          call find_order(mdle, norder)
c
c  .......determine the node orientation
          call find_orient(mdle, nedge_orient,nface_orient)
c
c  .......determine nodes coordinates
          call nodcor(mdle, xnod)
c
c  .......determine physical coordinates at xi_moved(1:3,idir)
 90       continue
          call invmap(NODES(mdle)%type,x_moved(1:3,idir),
     .                norder,nedge_orient,nface_orient,
     .                xnod, idec,xi_moved(1:3,idir))
          if (idec.ne.1) then
            write(*,*) 'compute_section: INCONSISTENCY'
            call pause
          endif
c
c  .......determine solution dof
          call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
c
c  .......determine physical coordinates and solution dof for two points
          iprint_soleval=1
          call soleval(NODES(mdle)%type,mdle,xi_moved(1:3,idir),
     .                 nedge_orient,nface_orient,norder,
     .                 xnod,zdofH,1,
     .                 x_moved2(1:3,idir),dxdxi,zsolH,zgradH)
c
c   ......checkpoint
          call norm(x_moved(1:3,idir)-x_moved2(1:3,idir),d)
          if (d.gt.GEOM_TOL) then
            write(*,*)'compute_section:INCONSISTENCY!'
            write(*,*)'x_moved',x_moved(1:3,idir)
            write(*,*)'x_moved2',x_moved2(1:3,idir)
            iprint=1
            iprint_invmap=iprint
            iprint_soleval=iprint
            go to 90
         endif
c AT THIS POINT
c         x (GCS) contains previous points
c         x_moved (GCS) contains ctual moved points
c         xi_moved (LCS) as above
c         mdle contains actual element midle node

c  .......store x_moved at x
          x(1:3,idir)=x_moved(1:3,idir)

c  .......compute the line coordinate
          call norm(x_moved(1:3,idir)-Xpoint(1:3), d)
          select case(idir)
          case(1)
            work(1,1,j) = d*(-1.d0)
            work(1,2,j) = work(1,1,j+1)
          case(2)
            work(1,2,j) = d
            work(1,1,j) = work(1,2,j-1)
          end select
c
          select case(idir)
          case(1)
            work(2,2,j) = work(2,1,j+1)
          case(2)
            work(2,1,j) = work(2,2,j-1)
          end select

c
          select case(NODES(mdle)%case)
c
c  .......elasticity
          case(2)
c
c  .........dilatation (divergence)
            zdil = zgradH(1,1) + zgradH(2,2) + zgradH(3,3)
            write(*,8050) zgradH(1,1),zgradH(2,2),zgradH(3,3)
 8050       format('zgradH(i,i) = ',3(2e12.5,2x))
c
c  .........get material data
            call getmat_elast(mdle, omeg,zEYoung,rni, rhos,ufactor)
            zmu  = zEYoung/2.d0/(1.d0+rni)
            zla  = zEYoung*rni/(1.d0+rni)/(1.d0-2.d0*rni)
            zpressure = - (2.d0/3.d0*zmu+zla)*zdil*ufactor
            ztraction = (2.d0*zmu*zgradH(1,1)+zla*zdil)*ufactor
            if (x(1,idir).gt.0.d0) then
              zuradial = zsolH(1)*ufactor
            else
              zuradial = -zsolH(1)*ufactor
            endif
c
            select case(ichoose)
            case(1)
              work(2,idir,j) = real(zpressure)
            case(2)
              work(2,idir,j) = dimag_part(zpressure)
            case(3)
              work(2,idir,j) = abs(zpressure)
            case(4)
              work(2,idir,j) = 20.d0*log10(
     .                      max(abs(zpressure),10.d0**(-5.d0/2.d0)))
            end select
c
c  .......acoustics
          case(1)
c
            call getmat_acous(mdle, ifPML,omeg,csound,rhof,ufactor)
            call incident_pressure(x(1:3,idir),zpinc,zdpinc)
            zpressure = zpinc(1) + zsolH(4)
            ztraction = -zpressure
            zpressure_grad(1:3) = zdpinc(1:3,1) + zgradH(4,1:3)
            if (x(1,idir).gt.0.d0) then
              zuradial = zpressure_grad(1)/rhof/omeg**2
            else
              zuradial = -zpressure_grad(1)/rhof/omeg**2
            endif
            select case(ichoose)
            case(1)
              work(2,idir,j) = real(zpressure)
              work(2,idir,j) = real(zpressure-zpinc(1))
            case(2)
              work(2,idir,j) = dimag_part(zpressure)
            case(3)
              work(2,idir,j) = abs(zpressure)
            case(4)
              work(2,idir,j) = 20.d0*log10(
     .                      max(abs(zpressure),10.d0**(-5.d0/2.d0)))
            end select
c
          case default
            write(*,*) 'comp_section: iel, ELEMS(iel)%physics = ',
     .                              iel, ELEMS(iel)%physics
            stop 1
          end select
ccc          write(*,*) 'j,x_moved2(1:3,idir),dreal_part(zpressure) = ',
ccc     .                j,x_moved2(1:3,idir),dreal_part(zpressure)



c
c  .......save number of elements in the direction
          jmax(idir) = abs(j)
c
          select case (idir)
          case(1)
            iplace=j-1
            work2(1,iplace)= -d
          case(2)
            iplace=j+1
            work2(1,iplace)= d
          end select
          work2(2,iplace)= dreal_part(zpressure)
          work2(3,iplace)= dimag_part(zpressure)
          work2(4,iplace)= abs(zpressure)
          work2(5,iplace)= dreal_part(zuradial)
          work2(6,iplace)= dimag_part(zuradial)
          work2(7,iplace)= abs(zuradial)
          work2(8,iplace) = dreal_part(ztraction)
          work2(9,iplace) = dimag_part(ztraction)
          work2(10,iplace) = abs(ztraction)
c
c  .....end of loop wrt moved points
        enddo
 100    continue
        call pause
c
c  ...end of loop wrt the directions
      enddo
c
c  ...total number of segments to draw is...
      Nrsegm = 1 + jmax(1) + jmax(2)
      if (Nrsegm.gt.Na) then
        write(*,*) 'comp_section: INCREASE Na, Nrsegm = ', Nrsegm
        return
      endif
c
c  ...reshuffle the data
      l=-jmax(1)-1
      do k=1,Nrsegm
        l=l+1
        Array(1,k) = work(1,1,l)
        Array(2,k) = work(1,2,l)
        Array(3,k) = work(2,1,l)
        Array(4,k) = work(2,2,l)
      enddo
c
      do j=-jmax(1)-1,jmax(2)+1
        write(NOUT,7004) work2(1:10,j)
 7004   format(f9.5,9e14.5)
      enddo
c
      deallocate(work,work2)
      write(*,*)'number of points',Nrsegm
c
      if (iprint.ge.1) then
        write(*,*) 'comp_section: DATA TO DRAW'
        call pause
        do k=1,Nrsegm
          write(*,7003) k,Array(1,k),Array(2,k),Array(3,k),Array(4,k)
 7003     format('   k,x1,x2,val1,val2 = ',i5,2x,2f8.3,2x,2e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine compute_section
c
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - find_xi
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 09
c
c   purpose            - given physical coordinates Xp of a point in
c                        space and geometry dof for an element,
c                        routine:
c                         1/ determines whether the point lies within
c                            the element;
c                        for a point within the element, routine
c                        determines then
c                         2/ master element coordinates of two
c                            intersection points of the line and
c                            the element boundary
c
c   arguments :
c     in:
c             Mdle     - element (middle node) number
c             Xp       - physical coordinates of a point
c             Vect     - a versor
c     out:
c             Idec     = 1 if the point is within the element
c                        0 otherwise
c             Xip1,Xip2  - master coordinates of the intersection
c                        points
c
c
c----------------------------------------------------------------------
c
      subroutine find_xi(Mdle,Xp,Vect,
     .                   Idec,Xip1,Xip2)
c
      use data_structure3D
      use control
#include "syscom.blk"
c
      dimension Xp(3), Vect(3)
c
c  ...element order, geometry dof
      dimension norder(19),xnod(3,MAXbrickH)
c
c  ...node orientations
      dimension nedge_orient(12),nface_orient(6)
c
      dimension Xip1(3),Xip2(3)
c
c  ...geometry
      dimension xi(3)
c
c  ...increment in xi
      dimension dxi(3)
c
c  ...work space
      dimension xw(3),xbeg(3),xend(3)
c
c
c-----------------------------------------------------------------------
c
      iprint=0 !2 for bisection
c
c  ...determine order of approximation
      call find_order(Mdle, norder)
c
c  ...determine the node orientation
      call find_orient(Mdle, nedge_orient,nface_orient)
c
c  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
c
c  ...check if the point is within the element
c      write(*,*)'call invmap:NODES(',Mdle,')%type',NODES(Mdle)%type
      call invmap(NODES(Mdle)%type,Xp,norder,nedge_orient,nface_orient,
     .            xnod, Idec,xi)
      if (Idec.eq.0) return
c
c  ...determine the element size(roughly)
      dmax = 0.d0
      do iv1=1,nvert(NODES(Mdle)%type)
        do iv2=1,nvert(NODES(Mdle)%type)
          call norm(xnod(1:3,iv2)-xnod(1:3,iv1), dist)
          dmax = max(dmax,dist)
        enddo
      enddo
c
c  ...find the entry and exit intersection points
      do iside=1,2
        select case(iside)
        case(1)
          xw(1:3) = Xp(1:3) - Vect(1:3)*dmax
        case(2)
          xw(1:3) = Xp(1:3) + Vect(1:3)*dmax
        end select
        call invmap(NODES(Mdle)%type,xw,norder,nedge_orient,
     .              nface_orient,xnod, idec1,xi)
        do while(idec1.eq.1)
          select case(iside)
          case(1)
            xw(1:3) = xw(1:3) - Vect(1:3)*dmax
          case(2)
            xw(1:3) = xw(1:3) + Vect(1:3)*dmax
          end select
          call invmap(NODES(Mdle)%type,xw,norder,nedge_orient,
     .                nface_orient,xnod, idec1,xi)
        enddo
        if(iprint.eq.2)then
          write(*,*)'find_xi: ',
     .      'found point in direction ',iside,' outside element',xw
          call pause
        endif
        xbeg(1:3) = xw(1:3); xend(1:3) = Xp(1:3)
        call norm(xend(1:3)-xbeg(1:3), d)
        do while(d.gt.GEOM_TOL)
          xw(1:3) = (xbeg(1:3)+xend(1:3))*.5d0
          if(iprint.eq.2)then
            write(*,*)'find_xi: bisection:xbeg',xbeg
            write(*,*)'find_xi: bisection:xend',xend
            write(*,*)'find_xi: bisection:xw',xw
            call pause
          endif
          call invmap(NODES(Mdle)%type,xw,norder,nedge_orient,
     .              nface_orient,xnod, idec1,xi)
          select case(idec1)
          case(0)
            xend(1:3) = xw(1:3)
            if(iprint.eq.2)write(*,*)'find_xi: bisection:xend=xw'
          case(1)
            xbeg(1:3) = xw(1:3)
            if(iprint.eq.2)write(*,*)'find_xi: bisection:xbeg=xw'
          end select
          call norm(xend(1:3)-xbeg(1:3), d)
        enddo
        select case(iside)
        case(1)
c  .......store entry point
          Xip1(1:3) = xw(1:3)
        case(2)
c  .......store exit point
          Xip2(1:3) = xw(1:3)
        end select
      enddo

      if(iprint.eq.1)then
        write(*,*)'find_xi: ',
     .    'point inside the element',Xp(1:3)
        write(*,*)'find_xi: ',
     .    'found two points on element faces',Xip1(1:3),Xip2(1:3)
        call pause
      endif
c
c
      end subroutine find_xi
