c----------------------------------------------------------------------
c
c   routine name       - finlimb
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 10
c
c   purpose            - routine finds 'numlev'- levels of solution values
c
c   arguments :
c     in:
c               Numlev - number of levels to plot solution values
c     out:
c               Solev  - limiting values for each level 
c
c----------------------------------------------------------------------
c
      subroutine finlimb(Numlev, Solev)
c
      use error
      use data_structure3D
      use graphmod
#include "syscom.blk"
c
      dimension Solev(NR_COLORS-10)
      dimension nedge_orient(12),nface_orient(6),norder(19)
c
c  ...geometry and actual dof
      dimension xnod(NDIMEN,MAXbrickH)

      allocatable :: zdofH(:,:),zdofE(:,:),zdofV(:,:),zdofQ(:,:)
c
c  ...master face coordinates of a point
      dimension t(2)
c
c  ...master and physical coordinates of a point, solution value
      dimension xi(3),xpoint(3)
      allocatable :: zsolH(:)
c
c
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'finlimb: Numlev = ',Numlev
        call pause
      endif
      allocate(zdofH(MAXEQNH,MAXbrickH),
     .         zdofE(MAXEQNE,MAXbrickE),
     .         zdofV(MAXEQNV,MAXbrickV),
     .         zdofQ(MAXEQNQ,MAXbrickQ),
     .         zsolH(MAXEQNH), stat=istat)
      if (istat.ne.SUCCESS) then
        call logic_error(ERR_ALLOC_FAILURE,
     .    __FILE__,__LINE__)
      endif
      
c
c  ...increment in master element coordinates
      dxi = DX
c
c  ...in search for biggest and smallest value
      solmax =-1.d10
      solmin = 1.d10
c
c  ...loop through elements      
      mdle=0
      do iel=1,NRELES
c
        call nelcon(mdle, mdle)
        call locate(mdle,IGINV,NRINVBL, loc)
        if (loc.gt.0) cycle
c
        call find_domain(mdle, ndom)
        if (NDOMAIN(ndom).eq.0) cycle
c
        call find_orient(mdle, nedge_orient,nface_orient)
        call find_order(mdle, norder)
        call nodcor(mdle, xnod)
        if (iprint.eq.1) then
          write(*,7002) mdle
 7002     format('finlimb: VERTEX COORDINATES FOR mdle = ',i5)
          do ivar=1,3
            write(*,7003) xnod(ivar,1:nvert(NODES(mdle)%type))
 7003       format(8(f8.5,2x))
          enddo
          call pause
        endif
c
        call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
        call celndof(NODES(mdle)%type,norder, 
     .               nrdofH,nrdofE,nrdofV,nrdofQ)
c
c  .....loop through element faces
        do if=1,nface(NODES(mdle)%type)
c
c  .......loop through points on the face
          nsub1=NRSUB
          do j=0,nsub1
            select case(face_type(NODES(mdle)%type,if))
            case('tria'); nsub2=NRSUB-j
            case('rect'); nsub2=NRSUB
            end select
            do i=0,nsub2
              t(1) = i*dxi; t(2) = j*dxi
              call compute_face(Numlev,mdle,if,
     .                          nedge_orient,nface_orient,norder,
     .                          xnod,zdofH,t,
     .                          xpoint,val)
              if (iprint.eq.1) then
                write(*,7001) mdle,if,i,j,xpoint,val
 7001           format('finlimb: mdle,if,i,j,xpoint,val = ',
     .                  i5,3i2,3f8.3,2x,e12.5)
                call pause
              endif
c
c  ...........update extremes  
              solmax = max(solmax,val)
              solmin = min(solmin,val)
            enddo
          enddo
c
c  .....end of loop through faces
        enddo
c
c  ...end of loop through elements
      enddo
c
c
      write(*,7010) solmin,solmax
 7010 format('finlimb: EXTREME VALUES = ',2e12.5)
c
      write(*,*) 'DO YOU WANT TO CHANGE THE RANGE OF COLORS ?'
      write(*,*) '0....NO'
      write(*,*) '1....USE COMMON BOUND FOR NEG AND POS VALUES'
      write(*,*) '2....SET UP YOUR OWN BOUNDS'
      read(*,*) idec
      select case(idec)
      case(1)
        solmax = max(abs(solmax),abs(solmin))
        solmin = -solmax
        write(*,7010) solmin,solmax
      case(2)
        write(*,*) 'finlinmb: GIVE LOWER AND UPPER BOUND'
        read(*,*) solmin,solmax
      end select
c
c  ...divide the whole range into levels
c
c  ...use the volume values....
      dsol = (solmax-solmin)/float(Numlev)
      Solev(1) = solmin
c
      do i=1,Numlev
        Solev(i+1) = Solev(1) +float(i)*dsol
      enddo
c
      deallocate(zdofH,zdofE,zdofV,zdofQ,zsolH,
     .  stat=istat)
      if (istat.ne.SUCCESS) then
        call logic_error(ERR_ALLOC_FAILURE,
     .    __FILE__,__LINE__)
      endif
c
c
      end
