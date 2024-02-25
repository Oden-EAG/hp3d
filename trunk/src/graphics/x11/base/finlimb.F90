#if HP3D_USE_X11

!> Purpose : routine finds 'numlev'- levels of solution values
!! @param[in]  Numlev - number of levels to plot solution values
!! @param[out] Solev  - limiting values for each level
#include "typedefs.h"
  subroutine finlimb(Numlev, Solev)
    use error
    use data_structure3D
    use graphmod
    implicit none
    ! ** Arguments
    !----------------------------------------
    integer,                          intent(in)  :: Numlev
    real(8), dimension(NR_COLORS-10), intent(out) :: Solev
    ! ** Locals
    !----------------------------------------
    integer :: mdle, ndom, nedge_orient(12),nface_orient(6),norder(19)
    integer :: nrdofH, nrdofE, nrdofV, nrdofQ
    ! geometry dof and master coordinates
    real(8), dimension(2) :: t
    real(8), dimension(3) :: xp
    real(8)               :: xnod(NDIMEN,MAXbrickH)
    integer               :: ntype
    VTYPE :: &
         zdofH(MAXEQNH,MAXbrickH), &
         zdofE(MAXEQNE,MAXbrickE), &
         zdofV(MAXEQNV,MAXbrickV), &
         zdofQ(MAXEQNQ,MAXbrickQ)

    real(8) :: dsol, dxi, solmax, solmin, val
    integer :: i, j, ivar, loc, iel, idec, iface, nsub
    !
#if HP3D_DEBUG
    integer :: iprint
    iprint=0
#endif
    !
    !----------------------------------------
    !
    t = 0.d0; xp = 0.d0; val = 0.d0
    !
#if HP3D_DEBUG
    if (iprint.eq.1) then
       write(*,*) 'finlimb: Numlev = ',Numlev
       call pause
    endif
#endif
    ! increment in master element coordinates
    dxi = DX

    ! search for biggest and smallest value
    solmax = -1.d10
    solmin =  1.d10

    mdle=0
    do iel=1,NRELES
       call nelcon(mdle, mdle)
       ntype = NODES(mdle)%ntype

       ! if it is not visible domain hide it
       call find_domain(mdle, ndom)
       if (NDOMAIN(ndom).eq.0) then
          cycle
       endif

       ! if the element is hidden list, hide it
       call locate(mdle,IGINV,NRINVBL, loc)
       if (loc.gt.0) then
          cycle
       endif

       call find_orient(mdle, nedge_orient,nface_orient)
       call find_order(mdle, norder)
       call nodcor(mdle, xnod)
#if HP3D_DEBUG
       if (iprint.eq.1) then
          write(*,7002) mdle
7002      format('finlimb: VERTEX COORDINATES FOR mdle = ',i5)
          do ivar=1,3
             write(*,7003) xnod(ivar,1:nvert(NODES(mdle)%ntype))
7003         format(8(f8.5,2x))
          enddo
          call pause
       endif
#endif
       call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
       call celndof(ntype,norder, &
                    nrdofH,nrdofE,nrdofV,nrdofQ)

       ! loop through element faces
       do iface=1,nface(ntype)

          do j=0,NRSUB
             select case(face_type(ntype,iface))
             case(TRIA);   nsub=NRSUB-j
             case(RECT);   nsub=NRSUB
             end select

             do i=0,nsub
                t(1) = i*dxi
                t(2) = j*dxi
                call compute_face( &
                     Numlev, &
                     mdle,iface,nedge_orient,nface_orient,norder, &
                     xnod, &
                     zdofH,zdofE,zdofV,zdofQ, &
                     t, xp,val)
#if HP3D_DEBUG
                if (iprint.eq.1) then
                   write(*,7001) mdle,iface,i,j,xp,val
7001               format('finlimb: mdle,iface,i,j,xp,val = ', &
                        i5,3i2,3f8.3,2x,e12.5)
                   call pause
                endif
#endif
                ! update extremes
                solmax = max(solmax,val)
                solmin = min(solmin,val)
             enddo
          enddo
       enddo
    enddo

    write(*,7010) solmin,solmax
7010 format('finlimb: EXTREME VALUES = ',2e12.5)

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
    !
    ! divide the whole range into levels
    ! use the volume values
    dsol = (solmax-solmin)/float(Numlev)
    Solev(1) = solmin

    do i=1,Numlev
       Solev(i+1) = Solev(1) +float(i)*dsol
    enddo
  end subroutine finlimb

#endif
