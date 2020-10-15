!----------------------------------------------------------------------------------------
!> Purpose : write scalar attribute to .h5 file
!!
!> @param[in ] Stype   - "Scalar"
!> @param[in ] Sname   - name for attribute's .h5 file
!> @param[in ] Snick   - attribute's nickname to be used in .h5 file
!> @param[in ] Scenter - "Node"
!> @param[in ] Domain    - number of subdomain to be dispalyed
!                         0 = domain 1 and 2
!                         1 = domain 1
!                         2 = domain 2
!> @param[in ] Scomp   - component of stress
!                         1 = sigma_rr
!                         2 = sigma_rt
!                         3 = sigma_tt
!                         4 = displacement
!> @param[out] Ic      - number of vertices of visualization object
!!
!> @date Mar 16
!----------------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine soln2vtk(Stype, Sname, Snick, Scenter, Domain, Scomp, Ic)
!
      use data_structure3D , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ , NODES,NRELES , &
                                    MAXbrickH,MAXbrickE,MAXbrickV,MAXbrickQ ,        &
                                    NRHVAR,NRQVAR
      use element_data
      use physics          , only : DTYPE,ADRES
      use upscale
      use isotropic_elast_material
!
      implicit none
      character(len=*), intent(in ) :: Stype
      character(len=*), intent(in ) :: Sname
      character(len=*), intent(in ) :: Snick
      character(len=*), intent(in ) :: Scenter
      integer,          intent(in ) :: Domain
      integer,          intent(in ) :: Scomp
      integer,          intent(out) :: Ic
!
      type(vis)      :: vis_obj
      character(len=4) :: type
      integer                             :: mdle, nflag, iel, iv, ndom
      real*8,  dimension(3)               :: xi,x
      integer, dimension(12)              :: nedge_orient
      integer, dimension(6)               :: nface_orient
      integer, dimension(19)              :: norder
      real*8,  dimension(3,MAXbrickH)     :: xnod
      VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
      VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
      VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
      VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
      real*8, dimension(3,3)              :: dxdxi
      VTYPE, dimension(  MAXEQNH  )       :: zsolH
      VTYPE, dimension(  MAXEQNH,3)       :: zgradH
      VTYPE, dimension(3,MAXEQNE  )       :: zsolE
      VTYPE, dimension(3,MAXEQNE  )       :: zcurlE
      VTYPE, dimension(3,MAXEQNV  )       :: zsolV
      VTYPE, dimension(  MAXEQNV  )       :: zdivV
      VTYPE, dimension(  MAXEQNQ  )       :: zsolQ
!
!  ...stiffness tensor
      real*8, dimension(3,3,3,3) :: C
!
!  ...stress and displacement
      real*8, dimension(3,3) :: sigma
      real*8, dimension(3)   :: u
!
!  ...miscellaneous
      real*8, dimension(3) :: val
      real*8  :: r2
      integer :: ibeg,iattr,icomp,isol,iload,ireal,i,j,m,n
!
!----------------------------------------------------------------------------------------
!
!     Step 0 : Clear vis
      if (Stype.eq."Scalar") call vis_scalar_clear
      if (Stype.eq."Vector") call vis_vector_clear
!
!     Step 1 : Write attribute of interest
      Ic=0
!
!     address of 1st component for the attribute
      ibeg=ADRES(1)
!
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
!
!       element info for computing solution
        call nodcor(         mdle, xnod)
        call solelm(         mdle, zdofH,zdofE,zdofV,zdofQ)
        call find_elem_nodes(mdle, norder,nedge_orient,nface_orient)
        call find_domain(    mdle, ndom)
!
!       SKIP IF WRONG DOMAIN
        if (Domain.ne.0) then
          if (ndom.ne.Domain) cycle
        endif
!
!       select appropriate visualization object
        type = NODES(mdle)%type
        vis_obj = vis_on_type(type)
!
!       loop over vertices of visualization object
        do iv=1,vis_obj%nr_vert
!
!         visualization point accounting for VLEVEL
          call get_vis_point(vis_obj, iv-1, xi)
!
!         compute element solution
          nflag=1
          call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                       x,dxdxi,zsolH,zgradH,zsolE,zcurlE,zsolV,zdivV,zsolQ                           )
!
!         compute radius^2
          r2 = x(2)**2+x(3)**2
!
!         compute stress and displacement in Cartesian coordinates
          ! if (ndom.eq.1) then
            !  DISPLACEMENT
            u = zsolH
            !  STRESS
            call getC(x, C)
            ! call getC(x,ndom, C)
            sigma=0.d0
            do i=1,3; do j=1,3; do m=1,3; do n=1,3
            sigma(i,j) = sigma(i,j)  &
                       + C(i,j,m,n)*zgradH(m,n)
            enddo; enddo; enddo; enddo
          ! else
          !   !  DISPLACEMENT
          !   u = zsolQ(1:3)
          !   !  STRESS
          !   !          ( sigma1   sigma4  sigma5 )
          !   !  sigma = ( sigma4   sigma2  sigma6 )
          !   !          ( sigma5   sigma6  sigma3 )
          !   sigma(1,1) = zsolQ(4)
          !   sigma(2,2) = zsolQ(5)
          !   sigma(3,3) = zsolQ(6)
          !   sigma(1,2) = zsolQ(7); sigma(2,1) = zsolQ(7) 
          !   sigma(1,3) = zsolQ(8); sigma(3,1) = zsolQ(8)
          !   sigma(2,3) = zsolQ(9); sigma(3,2) = zsolQ(9)
          ! endif
!
          select case(Scomp)
!
!         -- sigma_rr --
          case(1)
            val(1) =   sigma(2,2)*x(2)**2  &
                   + 2*sigma(2,3)*x(2)*x(3)  &
                   +   sigma(3,3)*x(3)**2
            val(1) = val(1)/r2
!
!         -- sigma_rt --
          case(2)
            val(1) = (sigma(3,3)-sigma(2,2))*x(2)*x(3)  &
                   + sigma(2,3)*(x(2)**2-x(3)**2)
            val(1) = val(1)/r2
!
!         -- sigma_tt --
          case(3)
            val(1) =   sigma(2,2)*x(3)**2  &
                   - 2*sigma(2,3)*x(2)*x(3)  &
                   +   sigma(3,3)*x(2)**2
            val(1) = val(1)/r2
!         -- displacement u --
          case(4)
            val = u
!
          endselect
!
!         add value to visualization object
          if (Stype.eq."Scalar") call vis_scalar_add_value(val(1))
          if (Stype.eq."Vector") call vis_vector_add_value(val)
!
        enddo
!
        Ic = Ic + vis_obj%nr_vert
!
      enddo
!
!     Step 2 : write to file and clear
      if (Stype.eq."Scalar") then
        call vis_scalar_write(Stype,len(Stype),Sname,  len(Sname), &
                              Snick,len(Snick),Scenter,len(Scenter))
        call vis_scalar_clear
      elseif (Stype.eq."Vector") then
        call vis_vector_write(Stype,len(Stype),Sname,  len(Sname), &
                              Snick,len(Snick),Scenter,len(Scenter))
        call vis_vector_clear
      endif
!
!
endsubroutine soln2vtk
