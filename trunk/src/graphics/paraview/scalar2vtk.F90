!----------------------------------------------------------------------------------------
!> Purpose : write scalar attribute to .h5 file
!!
!> @param[in ] Stype   - "Scalar" 
!> @param[in ] Sname   - name for attribute's .h5 file
!> @param[in ] Snick   - attribute's nickname to be used in .h5 file
!> @param[in ] Scenter - "Node"
!> @param[in ] Domain  - number of subdomain to be displayed (currently inactive)
!> @param[in ] Idx     - index identifying scalar attribute
!> @param[out] Ic      - number of vertices of visualization object
!!
!> @date Apr 2019
!----------------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine scalar2vtk(Stype,Sname,Snick,Scenter,Domain,Idx, Ic)
!
      use data_structure3D , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,NODES,NRELES, &
                                    MAXbrickH,MAXbrickE,MAXbrickV,MAXbrickQ,      &
                                    NRHVAR,NRQVAR
      use element_data
      use physics          , only : DTYPE,ADRES
      use upscale
!
      implicit none
!
      character(len=*), intent(in ) :: Stype
      character(len=*), intent(in ) :: Sname
      character(len=*), intent(in ) :: Snick
      character(len=*), intent(in ) :: Scenter
      integer,          intent(in ) :: Domain
      integer,          intent(in ) :: Idx
      integer,          intent(out) :: Ic
!
      type(vis)        :: vis_obj
      character(len=4) :: etype
      integer                             :: mdle, nflag, i, iv, nV
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
      real*8 :: val
      integer :: ibeg,iattr,icomp,isol,iload,ireal,ndom
!
      real*8, external :: dreal_part,dimag_part
!
!   ..OpenMP parallelization: auxiliary variables
      integer, dimension(NRELES) :: n_elem, n_vert_offset, n_elem_vert
      real*8, allocatable :: val_array(:)
!
!----------------------------------------------------------------------------------------
!
!     Step 0 : Clear vis scalar
      call vis_scalar_clear
!
!     Step 1 : Write attribute of interest
      Ic=0
!      
!     decode
      ireal = iabs(Idx)/Idx
      iload = iabs(Idx)/100
      iattr = iabs(Idx) - iload*100 ; iattr=iattr/10
      icomp = iabs(Idx) - iload*100 - iattr*10
!
!     address of 1st component for the attribute
      ibeg=ADRES(iattr)
!
!   ..create list of mdle node numbers,
!     and count number of visualization points
      mdle=0
      do i=1,NRELES
        call nelcon(mdle, mdle)
        n_elem(i) = mdle
        if (Domain.ne.0) then
          call find_domain(mdle, ndom)
          if (ndom.ne.Domain) cycle
        endif
        etype = NODES(mdle)%type
        vis_obj = vis_on_type(etype)
        n_vert_offset(i) = Ic
        n_elem_vert(i) = vis_obj%nr_vert
        Ic = Ic + vis_obj%nr_vert
      enddo
      allocate(val_array(Ic))
!
      nflag=1
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(mdle,ndom,etype,iv,nV,xi,               &
!$OMP         xnod,zdofH,zdofE,zdofV,zdofQ,           &
!$OMP         norder,nedge_orient,nface_orient,       &
!$OMP         x,dxdxi,zsolH,zgradH,zsolE,zcurlE,      &
!$OMP         zsolV,zdivV,zsolQ,isol,val)             &
!$OMP SCHEDULE(DYNAMIC)
      do i=1,NRELES
        mdle = n_elem(i)
!
!       SKIP IF WRONG DOMAIN
        if (Domain.ne.0) then
          call find_domain(mdle, ndom)
          if (ndom.ne.Domain) cycle
        endif
!
!       element info for computing solution
        call nodcor(         mdle, xnod)
        call solelm(         mdle, zdofH,zdofE,zdofV,zdofQ)
        call find_elem_nodes(mdle, norder,nedge_orient,nface_orient)
!
!       select appropriate visualization object  
        etype = NODES(mdle)%type
        nV = n_elem_vert(i)
!
!       loop over vertices of visualization object
        do iv=1,nV
!        
!         visualization point accounting for VLEVEL
          call get_vis_point(vis_on_type(etype),iv-1, xi)
!
!         compute element solution
          call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                       zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi,         &
                       zsolH,zgradH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
!
!         approximation space
          select case(DTYPE(iattr))
!
!         -- H1 --
          case('contin')
            isol = (iload-1)*NRHVAR + ibeg + icomp
!
!           REAL part
            if (ireal == 1) then
              val = dreal_part(zsolH(isol))
!
!           IMAGINARY part
            else
              val = dimag_part(zsolH(isol))
            endif
!
!         -- L2 --
          case('discon')
            isol = (iload-1)*NRQVAR + ibeg + icomp
!
!           REAL part
            if (ireal == 1) then
              val = dreal_part(zsolQ(isol))
!
!           IMAGINARY part
            else
              val = dimag_part(zsolQ(isol))
            endif  
          endselect
!
!         add value to visualization object
          val_array(n_vert_offset(i)+iv) = val
        enddo
!         
      enddo
!$OMP END PARALLEL DO
!
      do i=1,Ic
        call vis_scalar_add_value(val_array(i));
      enddo
      deallocate(val_array)
!
!     Step 2 : write to file
      call vis_scalar_write(Stype,len(Stype),Sname,  len(Sname), & 
                            Snick,len(Snick),Scenter,len(Scenter))
!
!     Step 3 : Clear 
      call vis_scalar_clear
!
!
end subroutine scalar2vtk
