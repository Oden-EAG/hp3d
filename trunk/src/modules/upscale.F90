!----------------------------------------------------------------------
!  module upscale
!
!> @brief   define data point for upscaling of VTK output
!> @date    Feb 2023
!----------------------------------------------------------------------
module upscale
!
   use node_types
   use paraview,     only : SECOND_ORDER_VIS, VIS_VTU
!
   implicit none
!
   save
!
!  define generic type for a visualization object
   type :: vis
!
!  ...number of vertices required by visualization level (0 - 3)
      integer :: NR_VERT
!
!  ...list of vertices' coordinates
      real(8), pointer :: VERTC(:,:) => null()
!
!  ...number of elements required by visualization level (0 - 3)
      integer :: NR_ELEM
!
!  ...list of elements' vertices (indexing starts at 0)
      integer, pointer :: ELEM(:,:) => null()
!
   endtype vis
!
   type(vis) :: TETR_VIS, PRIS_VIS, HEXA_VIS


contains


!> @brief return cell type (XDMF 2 or VTU: for vis object depending upon VIS_VTU Flag)
!> @param[in] Etype - Element type
!> @date Feb 2023
   integer function ivis_type(Etype)
!
      integer :: Etype
!
      if (SECOND_ORDER_VIS) then
         if(VIS_VTU) then
            select case(Etype)
               case(TETR,MDLN); ivis_type = 24
               case(PRIS,MDLP); ivis_type = 32
               case(BRIC,MDLB); ivis_type = 29
               case default
                  write(*,*) 'ivis_type'; stop
            end select
         else
            select case(Etype)
               case(TETR,MDLN); ivis_type = 38
               case(PRIS,MDLP); ivis_type = 41
               case(BRIC,MDLB); ivis_type = 50
               case default
                  write(*,*) 'ivis_type'; stop
            end select
         endif
      else
         if(VIS_VTU) then
            select case(Etype)
               case(TETR,MDLN); ivis_type = 10
               case(PRIS,MDLP); ivis_type = 13
               case(BRIC,MDLB); ivis_type = 12
               case default
                  write(*,*) 'ivis_type'; stop
            end select
         else
            select case(Etype)
               case(TETR,MDLN); ivis_type = 6
               case(PRIS,MDLP); ivis_type = 8
               case(BRIC,MDLB); ivis_type = 9
               case default
                  write(*,*) 'ivis_type'; stop
            end select
         endif
      endif
!
   end function ivis_type



!> @brief return number of points to describe vis object
!> @param[in] Etype - Element type
!> @date Feb 2023
   integer function nobj_conf(Etype)
!
      integer :: Etype
!
      if (SECOND_ORDER_VIS) then
         select case(Etype)
            case(TETR,MDLN); nobj_conf = 10
            case(PRIS,MDLP); nobj_conf = 18
            case(BRIC,MDLB); nobj_conf = 27
            case default
               write(*,*) 'nobj_conf: unexpected type. stop.'
               stop
         end select
      else
         select case(Etype)
            case(TETR,MDLN); nobj_conf = 4
            case(PRIS,MDLP); nobj_conf = 6
            case(BRIC,MDLB); nobj_conf = 8
            case default
               write(*,*) 'nobj_conf: unexpected type. stop.'
               stop
         end select
      endif
!
   end function nobj_conf



!> @brief return preloaded vis element configuration
!> @param[in] Etype - Element type
!> @date Feb 2023
   type(vis) function vis_on_type(Etype)
!
      integer :: Etype
!
      select case(Etype)
         case(TETR,MDLN); vis_on_type = TETR_VIS
         case(PRIS,MDLP); vis_on_type = PRIS_VIS
         case(BRIC,MDLB); vis_on_type = HEXA_VIS
      end select
!
   end function vis_on_type



!> @brief deallocate vis object (configuration)
!> @param[in,out] V
!> @date Feb 2023
   subroutine clear_vis(V)
!
      type(vis), intent(inout) :: V
!
      V%NR_VERT = 0
      if (associated(V%VERTC)) deallocate(V%VERTC)
!
      V%NR_ELEM = 0
      if (associated(V%ELEM )) deallocate(V%ELEM)
!
   end subroutine clear_vis



!> @brief set up vis element configuration
!> @param[in,out] V
!> @param[in]     Fp
!> @param[in]     Etype - Element type
!> @date Feb 2023
   subroutine load_vis(V,Fp,Etype)
!
      type(vis)       , intent(inout) :: V
      character(len=*), intent(in)    :: Fp
      integer         , intent(in)    :: Etype
!
      integer, parameter :: nin = 29
      integer            :: i, istat, ierr
!
      ierr = 0
!
!  ...clear_vis nullifies vis%elem and vis%vert pointers
      call clear_vis(V)
!
      open(unit=nin, file=Fp,form='formatted', &
           access='sequential',status='unknown')
!
      read(nin,*) V%NR_VERT
      allocate(V%VERTC(0:V%NR_VERT, 3), STAT=istat)
      V%VERTC = 0.d0
      ierr = ierr + istat
!
      do i=0,(V%NR_VERT-1)
         read(nin,*) V%VERTC(i,1), V%VERTC(i,2), V%VERTC(i,3)
      end do
!
      read(nin,*) V%NR_ELEM
      if (SECOND_ORDER_VIS) then
         allocate(V%ELEM(1:V%NR_ELEM, 0:27), STAT=istat)
      else
         allocate(V%ELEM(1:V%NR_ELEM, 0:8), STAT=istat)
      endif
      V%ELEM = 0
      ierr = ierr + istat
!
!
      if (SECOND_ORDER_VIS) then
         read(nin,*) V%ELEM(1,0)
         do i=1,V%NR_VERT
            read(nin,*) V%ELEM(1,i)
         enddo
      else
         do i=1,V%NR_ELEM
            select case (Etype)
            case(TETR,MDLN)
               read(nin,*) &
                  V%ELEM(i,0), &
                  V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3), V%ELEM(i,4)
            case(PRIS,MDLP)
               read(nin,*) &
                  V%ELEM(i,0), &
                  V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3),&
                  V%ELEM(i,4), V%ELEM(i,5), V%ELEM(i,6)
            case(BRIC,MDLB)
               read(nin,*) &
                  V%ELEM(i,0), &
                  V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3), V%ELEM(i,4), &
                  V%ELEM(i,5), V%ELEM(i,6), V%ELEM(i,7), V%ELEM(i,8)
            end select
         enddo
      end if
!
      close(nin)
      ierr = ierr + istat
      if (ierr.ne.0) then
         write(*,*) 'module upscale::Something wrong', ierr
      end if
!
   end subroutine load_vis



!> @brief return point coordinates of a point associated with a vis object
!> @param[in]  V
!> @param[in]  Idx
!> @param[out] Pt
!> @date Feb 2023
   subroutine get_vis_point(V,Idx, Pt)
!
      type(vis), intent(in) :: V
      integer  , intent(in) :: Idx
      real(8)  , intent(out):: Pt(3)
!
      Pt(1:3) = V%VERTC(Idx,1:3)
!
   end subroutine get_vis_point



!> @brief return list of vertices for vis object
!> @param[in]  V
!> @param[in]  Idx
!> @param[in]  Ioffs
!> @param[out] Iverl
!> @date Feb 2023
   subroutine get_vis_elem(V,Idx,Ioffs, Iverl)
!
      type(vis), intent(in)  :: V
      integer  , intent(in)  :: Idx, Ioffs
      integer  , intent(out) :: Iverl(27)

      if (SECOND_ORDER_VIS) then
         Iverl(1:27) = V%ELEM(Idx,1:27) + Ioffs
      else
         Iverl(1:8) = V%ELEM(Idx,1:8) + Ioffs
      endif
!
   end subroutine get_vis_elem



!> @brief return number of vis elements of the vis object
!> @param[in]  Etype  - Element type
!> @param[out] Nrelem - Number of vis elements for the type
!> @date Feb 2023
   subroutine get_vis_nrelem(Etype, Nrelem)
!
      integer, intent(in)  :: Etype
      integer, intent(out) :: Nrelem
!
      select case(Etype)
         case(TETR,MDLN); Nrelem = TETR_VIS%NR_ELEM
         case(PRIS,MDLP); Nrelem = PRIS_VIS%NR_ELEM
         case(BRIC,MDLB); Nrelem = HEXA_VIS%NR_ELEM
         case default
            write(*,*) 'get_vis_nrelem: unexpected type = ', S_Type(Etype)
            stop
      end select
!
   end subroutine get_vis_nrelem



!> @brief return list of vertices for vis object
!> @param[in]  Nout
!> @param[in]  V
!> @param[in]  Etype - Element type
!> @date Feb 2023
   subroutine disp_vis(Nout,V,Etype)
      type(vis), intent(in):: V
      integer  , intent(in):: Nout,Etype
!
      integer :: i
!
      write(Nout,*) S_Type(Etype), ' NR_VERT = ', V%NR_VERT
      do i=0,(V%NR_VERT-1)
         write(Nout,*) V%VERTC(i,1:3)
      end do
!
      write(Nout,*) S_Type(Etype), ' NR_ELEM = ', V%NR_ELEM
      do i=1,V%NR_ELEM
         select case (Etype)
         case(TETR,MDLN)
            write(Nout,*) &
               V%ELEM(i,0), &
               V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3), V%ELEM(i,4)
         case(PRIS,MDLP)
            write(Nout,*) &
               V%ELEM(i,0), &
               V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3),&
               V%ELEM(i,4), V%ELEM(i,5), V%ELEM(i,6)
         case(BRIC,MDLB)
            write(Nout,*) &
               V%ELEM(i,0), &
               V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3), V%ELEM(i,4), &
               V%ELEM(i,5), V%ELEM(i,6), V%ELEM(i,7), V%ELEM(i,8)
         end select
      end do
   end subroutine disp_vis
!

!> @brief returns number of points/nodes for VTU element type number 
!> @param[in] indx - VTU Element type number
!> @date March 2023
   integer function nobj_conf_VTU(indx)
!
      integer :: indx
!
      if (SECOND_ORDER_VIS) then
         select case(indx)
            case(24); nobj_conf_VTU = 10
            case(32); nobj_conf_VTU = 18
            case(29); nobj_conf_VTU = 27
            case default
               write(*,*) 'nobj_conf_VTU: unexpected type. stop.'
               stop
         end select
      else
         select case(indx)
            case(10); nobj_conf_VTU = 4
            case(13); nobj_conf_VTU = 6
            case(12); nobj_conf_VTU = 8
            case default
               write(*,*) 'nobj_conf_VTU: unexpected type. stop.',indx
               stop
         end select
      endif
!
   end function nobj_conf_VTU
end module upscale
