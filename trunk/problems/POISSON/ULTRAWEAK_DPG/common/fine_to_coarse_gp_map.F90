
!-----------------------------------------------------------------------
! routine name - fine_to_coarse_gp_map
!---------------------------------------------------------------------

! last revision: 7 Dec 2022

! arguments:

!in :
!   iel: index of the son
!   xi: gauss point correspoding to the son
!   xis: gauss point corresponding to the coarse element.
!   etype: element type

!-----------------------------------------------------------------------------

subroutine fine_to_coarse_gp_map(iel,xi,xis,etype)

implicit none

real(8), intent(in) :: xi(3)
integer, intent(in) :: iel
character(len=4), intent(in) :: etype

real(8), intent(out) :: xis(3)


select case(etype)
    case('mdlb')
        select case(iel)
            case(1)
                xis  = 0.5 * xi
            case(2)
                xis(1) = 0.5 + 0.5 * xi(1)
                xis(2) = 0.5 * xi(2)
                xis(3) = 0.5 * xi(3)
            case(3)
                xis(1) = 0.5 + 0.5 * xi(1)
                xis(2) = 0.5 + 0.5 * xi(2)
                xis(3) = 0.5 * xi(3)
            case(4)
                xis(1) = 0.5 * xi(1)
                xis(2) = 0.5 + 0.5 * xi(2)
                xis(3) = 0.5 * xi(3)
            case(5)
                xis(1) = 0.5 * xi(1)
                xis(2) = 0.5 * xi(2)
                xis(3) = 0.5 + 0.5 * xi(3)
            case(6)
                xis(1) = 0.5 + 0.5 * xi(1)
                xis(2) = 0.5 * xi(2)
                xis(3) = 0.5 + 0.5 * xi(3)
            case(7)
                xis(1) = 0.5 + 0.5 * xi(1)
                xis(2) = 0.5 + 0.5 * xi(2)
                xis(3) = 0.5 + 0.5 * xi(3)
            case(8)
                xis(1) = 0.5 * xi(1)
                xis(2) = 0.5 + 0.5 * xi(2)
                xis(3) = 0.5 + 0.5 * xi(3)
            case default
                write(*,*) "Wrong son"
        end select
    case default
        write(*,*) "not implemented yet"


end select






end subroutine fine_to_coarse_gp_map


subroutine  fine_to_subson_gp_map(etype,kref,sel,xi,xis)

implicit none
integer,    intent(in)  :: kref
integer,    intent(in)  :: sel
character(len=4), intent(in) :: etype

real(8), intent(in) :: xi(3)
real(8), intent(out) :: xis(3)


if(etype .eq. 'mdlb') then

    select case(kref)

        case(100)
            select case(sel)

                case(1,2)
                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,1.d0,0.5d0,0.5d0,xi,xis)
                case(4,3)
                    call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.d0,1.d0,0.5d0,0.5d0,xi,xis)
                case(5,6)
                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.5d0,1.d0,0.5d0,0.5d0,xi,xis)
                case(8,7)
                    call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.5d0,1.d0,0.5d0,0.5d0,xi,xis)

                case default
                    write(*,*) "wrong fine son"

            end select
        case(010)
            select case(sel)
                case(1,4)
                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,0.5d0,1.d0,0.5d0,xi,xis)
                case(2,3)
                    call singlefineson_to_subson_gp_map(0.5d0,0.d0,0.d0,0.5d0,1.d0,0.5d0,xi,xis)
                case(5,8)
                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.5d0,0.5d0,1.d0,0.5d0,xi,xis)
                case(6,7)
                    call singlefineson_to_subson_gp_map(0.5d0,0.d0,0.5d0,0.5d0,1.d0,0.5d0,xi,xis)
                case default
                    write(*,*) "wrong fine son"

            end select

        case(001)
            select case(sel)

                case(1,5)
                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,0.5d0,0.5d0,1.d0,xi,xis)
                case(2,6)
                    call singlefineson_to_subson_gp_map(0.5d0,0.d0,0.d0,0.5d0,0.5d0,1.d0,xi,xis)
                case(3,7)
                    call singlefineson_to_subson_gp_map(0.5d0,0.5d0,0.d0,0.5d0,0.5d0,1.d0,xi,xis)
                case(4,8)
                    call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.d0,0.5d0,0.5d0,1.d0,xi,xis)
                case default
                    write(*,*) "wrong fine son"

            end select

        case(110)
            select case(sel)
                case(1,2,3,4)

                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,1.d0,1.d0,0.5d0,xi,xis)

                case(5,6,7,8)

                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.5d0,1.d0,1.d0,0.5d0,xi,xis) 

            end select


        case(101)
            select case(sel)
                case(1,2,5,6)

                call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,1.d0,0.5d0,1.d0,xi,xis)

                case(4,3,8,7)

                call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.d0,1.d0,0.5d0,1.d0,xi,xis) 

            end select

        case(011)
            select case(sel)
                case(1,4,5,8)

                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,0.5d0,1.d0,1.d0,xi,xis)

                case(2,3,6,7)

                    call singlefineson_to_subson_gp_map(0.5d0,0.0d0,0.d0,0.5d0,1.d0,1.d0,xi,xis) 

            end select


        case(111)

             xis = xi



    end select


endif



end subroutine fine_to_subson_gp_map





subroutine singlefineson_to_subson_gp_map(ofx,ofy,ofz,fx,fy,fz,xi,xis)

    implicit none
    real(8), intent(in) :: ofx
    real(8), intent(in) :: ofy
    real(8), intent(in) :: ofz

    real(8), intent(in) :: fx
    real(8), intent(in) :: fy
    real(8), intent(in) :: fz

    real(8),dimension(3),   intent(in)    :: xi
    real(8),dimension(3),   intent(out)   :: xis

    xis(1) = ofx + fx * xi(1)
    xis(2) = ofy + fy * xi(2)
    xis(3) = ofz + fz * xi(3)


end subroutine singlefineson_to_subson_gp_map


subroutine subson_one_irregularity_map_isoref(etype,kref_intent,subson_impl,father_subson)


    implicit none

    character(len=4), intent(in) :: etype
    integer,    intent(in)  :: kref_intent !intented h-ref
    integer,    intent(in)  :: subson_impl  !isoref subson index when iso ref is implemented in place of intended h-refinement
    integer,    intent(out) :: father_subson ! subson which contains the isoref subson

    integer :: i,j

    if(etype .eq. "mdlb") then

        select  case(kref_intent)

            case(100)
                if((subson_impl .eq. 1) .or. (subson_impl .eq. 4) .or. (subson_impl .eq. 5) .or. (subson_impl .eq. 8)) then
                    father_subson = 1
                else
                    father_subson = 2
                endif

            case(010)
                if((subson_impl .eq. 1) .or. (subson_impl .eq. 2) .or. (subson_impl .eq. 5) .or. (subson_impl .eq. 6)) then
                    father_subson = 1
                else
                    father_subson = 2
                endif 

            case(001)
                if((subson_impl .eq. 1) .or. (subson_impl .eq. 2) .or. (subson_impl .eq. 3) .or. (subson_impl .eq. 4)) then
                    father_subson = 1
                else
                    father_subson = 2
                endif 

            case(110)
                if((subson_impl .eq. 1) .or. (subson_impl .eq. 5)) then
                    father_subson = 1
                else if((subson_impl .eq. 2) .or. (subson_impl .eq. 6)) then
                    father_subson = 2
                else if((subson_impl .eq. 3) .or. (subson_impl .eq. 7)) then
                    father_subson = 3
                else if((subson_impl .eq. 4) .or. (subson_impl .eq. 8)) then
                    father_subson = 4
                endif

            case(101)
                if((subson_impl .eq. 1) .or. (subson_impl .eq. 4)) then
                    father_subson = 1
                else if((subson_impl .eq. 2) .or. (subson_impl .eq. 3)) then
                    father_subson = 2
                else if((subson_impl .eq. 5) .or. (subson_impl .eq. 8)) then
                    father_subson = 3
                else if((subson_impl .eq. 6) .or. (subson_impl .eq. 7)) then
                    father_subson = 4
                endif
            
            case(011)
                if((subson_impl .eq. 1) .or. (subson_impl .eq. 2)) then
                    father_subson = 1
                else if((subson_impl .eq. 3) .or. (subson_impl .eq. 4)) then
                    father_subson = 2
                else if((subson_impl .eq. 5) .or. (subson_impl .eq. 6)) then
                    father_subson = 3
                else if((subson_impl .eq. 7) .or. (subson_impl .eq. 8)) then
                    father_subson = 4
                endif
            
            case(111)

                do i = 1,8

                    if(subson_impl .eq. i) then

                        father_subson = i
                    endif
                enddo

            case default
                write(*,*) " not a valid kref"
        end select



    endif




end subroutine subson_one_irregularity_map_isoref