!---------------------------------------------------------------------------------
subroutine generate_PRI_TET_PYR
!---------------------------------------------------------------------------------
! LATEST REVISION: Jul 09
!
! PURPOSE: routine generates new prisms, tets and pyramids
!---------------------------------------------------------------------------------
  use SPLIT_SURF
  IMPLICIT NONE
!---------------------------------------------------------------------------------
  integer :: ifig,nt,nr,ile,iv,iv1,iv2,iv3,ie,np,np1,np2,np3,flag
! FUNCTIONS
  integer :: mod3,mod4
!---------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'generate_PRI_TET_PYR: generating prisms, tets, pyramids...'
#endif
    flag = 0

! ..loop over figures to split
    do ifig = 1, nr_figs_to_split
      select case(figs_split(1,ifig))
! ......triangle
        case(1)
          nt = figs_split(2,ifig)
! ........determine how many points have been duplicated
          ile = 0
          do iv = 1,3
            np = TRIANGLES(nt)%VertNo(iv)
            if (new_point(np) .ne. np)  ile = ile + 1
          enddo
! ........select number of duplicated points
          select case(ile)
! ..........generate TETRA
            case(1)
              NRTETRA = NRTETRA + 1
              if (NRTETRA .gt. MAXTE) then
                write(*,*) 'generate_PRI_TET_PYR: INCREASE MAXTE'
                stop
              endif
              TETRAS(NRTETRA)%Type        = 'Linear'
              TETRAS(NRTETRA)%Domain      = NRDOMAIN
! ............loop over triangle vertices and identify duplicated point
              do iv = 1, 3
                np = TRIANGLES(nt)%VertNo(iv)
!  .............save local number of duplicated vertex
                if (new_point(np) .ne. np) then
                  iv1 = iv; exit
                endif
              enddo
!  ...........set up new tet
              do iv = 1,3
                iv2 = mod3(iv1 + iv - 1)
                TETRAS(NRTETRA)%VertNo(iv) = TRIANGLES(nt)%VertNo(iv2)
              enddo
              TETRAS(NRTETRA)%VertNo(4) = new_point(np)
!  ...........swap x and y-axis if necessary
              call check_orientation(3,NRTETRA)
#if I_PRINT >= 2
              write(*,*) 'generate_PRI_TET_PYR: GENERATED NEW TET = ',NRTETRA,TETRAS(NRTETRA)%VertNo(1:4)
#endif
! ..........generate PYRAMID
            case(2)
              NRPYRAM = NRPYRAM + 1
              if (NRPYRAM .gt. MAXPY) then
                write(*,*) 'generate_PRI_TET_PYR: INCREASE MAXPY'
                stop
              endif
              PYRAMIDS(NRPYRAM)%Type   = 'Linear'
              PYRAMIDS(NRPYRAM)%Domain = NRDOMAIN
              do ie = 1,3
                iv = ie;  iv1 = mod3(iv + 1)
                np  = TRIANGLES(nt)%VertNo(iv)
                np1 = TRIANGLES(nt)%VertNo(iv1)
                if ((new_point(np).ne.np).and.(new_point(np1).ne.np1))exit
              enddo
              PYRAMIDS(NRPYRAM)%VertNo(1) = np
              PYRAMIDS(NRPYRAM)%VertNo(2) = np1
              PYRAMIDS(NRPYRAM)%VertNo(3) = new_point(np1)
              PYRAMIDS(NRPYRAM)%VertNo(4) = new_point(np)
              iv2 = mod3(iv + 2)
              PYRAMIDS(NRPYRAM)%VertNo(5) = TRIANGLES(nt)%VertNo(iv2)
              call check_orientation(4,NRPYRAM)
! ..........generate a new prism
            case(3)
              NRPRISM = NRPRISM + 1
              if (NRPRISM .gt. MAXBT) then
                write(*,*) 'generate_PRI_TET_PYR: INCREASE MAXBT'
                stop
              endif
              PRISMS(NRPRISM)%Type   = 'Linear'
              PRISMS(NRPRISM)%Domain = NRDOMAIN
              do iv = 1, 3
                np = TRIANGLES(nt)%VertNo(iv)
                PRISMS(NRPRISM)%VertNo(iv)     = np
                PRISMS(NRPRISM)%VertNo(3 + iv) = new_point(np)
              enddo
              call check_orientation(1,NRPRISM)
            case default
              write(*,*) 'generate_PRI_TET_PYR: INCONSISTENCY 1, ile = ',ile
              stop
          end select
!
! ......split rectangle
        case(2)
          flag = 1
          nr = figs_split(2,ifig)
! ........determine how many points have been duplicated
          ile = 0
          do iv = 1, 4
            np = RECTANGLES(nr)%VertNo(iv)
            if (new_point(np) .ne. np)  ile = ile + 1
          enddo
!
          select case(ile)
!
! ........generate a new prism
          case(2)
            NRPRISM = NRPRISM + 1
            if (NRPRISM .gt. MAXBT) then
              write(*,*) 'generate_PRI_TET_PYR: INCREASE MAXBT'
              stop
            endif
            PRISMS(NRPRISM)%Type   = 'Linear'
            PRISMS(NRPRISM)%Domain = NRDOMAIN + 1
            do ie = 1, 4
              iv = ie;  iv1 = mod4(iv + 1)
              np  = RECTANGLES(nr)%VertNo(iv)
              np1 = RECTANGLES(nr)%VertNo(iv1)
              if ((new_point(np) .ne. np) .and. (new_point(np1) .eq. np1))  exit
            enddo
            PRISMS(NRPRISM)%VertNo(1) = np
            PRISMS(NRPRISM)%VertNo(2) = new_point(np)
            PRISMS(NRPRISM)%VertNo(3) = np1
            iv2 = mod4(iv + 2);  iv3 = mod4(iv + 3)
            np2 = RECTANGLES(nr)%VertNo(iv2)
            np3 = RECTANGLES(nr)%VertNo(iv3)
            PRISMS(NRPRISM)%VertNo(4) = np3
            PRISMS(NRPRISM)%VertNo(5) = new_point(np3)
            PRISMS(NRPRISM)%VertNo(6) = np2
            call check_orientation(1,NRPRISM)
#if I_PRINT >= 2
            write(*,7020) NRPRISM
 7020       format(' generate_PRI_TET_PYR: GENERATED NEW PRISM ',i6)
            write(*,7021) PRISMS(NRPRISM)%VertNo(1:6)
 7021       format(' WITH VERTEX POINTS =',6i6)
#endif
!
          case default
            write(*,*) 'generate_PRI_TET_PYR: INCONSISTENCY 2, ile = ',ile
            stop
        end select
      end select
!
! ..end of loop through figures to split
    enddo
    if (flag.eq.1) NRDOMAIN = NRDOMAIN + 1

#if I_PRINT >= 1
    write(*,*)'generate_PRI_TET_PYR: done!'
#endif
!
end subroutine generate_PRI_TET_PYR
