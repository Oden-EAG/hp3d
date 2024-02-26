!> @brief Saves common variables used in frontal solver
!> @date Feb 2024
module surfsc1
!
   implicit none
!
   save
!
   integer :: ISYM, NUMELM, IRESOL, NRHS, MA, IWRT,   &
              IASSEM, IERR, NNEGP, NPOSP, IPRSTR,     &
              IPRPIV, IPRDES, NFSOUT, NICMUL,         &
              IPFSDB, IPFSS1, IPFSST, IPFSLH, IPFSRH, &
              IPFSLF, IPFSRF, IPFSZD, IPFSLE, IPFSRE, &
              IPFSBF, IPFSBK, IPFSXX
!
#if HP3D_COMPLEX
   integer :: IDUMPWR
#endif
!
   contains
!
   integer function modr(I,J)
!
      integer, intent(in) :: I,J
!
      modr = I - (I/J)*J
!
   end function modr
!
end module surfsc1

! DESCRIPTION OF VARIABLES
!
!      COMMON /SURFC1/ ISYM, NUMELM, IRESOL, NRHS, MA, IWRT,
!     .                IASSEM, IERR, NNEGP, NPOSP, IPRSTR,
!     .                IPRPIV, IPRDES, NFSOUT, NICMUL,
!     .                IPFSDB, IPFSS1, IPFSST, IPFSLH, IPFSRH,
!     .                IPFSLF, IPFSRF, IPFSZD, IPFSLE, IPFSRE,
!     .                IPFSBF, IPFSBK, IPFSXX
!
! Frontal Solver Common Block
! ---------------------------
!  This common is used to pass control info between the calling program
!    and SURFS and between routines internal to SURFS
!  The variables in this common are set by users interfacing w/ SURFS.
!
! Important Notes to Users Interfacing w/ SURFS
! --------------------------------------------
! *) If you are using resolution, then take care in sub.solin2 not to
!    zero out or touch the element LHS workspace which is sent thru the
!    calling arguments when sub.solin2 is called with iflag=1
!    where;  iflag=1: get element rhs only ; iflag=2: get lhs and rhs
!
! *) A note on using resolution.  You must set ISYM to 3 or 4 so that
!    the appropriate eliminated lhs data is written to diskfile.
!    Then call surfs with IRESOL=0 at least once to obtain an eliminated lhs.
!    Subsequent calls may be made with IRESOL=1, using this elimnated lhs.
!
! *) SURFS uses scratch files for equation buffering, (when it has to)
!    These scratch files are opened with the following unit numbers:
!
!       No resolution:                       19, 18, 20  (ISYM=1,2)
!       Unsymetric w/ resolution capability: 21, 22, 23  (ISYM=3)
!       Symmetric w/ resolution capability:  24, 25, 26  (ISYM=4)
!
!    Note that seperate unit numbers are used for the resolution capability
!         so that one may solve different sets of equations and still be able
!         to perform a resolution of a previous set.
!         (ie: ISYM=3,IRESOL=0,... ISYM=2,ISYM=2,..., ISYM=3,IRESOL=1 )
!
!    Note that right now you cannot do two different ISYM=3,IRESOL=0 (or 4)
!          in a row since only one set of buffer tapes is available
!
!
!
! Description of Variables
! ------------------------
! NFSOUT = The output diskfile unit number for printed output (ie:6, in general)
!
! NICMUL = The multiplier to be used in packing and unpacking NDOF
!             from the destination vectors.
!             This number reflects the number of digits which you think
!             are necessary to hold the number
!             (ie: NICMUL=10 ==> NDOF=0-9,  NICMUL=100 ==> NDOF=0-99,  etc)
!
! ISYM  = 1 ; Symmetric coefficient matrix (without resolution capability)
!       = 2 ; Unsymmetric coefficient matrix (without resolution capability)
!       = 3 ; Unsymmetric coefficient matrix (with resolution capability)
!       = 4 ; Symmetric coefficient matrix (with resolution capability)
!
!
! NUMELM = Number of Elements
!
! IRESOL = 0 ; Complete Forward Elimination
!               (ie: both RHS & LHS)
!        = 1 ; Forward Elimination of RHS only
!               (ie: resolution)
!
! NRHS = Number of RHS to be processed
!         (NRHS may be 0 for IRESOL=0)
!
! MA = Length of array IA in SUB.SURFSP and the
!       length of the work array A in SUB.SURFSS
!       MA need not be the same in these two routines.
!
!wb **note: setting IWRT does not write the scratch buffers to a permanent
!           diskfile, this is the responsibilty of the calling program
!
! IWRT = 0 ; Minimum Number of buffer dumps.
!      = 1 ; The last buffer is written to tape.
!             This is essential if a problem is to be terminated
!             and later restarted using the resolution capabilty.
!
! IASSEM = 0 ; Element LHS & RHS  will be assembled into the
!               global arrays by the calling program
!        = 1 ; Element LHS & RHS will be assembled by SURFS
!
! IERR = 0 ; Successful execution
!      = 1 ; Not enough room allocated for array IA in SURFSP
!             (MA must be increased)
!      = 2 ; Zero pivot
!      = 3 ; Attempted resolution with ISYM=2
!      = 4 ; Attempted resolution with NRHS=0
!      = 5 ; Resolution with insufficient room in the B buffer
!             (NRHS must be decreased)
!      = 6 ; Insufficient room in the U buffer
!             (MA must be increased)
!
! NNEGP = Number of negative pivots
! NPOSP = Number of positive pivots
! **** Note: NNEGP,NPOSP are set to -1 for resolution runs, since they are
!            determined only during LHS elimination, and are not relevant
!
! Print Flags  (=0 ; No prints)
! -----------
! IPRSTR = 1 ; Storage allocations and processing times
! IPRPIV = 1 ; Pivots
! IPRDES = 1 ; Destination vectors
!
! DEBUG Print Flags
! =================
! **** Note: turning these flags on could produce a HUGE amount of output
!
! IPFSDB = 1 ; Turns on all of the flags below
!
! IPFSS1 = 1 ; Print the entering contents of this common block
! IPFSZD = 1 ; Print buffer dump info
! IPFSBF = 1 ; Print first few items in the buffer
!
! **** Note: the flags below may be flagged as: = 1 ; print the info for all elems,
!                                               = -nel ; print for this element only
!                                               = 0 ; no prints
!
! IPFSST = 1 ; Print storage variables
! IPFSLH = 1 ; Print element lhs
! IPFSRH = 1 ; Print element rhs
! IPFSLF = 1 ; Print assembled front lhs
! IPFSRF = 1 ; Print assembled front rhs
! IPFSLE = 1 ; Print eliminated lhs
! IPFSRE = 1 ; Print eliminated rhs
! IPFSBK = 1 ; Print backward elimination info
! IPFSXX = 1 : Print the solution sent to solout
