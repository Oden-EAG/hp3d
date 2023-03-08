!> @brief Saves common variables used in frontal solver
!> @date Mar 2023
module surfsc2
!
   save
!
   integer :: IB, IU, IL, IFB, IFU, IFL, MBUF, MW, MKF,  &
              MELEM, MFWR, MB, MDOF, MLDEST, NRHSF, NFN, &
              NFW, KFW, LFW, NDOFM, MFW, IDUMP, LENU,    &
              MBUFSV(2,2)
!
end module surfsc2

! DESCRIPTION OF VARIABLES
!
!      COMMON /SURFC2/ IB, IU, IL, IFB, IFU, IFL, MBUF, MW, MKF,
!     .                MELEM, MFWR, MB, MDOF, MLDEST, NRHSF, NFN,
!     .                NFW, KFW, LFW, NDOFM, MFW, IDUMP, LENU,
!     .                MBUFSV(2,2)
!
! Frontal Solver Common Block
! ---------------------------
!  This common is used to pass info between routines internal to SURFS
!
! Description of Variables
! ------------------------
! Workspace manipulation variables
! ---------------------------------
! BUFFER VARIABLES
! ----------------
! SURFS uses 3 buffers cut from a common workspace, into which the eliminated
!   equations are written.  The length of these buffers is based on available
!   workspace.  When the buffer fill they are dumped to a scratch diskfile to
!   make room for new equations.  The buffers are also dumped, full or not,
!   for resolution, cuz one cant count on the workspace coming back unaltered
!   when resolution is to begin
!
!   The buffers are:
!   U  :  The eliminated LHS buffer
!   L  :  The additional eliminated LHS buffer for unsymmetric w/ resolution
!   B  :  The eliminated RHS buffer
!
! IB    :  Current pointer into the RHS buffer (B)
!          Note: IB points at the Global RHS within the current buffer
!          (Also functions as the Length of the current RHS buffer)
!
! IFB   :  Record counter for the RHS buffer (B) diskfile
!
! IU    : Current pointer into the LHS buffer (U)
!          Note: IU points at a Global LHS equation within the current buffer
!          (Also functions as the Length of the current LHS buffer)
!        *NOTE: IU runs forward into the workspace: BUFF(1-->MBUF)
!
! IFU   :  Record counter for the LHS buffer (U) diskfile
!
! IL    : Current pointer into the LHS buffer (L) (Extra unsymmetric LHS buffer)
!          Note: IL points at a Global LHS equation within the current buffer
!          (Also functions as the Length of the current LHS buffer (L)
!        *NOTE: IL runs backward into the workspace BUFF(1<--MBUF)
!
! IFL   :  Record counter for the LHS buffer (L) diskfile
!
! IDUMP  : Flag indicating whether buffers have been dumped to diskfiles
!
! LENU   : Length of the U buffer when its all kept in core
!-----------------------------------------------------------------------------------
! SPACE CUT LENGTH VARIABLES
! --------------------------
! MELEM : Temp space required to hold the maximum elem LHS + RHS, during assembly
! MKF   : Temp space required to hold the LHS equations for the max frontwidth
! MBUF  : All space which is left over is allocated for buffering equations
! MFWR  : Temp space required to hold the LHS & RHS equations for the max frontwidth
! MW    : Total Temp space for front and individual elem temp storage
! MFWR  : Temp space required to hold the LHS & RHS for the max frontwidth
! MB    : Temp space which is available for buffering rhs equations
!
! NFN   : Total number of destination vectors in the model (from prefront)
! MFW   : Maximum DOF in the front at any given time (from prefront)
! MDOF  : Max nunber of DOF associatted w/ any elem (from prefront)
! MLDEST: Max nunber of DESTVECS associatted w/ any elem (from prefront)
!
! MBUFSV: Save the value of MBUF & MW to check against during resolution
! *** Note: we cannot let the size of MBUF shrink during resolution
!            because then the stuff that was buffered out will not fit back in
!         MBUFSV(1,i) for ISYM=3, MBUFSV(2,i) for ISYM=4
!
! FRONT LENGTH VARIABLES
! ----------------------
! NFW   : Number of DOF in the current front (while working on this elem)
!          This comes out of SUB.DEST as the true length
!          BUT the front is NOT allowed reduce, only to grow (for indexing reasons)
!          NFW is a dynamic counter for the current frontwidth
!
! KFW   : The length of the current front (used for indexing in unsymmetric storage)
!          *** the front is NOT allowed reduce, only to grow (for indexing reasons)
!          KFW only grows
!
! LFW   : FINAL number of DOF in the last front (while working on the previous elem)
!
! MISC. VARIABLES
! ---------------
! NDOFM : Total number of DOF for this elem  (from SUB.DEST)
! NRHSF : This is simply a copy of NRHS (from SURFSC1.BLK)
