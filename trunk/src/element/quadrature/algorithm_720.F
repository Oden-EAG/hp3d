C      ALGORITHM 720, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 19, NO. 3, SEPTEMBER, 1993, PP. 318-330.
C
C     This file summarizes the real(8) routines
C     in this directory. SINGLE PRECISION names are obtained
C     by replacing the leading D with S in all names.
C
C     DTEST1   Simple test driver for DCUTET.
C     DTEST2   More detailed test driver for DCUTET.
C              Sample output from a SUN 3/50 is included.
C     DCUTET   Main Integrator.DCUTET calls DCHTET and DADTET.
C     DCHTET   Checks the input to DCUTET.
C     DADTET   The adaptive integration routine. 
C              DADTET calls DTRTET and DRLTET.
C     DTRTET   Maintaines the heap of subtetrahedronds.
C     DRLTET   Computes estimates of integral and error over
C              subtetrahedrons.
C     DORTET   Computes fully symmetric sums of function evaluations.
      real(8) FUNCTION D1MACH(I)
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
C  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C
C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED integer(2)), YOU MAY NEED TO
C  CHANGE INTEGER TO integer(4) OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
      INTEGER SC
C
      real(8) DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES AND MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.
C     
       DATA SMALL(1),SMALL(2) /    1048576,          0 /
       DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
       DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
       DATA DIVER(1),DIVER(2) / 1018167296,          0 /
       DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /, SC/987/
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES AND 8087-BASED
C     MICROS, SUCH AS THE IBM PC AND AT&T 6300, IN WHICH THE LEAST
C     SIGNIFICANT BYTE IS STORED FIRST.
C
C      DATA SMALL(1),SMALL(2) /          0,    1048576 /
C      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /
C      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /
C      DATA DIVER(1),DIVER(2) /          0, 1018167296 /
C      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /, SC/987/
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1),SMALL(2) /    1048576,          0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /
C      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /
C      DATA DIVER(1),DIVER(2) /  873463808,          0 /
C      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA SMALL(1) / ZC00800000 /
C      DATA SMALL(2) / Z000000000 /
C
C      DATA LARGE(1) / ZDFFFFFFFF /
C      DATA LARGE(2) / ZFFFFFFFFF /
C
C      DATA RIGHT(1) / ZCC5800000 /
C      DATA RIGHT(2) / Z000000000 /
C
C      DATA DIVER(1) / ZCC6800000 /
C      DATA DIVER(2) / Z000000000 /
C
C      DATA LOG10(1) / ZD00E730E7 /
C      DATA LOG10(2) / ZC77800DC0 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O0000000000000000 /
C
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O0007777777777777 /
C
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O7770000000000000 /
C
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O7777777777777777 /
C
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN4 ON THE CDC 6000/7000 SERIES.
C
C      DATA SMALL(1) / 00564000000000000000B /
C      DATA SMALL(2) / 00000000000000000000B /
C
C      DATA LARGE(1) / 37757777777777777777B /
C      DATA LARGE(2) / 37157777777777777774B /
C
C      DATA RIGHT(1) / 15624000000000000000B /
C      DATA RIGHT(2) / 00000000000000000000B /
C
C      DATA DIVER(1) / 15634000000000000000B /
C      DATA DIVER(2) / 00000000000000000000B /
C
C      DATA LOG10(1) / 17164642023241175717B /
C      DATA LOG10(2) / 16367571421742254654B /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN5 ON THE CDC 6000/7000 SERIES.
C
C      DATA SMALL(1) / O"00564000000000000000" /
C      DATA SMALL(2) / O"00000000000000000000" /
C
C      DATA LARGE(1) / O"37757777777777777777" /
C      DATA LARGE(2) / O"37157777777777777774" /
C
C      DATA RIGHT(1) / O"15624000000000000000" /
C      DATA RIGHT(2) / O"00000000000000000000" /
C
C      DATA DIVER(1) / O"15634000000000000000" /
C      DATA DIVER(2) / O"00000000000000000000" /
C
C      DATA LOG10(1) / O"17164642023241175717" /
C      DATA LOG10(2) / O"16367571421742254654" /, SC/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1
C
C      DATA SMALL(1),SMALL(2) / '00100000'X, '00000000'X /
C      DATA LARGE(1),LARGE(2) / '7FFFFFFF'X, 'FFFFFFFF'X /
C      DATA RIGHT(1),RIGHT(2) / '3CC00000'X, '00000000'X /
C      DATA DIVER(1),DIVER(2) / '3CD00000'X, '00000000'X /
C      DATA LOG10(1),LOG10(2) / '3FF34413'X, '509F79FF'X /, SC/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA SMALL(1) / 201354000000000000000B /
C      DATA SMALL(2) / 000000000000000000000B /
C
C      DATA LARGE(1) / 577767777777777777777B /
C      DATA LARGE(2) / 000007777777777777776B /
C
C      DATA RIGHT(1) / 376434000000000000000B /
C      DATA RIGHT(2) / 000000000000000000000B /
C
C      DATA DIVER(1) / 376444000000000000000B /
C      DATA DIVER(2) / 000000000000000000000B /
C
C      DATA LOG10(1) / 377774642023241175717B /
C      DATA LOG10(2) / 000007571421742254654B /, SC/987/
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC DMACH(5)
C
C      DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/
C      DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/
C      DATA LOG10/40423K,42023K,50237K,74776K/, SC/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '37777577 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000333 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000334 /
C      DATA LOG10(1),LOG10(2) / '23210115, '10237777 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
C      DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C      DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
C      DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
C      DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /, SC/987/
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
C      DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
C      DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /
C      DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /
C      DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C      DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
C      DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C      DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
C      DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "047674776746 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
C
C      DATA SMALL(1),SMALL(2) / O00040000000, O00000000000 /
C      DATA LARGE(1),LARGE(2) / O17777777777, O37777777777 /
C      DATA RIGHT(1),RIGHT(2) / O04440000000, O00000000000 /
C      DATA DIVER(1),DIVER(2) / O04500000000, O00000000000 /
C      DATA LOG10(1),LOG10(2) / O07746420232, O20476747770 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /    128,      0 /
C      DATA SMALL(3),SMALL(4) /      0,      0 /
C
C      DATA LARGE(1),LARGE(2) /  32767,     -1 /
C      DATA LARGE(3),LARGE(4) /     -1,     -1 /
C
C      DATA RIGHT(1),RIGHT(2) /   9344,      0 /
C      DATA RIGHT(3),RIGHT(4) /      0,      0 /
C
C      DATA DIVER(1),DIVER(2) /   9472,      0 /
C      DATA DIVER(3),DIVER(4) /      0,      0 /
C
C      DATA LOG10(1),LOG10(2) /  16282,   8346 /
C      DATA LOG10(3),LOG10(4) / -31493, -12296 /, SC/987/
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA SMALL(3),SMALL(4) / O000000, O000000 /
C
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA LARGE(3),LARGE(4) / O177777, O177777 /
C
C      DATA RIGHT(1),RIGHT(2) / O022200, O000000 /
C      DATA RIGHT(3),RIGHT(4) / O000000, O000000 /
C
C      DATA DIVER(1),DIVER(2) / O022400, O000000 /
C      DATA DIVER(3),DIVER(4) / O000000, O000000 /
C
C      DATA LOG10(1),LOG10(2) / O037632, O020232 /
C      DATA LOG10(3),LOG10(4) / O102373, O147770 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
C     WITH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
C     SUPPLIED BY IGOR BRAY.
C
C      DATA SMALL(1),SMALL(2) / :10000000000, :00000100001 /
C      DATA LARGE(1),LARGE(2) / :17777777777, :37777677775 /
C      DATA RIGHT(1),RIGHT(2) / :10000000000, :00000000122 /
C      DATA DIVER(1),DIVER(2) / :10000000000, :00000000123 /
C      DATA LOG10(1),LOG10(2) / :11504046501, :07674600177 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000
C
C      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /
C      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /
C      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /
C      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /
C      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER
C
C      DATA SMALL(1),SMALL(2) /        128,           0 /
C      DATA LARGE(1),LARGE(2) /     -32769,          -1 /
C      DATA RIGHT(1),RIGHT(2) /       9344,           0 /
C      DATA DIVER(1),DIVER(2) /       9472,           0 /
C      DATA LOG10(1),LOG10(2) /  546979738,  -805796613 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER
C
C      DATA SMALL(1),SMALL(2) / Z00000080, Z00000000 /
C      DATA LARGE(1),LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C      DATA RIGHT(1),RIGHT(2) / Z00002480, Z00000000 /
C      DATA DIVER(1),DIVER(2) / Z00002500, Z00000000 /
C      DATA LOG10(1),LOG10(2) / Z209A3F9A, ZCFF884FB /, SC/987/
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2
C
C      DATA SMALL(1),SMALL(2) /       '80'X,        '0'X /
C      DATA LARGE(1),LARGE(2) / 'FFFF7FFF'X, 'FFFFFFFF'X /
C      DATA RIGHT(1),RIGHT(2) /     '2480'X,        '0'X /
C      DATA DIVER(1),DIVER(2) /     '2500'X,        '0'X /
C      DATA LOG10(1),LOG10(2) / '209A3F9A'X, 'CFF884FB'X /, SC/987/
C
C  ***  ISSUE STOP 779 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SC .NE. 987) STOP 779
      IF (I .LT. 1  .OR.  I .GT. 5) GOTO 999
      D1MACH = DMACH(I)
      RETURN
  999 WRITE(*,1999) I
 1999 FORMAT(' D1MACH - I OUT OF BOUNDS',I10)
      STOP
      END

      SUBROUTINE DADTET(NUMFUN,MDIV,VER,NUMTET,MINSUB,MAXSUB,FUNSUB,
     +                  EPSABS,EPSREL,LENVER,RESTAR,LENW,RESULT,ABSERR,
     +                  NEVAL,NSUB,IFAIL,VALUES,ERRORS,GREATE,WORK2,
     +                  WORK3,LIST,VACANT)
C***BEGIN PROLOGUE DADTET
C***REFER TO DCUTET
C***DESCRIPTION Computation of integrals over regions consisting of
C            a collection of tetrahedrons.
C
C            DADTET repeatedly
C            subdivides the tetrahedrons with greatest estimated  errors
C            and estimates the integrals and the errors over the new
C            new sub-tetrahedrons until the error request is met
C            or MAXPTS function evaluations have been used.
C
C            A 43 point integration rule
C            with all evaluation points inside the tetrahedron
C            is applied. The rule has polynomial degree 8.
C
C   ON ENTRY
C
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            MDIV is the number of tetrahedrons that are divided in
C            each subdivision step in DADTET.
C            MDIV is chosen default to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(8*MDIV,NPROC) = 0.
C     VER    Real array of dimension (3,4,LENVER).
C            VER(1,K,L), VER(2,K,L) and VER(3,K,L) are the x, y and z
C            coordinates respectively of vertex K in tetrahedron L.
C            On entry VER(*,*,L) must contain the vertices of the
C            NUMTET user specified tetrahedrons for L = 1,2,...,NUMTET.
C     NUMTET Integer.
C            The number of tetrahedrons.
C     MINSUB Integer.
C            The minimum allowed number of subregions.
C     MAXSUB Integer.
C            The maximum allowed number of subregions.
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (X,NUMFUN,FUNVLS)
C            Input parameters:
C              X(1)   The x-coordinate of the evaluation point.
C              X(2)   The y-coordinate of the evaluation point.
C              X(3)   The z-coordinate of the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C
C     LENVER Integer.
C            Defines the length of the array VER.
C
C            We let
C            MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS.
C            MAXSUB = 7*(MAXPTS-43*NUMTET)/(8*43) + NUMTET
C            LENVER should then be greater or equal to MAXSUB.
C
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1,
C            then we restart a previous attempt.
C            In this case the only parameters for DCUTET that may
C            be changed (with respect to the previous call of DCUTET)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C     LENW   Integer.
C            Length of the workspace WORK2.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute errors.
C     NEVAL  Integer.
C            Number of function evaluations used by DCUTET.
C     NSUB   Integer.
C            The number of tetrahedrons in the data structure.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit.
C
C              ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
C              function evaluations for all values of K,
C              1 <= K <= NUMFUN .
C
C            IFAIL = 1 if MAXSUB was too small for DADTET
C              to obtain the required accuracy. In this case DADTET
C              returns values of RESULT with estimated absolute
C              errors ABSERR.
C     VALUES Real array of dimension (NUMFUN,MAXSUB).
C            The estimated components of the integrals over the
C            subregions.
C     ERRORS Real array of dimension (NUMFUN,MAXSUB).
C            The estimated errors over the subregions.
C     GREATE Real array of dimension MAXSUB.
C            The greatest errors in each subregion.
C     WORK2  Real array of dimension 6*LENW.
C            Work array used in DRLTET.
C     WORK3  Real array of dimension LENW.
C            Work array used in DRLTET.
C     LIST   Integer array used in DTRTET of dimension LENVER.
C            Is a partially sorted list, where LIST(1) is the top
C            element in a heap of subregions.
C     VACANT Integer array of dimension MDIV
C            Used as intermediate storage of pointers.
C***ROUTINES CALLED DTRTET,DRLTET
C***END PROLOGUE DADTET
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NUMFUN,MDIV,NUMTET,MINSUB,MAXSUB,LENW,RESTAR,LENVER,NEVAL,
     +        NSUB,IFAIL,LIST(LENVER),VACANT(MDIV)
      real(8) VER(3,4,LENVER),EPSABS,EPSREL,RESULT(NUMFUN),
     +                 ABSERR(NUMFUN),VALUES(NUMFUN,MAXSUB),
     +                 ERRORS(NUMFUN,MAXSUB),GREATE(MAXSUB),WORK3(LENW),
     +                 WORK2(6*LENW)
C
C   Local variables.
C
C   SBRGNS is the number of stored subregions.
C   NDIV   The number of subregions to be divided in each main step.
C   POINTR Pointer to the position in the data structure where
C          the new subregions are to be stored.
C   G      is the homogeneous coordinates of the integration rule.
C   W      is the weights of the integration rule and the null rules.
C   TOP    is a pointer to the top element in the heap of subregions.
C   COUNT  is the number of times the termination criterium is
C          satisfied.
      INTEGER I,J,K,L,TOP,SIZE
      INTEGER SBRGNS,I1,I2,I3,I4,I5,I6,I7,I8
      INTEGER L1,K1
      INTEGER COUNT,GROUP,LEFT,NDIV,POINTR,INDEX,NUM
      real(8) VEROLD(3,4),VERNEW(3,6),SUMVAL(2),SUMERR(2)
C
C***FIRST EXECUTABLE STATEMENT DADTET
C
      COUNT = 0
      NUM = 43
      IF (RESTAR.EQ.1) THEN
          SBRGNS = NSUB
          GO TO 110
      END IF
C
C   Initiate RESULT and ABSERR.
C
      DO 10 J = 1,NUMFUN
          RESULT(J) = 0
          ABSERR(J) = 0
   10 CONTINUE
C
C   Apply DRLTET over the NUMTET tetrahedrons.
C   This loop may be run in parallel.
C
      DO 20 I = 1,NUMTET
          L1 = 1 + (I-1)*6*NUMFUN
          K1 = 1 + (I-1)*NUMFUN
          CALL DRLTET(VER(1,1,I),NUMFUN,FUNSUB,WORK2(L1),VALUES(1,I),
     +                ERRORS(1,I),GREATE(I),WORK3(K1))
   20 CONTINUE
      NEVAL = NUM*NUMTET
      SBRGNS = NUMTET
C
C   Add the computed values to RESULT and ABSERR.
C
      DO 40 I = 1,NUMTET
          DO 30 J = 1,NUMFUN
              RESULT(J) = RESULT(J) + VALUES(J,I)
              ABSERR(J) = ABSERR(J) + ERRORS(J,I)
   30     CONTINUE
   40 CONTINUE
C
C   Store results in heap.
C
      DO 50 I = 1,NUMTET
          INDEX = I
          CALL DTRTET(2,INDEX,GREATE,LIST,(INDEX))
   50 CONTINUE

C
C   We check for termination.
C
      IF (SBRGNS.LT.MINSUB) THEN
          GO TO 110
      END IF
      DO 60 J = 1,NUMFUN
          IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)) .AND.
     +        ABSERR(J).GT.EPSABS) THEN
              GO TO 110
          END IF
   60 CONTINUE
      IFAIL = 0
      COUNT = COUNT + 1
      GO TO 499
C
C***End initiation.
C
C***Begin loop while the error is too great,
C   and SBRGNS+7 is less than MAXSUB.
C
  110 IF (SBRGNS+7.LE.MAXSUB) THEN
C
C   If we are allowed to divide further,
C   prepare to apply basic rule over  the tetrahedrons produced
C   by dividing the
C   NDIV sub-tetrahedrons with greatest errors.
C   If MAXSUB and SBRGNS are great enough, NDIV = MDIV.
C
          NDIV = MAXSUB - SBRGNS
          NDIV = MIN(NDIV,MDIV,SBRGNS)
C
C   Divide each of the NDIV sub-tetrahedrons in eight new
C   sub-tetrahedrons, and compute integral and error over each.
C   When we pick a tetrahedron to divide in eight one of the
C   new sub-tetrahedron is stored in the original tetrahedrons'
C   position in the datastructure. Thus we get 7*NDIV new
C   elements in the datastructure after the subdivision. The size
C   of the datastructure before the sub-division
C   is stored in the variable SIZE, while SBRGNS is the size of
C   heap at any time.
          SIZE = SBRGNS
          DO 150 I = 1,NDIV
              POINTR = SIZE + 7* (NDIV+1-I)
C
C   Adjust RESULT and ABSERR. TOP is a pointer to the top of the heap.
C
              TOP = LIST(1)
              VACANT(I) = TOP
              DO 115 J = 1,NUMFUN
                  RESULT(J) = RESULT(J) - VALUES(J,TOP)
                  ABSERR(J) = ABSERR(J) - ERRORS(J,TOP)
  115         CONTINUE
C
C   Save the vertices.
C
              DO 130 L = 1,3
                  DO 135 K = 1,4
                      VEROLD(L,K) = VER(L,K,TOP)
  135             CONTINUE
  130         CONTINUE
C
C   Adjust the heap.
C
              CALL DTRTET(1,SBRGNS,GREATE,LIST,K)
C
C   Compute the eight new tetrahedrons. Store 7 in consecutive
C   positions and 1 in the vacant position.
C
              I1 = TOP
              I2 = POINTR - 6
              I3 = POINTR - 5
              I4 = POINTR - 4
              I5 = POINTR - 3
              I6 = POINTR - 2
              I7 = POINTR - 1
              I8 = POINTR
C
C   Compute the 6 division points.
C
              DO 140 L = 1,3
                  VERNEW(L,1) = (VEROLD(L,1)+VEROLD(L,2))/2
                  VERNEW(L,2) = (VEROLD(L,2)+VEROLD(L,3))/2
                  VERNEW(L,3) = (VEROLD(L,1)+VEROLD(L,3))/2
                  VERNEW(L,4) = (VEROLD(L,1)+VEROLD(L,4))/2
                  VERNEW(L,5) = (VEROLD(L,2)+VEROLD(L,4))/2
                  VERNEW(L,6) = (VEROLD(L,3)+VEROLD(L,4))/2
  140         CONTINUE
C
              DO 145 L = 1,3
                  VER(L,1,I1) = VEROLD(L,1)
                  VER(L,2,I1) = VERNEW(L,1)
                  VER(L,3,I1) = VERNEW(L,3)
                  VER(L,4,I1) = VERNEW(L,4)
                  VER(L,1,I2) = VEROLD(L,2)
                  VER(L,2,I2) = VERNEW(L,1)
                  VER(L,3,I2) = VERNEW(L,2)
                  VER(L,4,I2) = VERNEW(L,5)
                  VER(L,1,I3) = VEROLD(L,3)
                  VER(L,2,I3) = VERNEW(L,2)
                  VER(L,3,I3) = VERNEW(L,3)
                  VER(L,4,I3) = VERNEW(L,6)
                  VER(L,1,I4) = VEROLD(L,4)
                  VER(L,2,I4) = VERNEW(L,4)
                  VER(L,3,I4) = VERNEW(L,5)
                  VER(L,4,I4) = VERNEW(L,6)
                  VER(L,1,I5) = VERNEW(L,1)
                  VER(L,2,I5) = VERNEW(L,3)
                  VER(L,3,I5) = VERNEW(L,4)
                  VER(L,4,I5) = VERNEW(L,5)
                  VER(L,1,I6) = VERNEW(L,1)
                  VER(L,2,I6) = VERNEW(L,2)
                  VER(L,3,I6) = VERNEW(L,3)
                  VER(L,4,I6) = VERNEW(L,5)
                  VER(L,1,I7) = VERNEW(L,2)
                  VER(L,2,I7) = VERNEW(L,3)
                  VER(L,3,I7) = VERNEW(L,5)
                  VER(L,4,I7) = VERNEW(L,6)
                  VER(L,1,I8) = VERNEW(L,3)
                  VER(L,2,I8) = VERNEW(L,4)
                  VER(L,3,I8) = VERNEW(L,5)
                  VER(L,4,I8) = VERNEW(L,6)
  145         CONTINUE
  150     CONTINUE
C
C   Apply basic rule.
C   This loop may be run in parallel.
C
          DO 200 I = 1,8*NDIV
              IF (I.LE.NDIV) THEN
                  INDEX = VACANT(I)
              ELSE
                  INDEX = SBRGNS + I
              END IF
              L1 = 1 + (I-1)*6*NUMFUN
              K1 = 1 + (I-1)*NUMFUN
              CALL DRLTET(VER(1,1,INDEX),NUMFUN,FUNSUB,WORK2(L1),
     +                    VALUES(1,INDEX),ERRORS(1,INDEX),GREATE(INDEX),
     +                    WORK3(K1))
  200     CONTINUE
          NEVAL = NEVAL + 8*NDIV*NUM
C
C   Add new contributions to RESULT and ABSERR.
C
          DO 220 I = 1,8*NDIV
              IF (I.LE.NDIV) THEN
                  INDEX = VACANT(I)
              ELSE
                  INDEX = SBRGNS + I
              END IF
              DO 210 J = 1,NUMFUN
                  RESULT(J) = RESULT(J) + VALUES(J,INDEX)
                  ABSERR(J) = ABSERR(J) + ERRORS(J,INDEX)
  210         CONTINUE
  220     CONTINUE
C
C   Store results in heap.
C
          DO 250 I = 1,8*NDIV
              IF (I.LE.NDIV) THEN
                  INDEX = VACANT(I)
              ELSE
                  INDEX = SBRGNS + I
              END IF
              J = SBRGNS + I
              CALL DTRTET(2,J,GREATE,LIST,INDEX)
  250     CONTINUE
          SBRGNS = SBRGNS + 8*NDIV
C
C   Check for termination.
C
          IF (SBRGNS.LT.MINSUB) THEN
              GO TO 110
          END IF
          DO 255 J = 1,NUMFUN
              IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)) .AND.
     +            ABSERR(J).GT.EPSABS) THEN
                  GO TO 110
              END IF
  255     CONTINUE
          IFAIL = 0
          COUNT = COUNT + 1
C
C   Else we did not succeed with the
C   given value of MAXSUB.
C
      ELSE
          IFAIL = 1
      END IF
C
C   Compute more accurate values of RESULT and ABSERR. This is done
C   to reduce the effect of roundoff on final results. Large
C   intermediate sums in the computation may course large, unnecessary
C   roundoff errors. Thus recomputing the sums of errors and estimates
C   and in addition grouping the sums in three groups should remove this
C   problem.
C
  499 CONTINUE
      DO 550 J = 1,NUMFUN
          RESULT(J) = 0
          ABSERR(J) = 0
          GROUP = SBRGNS** (1.d0/3.d0)
          IF (GROUP**3.NE.SBRGNS) THEN
              GROUP = GROUP + 1
          END IF
          LEFT = SBRGNS
          DO 520 I2 = 0,SBRGNS - 1,GROUP**2
              SUMVAL(2) = 0
              SUMERR(2) = 0
              DO 510 I1 = 0,MIN(LEFT,GROUP**2) - 1,GROUP
                  SUMVAL(1) = 0
                  SUMERR(1) = 0
                  DO 500 I = 1 + I1 + I2,MIN(LEFT,GROUP) + I1 + I2
                      SUMVAL(1) = SUMVAL(1) + VALUES(J,I)
                      SUMERR(1) = SUMERR(1) + ERRORS(J,I)
  500             CONTINUE
                  LEFT = LEFT - MIN(LEFT,GROUP)
                  SUMVAL(2) = SUMVAL(2) + SUMVAL(1)
                  SUMERR(2) = SUMERR(2) + SUMERR(1)
  510         CONTINUE
              RESULT(J) = RESULT(J) + SUMVAL(2)
              ABSERR(J) = ABSERR(J) + SUMERR(2)
  520     CONTINUE
  550 CONTINUE
C
C   Check again for termination
C
      IF (IFAIL.EQ.0) THEN
          IF (COUNT.LE.1) THEN
              DO 560 J = 1,NUMFUN
                  IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)) .AND.
     +                ABSERR(J).GT.EPSABS) THEN
                      GO TO 110
                  END IF
  560         CONTINUE
          END IF
      END IF
C
      NSUB = SBRGNS
      RETURN
C
C***END DADTET
C
      END
      SUBROUTINE DCHTET(NUMFUN,MDIV,VER,NUMTET,MINPTS,MAXPTS,EPSABS,
     +                  EPSREL,LENVER,NW,RESTAR,MAXSUB,MINSUB,IFAIL)
C***BEGIN PROLOGUE DCHTET
C***REFER TO DCUTET
C***PURPOSE  DCHTET checks the validity of the
C            input parameters to DCUTET.
C***DESCRIPTION
C            DCHTET computes MAXSUB, MINSUB and IFAIL as
C            functions of the input parameters to DCUTET,
C            and checks the validity of the input parameters to DCUTET.
C
C   ON ENTRY
C
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            MDIV is the number of tetrahedrons that are divided in
C            each subdivision step in DADTET.
C            MDIV is chosen default to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(8*MDIV,NPROC) = 0.
C     VER    Real array of dimension (3,4,LENVER).
C            VER(1,K,L), VER(2,K,L) and VER(3,K,L) are the x, y and z
C            coordinates respectively of vertex K in tetrahedron L.
C            On entry VER(*,*,L) must contain the vertices of the
C            NUMTET user specified tetrahedrons for L = 1,2,...,NUMTET.
C     NUMTET Integer.
C            The number of tetrahedrons.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            MAXPTS >= NUM*NUMTET and MAXPTS >= MINPTS
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     LENVER Integer.
C            Defines the length of the array VER.
C
C            We let
C            MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS.
C            MAXSUB = 7*(MAXPTS-NUM*NUMTET)/(8*NUM) + NUMTET
C            LENVER should then be greater or equal to MAXSUB.
C
C     NW     Integer.
C            Defines the length of the working array WORK.
C
C            We let
C            MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS.
C            MAXSUB = 7*(MAXPTS-NUM*NUMTET)/(8*NUM) + NUMTET
C            NW should then be greater or equal to
C            MAXSUB*(2*NUMFUN+1) + 6*MAX(8*MDIV,NUMTET)*NUMFUN + 1
C
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1,
C            then we restart a previous attempt.
C            In this case the only parameters for DCUTET that may
C            be changed (with respect to the previous call of DCUTET)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C
C
C   ON RETURN
C
C     MAXSUB Integer.
C            The maximum allowed number of subregions for the
C            given values of MAXPTS.
C     MINSUB Integer.
C            The minimum allowed number of subregions for the given
C            values of MINPTS.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit.
C            IFAIL = 2 if NUMFUN is less than 1.
C            IFAIL = 3 if the volume of one of the initially given
C                      tetrahedrons is zero.
C            IFAIL = 4 if MAXPTS is less than NUM*NUMTET.
C            IFAIL = 5 if MAXPTS is less than MINPTS.
C            IFAIL = 6 if EPSABS <= 0 and EPSREL <= 0.
C            IFAIL = 7 if LENVER is less than MAXSUB.
C            IFAIL = 8 if NW is less than MAXSUB*(2*NUMFUN+1) +
C                      NUMFUN*6*MAX(8*MDIV,NUMTET) + 1.
C            IFAIL = 9 if illegal value of RESTAR.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE DCHTET
C
C   Global variables.
C
      INTEGER NUMFUN,MDIV,NUMTET,MINPTS,MAXPTS,NW,MAXSUB,MINSUB,RESTAR,
     +        LENVER,IFAIL
      real(8) VER(3,4,NUMTET),EPSABS,EPSREL
C
C   Local variables. 
C   NUM is set equal to the number of evaluation points used
C   by the cubature rule over each tetrahedron.
C
      INTEGER LIMIT,J,NUM
      real(8) VOLUME
C
C***FIRST EXECUTABLE STATEMENT DCHTET
C
      IFAIL = 0
      NUM = 43
C
C   Compute MAXSUB and MINSUB.
C
      MAXSUB = 7* (MAXPTS-NUM*NUMTET)/ (8*NUM) + NUMTET
      MINSUB = 7* (MINPTS-NUM*NUMTET)/ (8*NUM) + NUMTET
      IF (MOD(7* (MINPTS-NUM*NUMTET),8*NUM).GT.0) THEN
          MINSUB = MINSUB + 1
      END IF
      MINSUB = MAX(NUMTET,MINSUB)
C
C   Check on positive NUMFUN.
C
      IF (NUMFUN.LT.1) THEN
          IFAIL = 2
          RETURN
      END IF
C
C   Check on legal volume of each subregion of integration.
C
      DO 10 J = 1,NUMTET
           VOLUME = (VER(1,2,J)-VER(1,1,J))*
     +     ((VER(2,3,J)-VER(2,1,J))*(VER(3,4,J)-VER(3,1,J)) -
     +      (VER(2,4,J)-VER(2,1,J))*(VER(3,3,J)-VER(3,1,J))) -
     +              (VER(2,2,J)-VER(2,1,J))*
     +     ((VER(1,3,J)-VER(1,1,J))*(VER(3,4,J)-VER(3,1,J)) - 
     +      (VER(1,4,J)-VER(1,1,J))*(VER(3,3,J)-VER(3,1,J))) +
     +              (VER(3,2,J)-VER(3,1,J))* 
     +     ((VER(1,3,J)-VER(1,1,J))*(VER(2,4,J)-VER(2,1,J)) -  
     +      (VER(1,4,J)-VER(1,1,J))*(VER(2,3,J)-VER(2,1,J)))
           VOLUME = VOLUME/6
          IF (VOLUME.EQ.0) THEN
              IFAIL = 3
              RETURN
          END IF
   10 CONTINUE
C
C   Check on MAXPTS >= NUM*NUMTET.
C
      IF (MAXPTS.LT.NUM*NUMTET) THEN
          IFAIL = 4
          RETURN
      END IF
C
C   Check on MAXPTS >= MINPTS.
C
      IF (MAXPTS.LT.MINPTS) THEN
          IFAIL = 5
          RETURN
      END IF
C
C   Check on legal accuracy requests.
C
      IF (EPSABS.LE.0 .AND. EPSREL.LE.0) THEN
          IFAIL = 6
          RETURN
      END IF
C
C   Check on big enough LENVER.
C
      IF (LENVER.LT.MAXSUB) THEN
          IFAIL = 7
          RETURN
      END IF
C
C   Check on big enough NW.
C
      LIMIT = MAXSUB* (2*NUMFUN+1) + 7*MAX(8*MDIV,NUMTET)*NUMFUN + 1
      IF (NW.LT.LIMIT) THEN
          IFAIL = 8
          RETURN
      END IF
C
C    Check on legal RESTAR.
C
      IF (RESTAR.NE.0 .AND. RESTAR.NE.1) THEN
          IFAIL = 9
      END IF
      RETURN
C
C***END DCHTET
C
      END
      SUBROUTINE DCUTET(FUNSUB,NUMFUN,VER,NUMTET,MINPTS,MAXPTS,EPSABS,
     +                  EPSREL,LENVER,NW,RESTAR,RESULT,ABSERR,NEVAL,
     +                  IFAIL,WORK,IWORK)
C***BEGIN PROLOGUE DCUTET
C***DATE WRITTEN   900220   (YYMMDD)
C***REVISION DATE  900307   (YYMMDD)
C***CATEGORY NO. H2B2A1/H2B2A2
C***KEYWORDS QUADRATURE, THREE DIMENSIONAL, ADAPTIVE, CUBATURE
C***AUTHORS
C            Jarle Berntsen, The Computing Centre,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544055
C            Email..  jarle@eik.ii.uib.no
C
C            Ronald Cools, Dept. of Computer Science,
C            Katholieke Universiteit Leuven, Celestijnenlaan 200A,
C            B-3030 Heverlee, Belgium
C            Phone..  32-16-201015 (3562)
C            Email..  ronald@cs.kuleuven.ac.be
C
C            Terje O. Espelid, Department of Informatics,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544180
C            Email..  terje@eik.ii.uib.no
C
C***PURPOSE  To compute the three-dimensional integral over a region
C            consisting of NUMTET tetrahedrons.
C***DESCRIPTION
C
C            The routine calculates an approximation to a given
C            vector of definite integrals
C
C            I  I  I (F ,F ,...,F )      DX(3)DX(2)DX(1),
C                      1  2      NUMFUN
C
C            where the region of integration is a collection of
C            NUMTET tetrahedrons and
C            where F = F (X(1),X(2),X(3)), J = 1,2,...,NUMFUN.
C                   J   J
C
C
C            A globally adaptive strategy is applied in order to
C            compute approximations RESULT(K)
C            hopefully satisfying, for each component of I, the
C            following claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C
C            DCUTET is a driver for the integration routine
C            DADTET.
C
C            DADTET repeatedly
C            subdivides the tetrahedrons with greatest estimated  errors
C            and estimates the integrals and the
C            errors over the new subtetrahedrons
C            until the error request
C            is met or MAXPTS function evaluations have been used.
C
C            A 43 point integration rule
C            with all evaluation points inside the tetrahedron
C            is applied. The rule has polynomial degree 8.
C
C            If the values of the input parameters EPSABS
C            or EPSREL are selected great enough,
C            an integration rule is applied over each tetrahedron and
C            the results are added up to give the approximations
C            RESULT(K). No further subdivision
C            of the tetrahedrons will then be applied.
C
C            When DCUTET computes estimates to a vector of
C            integrals, all components of the vector are given
C            the same treatment. That is, I(F ) and I(F ) for
C                                            J         K
C            J not equal to K, are estimated with the same
C            subdivision of the region of integration.
C            For integrals with enough similarity, we may save
C            time by applying DCUTET to all integrands in one call.
C            For integrals that varies continuously as functions of
C            some parameter, the estimates produced by DCUTET will
C            also vary continuously when the same subdivision is
C            applied to all components. This will generally not be
C            the case when the different components are given
C            separate treatment.
C
C            On the other hand this feature should be used with
C            caution when the different components of the integrals
C            require clearly different subdivisions.
C
C   ON ENTRY
C
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (X,NUMFUN,FUNVLS)
C            Input parameters:
C              X(1)   The x-coordinate of the evaluation point.
C              X(2)   The y-coordinate of the evaluation point.
C              X(3)   The z-coordinate of the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     NUMFUN Integer.
C            Number of components of the integral.
C     VER    Real array of dimension (3,4,LENVER).
C            VER(1,K,L), VER(2,K,L) and VER(3,K,L) are the x, y and z
C            coordinates respectively of vertex K in tetrahedron L.
C            On entry VER(*,*,L) must contain the vertices of the
C            NUMTET user specified tetrahedrons for L = 1,2,...,NUMTET.
C     NUMTET Integer.
C            The number of tetrahedrons.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is 43.
C
C            MAXPTS >= 43*NUMTET and MAXPTS >= MINPTS
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     LENVER Integer.
C            Defines the length of the array VER.
C
C            We let
C            MAXSUB denote the maximum allowed number of subregions
C            for the given value of MAXPTS.
C            MAXSUB = 7*(MAXPTS-43*NUMTET)/(8*43) + NUMTET
C            LENVER should be greater or equal to MAXSUB.
C
C     NW     Integer.
C            Defines the length of the working array WORK.
C
C            We let
C            MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS.
C            MAXSUB = 7*(MAXPTS-43*NUMTET)/(8*43) + NUMTET
C            NW should then be greater or equal to
C            MAXSUB*(2*NUMFUN+1) + MAX(8*6*MDIV,6*NUMTET)*NUMFUN + 1
C            where MDIV is the number of tetrahedrons that are divided
C            in each subdivision step.
C            MDIV is default set internally in DCUTET equal to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(8*MDIV,NPROC) = 0.
C
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1,
C            then we restart a previous attempt.
C            In this case the only parameters for DCUTET that may
C            be changed (with respect to the previous call of DCUTET)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute errors.
C     NEVAL  Integer.
C            Number of function evaluations used by DCUTET.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit.
C
C              ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
C              function evaluations for all values of K,
C              1 <= K <= NUMFUN .
C
C            IFAIL = 1 if MAXPTS was too small for DCUTET
C              to obtain the required accuracy. In this case DCUTET
C              returns values of RESULT with estimated absolute
C              errors ABSERR.
C            IFAIL = 2 if NUMFUN is less than 1.
C            IFAIL = 3 if the volume of one of the initially given
C                      tetrahedrons is zero.
C            IFAIL = 4 if MAXPTS is less than 43*NUMTET.
C            IFAIL = 5 if MAXPTS is less than MINPTS.
C            IFAIL = 6 if EPSABS <= 0 and EPSREL <= 0.
C            IFAIL = 7 if LENVER is less than MAXSUB.
C            IFAIL = 8 if NW is less than MAXSUB*(2*NUMFUN+1) +
C                      NUMFUN*MAX(8*6*MDIV,6*NUMTET) + 1.
C            IFAIL = 9 if unlegal RESTAR.
C     VER    Real array of dimension (3,4,LENVER).
C            On exit VER
C            contains the vertices of the tetrahedrons used to produce
C            the approximations to the integrals.
C     WORK   Real array of dimension NW.
C            Used as working storage.
C            WORK(NW) = NSUB, the number of subregions in the data
C            structure.
C            Let WRKSUB=(NW-1-NUMFUN*7*MAX(8*MDIV,NUMTET))/
C                        (2*NUMFUN+1).
C            WORK(1),...,WORK(NUMFUN*WRKSUB) contain
C              the estimated components of the integrals over the
C              subregions.
C            WORK(NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB) contain
C              the estimated errors over the subregions.
C            WORK(2*NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*
C              WRKSUB+WRKSUB) contain the greatest errors
C              in each subregion.
C            WORK((2*NUMFUN+1)*WRKSUB),...,WORK(NW - 1) is used as
C              temporary storage in DADTET.
C     IWORK  Integer array of dimension LENVER + MDIV.
C            Used as working storage.
C
C           SAMPLE PROGRAM
C   The following program will approximate the integral of
C                     exp(x*x+y*y+z*z)
C   over the tetrahedron (0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.)
C   and print out the values of the estimated integral, the estimated
C   error, the number of function evaluations and IFAIL.
C     PROGRAM DTEST1
C     real(8) VER(3,4,50),EPSABS,EPSREL,RESULT(1),ABSERR(1),
C    +                 WORK(1000)
C     INTEGER NUMFUN,NUMTET,MINPTS,MAXPTS,LENVER,NW,RESTAR,NEVAL,
C    +        IFAIL,IWORK(51)
C     EXTERNAL F
C     VER(1,1,1) = 0.D0
C     VER(1,2,1) = 1.D0
C     VER(1,3,1) = 0.D0
C     VER(1,4,1) = 0.D0
C     VER(2,1,1) = 0.D0
C     VER(2,2,1) = 0.D0
C     VER(2,3,1) = 1.D0
C     VER(2,4,1) = 0.D0
C     VER(3,1,1) = 0.D0
C     VER(3,2,1) = 0.D0
C     VER(3,3,1) = 0.D0
C     VER(3,4,1) = 1.D0
C     EPSABS = 0.D0
C     EPSREL = 1.D-5
C     NUMFUN = 1
C     NUMTET = 1
C     MINPTS = 0
C     MAXPTS = 1000
C     LENVER = 50
C     NW = 1000
C     RESTAR = 0
C     CALL DCUTET(F,NUMFUN,VER,NUMTET,MINPTS,MAXPTS,EPSABS,
C    +            EPSREL,LENVER,NW,RESTAR,RESULT,ABSERR,NEVAL,
C    +            IFAIL,WORK,IWORK)
C     WRITE(*,*)'RESULT=',RESULT(1)
C     WRITE(*,*)'ERROR =',ABSERR(1)
C     WRITE(*,*)'NEVAL =',NEVAL
C     WRITE(*,*)'IFAIL =',IFAIL
C     END
C     SUBROUTINE F(X,NUMFUN,FUNVLS)
C     INTEGER NUMFUN
C     real(8) X(3),FUNVLS(NUMFUN)
C     FUNVLS(1) = EXP(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
C     RETURN
C     END
C
C   Output produced on a SUN 3/50
C
C       RESULT=   0.22779999698986
C       ERROR =    1.9752829160298D-06
C       NEVAL =  43
C       IFAIL =  0
C
C***LONG DESCRIPTION
C
C   The information for each tetrahedron is contained in the
C   data structures VER,  WORK and IWORK.
C   VER contains the coordinates of the tetrahedrons.
C   When passed on to DADTET, WORK is split into four
C   arrays VALUES, ERRORS, GREATE and WORK2.
C   VALUES contains the estimated values of the integrals.
C   ERRORS contains the estimated errors.
C   GREATE contains the greatest estimated error for each tetrahedron.
C   The data structures for the tetrahedrons are in DADTET organized
C   as a heap and the size of GREATE(I) defines the position of
C   tetrahedron I in the heap. The heap is maintained by the program
C   DTRTET and we use a partially ordered list of pointers, saved in
C   IWORK. When passed on to DADTET, IWORK is split into two
C   arrays LIST and VACANT. LIST is a partially ordered list
C   such that GREATE(LIST(1)) is the maximum error estimate
C   for all subtetrahedrons in our datastructure.
C   VACANT is a work space needed in th updating process of the list.
C
C   The subroutine DADTET is written for efficient execution on shared
C   memory parallel computers. On a computer with NPROC processors we
C   will in each subdivision step divide MDIV tetrahedrons, where MDIV
C   is chosen such that MOD(8*MDIV,NPROC) = 0, in totally 8*MDIV new
C   tetrahedrons.  Each processor will then compute estimates of the
C   integrals and errors over 8*MDIV/NPROC tetrahedrons in each
C   subdivision step.
C   The subroutine for estimating the integral and the error over
C   each subregion, DRLTET, uses WORK2 as a work array.
C   We must make sure that each processor writes its results to
C   separate parts of the memory, and therefore the sizes of WORK and
C   WORK2 are functions of MDIV.
C   In order to achieve parallel processing of tetrahedrons, compiler
C   directives should be placed in front of the DO 20 and the DO 200
C   loops in DADTET on machines like Alliant and CRAY.
C
C***REFERENCES
C   J.Berntsen, R. Cools and T.O.Espelid, An Algorithm for Adaptive
C   Cubature Over a Collection of Tetrahedrons, to be published.
C
C***ROUTINES CALLED DCHTET,DADTET
C***END PROLOGUE DCUTET
C
C   Parameters
C
C   MDIV   Integer.
C          MDIV is the number of tetrahedrons that are divided in
C          each subdivision step in DADTET.
C          MDIV is chosen default to 1.
C          For efficient execution on parallel computers
C          with NPROC processors MDIV should be set equal to
C          the smallest integer such that MOD(8*MDIV,NPROC) = 0.
C
      INTEGER MDIV
      PARAMETER (MDIV=1)
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NUMFUN,NUMTET,MINPTS,MAXPTS,LENVER,NW,RESTAR,NEVAL,IFAIL,
     +        IWORK(LENVER+MDIV)
      real(8) VER(3,4,LENVER),EPSABS,EPSREL,RESULT(NUMFUN),
     +                 ABSERR(NUMFUN),WORK(NW)
C
C   Local variables.
C
C   MAXSUB Integer.
C          The maximum allowed number of subdivisions
C          for the given values of MAXPTS.
C   MINSUB Integer.
C          The minimum allowed number of subregions for the given
C          values of MINPTS.
C   WRKSUB Integer.
C          The maximum allowed number of subregions as a function
C          of NW, NUMFUN and MDIV. This determines the length
C          of the main work arrays.
C
      INTEGER MAXSUB,MINSUB,NSUB,LENW
      INTEGER WRKSUB,I1,I2,I3,I4,I5
C
C***FIRST EXECUTABLE STATEMENT DCUTET
C
C   Call DCHTET to compute MAXSUB and MINSUB,
C   and to check the input parameters.
C
      CALL DCHTET(NUMFUN,MDIV,VER,NUMTET,MINPTS,MAXPTS,EPSABS,EPSREL,
     +            LENVER,NW,RESTAR,MAXSUB,MINSUB,IFAIL)
C
C   All tests on input data passed?
C
      IF (IFAIL.NE.0) THEN
          RETURN
      END IF
C
C   Split up the work space.
C
      WRKSUB = (NW-1-NUMFUN*7*MAX(8*MDIV,NUMTET))/ (2*NUMFUN+1)
      I1 = 1
      I2 = I1 + WRKSUB*NUMFUN
      I3 = I2 + WRKSUB*NUMFUN
      I4 = I3 + WRKSUB
      I5 = I4 + 6*MAX(8*MDIV,NUMTET)*NUMFUN
C
C   On restart runs the number of subregions from the
C   previous call is assigned to NSUB.
C
      IF (RESTAR.EQ.1) THEN
          NSUB = WORK(NW)
      END IF
C
C   Compute the size of the temporary work space needed in DADTET.
C
      LENW = MAX(8*MDIV,NUMTET)*NUMFUN
      CALL DADTET(NUMFUN,MDIV,VER,NUMTET,MINSUB,MAXSUB,FUNSUB,EPSABS,
     +            EPSREL,LENVER,RESTAR,LENW,RESULT,ABSERR,NEVAL,NSUB,
     +            IFAIL,WORK(I1),WORK(I2),WORK(I3),WORK(I4),WORK(I5),
     +            IWORK(1),IWORK(1+LENVER))
      WORK(NW) = NSUB
      RETURN
C
C***END DCUTET
C
      END
      SUBROUTINE DORTET(TYPE,GENER,VER,NUMFUN,FUNSUB,SUMVAL,WORK)
C***BEGIN PROLOGUE DORTET
C***PURPOSE  To compute the sum of function values over all points
C            of an orbit.
C***DESCRIPTION
C   ON ENTRY
C
C   TYPE   Integer
C          The type of the orbit.
C   GENER  Integer array of dimension (3).
C          The generator for the orbit in homogeneous coordinates.
C   VER    Real array of dimension (3,4).
C          The coordinates of the vertices of the tetrahedron.
C          vertex i -> ( ver(1,i),ver(2,i),ver(3,i) )
C   NUMFUN Integer.
C          Number of components of the vector integrand.
C   FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (X,NUMFUN,FUNVLS)
C            Input parameters:
C              X(1)   The x-coordinate of the evaluation point.
C              X(2)   The y-coordinate of the evaluation point.
C              X(3)   The z-coordinate of the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of the vector integrand.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C   WORK   Real array of dimension (NUMFUN).
C          A work array
C
C   ON RETURN
C
C   SUMVAL Real array of dimension (NUMFUN).
C          The sum of function values over all points
C          of the given orbit.
C
C***END PROLOGUE  DORTET
C
C  Global variables
C
      INTEGER NUMFUN,TYPE
      real(8) GENER(3),SUMVAL(NUMFUN),VER(3,4),WORK(NUMFUN)
      EXTERNAL FUNSUB
C
C  Local variables
C
      INTEGER I,J,NUMBER
      real(8) X(3,12),Z1,Z2,Z3
C
C  The array 'gener' contains the necessary information
C   for the generator
C
C***FIRST EXECUTABLE STATEMENT  DORTET
      GO TO (100,110,120,130),TYPE + 1
C
C  Generator with homogeneous coordinates (1/4,1/4,1/4;1/4)
C
  100 NUMBER = 1
      X(1,1) = (VER(1,1)+VER(1,2)+VER(1,3)+VER(1,4))/4
      X(2,1) = (VER(2,1)+VER(2,2)+VER(2,3)+VER(2,4))/4
      X(3,1) = (VER(3,1)+VER(3,2)+VER(3,3)+VER(3,4))/4
      GO TO 200
C
C  Generator with homogeneous coordinates (z1,z2,z2;z2)
C
  110 NUMBER = 4
      Z1 = GENER(1)
      Z2 = GENER(2)
      DO 115 J = 1,3
          X(J,1) = Z1*VER(J,1) + Z2* (VER(J,2)+VER(J,3)+VER(J,4))
          X(J,2) = Z1*VER(J,2) + Z2* (VER(J,1)+VER(J,3)+VER(J,4))
          X(J,3) = Z1*VER(J,3) + Z2* (VER(J,2)+VER(J,1)+VER(J,4))
          X(J,4) = Z1*VER(J,4) + Z2* (VER(J,2)+VER(J,3)+VER(J,1))
  115 CONTINUE
      GO TO 200
C
C  Generator with homogeneous coordinates (z1,z1,z2;z2)
C
  120 NUMBER = 6
      Z1 = GENER(1)
      Z2 = GENER(2)
      DO 125 J = 1,3
          X(J,1) = Z1* (VER(J,1)+VER(J,2)) + Z2* (VER(J,3)+VER(J,4))
          X(J,2) = Z1* (VER(J,1)+VER(J,3)) + Z2* (VER(J,2)+VER(J,4))
          X(J,3) = Z1* (VER(J,1)+VER(J,4)) + Z2* (VER(J,3)+VER(J,2))
          X(J,4) = Z1* (VER(J,2)+VER(J,3)) + Z2* (VER(J,1)+VER(J,4))
          X(J,5) = Z1* (VER(J,2)+VER(J,4)) + Z2* (VER(J,1)+VER(J,3))
          X(J,6) = Z1* (VER(J,3)+VER(J,4)) + Z2* (VER(J,1)+VER(J,2))
  125 CONTINUE
      GO TO 200
C
C  Generator with homogeneous coordinates (z1,z2,z3;z3)
C
  130 NUMBER = 12
      Z1 = GENER(1)
      Z2 = GENER(2)
      Z3 = GENER(3)
      DO 135 J = 1,3
          X(J,1) = Z1*VER(J,1) + Z2*VER(J,2) + Z3* (VER(J,3)+VER(J,4))
          X(J,2) = Z1*VER(J,1) + Z2*VER(J,3) + Z3* (VER(J,2)+VER(J,4))
          X(J,3) = Z1*VER(J,1) + Z2*VER(J,4) + Z3* (VER(J,2)+VER(J,3))
          X(J,4) = Z1*VER(J,2) + Z2*VER(J,1) + Z3* (VER(J,3)+VER(J,4))
          X(J,5) = Z1*VER(J,2) + Z2*VER(J,3) + Z3* (VER(J,1)+VER(J,4))
          X(J,6) = Z1*VER(J,2) + Z2*VER(J,4) + Z3* (VER(J,1)+VER(J,3))
          X(J,7) = Z1*VER(J,3) + Z2*VER(J,1) + Z3* (VER(J,2)+VER(J,4))
          X(J,8) = Z1*VER(J,3) + Z2*VER(J,2) + Z3* (VER(J,1)+VER(J,4))
          X(J,9) = Z1*VER(J,3) + Z2*VER(J,4) + Z3* (VER(J,1)+VER(J,2))
          X(J,10) = Z1*VER(J,4) + Z2*VER(J,1) + Z3* (VER(J,2)+VER(J,3))
          X(J,11) = Z1*VER(J,4) + Z2*VER(J,2) + Z3* (VER(J,1)+VER(J,3))
          X(J,12) = Z1*VER(J,4) + Z2*VER(J,3) + Z3* (VER(J,1)+VER(J,2))
  135 CONTINUE
  200 CONTINUE
      CALL FUNSUB(X(1,1),NUMFUN,SUMVAL)
      DO 220 J = 2,NUMBER
          CALL FUNSUB(X(1,J),NUMFUN,WORK)
          DO 210 I = 1,NUMFUN
              SUMVAL(I) = SUMVAL(I) + WORK(I)
  210     CONTINUE
  220 CONTINUE
C
C***END DORTET
C
      RETURN
      END
      SUBROUTINE DRLTET(VER,NUMFUN,FUNSUB,NULL,BASVAL,RGNERR,GREATE,
     +                  SUMVAL)
C***BEGIN PROLOGUE DRLTET
C***REFER TO DCUTET
C***PURPOSE  To compute basic integration rule values and
C            corresponding error estimates.
C***DESCRIPTION DRLTET computes basic integration rule values
C            for a vector of integrands over a tetrahedron.
C            DRLTET also computes estimates for the errors by
C            using several null rule approximations.
C   ON ENTRY
C
C   VER    Real array of dimension (3,4).
C          The coordinates of the vertices of the tetrahedron.
C          vertex i -> ( ver(1,i),ver(2,i),ver(3,i) )
C   NUMFUN Integer.
C          Number of components of the vector integrand.
C   FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (X,NUMFUN,FUNVLS)
C            Input parameters:
C              X(1)   The x-coordinate of the evaluation point.
C              X(2)   The y-coordinate of the evaluation point.
C              X(3)   The z-coordinate of the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of the vector integrand.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C   NULL   Real array of dimension (NUMFUN,6).
C          A work array.
C   SUMVAL Real array of dimension (NUMFUN).
C          A work array
C
C   ON RETURN
C
C   BASVAL Real array of dimension NUMFUN.
C          The values for the basic rule for each component
C          of the integrand.
C   RGNERR Real array of dimension NUMFUN.
C          The error estimates for each component of the integrand.
C   GREATE Real.
C          The greatest error component of the integrand.
C
C***REFERENCES Beckers M., A. Haegemans,
C              The construction of cubature formula for the tetrahedron
C              Report TW128, K.U. Leuven (1990).
C***ROUTINES CALLED
C                   D1MACH,DORTET,FUNSUB
C***END PROLOGUE DRLTET
C
C   Parameters
C
C   ORBITS Integer
C          The number of orbits of the cubature formula and null rules
C   CRIVAL Real
C          The decision to choose the optimistic part of the error
C          estimator is based on CRIVAL
C   FACOPT Real
C          FACOPT is the safety coefficient used in the optimistic part
C          of the error estimator.
C   FACMED Real
C          FACMED is the safety coefficient used in the non-optimistic
C          part of the error estimator. FACMED is related to CRIVAL 
C          and FACOPT.
C
      INTEGER ORBITS
      real(8) CRIVAL,FACOPT,FACMED
      PARAMETER (ORBITS=7)
      PARAMETER (CRIVAL=0.5,FACMED=5,FACOPT=FACMED/CRIVAL)
C
C   Global variables.
C
      INTEGER NUMFUN
      real(8) VER(3,4),BASVAL(NUMFUN),RGNERR(NUMFUN),GREATE,
     +                 SUMVAL(NUMFUN),NULL(NUMFUN,6)
      EXTERNAL FUNSUB
C
C   Local variables.
C
C   K      Integer array of dimension (0:3) that contains the structure
C          parameters. K(I) = number of orbits of type I.
C   TYPE1  Real array of dimension (K(1)).
C          Contains the first homogeneous coordinate of the generators
C          of type 1
C   TYPE2  Real array of dimension (K(2)).
C          Contains the first homogeneous coordinate of the generators
C          of type 2
C   TYPE3  Real array of dimension (2,K(2)).
C          Contains the first two homogeneous coordinates of
C          the generators of type 3.
C   WEIGHT Real array of dimension (9,ORBITS).
C          The weights of the cubature formula and the null rules.
C          WEIGHT(1,1) ,..., WEIGHT(1,ORBITS) are the weights of the
C                cubature formula
C          WEIGHT(I,1) ,..., WEIGHT(I,ORBITS) for I > 1, are the weights
C                of the null rules
C
C
      INTEGER TYPE,NR,K(0:3),I,J,P
      real(8) Z(3),TYPE1(3),TYPE2(1),TYPE3(2,2),NOISE,D1MACH,
     +                 WEIGHT(7,ORBITS),VOLUME,DEG4,DEG3,DEG1,R2,R1,R,
     +                 TRES
C
C  Cubature formula of degree 8 with 43 points
C
C  Structure parameters
C
      DATA (K(I),I=0,3)/1,3,1,2/
C
C  Information for the generators
C
      DATA (TYPE1(I),I=1,3)/0.379510205167980387748057300876D+00,
     +     0.753689235068359830728182577696D+00,
     +     0.982654148484406008240470085259D+00/
      DATA (TYPE2(I),I=1,1)/0.449467259981105775574375471447D+00/
      DATA ((TYPE3(I,J),I=1,2),J=1,2)/
     +     0.506227344977843677082264893876D+00,
     +     0.356395827885340437169173969841D-01,
     +     0.736298458958971696943019005441D+00,
     +     0.190486041934633455699433285302D+00/
C
C  Weights of the cubature formula
C
      DATA (WEIGHT(1,I),I=1,ORBITS)/-
     +     0.123001131951839495043519102752D+00,
     +     0.855018349372014074906384482699D-01,
     +     0.118021998788034059253768205083D-01,
     +     0.101900465455732427902646736855D-02,
     +     0.274781029468036908044610867719D-01,
     +     0.342269148520915110408153517904D-01,
     +     0.128431148469725555789001180031D-01/
C
C  Weights of the null rule of degree 5
C
      DATA (WEIGHT(2,I),I=1,ORBITS)/0.211921237628032658308230999090D+00
     +     ,-0.660207516445726284649283745987D-01,
     +     0.225058824086711710443385047042D-01,
     +     -0.375962972067425589765730699401D-03,
     +     0.710066020561055159657284834784D-02,
     +     0.156515256061747694921427149028D-02,
     +     -0.814530839643584660306807872526D-02/
C
C  Weights of null rule of degree 4
C
      DATA (WEIGHT(3,I),I=1,ORBITS)/-
     +     0.508105488137100551376844924797D-01,
     +     0.104596681151665328209751420525D-01,
     +     0.927471438532788763594989973184D-01,
     +     0.210489990008917994323967321174D-02,
     +     0.379184172251962722213408547663D-01,
     +     -0.111747242913563605790923001557D-01,
     +     -0.386541758762774673113423570465D-01/
C
C  Weights of first null rule of degree 3
C
      DATA (WEIGHT(4,I),I=1,ORBITS)/-
     +     0.775992773232808462404390159802D-01,
     +     -0.527453289659022924847298408064D-01,
     +     0.145876238555932704488677626554D-01,
     +     0.739374873393616192857532718429D-02,
     +     -0.374618791364332892611678523428D-01,
     +     0.538502846550653076078817013885D-01,
     +     -0.183980865177843057548322735665D-01/
C
C  Weights of second null rule of degree 3
C
      DATA (WEIGHT(5,I),I=1,ORBITS)/0.181767621501470154602720474731D-01
     +     ,0.179938831310058580533178529022D-01,
     +     0.713210362750414891598257378898D-01,
     +     -0.443935688958258805893448212636D-01,
     +     -0.657639036547720234169662790056D-01,
     +     -0.101551807522541414699808460583D-01,
     +     0.265486188970540796821750584204D-01/
C
C  Weights of null rule of degree 2
C
      DATA (WEIGHT(6,I),I=1,ORBITS)/-
     +     0.867629853722843888927184699428D-01,
     +     -0.715881271235661902772072127812D-01,
     +     0.886720767790426261677273459523D-02,
     +     -0.577885573028655167063092577589D-01,
     +     0.430310167581202031805055255554D-01,
     +     -0.606467834856775537069463817445D-02,
     +     0.319492443333738343104163265406D-01/
C
C   Weights of null rule of degree 1
C
      DATA (WEIGHT(7,I),I=1,ORBITS)/0.510374015624925451319499382594D-01
     +     ,0.463998830432033721597269299429D-01,
     +     -0.191086148397852799983451475821D-01,
     +     -0.973768821003670776204287367278D-01,
     +     0.180352562073914141268335496511D-01,
     +     0.277129527093489643801598303110D-01,
     +     -0.176218263109360550515567818653D-01/
C
C  The number of points used by the cubature formula is
C     NUM    = K(0) + 4*K(1) + 6*K(2) + 12*K(3) = 43
C
C***FIRST EXECUTABLE STATEMENT DRLTRI
      TRES = 50*D1MACH(4)
C
C  Compute the volume of the tetrahedron
C
      VOLUME = (VER(1,2)-VER(1,1))* ((VER(2,3)-VER(2,1))*
     +         (VER(3,4)-VER(3,1))- (VER(2,4)-VER(2,1))*
     +         (VER(3,3)-VER(3,1))) - (VER(2,2)-VER(2,1))*
     +         ((VER(1,3)-VER(1,1))* (VER(3,4)-VER(3,1))-
     +         (VER(1,4)-VER(1,1))* (VER(3,3)-VER(3,1))) +
     +         (VER(3,2)-VER(3,1))* ((VER(1,3)-VER(1,1))*
     +         (VER(2,4)-VER(2,1))- (VER(1,4)-VER(1,1))*
     +         (VER(2,3)-VER(2,1)))
      VOLUME = ABS(VOLUME)/6
C
C  Initialise BASVAL and NULL
C
      DO 10 J = 1,NUMFUN
          BASVAL(J) = 0.d0
          DO 20 I = 1,6
              NULL(J,I) = 0.d0
   20     CONTINUE
   10 CONTINUE
      P = 1
C
C  Compute contributions from orbits with 1, 4, 6 and 12  points
C
      DO 4 TYPE = 0,3
          DO 3 NR = 1,K(TYPE)
              IF (TYPE.EQ.1) THEN
C                 Generator ( z(1) , z(2), z(2) ; z(2) )
                  Z(1) = TYPE1(NR)
                  Z(2) = (1-Z(1))/3
              ELSE IF (TYPE.EQ.2) THEN
C                 Generator ( z(1) , z(1), z(2) ; z(2) )
                  Z(1) = TYPE2(NR)
                  Z(2) = (1-2*Z(1))/2
              ELSE
C                 Generator ( z(1) , z(2), z(3) ; z(3) )
                  Z(1) = TYPE3(1,NR)
                  Z(2) = TYPE3(2,NR)
                  Z(3) = (1-Z(1)-Z(2))/2
              END IF
C                                RGNERR is work array for DORTHE
              CALL DORTET(TYPE,Z,VER,NUMFUN,FUNSUB,SUMVAL,RGNERR)
              DO 2 J = 1,NUMFUN
                  BASVAL(J) = BASVAL(J) + WEIGHT(1,P)*SUMVAL(J)
                  DO 1 I = 1,6
                      NULL(J,I) = NULL(J,I) + WEIGHT(I+1,P)*SUMVAL(J)
    1             CONTINUE
    2         CONTINUE
              P = P + 1
    3     CONTINUE
    4 CONTINUE
C
C  Compute error estimates
C
      GREATE = 0
      DO 30 J = 1,NUMFUN
          NOISE = ABS(BASVAL(J))*TRES
          DEG4 = SQRT(NULL(J,1)**2+NULL(J,2)**2)
          DEG3 = SQRT(NULL(J,3)**2+NULL(J,4)**2)
          IF (DEG4.LE.NOISE) THEN
              RGNERR(J) = NOISE
          ELSE
              DEG1 = SQRT(NULL(J,5)**2+NULL(J,6)**2)
              IF (DEG3.NE.0) THEN
                  R1 = (DEG4/DEG3)**2
              ELSE
                  R1 = 1
              END IF
              IF (DEG1.NE.0) THEN
                  R2 = DEG3/DEG1
              ELSE
                  R2 = 1
              END IF
              R = MAX(R1,R2)
              IF (R.GE.CRIVAL) THEN
                  RGNERR(J) = FACMED*R*DEG4
              ELSE
                  RGNERR(J) = FACOPT* (R**2)*DEG4
              END IF
              RGNERR(J) = MAX(NOISE,RGNERR(J))
          END IF
          RGNERR(J) = VOLUME*RGNERR(J)
          BASVAL(J) = VOLUME*BASVAL(J)
          GREATE = MAX(GREATE,RGNERR(J))
   30 CONTINUE
C
C***END DRLTET
C
      RETURN
      END
      PROGRAM DTEST1
      real(8) VER(3,4,50),EPSABS,EPSREL,RESULT(1),ABSERR(1),
     +                 WORK(1000)
      INTEGER NUMFUN,NUMTET,MINPTS,MAXPTS,LENVER,NW,RESTAR,NEVAL,
     +        IFAIL,IWORK(51)
      EXTERNAL F
      VER(1,1,1) = 0.D0
      VER(1,2,1) = 1.D0
      VER(1,3,1) = 0.D0
      VER(1,4,1) = 0.D0
      VER(2,1,1) = 0.D0
      VER(2,2,1) = 0.D0
      VER(2,3,1) = 1.D0
      VER(2,4,1) = 0.D0
      VER(3,1,1) = 0.D0
      VER(3,2,1) = 0.D0
      VER(3,3,1) = 0.D0
      VER(3,4,1) = 1.D0
      EPSABS = 0.D0
      EPSREL = 1.D-5
      NUMFUN = 1
      NUMTET = 1
      MINPTS = 0
      MAXPTS = 1000
      LENVER = 50
      NW = 1000
      RESTAR = 0
      CALL DCUTET(F,NUMFUN,VER,NUMTET,MINPTS,MAXPTS,EPSABS,
     +            EPSREL,LENVER,NW,RESTAR,RESULT,ABSERR,NEVAL,
     +            IFAIL,WORK,IWORK)
      WRITE(*,*)'RESULT=',RESULT(1)
      WRITE(*,*)'ERROR =',ABSERR(1)
      WRITE(*,*)'NEVAL =',NEVAL
      WRITE(*,*)'IFAIL =',IFAIL
      END
      SUBROUTINE F(X,NUMFUN,FUNVLS)
      INTEGER NUMFUN
      real(8) X(3),FUNVLS(NUMFUN)
      FUNVLS(1) = EXP(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
      RETURN
      END
      PROGRAM DTEST2
C
C   DTEST2 checks that DCUTET integrates to machine precision
C   all monomials of degree less than or equal to 8.
C   DTEST2 checks that DCUTET integrates to machine precision
C   all invariant polynomials that were used to construct the
C   cubature formula of degree 8 used by DRLTET.
C   DTEST2 checks that the restart feature of DCUTET works and
C   that the routine works with vectors of functions.
C
C   Output produced on a SUN 3/50
C      
C     INTEGRATE MONOMIALS OF DEGREE 0
C          SUBROUTINE CALLS =     43, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.11E-13   1.0000000000000000
C      
C     INTEGRATE MONOMIALS OF DEGREE 1
C          SUBROUTINE CALLS =     43, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.11E-13   1.0000000000000000
C             2       0.00E+00       0.11E-13   1.0000000000000000
C             3       0.22E-15       0.11E-13   1.0000000000000002
C      
C     INTEGRATE MONOMIALS OF DEGREE 2
C          SUBROUTINE CALLS =     43, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.11E-13   1.0000000000000000
C             2       0.00E+00       0.11E-13   1.0000000000000000
C             3       0.00E+00       0.11E-13   1.0000000000000000
C             4       0.00E+00       0.11E-13   1.0000000000000000
C             5       0.00E+00       0.11E-13   1.0000000000000000
C             6       0.00E+00       0.11E-13   1.0000000000000000
C      
C     INTEGRATE MONOMIALS OF DEGREE 3
C          SUBROUTINE CALLS =     43, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.11E-13   1.0000000000000000
C             2       0.00E+00       0.11E-13   1.0000000000000000
C             3       0.22E-15       0.11E-13   0.9999999999999998
C             4       0.00E+00       0.11E-13   1.0000000000000000
C             5       0.00E+00       0.11E-13   1.0000000000000000
C             6       0.00E+00       0.11E-13   1.0000000000000000
C             7       0.00E+00       0.11E-13   1.0000000000000000
C             8       0.22E-15       0.11E-13   0.9999999999999998
C             9       0.00E+00       0.11E-13   1.0000000000000000
C            10       0.00E+00       0.11E-13   1.0000000000000000
C      
C     INTEGRATE MONOMIALS OF DEGREE 4
C          SUBROUTINE CALLS =     43, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.11E-13   1.0000000000000000
C             2       0.22E-15       0.11E-13   0.9999999999999998
C             3       0.00E+00       0.11E-13   1.0000000000000000
C             4       0.00E+00       0.11E-13   1.0000000000000000
C             5       0.00E+00       0.11E-13   1.0000000000000000
C             6       0.00E+00       0.11E-13   1.0000000000000000
C             7       0.00E+00       0.11E-13   1.0000000000000000
C             8       0.00E+00       0.11E-13   1.0000000000000000
C             9       0.22E-15       0.11E-13   0.9999999999999998
C            10       0.00E+00       0.11E-13   1.0000000000000000
C            11       0.00E+00       0.11E-13   1.0000000000000000
C            12       0.00E+00       0.11E-13   1.0000000000000000
C            13       0.22E-15       0.11E-13   0.9999999999999998
C            14       0.22E-15       0.11E-13   0.9999999999999998
C            15       0.00E+00       0.11E-13   1.0000000000000000
C      
C     INTEGRATE MONOMIALS OF DEGREE 5
C          SUBROUTINE CALLS =  12771, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.11E-15       0.33E-06   0.9999999999999999
C             2       0.22E-15       0.62E-05   0.9999999999999998
C             3       0.22E-15       0.20E-05   0.9999999999999998
C             4       0.22E-15       0.46E-05   0.9999999999999998
C             5       0.11E-15       0.58E-05   0.9999999999999999
C             6       0.00E+00       0.31E-06   1.0000000000000000
C             7       0.22E-15       0.26E-06   0.9999999999999998
C             8       0.22E-15       0.64E-05   0.9999999999999998
C             9       0.00E+00       0.51E-05   1.0000000000000000
C            10       0.22E-15       0.89E-05   0.9999999999999998
C            11       0.11E-15       0.43E-05   0.9999999999999999
C            12       0.11E-15       0.94E-07   0.9999999999999999
C            13       0.00E+00       0.31E-06   1.0000000000000000
C            14       0.00E+00       0.64E-05   1.0000000000000000
C            15       0.22E-15       0.61E-05   0.9999999999999998
C            16       0.22E-15       0.75E-07   0.9999999999999998
C            17       0.22E-15       0.60E-05   0.9999999999999998
C            18       0.22E-15       0.20E-05   0.9999999999999998
C            19       0.22E-15       0.26E-06   0.9999999999999998
C            20       0.11E-15       0.45E-05   0.9999999999999999
C            21       0.11E-15       0.28E-06   0.9999999999999999
C      
C     INTEGRATE MONOMIALS OF DEGREE 6
C          SUBROUTINE CALLS =  16899, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.55E-06   1.0000000000000000
C             2       0.00E+00       0.25E-05   1.0000000000000000
C             3       0.22E-15       0.61E-05   0.9999999999999998
C             4       0.11E-15       0.24E-05   0.9999999999999999
C             5       0.22E-15       0.95E-05   0.9999999999999998
C             6       0.00E+00       0.30E-05   1.0000000000000000
C             7       0.00E+00       0.46E-06   1.0000000000000000
C             8       0.00E+00       0.26E-05   1.0000000000000000
C             9       0.22E-15       0.55E-05   0.9999999999999998
C            10       0.11E-15       0.45E-05   0.9999999999999999
C            11       0.22E-15       0.54E-05   0.9999999999999998
C            12       0.22E-15       0.48E-05   0.9999999999999998
C            13       0.00E+00       0.36E-05   1.0000000000000000
C            14       0.11E-15       0.16E-05   0.9999999999999999
C            15       0.11E-15       0.56E-05   0.9999999999999999
C            16       0.22E-15       0.99E-05   0.9999999999999998
C            17       0.11E-15       0.59E-05   0.9999999999999999
C            18       0.22E-15       0.10E-04   0.9999999999999998
C            19       0.00E+00       0.23E-05   1.0000000000000000
C            20       0.11E-15       0.55E-05   0.9999999999999999
C            21       0.22E-15       0.53E-05   0.9999999999999998
C            22       0.11E-15       0.26E-05   0.9999999999999999
C            23       0.22E-15       0.15E-05   0.9999999999999998
C            24       0.22E-15       0.63E-05   0.9999999999999998
C            25       0.11E-15       0.62E-05   0.9999999999999999
C            26       0.00E+00       0.23E-05   1.0000000000000000
C            27       0.00E+00       0.26E-05   1.0000000000000000
C            28       0.00E+00       0.56E-06   1.0000000000000000
C      
C     INTEGRATE MONOMIALS OF DEGREE 7
C          SUBROUTINE CALLS =  23091, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.11E-15       0.24E-06   0.9999999999999999
C             2       0.00E+00       0.76E-06   1.0000000000000000
C             3       0.11E-15       0.20E-05   0.9999999999999999
C             4       0.00E+00       0.76E-06   1.0000000000000000
C             5       0.00E+00       0.76E-06   1.0000000000000000
C             6       0.00E+00       0.39E-05   1.0000000000000000
C             7       0.00E+00       0.38E-05   1.0000000000000000
C             8       0.00E+00       0.32E-06   1.0000000000000000
C             9       0.00E+00       0.93E-06   1.0000000000000000
C            10       0.22E-15       0.18E-05   0.9999999999999998
C            11       0.11E-15       0.30E-05   0.9999999999999999
C            12       0.00E+00       0.49E-05   1.0000000000000000
C            13       0.11E-15       0.49E-05   0.9999999999999999
C            14       0.11E-15       0.45E-05   0.9999999999999999
C            15       0.00E+00       0.82E-05   1.0000000000000000
C            16       0.22E-15       0.13E-05   0.9999999999999998
C            17       0.22E-15       0.44E-05   0.9999999999999998
C            18       0.11E-15       0.88E-05   0.9999999999999999
C            19       0.00E+00       0.69E-05   1.0000000000000000
C            20       0.00E+00       0.81E-05   1.0000000000000000
C            21       0.11E-15       0.80E-05   0.9999999999999999
C            22       0.00E+00       0.15E-05   1.0000000000000000
C            23       0.00E+00       0.76E-05   1.0000000000000000
C            24       0.11E-15       0.85E-05   0.9999999999999999
C            25       0.00E+00       0.77E-05   1.0000000000000000
C            26       0.00E+00       0.56E-05   1.0000000000000000
C            27       0.00E+00       0.14E-05   1.0000000000000000
C            28       0.00E+00       0.40E-05   1.0000000000000000
C            29       0.22E-15       0.30E-05   0.9999999999999998
C            30       0.00E+00       0.13E-05   1.0000000000000000
C            31       0.11E-15       0.19E-05   0.9999999999999999
C            32       0.22E-15       0.22E-05   0.9999999999999998
C            33       0.22E-15       0.93E-06   0.9999999999999998
C            34       0.00E+00       0.13E-05   1.0000000000000000
C            35       0.00E+00       0.98E-06   1.0000000000000000
C            36       0.00E+00       0.30E-06   1.0000000000000000
C      
C     INTEGRATE MONOMIALS OF DEGREE 8
C          SUBROUTINE CALLS =  39947, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.89E-07   1.0000000000000000
C             2       0.00E+00       0.11E-05   1.0000000000000000
C             3       0.00E+00       0.12E-05   1.0000000000000000
C             4       0.22E-15       0.15E-06   0.9999999999999998
C             5       0.22E-15       0.70E-06   0.9999999999999998
C             6       0.00E+00       0.11E-06   1.0000000000000000
C             7       0.00E+00       0.93E-06   1.0000000000000000
C             8       0.00E+00       0.62E-06   1.0000000000000000
C             9       0.00E+00       0.89E-07   1.0000000000000000
C            10       0.00E+00       0.31E-06   1.0000000000000000
C            11       0.00E+00       0.17E-05   1.0000000000000000
C            12       0.00E+00       0.39E-05   1.0000000000000000
C            13       0.00E+00       0.43E-05   1.0000000000000000
C            14       0.22E-15       0.36E-05   0.9999999999999998
C            15       0.00E+00       0.27E-05   1.0000000000000000
C            16       0.00E+00       0.14E-05   1.0000000000000000
C            17       0.00E+00       0.74E-06   1.0000000000000000
C            18       0.00E+00       0.36E-06   1.0000000000000000
C            19       0.00E+00       0.21E-05   1.0000000000000000
C            20       0.11E-15       0.90E-05   0.9999999999999999
C            21       0.00E+00       0.39E-05   1.0000000000000000
C            22       0.22E-15       0.54E-05   0.9999999999999998
C            23       0.00E+00       0.49E-05   1.0000000000000000
C            24       0.00E+00       0.23E-05   1.0000000000000000
C            25       0.00E+00       0.17E-06   1.0000000000000000
C            26       0.00E+00       0.55E-05   1.0000000000000000
C            27       0.00E+00       0.82E-05   1.0000000000000000
C            28       0.00E+00       0.41E-05   1.0000000000000000
C            29       0.00E+00       0.33E-05   1.0000000000000000
C            30       0.00E+00       0.20E-06   1.0000000000000000
C            31       0.00E+00       0.65E-06   1.0000000000000000
C            32       0.00E+00       0.60E-05   1.0000000000000000
C            33       0.22E-15       0.95E-05   0.9999999999999998
C            34       0.00E+00       0.30E-05   1.0000000000000000
C            35       0.00E+00       0.76E-06   1.0000000000000000
C            36       0.00E+00       0.23E-06   1.0000000000000000
C            37       0.00E+00       0.23E-05   1.0000000000000000
C            38       0.00E+00       0.45E-05   1.0000000000000000
C            39       0.00E+00       0.16E-06   1.0000000000000000
C            40       0.22E-15       0.67E-06   0.9999999999999998
C            41       0.00E+00       0.14E-05   1.0000000000000000
C            42       0.00E+00       0.76E-06   1.0000000000000000
C            43       0.00E+00       0.39E-06   1.0000000000000000
C            44       0.00E+00       0.10E-05   1.0000000000000000
C            45       0.00E+00       0.11E-06   1.0000000000000000
C      
C     INTEGRATE INVARIANT POLYNOMIALS
C          SUBROUTINE CALLS =   1763, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.11E-13   1.0000000000000000
C             2       0.00E+00       0.11E-13   1.0000000000000000
C             3       0.00E+00       0.11E-13   1.0000000000000000
C             4       0.11E-15       0.11E-13   0.9999999999999999
C             5       0.22E-15       0.11E-13   1.0000000000000002
C             6       0.11E-15       0.35E-07   0.9999999999999999
C             7       0.11E-15       0.37E-03   0.9999999999999999
C             8       0.00E+00       0.63E-02   1.0000000000000000
C             9       0.00E+00       0.32E-01   1.0000000000000000
C            10       0.22E-15       0.33E-05   1.0000000000000002
C            11       0.22E-15       0.35E-04   1.0000000000000002
C            12       0.11E-15       0.92E-03   0.9999999999999999
C            13       0.44E-15       0.20E-01   0.9999999999999996
C            14       0.44E-15       0.35E-01   1.0000000000000004
C            15       0.00E+00       0.11E-01   1.0000000000000000
C          SUBROUTINE CALLS =   2795, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.11E-13   1.0000000000000000
C             2       0.22E-15       0.11E-13   1.0000000000000002
C             3       0.00E+00       0.13E-13   1.0000000000000000
C             4       0.11E-15       0.11E-13   0.9999999999999999
C             5       0.22E-15       0.11E-13   1.0000000000000002
C             6       0.11E-15       0.47E-05   0.9999999999999999
C             7       0.22E-15       0.10E-03   0.9999999999999998
C             8       0.00E+00       0.16E-02   1.0000000000000000
C             9       0.00E+00       0.80E-02   1.0000000000000000
C            10       0.00E+00       0.22E-04   1.0000000000000000
C            11       0.22E-15       0.48E-03   1.0000000000000002
C            12       0.00E+00       0.25E-03   1.0000000000000000
C            13       0.44E-15       0.52E-02   0.9999999999999996
C            14       0.22E-15       0.91E-02   1.0000000000000002
C            15       0.00E+00       0.33E-02   1.0000000000000000
C          SUBROUTINE CALLS =   3139, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.11E-13   1.0000000000000000
C             2       0.22E-15       0.11E-13   1.0000000000000002
C             3       0.00E+00       0.13E-13   1.0000000000000000
C             4       0.11E-15       0.11E-13   0.9999999999999999
C             5       0.00E+00       0.11E-13   1.0000000000000000
C             6       0.11E-15       0.57E-05   0.9999999999999999
C             7       0.22E-15       0.73E-05   0.9999999999999998
C             8       0.11E-15       0.38E-04   0.9999999999999999
C             9       0.00E+00       0.61E-04   1.0000000000000000
C            10       0.00E+00       0.27E-04   1.0000000000000000
C            11       0.22E-15       0.52E-03   1.0000000000000002
C            12       0.22E-15       0.23E-04   0.9999999999999998
C            13       0.56E-15       0.20E-03   0.9999999999999994
C            14       0.44E-15       0.42E-03   1.0000000000000004
C            15       0.00E+00       0.64E-03   1.0000000000000000
C          SUBROUTINE CALLS =  17587, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.11E-13   1.0000000000000000
C             2       0.11E-15       0.11E-13   0.9999999999999999
C             3       0.00E+00       0.12E-13   1.0000000000000000
C             4       0.00E+00       0.11E-13   1.0000000000000000
C             5       0.00E+00       0.11E-13   1.0000000000000000
C             6       0.00E+00       0.47E-06   1.0000000000000000
C             7       0.00E+00       0.27E-06   1.0000000000000000
C             8       0.22E-15       0.39E-05   1.0000000000000002
C             9       0.00E+00       0.89E-05   1.0000000000000000
C            10       0.11E-15       0.40E-05   0.9999999999999999
C            11       0.11E-15       0.17E-04   0.9999999999999999
C            12       0.00E+00       0.57E-05   1.0000000000000000
C            13       0.44E-15       0.60E-04   1.0000000000000004
C            14       0.33E-15       0.59E-04   0.9999999999999997
C            15       0.22E-15       0.97E-04   0.9999999999999998
C          SUBROUTINE CALLS =  26187, IFAIL =  0
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.11E-13   1.0000000000000000
C             2       0.22E-15       0.11E-13   0.9999999999999998
C             3       0.00E+00       0.12E-13   1.0000000000000000
C             4       0.11E-15       0.11E-13   0.9999999999999999
C             5       0.22E-15       0.11E-13   0.9999999999999998
C             6       0.22E-15       0.44E-06   0.9999999999999998
C             7       0.11E-15       0.63E-07   0.9999999999999999
C             8       0.22E-15       0.58E-06   0.9999999999999998
C             9       0.00E+00       0.47E-05   1.0000000000000000
C            10       0.11E-15       0.16E-05   0.9999999999999999
C            11       0.22E-15       0.91E-05   0.9999999999999998
C            12       0.11E-15       0.24E-06   0.9999999999999999
C            13       0.22E-15       0.26E-05   1.0000000000000002
C            14       0.56E-15       0.25E-05   0.9999999999999994
C            15       0.22E-15       0.59E-05   0.9999999999999998
C          SUBROUTINE CALLS =  49923, IFAIL =  1
C            N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C             1       0.00E+00       0.11E-13   1.0000000000000000
C             2       0.11E-15       0.11E-13   0.9999999999999999
C             3       0.00E+00       0.12E-13   1.0000000000000000
C             4       0.11E-15       0.11E-13   0.9999999999999999
C             5       0.11E-15       0.11E-13   0.9999999999999999
C             6       0.22E-15       0.12E-06   1.0000000000000002
C             7       0.00E+00       0.30E-07   1.0000000000000000
C             8       0.11E-15       0.24E-06   0.9999999999999999
C             9       0.11E-15       0.45E-06   0.9999999999999999
C            10       0.11E-15       0.32E-06   0.9999999999999999
C            11       0.22E-15       0.14E-05   0.9999999999999998
C            12       0.11E-15       0.13E-06   0.9999999999999999
C            13       0.22E-15       0.14E-05   1.0000000000000002
C            14       0.67E-15       0.12E-05   0.9999999999999993
C            15       0.22E-15       0.18E-05   0.9999999999999998
C
C   End of output
C
      INTEGER LENVER,NW,DEGREE
      PARAMETER (LENVER=1300,NW=200000)
      real(8) VER(3,4,LENVER),EPSREL,WORK(NW),FACTOR
      INTEGER NUMFUN,MAXPTS,RESTAR,IWORK(LENVER+1),I
      COMMON /INFOFC/DEGREE
      EXTERNAL FMONOM,FINVAR
C
C  Test for integrating monomials of degree <= 8
C
      DO 10 DEGREE = 0,8
          WRITE (*,90) DEGREE
   90     FORMAT (' ',/'INTEGRATE MONOMIALS OF DEGREE ',I1)
          VER(1,1,1) = 0.D0
          VER(1,2,1) = 1.D0
          VER(1,3,1) = 0.D0
          VER(1,4,1) = 0.D0
          VER(2,1,1) = 0.D0
          VER(2,2,1) = 0.D0
          VER(2,3,1) = 1.D0
          VER(2,4,1) = 0.D0
          VER(3,1,1) = 0.D0
          VER(3,2,1) = 0.D0
          VER(3,3,1) = 0.D0
          VER(3,4,1) = 1.D0
          EPSREL = 1.D-5
          NUMFUN = FACTOR(DEGREE+2)/ (FACTOR(DEGREE)*2)
          MAXPTS = 50000
          RESTAR = 0
          CALL ATEST(FMONOM,NUMFUN,VER,1,MAXPTS,EPSREL,LENVER,NW,RESTAR,
     +               WORK,IWORK)
   10 CONTINUE
c
c Test restart feature on the special basis of invariant polynomials
c
      WRITE (*,91)
   91 FORMAT (' ',/'INTEGRATE INVARIANT POLYNOMIALS')
      EPSREL = 1
      NUMFUN = 15
      MAXPTS = 50000
      RESTAR = 0
      VER(1,1,1) = 0.D0
      VER(1,2,1) = SQRT(2.d0)/3.d0
      VER(1,3,1) = SQRT(2.d0)/3.d0
      VER(1,4,1) = -2.d0*SQRT(2.d0)/3.d0
      VER(2,1,1) = 0.D0
      VER(2,2,1) = SQRT(6.d0)/3.d0
      VER(2,3,1) = -SQRT(6.d0)/3.d0
      VER(2,4,1) = 0.D0
      VER(3,1,1) = 1.D0
      VER(3,2,1) = -1.d0/3.d0
      VER(3,3,1) = -1.d0/3.d0
      VER(3,4,1) = -1.D0/3.d0
      DO 20 I = 1,6
          EPSREL = EPSREL/10
          CALL ATEST(FINVAR,NUMFUN,VER,1,MAXPTS,EPSREL,LENVER,NW,RESTAR,
     +               WORK,IWORK)
          RESTAR = 1
   20 CONTINUE
      STOP

      END
      SUBROUTINE ATEST(F,NUMFUN,VER,NUMTET,MAXPTS,EPSREL,LENVER,NW,
     +                 RESTAR,WORK,IWORK)
      INTEGER NUMFUN,NUMTET,MAXPTS,LENVER,NW,RESTAR,IWORK(LENVER+1),
     +        NEVAL,IFAIL,N
      real(8) VER(3,4,LENVER),EPSREL,WORK(NW),RESULT(165),
     +                 ABSERR(165)
      EXTERNAL F

      CALL DCUTET(F,NUMFUN,VER,NUMTET,0,MAXPTS,0.d0,EPSREL,LENVER,NW,
     +            RESTAR,RESULT,ABSERR,NEVAL,IFAIL,WORK,IWORK)
      WRITE (*,90) NEVAL,IFAIL
   90 FORMAT (5X,'SUBROUTINE CALLS = ',I6,', IFAIL = ',I2)
      WRITE (*,91)
   91 FORMAT (7X,'N',3X,'ABSOLUTE ERROR',2X,'ESTIMATED ERROR',4X,
     +       'INTEGRAL')
      DO 10 N = 1,NUMFUN
          WRITE (*,92) N,ABS(RESULT(N)-1),ABSERR(N),RESULT(N)
   92     FORMAT (6X,I3,2E15.2,F21.16)
   10 CONTINUE
      END
      SUBROUTINE FMONOM(X,NUMFUN,VALUES)
C
C   Evaluate the monomials of degree DEGREE
C
      INTEGER NUMFUN,COUNT,DEGREE,I,J
      real(8) X(3),VALUES(NUMFUN),MOMENT
      COMMON /INFOFC/DEGREE

      COUNT = 1
      DO 20 J = 0,DEGREE
          DO 10 I = 0,DEGREE - J
              VALUES(COUNT) = X(1)**J*X(2)**I*X(3)** (DEGREE-J-I)
              VALUES(COUNT) = VALUES(COUNT)/MOMENT(J,I,DEGREE-J-I)
              COUNT = COUNT + 1
              IF (COUNT.GT.NUMFUN) RETURN
   10     CONTINUE
   20 CONTINUE
      RETURN
      END
      real(8) FUNCTION FACTOR(L)
C
C  Computes the factorial of L
C
      INTEGER L,K

      FACTOR = 1
      DO 10 K = 2,L
          FACTOR = FACTOR*K
   10 CONTINUE
      RETURN
      END
      real(8) FUNCTION MOMENT(I,J,K)
C
C  Computes the integral of x(1)**I * x(2)**J * x(3)**K
C
      INTEGER I,J,K
      real(8) FACTOR

      MOMENT = FACTOR(I)*FACTOR(J)*FACTOR(K)/FACTOR(3+I+J+K)
      RETURN
      END
      SUBROUTINE FINVAR(X,NUMFUN,VALUES)
C
C  Evaluates the invariant polynomials used to construct
C  a cubature formula.
C
      INTEGER NUMFUN,I
      real(8) X(3),VALUES(NUMFUN),P2,P3,P4,INVAR(15)
c
c  Basic invariant polynomials
c
      P2 = X(1)**2 + X(2)**2 + X(3)**2
      P3 = - (SQRT(2.d0)* (X(1)**2-3.*X(2)**2)*X(1)+3.*X(1)**2*X(3)+
     +     3.*X(2)**2*X(3)-2.*X(3)**3)/2.
      P4 = (4.*SQRT(2.d0)* (X(1)**2-3.*X(2)**2)*X(1)*X(3)+
     +     6.* (2.*X(2)**2+X(3)**2)*X(1)**2+6.*X(1)**4+6.*X(2)**4+
     +     6.*X(2)**2*X(3)**2+7.*X(3)**4)/7.
c
c  Basis used to construct cubature formula
c
      INVAR(1) = 1
      INVAR(2) = 5*P2
      INVAR(3) = 15*P3
      INVAR(4) = 315*P2**2/19
      INVAR(5) = (945.d0/2.d0)* (7.d0/24.d0)* (P2**2-P4)
      INVAR(6) = 35*P2*P3
      INVAR(7) = (1701.d0/41.d0)*P2**3
      INVAR(8) = (3645.d0/2.d0)* (7.d0/24.d0)* (P2**3-P2*P4)
      INVAR(9) = (-8505.d0/8.d0)* (3*P2**3+4*P3**2-7*P2*P4)/12
      INVAR(10) = (18900.d0/67.d0)*P2**2*P3/4
      INVAR(11) = (14175.d0/8.d0)* (7.d0/6.d0)* (P2**2-P4)*P3
      INVAR(12) = (623700.d0/1829.d0)*P2**4/4
      INVAR(13) = (155925.d0/26.d0)* (7.d0/24.d0)* (P2**4-P2**2*P4)
      INVAR(14) = (-467775.d0/152.d0)* (1.d0/12.d0)*
     +            (3*P2**4+4*P2*P3**2-7*P2**2*P4)
      INVAR(15) = (155925.d0/232.d0)* (7.d0/12.d0)*
     +            (7*P4**2-4*P2*P3**2-3*P2**2*P4)
      DO 10 I = 1,NUMFUN
          VALUES(I) = (9*SQRT(3.d0)/8)*INVAR(I)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE DTRTET(DVFLAG,SBRGNS,GREATE,LIST,NEW)
C***BEGIN PROLOGUE DTRTET
C***REFER TO DCUTET
C***PURPOSE  DTRTET maintains a heap of subregions.
C***DESCRIPTION DTRTET maintains a heap of subregions.
C            The subregions are stored in a partially sorted
C            binary tree, ordered according to the size of the
C            greatest error estimates of each subregion(GREATE).
C            The subregion with greatest error estimate is in the
C            first position of the heap.
C
C   PARAMETERS
C
C     DVFLAG Integer.
C            If DVFLAG = 1, we remove the subregion with
C            greatest error from the heap.
C            If DVFLAG = 2, we insert a new subregion in the heap.
C     SBRGNS Integer.
C            Number of subregions in the heap.
C     GREATE Real array of dimension SBRGNS.
C            Used to store the greatest estimated errors in
C            all subregions.
C     LIST   Integer array of dimension SBRGNS.
C            Used as a partially ordered list of pointers to the
C            different subregions. This list is a heap where the
C            element on top of the list is the subregion with the
C            greatest error estimate.
C     NEW    Integer.
C            Index to the new region to be inserted in the heap.
C***ROUTINES CALLED-NONE
C***END PROLOGUE DTRTET
C
C   Global variables.
C
      INTEGER DVFLAG,NEW,SBRGNS,LIST(*)
      real(8) GREATE(*)
C
C   Local variables.
C
C   GREAT  is used as intermediate storage for the greatest error of a
C          subregion.
C   SUBRGN Position of child/parent subregion in the heap.
C   SUBTMP Position of parent/child subregion in the heap.
      INTEGER SUBRGN,SUBTMP
      real(8) GREAT
C
C***FIRST EXECUTABLE STATEMENT DTRTET
C
C    If DVFLAG = 1, we will reduce the partial ordered list by the
C    element with greatest estimated error. Thus the element in
C    in the heap with index LIST(1) is vacant and can be used later.
C    Reducing the heap by one element implies that the last element
C    should be re-positioned.
C
      IF (DVFLAG.EQ.1) THEN
          GREAT = GREATE(LIST(SBRGNS))
          SBRGNS = SBRGNS - 1
          SUBRGN = 1
   20     SUBTMP = 2*SUBRGN
          IF (SUBTMP.LE.SBRGNS) THEN
              IF (SUBTMP.NE.SBRGNS) THEN
C
C   Find max. of left and right child.
C
                  IF (GREATE(LIST(SUBTMP)).LT.
     +                GREATE(LIST(SUBTMP+1))) THEN
                      SUBTMP = SUBTMP + 1
                  END IF
              END IF
C
C   Compare max.child with parent.
C   If parent is max., then done.
C
              IF (GREAT.LT.GREATE(LIST(SUBTMP))) THEN
C
C   Move the pointer at position subtmp up the heap.
C
                  LIST(SUBRGN) = LIST(SUBTMP)
                  SUBRGN = SUBTMP
                  GO TO 20
              END IF
          END IF
C
C   Update the pointer.
C
          IF (SBRGNS.GT.0) THEN
              LIST(SUBRGN) = LIST(SBRGNS+1)
          END IF
      ELSE IF (DVFLAG.EQ.2) THEN
C
C   If DVFLAG = 2, find the position for the NEW region in the heap.
C
          GREAT = GREATE(NEW)
          SUBRGN = SBRGNS
   40     SUBTMP = SUBRGN/2
          IF (SUBTMP.GE.1) THEN
C
C   Compare max.child with parent.
C   If parent is max, then done.
C
              IF (GREAT.GT.GREATE(LIST(SUBTMP))) THEN
C
C   Move the pointer at position subtmp down the heap.
C
                  LIST(SUBRGN) = LIST(SUBTMP)
                  SUBRGN = SUBTMP
                  GO TO 40
              END IF
          END IF
C
C    Set the pointer to the new region in the heap.
C
          LIST(SUBRGN) = NEW
      END IF
C
C***END DTRTET
C
      RETURN
      END
      REAL FUNCTION R1MACH(I)
C
C  SINGLE-PRECISION MACHINE CONSTANTS
C
C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  R1MACH(5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C
C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED integer(2)), YOU MAY NEED TO
C  CHANGE INTEGER TO integer(4) OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC
C
      REAL RMACH(5)
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
       DATA SMALL(1) /     8388608 /
       DATA LARGE(1) /  2139095039 /
       DATA RIGHT(1) /   864026624 /
       DATA DIVER(1) /   872415232 /
       DATA LOG10(1) /  1050288283 /, SC/987/
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1) /    1048576 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  990904320 /
C      DATA DIVER(1) / 1007681536 /
C      DATA LOG10(1) / 1091781651 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA RMACH(1) / Z400800000 /
C      DATA RMACH(2) / Z5FFFFFFFF /
C      DATA RMACH(3) / Z4E9800000 /
C      DATA RMACH(4) / Z4EA800000 /
C      DATA RMACH(5) / Z500E730E8 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.
C
C      DATA RMACH(1) / O1771000000000000 /
C      DATA RMACH(2) / O0777777777777777 /
C      DATA RMACH(3) / O1311000000000000 /
C      DATA RMACH(4) / O1301000000000000 /
C      DATA RMACH(5) / O1157163034761675 /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN4 ON THE CDC 6000/7000 SERIES.
C
C      DATA RMACH(1) / 00564000000000000000B /
C      DATA RMACH(2) / 37767777777777777776B /
C      DATA RMACH(3) / 16414000000000000000B /
C      DATA RMACH(4) / 16424000000000000000B /
C      DATA RMACH(5) / 17164642023241175720B /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN5 ON THE CDC 6000/7000 SERIES.
C
C      DATA RMACH(1) / O"00564000000000000000" /
C      DATA RMACH(2) / O"37767777777777777776" /
C      DATA RMACH(3) / O"16414000000000000000" /
C      DATA RMACH(4) / O"16424000000000000000" /
C      DATA RMACH(5) / O"17164642023241175720" /, SC/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA RMACH(1) / '00800000'X /
C      DATA RMACH(2) / '7FFFFFFF'X /
C      DATA RMACH(3) / '34800000'X /
C      DATA RMACH(4) / '35000000'X /
C      DATA RMACH(5) / '3F9A209B'X /, SC/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA RMACH(1) / 200034000000000000000B /
C      DATA RMACH(2) / 577767777777777777776B /
C      DATA RMACH(3) / 377224000000000000000B /
C      DATA RMACH(4) / 377234000000000000000B /
C      DATA RMACH(5) / 377774642023241175720B /, SC/987/
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC RMACH(5)
C
C      DATA SMALL/20K,0/,LARGE/77777K,177777K/
C      DATA RIGHT/35420K,0/,DIVER/36020K,0/
C      DATA LOG10/40423K,42023K/, SC/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '00000177 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000353 /
C      DATA LOG10(1),LOG10(2) / '23210115, '00000377 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA RMACH(1) / O402400000000 /
C      DATA RMACH(2) / O376777777777 /
C      DATA RMACH(3) / O714400000000 /
C      DATA RMACH(4) / O716400000000 /
C      DATA RMACH(5) / O776464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA RMACH(1) / Z00100000 /
C      DATA RMACH(2) / Z7FFFFFFF /
C      DATA RMACH(3) / Z3B100000 /
C      DATA RMACH(4) / Z3C100000 /
C      DATA RMACH(5) / Z41134413 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA RMACH(1) / Z'00100000' /
C      DATA RMACH(2) / Z'7EFFFFFF' /
C      DATA RMACH(3) / Z'3B100000' /
C      DATA RMACH(4) / Z'3C100000' /
C      DATA RMACH(5) / Z'41134413' /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C
C      DATA RMACH(1) / "000400000000 /
C      DATA RMACH(2) / "377777777777 /
C      DATA RMACH(3) / "146400000000 /
C      DATA RMACH(4) / "147400000000 /
C      DATA RMACH(5) / "177464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1) /    8388608 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  880803840 /
C      DATA DIVER(1) /  889192448 /
C      DATA LOG10(1) / 1067065499 /, SC/987/
C
C      DATA RMACH(1) / O00040000000 /
C      DATA RMACH(2) / O17777777777 /
C      DATA RMACH(3) / O06440000000 /
C      DATA RMACH(4) / O06500000000 /
C      DATA RMACH(5) / O07746420233 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /   128,     0 /
C      DATA LARGE(1),LARGE(2) / 32767,    -1 /
C      DATA RIGHT(1),RIGHT(2) / 13440,     0 /
C      DATA DIVER(1),DIVER(2) / 13568,     0 /
C      DATA LOG10(1),LOG10(2) / 16282,  8347 /, SC/987/
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA RIGHT(1),RIGHT(2) / O032200, O000000 /
C      DATA DIVER(1),DIVER(2) / O032400, O000000 /
C      DATA LOG10(1),LOG10(2) / O037632, O020233 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA SMALL(1) / $00800000 /
C      DATA LARGE(1) / $7F7FFFFF /
C      DATA RIGHT(1) / $33800000 /
C      DATA DIVER(1) / $34000000 /
C      DATA LOG10(1) / $3E9A209B /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA RMACH(1) / O000400000000 /
C      DATA RMACH(2) / O377777777777 /
C      DATA RMACH(3) / O146400000000 /
C      DATA RMACH(4) / O147400000000 /
C      DATA RMACH(5) / O177464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER.
C
C      DATA SMALL(1) /       128 /
C      DATA LARGE(1) /    -32769 /
C      DATA RIGHT(1) /     13440 /
C      DATA DIVER(1) /     13568 /
C      DATA LOG10(1) / 547045274 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER.
C
C      DATA RMACH(1) / Z00000080 /
C      DATA RMACH(2) / ZFFFF7FFF /
C      DATA RMACH(3) / Z00003480 /
C      DATA RMACH(4) / Z00003500 /
C      DATA RMACH(5) / Z209B3F9A /, SC/987/
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2.
C
C      DATA RMACH(1) /       '80'X /
C      DATA RMACH(2) / 'FFFF7FFF'X /
C      DATA RMACH(3) /     '3480'X /
C      DATA RMACH(4) /     '3500'X /
C      DATA RMACH(5) / '209B3F9A'X /, SC/987/
C
C  ***  ISSUE STOP 778 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SC .NE. 987) STOP 778
      IF (I .LT. 1  .OR.  I .GT. 5) GOTO 999
      R1MACH = RMACH(I)
      RETURN
  999 WRITE(*,1999) I
 1999 FORMAT(' R1MACH - I OUT OF BOUNDS',I10)
      STOP
      END

      SUBROUTINE SADTET(NUMFUN,MDIV,VER,NUMTET,MINSUB,MAXSUB,FUNSUB,
     +                  EPSABS,EPSREL,LENVER,RESTAR,LENW,RESULT,ABSERR,
     +                  NEVAL,NSUB,IFAIL,VALUES,ERRORS,GREATE,WORK2,
     +                  WORK3,LIST,VACANT)
C***BEGIN PROLOGUE SADTET
C***REFER TO SCUTET
C***DESCRIPTION Computation of integrals over regions consisting of
C            a collection of tetrahedrons.
C
C            SADTET repeatedly
C            subdivides the tetrahedrons with greatest estimated  errors
C            and estimates the integrals and the errors over the new
C            new sub-tetrahedrons until the error request is met
C            or MAXPTS function evaluations have been used.
C
C            A 43 point integration rule
C            with all evaluation points inside the tetrahedron
C            is applied. The rule has polynomial degree 8.
C
C   ON ENTRY
C
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            MDIV is the number of tetrahedrons that are divided in
C            each subdivision step in SADTET.
C            MDIV is chosen default to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(8*MDIV,NPROC) = 0.
C     VER    Real array of dimension (3,4,LENVER).
C            VER(1,K,L), VER(2,K,L) and VER(3,K,L) are the x, y and z
C            coordinates respectively of vertex K in tetrahedron L.
C            On entry VER(*,*,L) must contain the vertices of the
C            NUMTET user specified tetrahedrons for L = 1,2,...,NUMTET.
C     NUMTET Integer.
C            The number of tetrahedrons.
C     MINSUB Integer.
C            The minimum allowed number of subregions.
C     MAXSUB Integer.
C            The maximum allowed number of subregions.
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (X,NUMFUN,FUNVLS)
C            Input parameters:
C              X(1)   The x-coordinate of the evaluation point.
C              X(2)   The y-coordinate of the evaluation point.
C              X(3)   The z-coordinate of the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C
C     LENVER Integer.
C            Defines the length of the array VER.
C
C            We let
C            MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS.
C            MAXSUB = 7*(MAXPTS-43*NUMTET)/(8*43) + NUMTET
C            LENVER should then be greater or equal to MAXSUB.
C
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1,
C            then we restart a previous attempt.
C            In this case the only parameters for DCUTET that may
C            be changed (with respect to the previous call of DCUTET)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C     LENW   Integer.
C            Length of the workspace WORK2.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute errors.
C     NEVAL  Integer.
C            Number of function evaluations used by DCUTET.
C     NSUB   Integer.
C            The number of tetrahedrons in the data structure.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit.
C
C              ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
C              function evaluations for all values of K,
C              1 <= K <= NUMFUN .
C
C            IFAIL = 1 if MAXSUB was too small for SADTET
C              to obtain the required accuracy. In this case SADTET
C              returns values of RESULT with estimated absolute
C              errors ABSERR.
C     VALUES Real array of dimension (NUMFUN,MAXSUB).
C            The estimated components of the integrals over the
C            subregions.
C     ERRORS Real array of dimension (NUMFUN,MAXSUB).
C            The estimated errors over the subregions.
C     GREATE Real array of dimension MAXSUB.
C            The greatest errors in each subregion.
C     WORK2  Real array of dimension 6*LENW.
C            Work array used in SRLTET.
C     WORK3  Real array of dimension LENW.
C            Work array used in SRLTET.
C     LIST   Integer array used in STRTET of dimension LENVER.
C            Is a partially sorted list, where LIST(1) is the top
C            element in a heap of subregions.
C     VACANT Integer array of dimension MDIV
C            Used as intermediate storage of pointers.
C***ROUTINES CALLED STRTET,SRLTET
C***END PROLOGUE SADTET
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NUMFUN,MDIV,NUMTET,MINSUB,MAXSUB,LENW,RESTAR,LENVER,NEVAL,
     +        NSUB,IFAIL,LIST(LENVER),VACANT(MDIV)
      REAL VER(3,4,LENVER),EPSABS,EPSREL,RESULT(NUMFUN),ABSERR(NUMFUN),
     +     VALUES(NUMFUN,MAXSUB),ERRORS(NUMFUN,MAXSUB),GREATE(MAXSUB),
     +     WORK3(LENW),WORK2(6*LENW)
C
C   Local variables.
C
C   SBRGNS is the number of stored subregions.
C   NDIV   The number of subregions to be divided in each main step.
C   POINTR Pointer to the position in the data structure where
C          the new subregions are to be stored.
C   G      is the homogeneous coordinates of the integration rule.
C   W      is the weights of the integration rule and the null rules.
C   TOP    is a pointer to the top element in the heap of subregions.
C   COUNT  is the number of times the termination criterium is
C          satisfied.
      INTEGER I,J,K,L,TOP,SIZE
      INTEGER SBRGNS,I1,I2,I3,I4,I5,I6,I7,I8
      INTEGER L1,K1
      INTEGER COUNT,GROUP,LEFT,NDIV,POINTR,INDEX,NUM
      REAL VEROLD(3,4),VERNEW(3,6),SUMVAL(2),SUMERR(2)
C
C***FIRST EXECUTABLE STATEMENT SADTET
C
      COUNT = 0
      NUM = 43
      IF (RESTAR.EQ.1) THEN
          SBRGNS = NSUB
          GO TO 110
      END IF
C
C   Initiate RESULT and ABSERR.
C
      DO 10 J = 1,NUMFUN
          RESULT(J) = 0
          ABSERR(J) = 0
   10 CONTINUE
C
C   Apply SRLTET over the NUMTET tetrahedrons.
C   This loop may be run in parallel.
C
      DO 20 I = 1,NUMTET
          L1 = 1 + (I-1)*6*NUMFUN
          K1 = 1 + (I-1)*NUMFUN
          CALL SRLTET(VER(1,1,I),NUMFUN,FUNSUB,WORK2(L1),VALUES(1,I),
     +                ERRORS(1,I),GREATE(I),WORK3(K1))
   20 CONTINUE
      NEVAL = NUM*NUMTET
      SBRGNS = NUMTET
C
C   Add the computed values to RESULT and ABSERR.
C
      DO 40 I = 1,NUMTET
          DO 30 J = 1,NUMFUN
              RESULT(J) = RESULT(J) + VALUES(J,I)
              ABSERR(J) = ABSERR(J) + ERRORS(J,I)
   30     CONTINUE
   40 CONTINUE
C
C   Store results in heap.
C
      DO 50 I = 1,NUMTET
          INDEX = I
          CALL STRTET(2,INDEX,GREATE,LIST, (INDEX))
   50 CONTINUE

C
C   We check for termination.
C
      IF (SBRGNS.LT.MINSUB) THEN
          GO TO 110
      END IF
      DO 60 J = 1,NUMFUN
          IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)) .AND.
     +        ABSERR(J).GT.EPSABS) THEN
              GO TO 110
          END IF
   60 CONTINUE
      IFAIL = 0
      COUNT = COUNT + 1
      GO TO 499
C
C***End initiation.
C
C***Begin loop while the error is too great,
C   and SBRGNS+7 is less than MAXSUB.
C
  110 IF (SBRGNS+7.LE.MAXSUB) THEN
C
C   If we are allowed to divide further,
C   prepare to apply basic rule over  the tetrahedrons produced
C   by dividing the
C   NDIV sub-tetrahedrons with greatest errors.
C   If MAXSUB and SBRGNS are great enough, NDIV = MDIV.
C
          NDIV = MAXSUB - SBRGNS
          NDIV = MIN(NDIV,MDIV,SBRGNS)
C
C   Divide each of the NDIV sub-tetrahedrons in eight new
C   sub-tetrahedrons, and compute integral and error over each.
C   When we pick a tetrahedron to divide in eight one of the
C   new sub-tetrahedron is stored in the original tetrahedrons'
C   position in the datastructure. Thus we get 7*NDIV new
C   elements in the datastructure after the subdivision. The size
C   of the datastructure before the sub-division
C   is stored in the variable SIZE, while SBRGNS is the size of
C   heap at any time.
          SIZE = SBRGNS
          DO 150 I = 1,NDIV
              POINTR = SIZE + 7* (NDIV+1-I)
C
C   Adjust RESULT and ABSERR. TOP is a pointer to the top of the heap.
C
              TOP = LIST(1)
              VACANT(I) = TOP
              DO 115 J = 1,NUMFUN
                  RESULT(J) = RESULT(J) - VALUES(J,TOP)
                  ABSERR(J) = ABSERR(J) - ERRORS(J,TOP)
  115         CONTINUE
C
C   Save the vertices.
C
              DO 130 L = 1,3
                  DO 135 K = 1,4
                      VEROLD(L,K) = VER(L,K,TOP)
  135             CONTINUE
  130         CONTINUE
C
C   Adjust the heap.
C
              CALL STRTET(1,SBRGNS,GREATE,LIST,K)
C
C   Compute the eight new tetrahedrons. Store 7 in consecutive
C   positions and 1 in the vacant position.
C
              I1 = TOP
              I2 = POINTR - 6
              I3 = POINTR - 5
              I4 = POINTR - 4
              I5 = POINTR - 3
              I6 = POINTR - 2
              I7 = POINTR - 1
              I8 = POINTR
C
C   Compute the 6 division points.
C
              DO 140 L = 1,3
                  VERNEW(L,1) = (VEROLD(L,1)+VEROLD(L,2))/2
                  VERNEW(L,2) = (VEROLD(L,2)+VEROLD(L,3))/2
                  VERNEW(L,3) = (VEROLD(L,1)+VEROLD(L,3))/2
                  VERNEW(L,4) = (VEROLD(L,1)+VEROLD(L,4))/2
                  VERNEW(L,5) = (VEROLD(L,2)+VEROLD(L,4))/2
                  VERNEW(L,6) = (VEROLD(L,3)+VEROLD(L,4))/2
  140         CONTINUE
C
              DO 145 L = 1,3
                  VER(L,1,I1) = VEROLD(L,1)
                  VER(L,2,I1) = VERNEW(L,1)
                  VER(L,3,I1) = VERNEW(L,3)
                  VER(L,4,I1) = VERNEW(L,4)
                  VER(L,1,I2) = VEROLD(L,2)
                  VER(L,2,I2) = VERNEW(L,1)
                  VER(L,3,I2) = VERNEW(L,2)
                  VER(L,4,I2) = VERNEW(L,5)
                  VER(L,1,I3) = VEROLD(L,3)
                  VER(L,2,I3) = VERNEW(L,2)
                  VER(L,3,I3) = VERNEW(L,3)
                  VER(L,4,I3) = VERNEW(L,6)
                  VER(L,1,I4) = VEROLD(L,4)
                  VER(L,2,I4) = VERNEW(L,4)
                  VER(L,3,I4) = VERNEW(L,5)
                  VER(L,4,I4) = VERNEW(L,6)
                  VER(L,1,I5) = VERNEW(L,1)
                  VER(L,2,I5) = VERNEW(L,3)
                  VER(L,3,I5) = VERNEW(L,4)
                  VER(L,4,I5) = VERNEW(L,5)
                  VER(L,1,I6) = VERNEW(L,1)
                  VER(L,2,I6) = VERNEW(L,2)
                  VER(L,3,I6) = VERNEW(L,3)
                  VER(L,4,I6) = VERNEW(L,5)
                  VER(L,1,I7) = VERNEW(L,2)
                  VER(L,2,I7) = VERNEW(L,3)
                  VER(L,3,I7) = VERNEW(L,5)
                  VER(L,4,I7) = VERNEW(L,6)
                  VER(L,1,I8) = VERNEW(L,3)
                  VER(L,2,I8) = VERNEW(L,4)
                  VER(L,3,I8) = VERNEW(L,5)
                  VER(L,4,I8) = VERNEW(L,6)
  145         CONTINUE
  150     CONTINUE
C
C   Apply basic rule.
C   This loop may be run in parallel.
C
          DO 200 I = 1,8*NDIV
              IF (I.LE.NDIV) THEN
                  INDEX = VACANT(I)
              ELSE
                  INDEX = SBRGNS + I
              END IF
              L1 = 1 + (I-1)*6*NUMFUN
              K1 = 1 + (I-1)*NUMFUN
              CALL SRLTET(VER(1,1,INDEX),NUMFUN,FUNSUB,WORK2(L1),
     +                    VALUES(1,INDEX),ERRORS(1,INDEX),GREATE(INDEX),
     +                    WORK3(K1))
  200     CONTINUE
          NEVAL = NEVAL + 8*NDIV*NUM
C
C   Add new contributions to RESULT and ABSERR.
C
          DO 220 I = 1,8*NDIV
              IF (I.LE.NDIV) THEN
                  INDEX = VACANT(I)
              ELSE
                  INDEX = SBRGNS + I
              END IF
              DO 210 J = 1,NUMFUN
                  RESULT(J) = RESULT(J) + VALUES(J,INDEX)
                  ABSERR(J) = ABSERR(J) + ERRORS(J,INDEX)
  210         CONTINUE
  220     CONTINUE
C
C   Store results in heap.
C
          DO 250 I = 1,8*NDIV
              IF (I.LE.NDIV) THEN
                  INDEX = VACANT(I)
              ELSE
                  INDEX = SBRGNS + I
              END IF
              J = SBRGNS + I
              CALL STRTET(2,J,GREATE,LIST,INDEX)
  250     CONTINUE
          SBRGNS = SBRGNS + 8*NDIV
C
C   Check for termination.
C
          IF (SBRGNS.LT.MINSUB) THEN
              GO TO 110
          END IF
          DO 255 J = 1,NUMFUN
              IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)) .AND.
     +            ABSERR(J).GT.EPSABS) THEN
                  GO TO 110
              END IF
  255     CONTINUE
          IFAIL = 0
          COUNT = COUNT + 1
C
C   Else we did not succeed with the
C   given value of MAXSUB.
C
      ELSE
          IFAIL = 1
      END IF
C
C   Compute more accurate values of RESULT and ABSERR. This is done
C   to reduce the effect of roundoff on final results. Large
C   intermediate sums in the computation may course large, unnecessary
C   roundoff errors. Thus recomputing the sums of errors and estimates
C   and in addition grouping the sums in three groups should remove this
C   problem.
C
  499 CONTINUE
      DO 550 J = 1,NUMFUN
          RESULT(J) = 0
          ABSERR(J) = 0
          GROUP = SBRGNS** (1./3.)
          IF (GROUP**3.NE.SBRGNS) THEN
              GROUP = GROUP + 1
          END IF
          LEFT = SBRGNS
          DO 520 I2 = 0,SBRGNS - 1,GROUP**2
              SUMVAL(2) = 0
              SUMERR(2) = 0
              DO 510 I1 = 0,MIN(LEFT,GROUP**2) - 1,GROUP
                  SUMVAL(1) = 0
                  SUMERR(1) = 0
                  DO 500 I = 1 + I1 + I2,MIN(LEFT,GROUP) + I1 + I2
                      SUMVAL(1) = SUMVAL(1) + VALUES(J,I)
                      SUMERR(1) = SUMERR(1) + ERRORS(J,I)
  500             CONTINUE
                  LEFT = LEFT - MIN(LEFT,GROUP)
                  SUMVAL(2) = SUMVAL(2) + SUMVAL(1)
                  SUMERR(2) = SUMERR(2) + SUMERR(1)
  510         CONTINUE
              RESULT(J) = RESULT(J) + SUMVAL(2)
              ABSERR(J) = ABSERR(J) + SUMERR(2)
  520     CONTINUE
  550 CONTINUE
C
C   Check again for termination
C
      IF (IFAIL.EQ.0) THEN
          IF (COUNT.LE.1) THEN
              DO 560 J = 1,NUMFUN
                  IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)) .AND.
     +                ABSERR(J).GT.EPSABS) THEN
                      GO TO 110
                  END IF
  560         CONTINUE
          END IF
      END IF
C
      NSUB = SBRGNS
      RETURN
C
C***END SADTET
C
      END
      SUBROUTINE SCHTET(NUMFUN,MDIV,VER,NUMTET,MINPTS,MAXPTS,EPSABS,
     +                  EPSREL,LENVER,NW,RESTAR,MAXSUB,MINSUB,IFAIL)
C***BEGIN PROLOGUE SCHTET
C***REFER TO SCUTET
C***PURPOSE  SCHTET checks the validity of the
C            input parameters to SCUTET.
C***DESCRIPTION
C            SCHTET computes MAXSUB, MINSUB and IFAIL as
C            functions of the input parameters to SCUTET,
C            and checks the validity of the input parameters to SCUTET.
C
C   ON ENTRY
C
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            MDIV is the number of tetrahedrons that are divided in
C            each subdivision step in SADTET.
C            MDIV is chosen default to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(8*MDIV,NPROC) = 0.
C     VER    Real array of dimension (3,4,LENVER).
C            VER(1,K,L), VER(2,K,L) and VER(3,K,L) are the x, y and z
C            coordinates respectively of vertex K in tetrahedron L.
C            On entry VER(*,*,L) must contain the vertices of the
C            NUMTET user specified tetrahedrons for L = 1,2,...,NUMTET.
C     NUMTET Integer.
C            The number of tetrahedrons.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            MAXPTS >= NUM*NUMTET and MAXPTS >= MINPTS
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     LENVER Integer.
C            Defines the length of the array VER.
C
C            We let
C            MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS.
C            MAXSUB = 7*(MAXPTS-NUM*NUMTET)/(8*NUM) + NUMTET
C            LENVER should then be greater or equal to MAXSUB.
C
C     NW     Integer.
C            Defines the length of the working array WORK.
C
C            We let
C            MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS.
C            MAXSUB = 7*(MAXPTS-NUM*NUMTET)/(8*NUM) + NUMTET
C            NW should then be greater or equal to
C            MAXSUB*(2*NUMFUN+1) + 6*MAX(8*MDIV,NUMTET)*NUMFUN + 1
C
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1,
C            then we restart a previous attempt.
C            In this case the only parameters for SCUTET that may
C            be changed (with respect to the previous call of SCUTET)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C
C
C   ON RETURN
C
C     MAXSUB Integer.
C            The maximum allowed number of subregions for the
C            given values of MAXPTS.
C     MINSUB Integer.
C            The minimum allowed number of subregions for the given
C            values of MINPTS.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit.
C            IFAIL = 2 if NUMFUN is less than 1.
C            IFAIL = 3 if the volume of one of the initially given
C                      tetrahedrons is zero.
C            IFAIL = 4 if MAXPTS is less than NUM*NUMTET.
C            IFAIL = 5 if MAXPTS is less than MINPTS.
C            IFAIL = 6 if EPSABS <= 0 and EPSREL <= 0.
C            IFAIL = 7 if LENVER is less than MAXSUB.
C            IFAIL = 8 if NW is less than MAXSUB*(2*NUMFUN+1) +
C                      NUMFUN*6*MAX(8*MDIV,NUMTET) + 1.
C            IFAIL = 9 if illegal value of RESTAR.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE SCHTET
C
C   Global variables.
C
      INTEGER NUMFUN,MDIV,NUMTET,MINPTS,MAXPTS,NW,MAXSUB,MINSUB,RESTAR,
     +        LENVER,IFAIL
      REAL VER(3,4,NUMTET),EPSABS,EPSREL
C
C   Local variables.
C   NUM is set equal to the number of evaluation points used
C   by the cubature rule over each tetrahedron.
C
      INTEGER LIMIT,J,NUM
      REAL VOLUME
C
C***FIRST EXECUTABLE STATEMENT SCHTET
C
      IFAIL = 0
      NUM = 43
C
C   Compute MAXSUB and MINSUB.
C
      MAXSUB = 7* (MAXPTS-NUM*NUMTET)/ (8*NUM) + NUMTET
      MINSUB = 7* (MINPTS-NUM*NUMTET)/ (8*NUM) + NUMTET
      IF (MOD(7* (MINPTS-NUM*NUMTET),8*NUM).GT.0) THEN
          MINSUB = MINSUB + 1
      END IF
      MINSUB = MAX(NUMTET,MINSUB)
C
C   Check on positive NUMFUN.
C
      IF (NUMFUN.LT.1) THEN
          IFAIL = 2
          RETURN
      END IF
C
C   Check on legal volume of each subregion of integration.
C
      DO 10 J = 1,NUMTET
          VOLUME = (VER(1,2,J)-VER(1,1,J))*
     +             ((VER(2,3,J)-VER(2,1,J))* (VER(3,4,J)-VER(3,1,J))-
     +             (VER(2,4,J)-VER(2,1,J))* (VER(3,3,J)-VER(3,1,J))) -
     +             (VER(2,2,J)-VER(2,1,J))* ((VER(1,3,J)-VER(1,1,J))*
     +             (VER(3,4,J)-VER(3,1,J))- (VER(1,4,J)-VER(1,1,J))*
     +             (VER(3,3,J)-VER(3,1,J))) +
     +             (VER(3,2,J)-VER(3,1,J))* ((VER(1,3,J)-VER(1,1,J))*
     +             (VER(2,4,J)-VER(2,1,J))- (VER(1,4,J)-VER(1,1,J))*
     +             (VER(2,3,J)-VER(2,1,J)))
          VOLUME = VOLUME/6
          IF (VOLUME.EQ.0) THEN
              IFAIL = 3
              RETURN
          END IF
   10 CONTINUE
C
C   Check on MAXPTS >= NUM*NUMTET.
C
      IF (MAXPTS.LT.NUM*NUMTET) THEN
          IFAIL = 4
          RETURN
      END IF
C
C   Check on MAXPTS >= MINPTS.
C
      IF (MAXPTS.LT.MINPTS) THEN
          IFAIL = 5
          RETURN
      END IF
C
C   Check on legal accuracy requests.
C
      IF (EPSABS.LE.0 .AND. EPSREL.LE.0) THEN
          IFAIL = 6
          RETURN
      END IF
C
C   Check on big enough LENVER.
C
      IF (LENVER.LT.MAXSUB) THEN
          IFAIL = 7
          RETURN
      END IF
C
C   Check on big enough NW.
C
      LIMIT = MAXSUB* (2*NUMFUN+1) + 7*MAX(8*MDIV,NUMTET)*NUMFUN + 1
      IF (NW.LT.LIMIT) THEN
          IFAIL = 8
          RETURN
      END IF
C
C    Check on legal RESTAR.
C
      IF (RESTAR.NE.0 .AND. RESTAR.NE.1) THEN
          IFAIL = 9
      END IF
      RETURN
C
C***END SCHTET
C
      END
      SUBROUTINE SCUTET(FUNSUB,NUMFUN,VER,NUMTET,MINPTS,MAXPTS,EPSABS,
     +                  EPSREL,LENVER,NW,RESTAR,RESULT,ABSERR,NEVAL,
     +                  IFAIL,WORK,IWORK)
C***BEGIN PROLOGUE SCUTET
C***DATE WRITTEN   900220   (YYMMDD)
C***REVISION DATE  900307   (YYMMDD)
C***CATEGORY NO. H2B2A1/H2B2A2
C***KEYWORDS QUADRATURE, THREE DIMENSIONAL, ADAPTIVE, CUBATURE
C***AUTHORS
C            Jarle Berntsen, The Computing Centre,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544055
C            Email..  jarle@eik.ii.uib.no
C
C            Ronald Cools, Dept. of Computer Science,
C            Katholieke Universiteit Leuven, Celestijnenlaan 200A,
C            B-3030 Heverlee, Belgium
C            Phone..  32-16-201015 (3562)
C            Email..  ronald@cs.kuleuven.ac.be
C
C            Terje O. Espelid, Department of Informatics,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544180
C            Email..  terje@eik.ii.uib.no
C
C***PURPOSE  To compute the three-dimensional integral over a region
C            consisting of NUMTET tetrahedrons.
C***DESCRIPTION
C
C            The routine calculates an approximation to a given
C            vector of definite integrals
C
C            I  I  I (F ,F ,...,F )      DX(3)DX(2)DX(1),
C                      1  2      NUMFUN
C
C            where the region of integration is a collection of
C            NUMTET tetrahedrons and
C            where F = F (X(1),X(2),X(3)), J = 1,2,...,NUMFUN.
C                   J   J
C
C
C            A globally adaptive strategy is applied in order to
C            compute approximations RESULT(K)
C            hopefully satisfying, for each component of I, the
C            following claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C
C            SCUTET is a driver for the integration routine
C            SADTET.
C
C            SADTET repeatedly
C            subdivides the tetrahedrons with greatest estimated  errors
C            and estimates the integrals and the
C            errors over the new subtetrahedrons
C            until the error request
C            is met or MAXPTS function evaluations have been used.
C
C            A 43 point integration rule
C            with all evaluation points inside the tetrahedron
C            is applied. The rule has polynomial degree 8.
C
C            If the values of the input parameters EPSABS
C            or EPSREL are selected great enough,
C            an integration rule is applied over each tetrahedron and
C            the results are added up to give the approximations
C            RESULT(K). No further subdivision
C            of the tetrahedrons will then be applied.
C
C            When SCUTET computes estimates to a vector of
C            integrals, all components of the vector are given
C            the same treatment. That is, I(F ) and I(F ) for
C                                            J         K
C            J not equal to K, are estimated with the same
C            subdivision of the region of integration.
C            For integrals with enough similarity, we may save
C            time by applying SCUTET to all integrands in one call.
C            For integrals that varies continuously as functions of
C            some parameter, the estimates produced by SCUTET will
C            also vary continuously when the same subdivision is
C            applied to all components. This will generally not be
C            the case when the different components are given
C            separate treatment.
C
C            On the other hand this feature should be used with
C            caution when the different components of the integrals
C            require clearly different subdivisions.
C
C   ON ENTRY
C
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (X,NUMFUN,FUNVLS)
C            Input parameters:
C              X(1)   The x-coordinate of the evaluation point.
C              X(2)   The y-coordinate of the evaluation point.
C              X(3)   The z-coordinate of the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     NUMFUN Integer.
C            Number of components of the integral.
C     VER    Real array of dimension (3,4,LENVER).
C            VER(1,K,L), VER(2,K,L) and VER(3,K,L) are the x, y and z
C            coordinates respectively of vertex K in tetrahedron L.
C            On entry VER(*,*,L) must contain the vertices of the
C            NUMTET user specified tetrahedrons for L = 1,2,...,NUMTET.
C     NUMTET Integer.
C            The number of tetrahedrons.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is 43.
C
C            MAXPTS >= 43*NUMTET and MAXPTS >= MINPTS
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     LENVER Integer.
C            Defines the length of the array VER.
C
C            We let
C            MAXSUB denote the maximum allowed number of subregions
C            for the given value of MAXPTS.
C            MAXSUB = 7*(MAXPTS-43*NUMTET)/(8*43) + NUMTET
C            LENVER should be greater or equal to MAXSUB.
C
C     NW     Integer.
C            Defines the length of the working array WORK.
C
C            We let
C            MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS.
C            MAXSUB = 7*(MAXPTS-43*NUMTET)/(8*43) + NUMTET
C            NW should then be greater or equal to
C            MAXSUB*(2*NUMFUN+1) + MAX(8*6*MDIV,6*NUMTET)*NUMFUN + 1
C            where MDIV is the number of tetrahedrons that are divided
C            in each subdivision step.
C            MDIV is default set internally in SCUTET equal to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(8*MDIV,NPROC) = 0.
C
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1,
C            then we restart a previous attempt.
C            In this case the only parameters for SCUTET that may
C            be changed (with respect to the previous call of SCUTET)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute errors.
C     NEVAL  Integer.
C            Number of function evaluations used by SCUTET.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit.
C
C              ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
C              function evaluations for all values of K,
C              1 <= K <= NUMFUN .
C
C            IFAIL = 1 if MAXPTS was too small for SCUTET
C              to obtain the required accuracy. In this case SCUTET
C              returns values of RESULT with estimated absolute
C              errors ABSERR.
C            IFAIL = 2 if NUMFUN is less than 1.
C            IFAIL = 3 if the volume of one of the initially given
C                      tetrahedrons is zero.
C            IFAIL = 4 if MAXPTS is less than 43*NUMTET.
C            IFAIL = 5 if MAXPTS is less than MINPTS.
C            IFAIL = 6 if EPSABS <= 0 and EPSREL <= 0.
C            IFAIL = 7 if LENVER is less than MAXSUB.
C            IFAIL = 8 if NW is less than MAXSUB*(2*NUMFUN+1) +
C                      NUMFUN*MAX(8*6*MDIV,6*NUMTET) + 1.
C            IFAIL = 9 if unlegal RESTAR.
C     VER    Real array of dimension (3,4,LENVER).
C            On exit VER
C            contains the vertices of the tetrahedrons used to produce
C            the approximations to the integrals.
C     WORK   Real array of dimension NW.
C            Used as working storage.
C            WORK(NW) = NSUB, the number of subregions in the data
C            structure.
C            Let WRKSUB=(NW-1-NUMFUN*7*MAX(8*MDIV,NUMTET))/
C                        (2*NUMFUN+1).
C            WORK(1),...,WORK(NUMFUN*WRKSUB) contain
C              the estimated components of the integrals over the
C              subregions.
C            WORK(NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB) contain
C              the estimated errors over the subregions.
C            WORK(2*NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*
C              WRKSUB+WRKSUB) contain the greatest errors
C              in each subregion.
C            WORK((2*NUMFUN+1)*WRKSUB),...,WORK(NW - 1) is used as
C              temporary storage in SADTET.
C     IWORK  Integer array of dimension LENVER + MDIV.
C            Used as working storage.
C
C           SAMPLE PROGRAM
C   The following program will approximate the integral of
C                     exp(x*x+y*y+z*z)
C   over the tetrahedron (0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.)
C   and print out the values of the estimated integral, the estimated
C   error, the number of function evaluations and IFAIL.
C     PROGRAM STEST1
C     real(8) VER(3,4,50),EPSABS,EPSREL,RESULT(1),ABSERR(1),
C    +                 WORK(1000)
C     INTEGER NUMFUN,NUMTET,MINPTS,MAXPTS,LENVER,NW,RESTAR,NEVAL,
C    +        IFAIL,IWORK(51)
C     EXTERNAL F
C     VER(1,1,1) = 0.D0
C     VER(1,2,1) = 1.D0
C     VER(1,3,1) = 0.D0
C     VER(1,4,1) = 0.D0
C     VER(2,1,1) = 0.D0
C     VER(2,2,1) = 0.D0
C     VER(2,3,1) = 1.D0
C     VER(2,4,1) = 0.D0
C     VER(3,1,1) = 0.D0
C     VER(3,2,1) = 0.D0
C     VER(3,3,1) = 0.D0
C     VER(3,4,1) = 1.D0
C     EPSABS = 0.D0
C     EPSREL = 1.D-5
C     NUMFUN = 1
C     NUMTET = 1
C     MINPTS = 0
C     MAXPTS = 1000
C     LENVER = 50
C     NW = 1000
C     RESTAR = 0
C     CALL SCUTET(F,NUMFUN,VER,NUMTET,MINPTS,MAXPTS,EPSABS,
C    +            EPSREL,LENVER,NW,RESTAR,RESULT,ABSERR,NEVAL,
C    +            IFAIL,WORK,IWORK)
C     WRITE(*,*)'RESULT=',RESULT(1)
C     WRITE(*,*)'ERROR =',ABSERR(1)
C     WRITE(*,*)'NEVAL =',NEVAL
C     WRITE(*,*)'IFAIL =',IFAIL
C     END
C     SUBROUTINE F(X,NUMFUN,FUNVLS)
C     INTEGER NUMFUN
C     real(8) X(3),FUNVLS(NUMFUN)
C     FUNVLS(1) = EXP(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
C     RETURN
C     END
C
C   Output produced on a SUN 3/50
C
C       RESULT=   0.227800
C       ERROR =    1.97534E-06
C       NEVAL =  43
C       IFAIL =  0
C
C***LONG DESCRIPTION
C
C   The information for each tetrahedron is contained in the
C   data structures VER,  WORK and IWORK.
C   VER contains the coordinates of the tetrahedrons.
C   When passed on to SADTET, WORK is split into four
C   arrays VALUES, ERRORS, GREATE and WORK2.
C   VALUES contains the estimated values of the integrals.
C   ERRORS contains the estimated errors.
C   GREATE contains the greatest estimated error for each tetrahedron.
C   The data structures for the tetrahedrons are in SADTET organized
C   as a heap and the size of GREATE(I) defines the position of
C   tetrahedron I in the heap. The heap is maintained by the program
C   STRTET and we use a partially ordered list of pointers, saved in
C   IWORK. When passed on to SADTET, IWORK is split into two
C   arrays LIST and VACANT. LIST is a partially ordered list
C   such that GREATE(LIST(1)) is the maximum error estimate
C   for all subtetrahedrons in our datastructure.
C   VACANT is a work space needed in th updating process of the list.
C
C   The subroutine SADTET is written for efficient execution on shared
C   memory parallel computers. On a computer with NPROC processors we
C   will in each subdivision step divide MDIV tetrahedrons, where MDIV
C   is chosen such that MOD(8*MDIV,NPROC) = 0, in totally 8*MDIV new
C   tetrahedrons.  Each processor will then compute estimates of the
C   integrals and errors over 8*MDIV/NPROC tetrahedrons in each
C   subdivision step.
C   The subroutine for estimating the integral and the error over
C   each subregion, SRLTET, uses WORK2 as a work array.
C   We must make sure that each processor writes its results to
C   separate parts of the memory, and therefore the sizes of WORK and
C   WORK2 are functions of MDIV.
C   In order to achieve parallel processing of tetrahedrons, compiler
C   directives should be placed in front of the DO 20 and the DO 200
C   loops in SADTET on machines like Alliant and CRAY.
C
C***REFERENCES
C   J.Berntsen, R. Cools and T.O.Espelid, An Algorithm for Adaptive
C   Cubature Over a Collection of Tetrahedrons, to be published.
C
C***ROUTINES CALLED SCHTET,SADTET
C***END PROLOGUE SCUTET
C
C   Parameters
C
C   MDIV   Integer.
C          MDIV is the number of tetrahedrons that are divided in
C          each subdivision step in SADTET.
C          MDIV is chosen default to 1.
C          For efficient execution on parallel computers
C          with NPROC processors MDIV should be set equal to
C          the smallest integer such that MOD(8*MDIV,NPROC) = 0.
C
      INTEGER MDIV
      PARAMETER (MDIV=1)
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NUMFUN,NUMTET,MINPTS,MAXPTS,LENVER,NW,RESTAR,NEVAL,IFAIL,
     +        IWORK(LENVER+MDIV)
      REAL VER(3,4,LENVER),EPSABS,EPSREL,RESULT(NUMFUN),ABSERR(NUMFUN),
     +     WORK(NW)
C
C   Local variables.
C
C   MAXSUB Integer.
C          The maximum allowed number of subdivisions
C          for the given values of MAXPTS.
C   MINSUB Integer.
C          The minimum allowed number of subregions for the given
C          values of MINPTS.
C   WRKSUB Integer.
C          The maximum allowed number of subregions as a function
C          of NW, NUMFUN and MDIV. This determines the length
C          of the main work arrays.
C
      INTEGER MAXSUB,MINSUB,NSUB,LENW
      INTEGER WRKSUB,I1,I2,I3,I4,I5
C
C***FIRST EXECUTABLE STATEMENT SCUTET
C
C   Call SCHTET to compute MAXSUB and MINSUB,
C   and to check the input parameters.
C
      CALL SCHTET(NUMFUN,MDIV,VER,NUMTET,MINPTS,MAXPTS,EPSABS,EPSREL,
     +            LENVER,NW,RESTAR,MAXSUB,MINSUB,IFAIL)
C
C   All tests on input data passed?
C
      IF (IFAIL.NE.0) THEN
          RETURN
      END IF
C
C   Split up the work space.
C
      WRKSUB = (NW-1-NUMFUN*7*MAX(8*MDIV,NUMTET))/ (2*NUMFUN+1)
      I1 = 1
      I2 = I1 + WRKSUB*NUMFUN
      I3 = I2 + WRKSUB*NUMFUN
      I4 = I3 + WRKSUB
      I5 = I4 + 6*MAX(8*MDIV,NUMTET)*NUMFUN
C
C   On restart runs the number of subregions from the
C   previous call is assigned to NSUB.
C
      IF (RESTAR.EQ.1) THEN
          NSUB = WORK(NW)
      END IF
C
C   Compute the size of the temporary work space needed in SADTET.
C
      LENW = MAX(8*MDIV,NUMTET)*NUMFUN
      CALL SADTET(NUMFUN,MDIV,VER,NUMTET,MINSUB,MAXSUB,FUNSUB,EPSABS,
     +            EPSREL,LENVER,RESTAR,LENW,RESULT,ABSERR,NEVAL,NSUB,
     +            IFAIL,WORK(I1),WORK(I2),WORK(I3),WORK(I4),WORK(I5),
     +            IWORK(1),IWORK(1+LENVER))
      WORK(NW) = NSUB
      RETURN
C
C***END SCUTET
C
      END
      SUBROUTINE SORTET(TYPE,GENER,VER,NUMFUN,FUNSUB,SUMVAL,WORK)
C***BEGIN PROLOGUE SORTET
C***PURPOSE  To compute the sum of function values over all points
C            of an orbit.
C***DESCRIPTION
C   ON ENTRY
C
C   TYPE   Integer
C          The type of the orbit.
C   GENER  Integer array of dimension (3).
C          The generator for the orbit in homogeneous coordinates.
C   VER    Real array of dimension (3,4).
C          The coordinates of the vertices of the tetrahedron.
C          vertex i -> ( ver(1,i),ver(2,i),ver(3,i) )
C   NUMFUN Integer.
C          Number of components of the vector integrand.
C   FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (X,NUMFUN,FUNVLS)
C            Input parameters:
C              X(1)   The x-coordinate of the evaluation point.
C              X(2)   The y-coordinate of the evaluation point.
C              X(3)   The z-coordinate of the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of the vector integrand.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C   WORK   Real array of dimension (NUMFUN).
C          A work array
C
C   ON RETURN
C
C   SUMVAL Real array of dimension (NUMFUN).
C          The sum of function values over all points
C          of the given orbit.
C
C***END PROLOGUE  SORTET
C
C  Global variables
C
      INTEGER NUMFUN,TYPE
      REAL GENER(3),SUMVAL(NUMFUN),VER(3,4),WORK(NUMFUN)
      EXTERNAL FUNSUB
C
C  Local variables
C
      INTEGER I,J,NUMBER
      REAL X(3,12),Z1,Z2,Z3
C
C  The array 'gener' contains the necessary information
C   for the generator
C
C***FIRST EXECUTABLE STATEMENT  SORTET
      GO TO (100,110,120,130) TYPE + 1
C
C  Generator with homogeneous coordinates (1/4,1/4,1/4;1/4)
C
  100 NUMBER = 1
      X(1,1) = (VER(1,1)+VER(1,2)+VER(1,3)+VER(1,4))/4
      X(2,1) = (VER(2,1)+VER(2,2)+VER(2,3)+VER(2,4))/4
      X(3,1) = (VER(3,1)+VER(3,2)+VER(3,3)+VER(3,4))/4
      GO TO 200
C
C  Generator with homogeneous coordinates (z1,z2,z2;z2)
C
  110 NUMBER = 4
      Z1 = GENER(1)
      Z2 = GENER(2)
      DO 115 J = 1,3
          X(J,1) = Z1*VER(J,1) + Z2* (VER(J,2)+VER(J,3)+VER(J,4))
          X(J,2) = Z1*VER(J,2) + Z2* (VER(J,1)+VER(J,3)+VER(J,4))
          X(J,3) = Z1*VER(J,3) + Z2* (VER(J,2)+VER(J,1)+VER(J,4))
          X(J,4) = Z1*VER(J,4) + Z2* (VER(J,2)+VER(J,3)+VER(J,1))
  115 CONTINUE
      GO TO 200
C
C  Generator with homogeneous coordinates (z1,z1,z2;z2)
C
  120 NUMBER = 6
      Z1 = GENER(1)
      Z2 = GENER(2)
      DO 125 J = 1,3
          X(J,1) = Z1* (VER(J,1)+VER(J,2)) + Z2* (VER(J,3)+VER(J,4))
          X(J,2) = Z1* (VER(J,1)+VER(J,3)) + Z2* (VER(J,2)+VER(J,4))
          X(J,3) = Z1* (VER(J,1)+VER(J,4)) + Z2* (VER(J,3)+VER(J,2))
          X(J,4) = Z1* (VER(J,2)+VER(J,3)) + Z2* (VER(J,1)+VER(J,4))
          X(J,5) = Z1* (VER(J,2)+VER(J,4)) + Z2* (VER(J,1)+VER(J,3))
          X(J,6) = Z1* (VER(J,3)+VER(J,4)) + Z2* (VER(J,1)+VER(J,2))
  125 CONTINUE
      GO TO 200
C
C  Generator with homogeneous coordinates (z1,z2,z3;z3)
C
  130 NUMBER = 12
      Z1 = GENER(1)
      Z2 = GENER(2)
      Z3 = GENER(3)
      DO 135 J = 1,3
          X(J,1) = Z1*VER(J,1) + Z2*VER(J,2) + Z3* (VER(J,3)+VER(J,4))
          X(J,2) = Z1*VER(J,1) + Z2*VER(J,3) + Z3* (VER(J,2)+VER(J,4))
          X(J,3) = Z1*VER(J,1) + Z2*VER(J,4) + Z3* (VER(J,2)+VER(J,3))
          X(J,4) = Z1*VER(J,2) + Z2*VER(J,1) + Z3* (VER(J,3)+VER(J,4))
          X(J,5) = Z1*VER(J,2) + Z2*VER(J,3) + Z3* (VER(J,1)+VER(J,4))
          X(J,6) = Z1*VER(J,2) + Z2*VER(J,4) + Z3* (VER(J,1)+VER(J,3))
          X(J,7) = Z1*VER(J,3) + Z2*VER(J,1) + Z3* (VER(J,2)+VER(J,4))
          X(J,8) = Z1*VER(J,3) + Z2*VER(J,2) + Z3* (VER(J,1)+VER(J,4))
          X(J,9) = Z1*VER(J,3) + Z2*VER(J,4) + Z3* (VER(J,1)+VER(J,2))
          X(J,10) = Z1*VER(J,4) + Z2*VER(J,1) + Z3* (VER(J,2)+VER(J,3))
          X(J,11) = Z1*VER(J,4) + Z2*VER(J,2) + Z3* (VER(J,1)+VER(J,3))
          X(J,12) = Z1*VER(J,4) + Z2*VER(J,3) + Z3* (VER(J,1)+VER(J,2))
  135 CONTINUE
  200 CONTINUE
      CALL FUNSUB(X(1,1),NUMFUN,SUMVAL)
      DO 220 J = 2,NUMBER
          CALL FUNSUB(X(1,J),NUMFUN,WORK)
          DO 210 I = 1,NUMFUN
              SUMVAL(I) = SUMVAL(I) + WORK(I)
  210     CONTINUE
  220 CONTINUE
C
C***END SORTET
C
      RETURN
      END
      SUBROUTINE SRLTET(VER,NUMFUN,FUNSUB,NULL,BASVAL,RGNERR,GREATE,
     +                  SUMVAL)
C***BEGIN PROLOGUE SRLTET
C***REFER TO DCUTET
C***PURPOSE  To compute basic integration rule values and
C            corresponding error estimates.
C***DESCRIPTION SRLTET computes basic integration rule values
C            for a vector of integrands over a tetrahedron.
C            SRLTET also computes estimates for the errors by
C            using several null rule approximations.
C   ON ENTRY
C
C   VER    Real array of dimension (3,4).
C          The coordinates of the vertices of the tetrahedron.
C          vertex i -> ( ver(1,i),ver(2,i),ver(3,i) )
C   NUMFUN Integer.
C          Number of components of the vector integrand.
C   FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (X,NUMFUN,FUNVLS)
C            Input parameters:
C              X(1)   The x-coordinate of the evaluation point.
C              X(2)   The y-coordinate of the evaluation point.
C              X(3)   The z-coordinate of the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of the vector integrand.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C   NULL   Real array of dimension (NUMFUN,6).
C          A work array.
C   SUMVAL Real array of dimension (NUMFUN).
C          A work array
C
C   ON RETURN
C
C   BASVAL Real array of dimension NUMFUN.
C          The values for the basic rule for each component
C          of the integrand.
C   RGNERR Real array of dimension NUMFUN.
C          The error estimates for each component of the integrand.
C   GREATE Real.
C          The greatest error component of the integrand.
C
C***REFERENCES Beckers M., A. Haegemans,
C              The construction of cubature formula for the tetrahedron
C              Report TW128, K.U. Leuven (1990).
C***ROUTINES CALLED
C                   R1MACH,SORTET,FUNSUB
C***END PROLOGUE SRLTET
C
C   Parameters
C
C   ORBITS Integer
C          The number of orbits of the cubature formula and null rules
C   CRIVAL Real
C          The decision to choose the optimistic part of the error
C          estimator is based on CRIVAL
C   FACOPT Real
C          FACOPT is the safety coefficient used in the optimistic part
C          of the error estimator.
C   FACMED Real
C          FACMED is the safety coefficient used in the non-optimistic
C          part of the error estimator. FACMED is related to CRIVAL
C          and FACOPT.
C
      INTEGER ORBITS
      REAL CRIVAL,FACOPT,FACMED
      PARAMETER (ORBITS=7)
      PARAMETER (CRIVAL=0.5,FACMED=5,FACOPT=FACMED/CRIVAL)
C
C   Global variables.
C
      INTEGER NUMFUN
      REAL VER(3,4),BASVAL(NUMFUN),RGNERR(NUMFUN),GREATE,SUMVAL(NUMFUN),
     +     NULL(NUMFUN,6)
      EXTERNAL FUNSUB
C
C   Local variables.
C
C   K      Integer array of dimension (0:3) that contains the structure
C          parameters. K(I) = number of orbits of type I.
C   TYPE1  Real array of dimension (K(1)).
C          Contains the first homogeneous coordinate of the generators
C          of type 1
C   TYPE2  Real array of dimension (K(2)).
C          Contains the first homogeneous coordinate of the generators
C          of type 2
C   TYPE3  Real array of dimension (2,K(2)).
C          Contains the first two homogeneous coordinates of
C          the generators of type 3.
C   WEIGHT Real array of dimension (9,ORBITS).
C          The weights of the cubature formula and the null rules.
C          WEIGHT(1,1) ,..., WEIGHT(1,ORBITS) are the weights of the
C                cubature formula
C          WEIGHT(I,1) ,..., WEIGHT(I,ORBITS) for I > 1, are the weights
C                of the null rules
C
C
      INTEGER TYPE,NR,K(0:3),I,J,P
      REAL Z(3),TYPE1(3),TYPE2(1),TYPE3(2,2),NOISE,R1MACH,
     +     WEIGHT(7,ORBITS),VOLUME,DEG4,DEG3,DEG1,R2,R1,R,TRES
C
C  Cubature formula of degree 8 with 43 points
C
C  Structure parameters
C
      DATA (K(I),I=0,3)/1,3,1,2/
C
C  Information for the generators
C
      DATA (TYPE1(I),I=1,3)/0.379510205167980387748057300876E+00,
     +     0.753689235068359830728182577696E+00,
     +     0.982654148484406008240470085259E+00/
      DATA (TYPE2(I),I=1,1)/0.449467259981105775574375471447E+00/
      DATA ((TYPE3(I,J),I=1,2),J=1,2)/
     +     0.506227344977843677082264893876E+00,
     +     0.356395827885340437169173969841E-01,
     +     0.736298458958971696943019005441E+00,
     +     0.190486041934633455699433285302E+00/
C
C  Weights of the cubature formula
C
      DATA (WEIGHT(1,I),I=1,ORBITS)/-
     +     0.123001131951839495043519102752E+00,
     +     0.855018349372014074906384482699E-01,
     +     0.118021998788034059253768205083E-01,
     +     0.101900465455732427902646736855E-02,
     +     0.274781029468036908044610867719E-01,
     +     0.342269148520915110408153517904E-01,
     +     0.128431148469725555789001180031E-01/
C
C  Weights of the null rule of degree 5
C
      DATA (WEIGHT(2,I),I=1,ORBITS)/0.211921237628032658308230999090E+00
     +     ,-0.660207516445726284649283745987E-01,
     +     0.225058824086711710443385047042E-01,
     +     -0.375962972067425589765730699401E-03,
     +     0.710066020561055159657284834784E-02,
     +     0.156515256061747694921427149028E-02,
     +     -0.814530839643584660306807872526E-02/
C
C  Weights of null rule of degree 4
C
      DATA (WEIGHT(3,I),I=1,ORBITS)/-
     +     0.508105488137100551376844924797E-01,
     +     0.104596681151665328209751420525E-01,
     +     0.927471438532788763594989973184E-01,
     +     0.210489990008917994323967321174E-02,
     +     0.379184172251962722213408547663E-01,
     +     -0.111747242913563605790923001557E-01,
     +     -0.386541758762774673113423570465E-01/
C
C  Weights of first null rule of degree 3
C
      DATA (WEIGHT(4,I),I=1,ORBITS)/-
     +     0.775992773232808462404390159802E-01,
     +     -0.527453289659022924847298408064E-01,
     +     0.145876238555932704488677626554E-01,
     +     0.739374873393616192857532718429E-02,
     +     -0.374618791364332892611678523428E-01,
     +     0.538502846550653076078817013885E-01,
     +     -0.183980865177843057548322735665E-01/
C
C  Weights of second null rule of degree 3
C
      DATA (WEIGHT(5,I),I=1,ORBITS)/0.181767621501470154602720474731E-01
     +     ,0.179938831310058580533178529022E-01,
     +     0.713210362750414891598257378898E-01,
     +     -0.443935688958258805893448212636E-01,
     +     -0.657639036547720234169662790056E-01,
     +     -0.101551807522541414699808460583E-01,
     +     0.265486188970540796821750584204E-01/
C
C  Weights of null rule of degree 2
C
      DATA (WEIGHT(6,I),I=1,ORBITS)/-
     +     0.867629853722843888927184699428E-01,
     +     -0.715881271235661902772072127812E-01,
     +     0.886720767790426261677273459523E-02,
     +     -0.577885573028655167063092577589E-01,
     +     0.430310167581202031805055255554E-01,
     +     -0.606467834856775537069463817445E-02,
     +     0.319492443333738343104163265406E-01/
C
C   Weights of null rule of degree 1
C
      DATA (WEIGHT(7,I),I=1,ORBITS)/0.510374015624925451319499382594E-01
     +     ,0.463998830432033721597269299429E-01,
     +     -0.191086148397852799983451475821E-01,
     +     -0.973768821003670776204287367278E-01,
     +     0.180352562073914141268335496511E-01,
     +     0.277129527093489643801598303110E-01,
     +     -0.176218263109360550515567818653E-01/
C
C  The number of points used by the cubature formula is
C     NUM    = K(0) + 4*K(1) + 6*K(2) + 12*K(3) = 43
C
C***FIRST EXECUTABLE STATEMENT DRLTRI
      TRES = 50*R1MACH(4)
C
C  Compute the volume of the tetrahedron
C
      VOLUME = (VER(1,2)-VER(1,1))* ((VER(2,3)-VER(2,1))*
     +         (VER(3,4)-VER(3,1))- (VER(2,4)-VER(2,1))*
     +         (VER(3,3)-VER(3,1))) - (VER(2,2)-VER(2,1))*
     +         ((VER(1,3)-VER(1,1))* (VER(3,4)-VER(3,1))-
     +         (VER(1,4)-VER(1,1))* (VER(3,3)-VER(3,1))) +
     +         (VER(3,2)-VER(3,1))* ((VER(1,3)-VER(1,1))*
     +         (VER(2,4)-VER(2,1))- (VER(1,4)-VER(1,1))*
     +         (VER(2,3)-VER(2,1)))
      VOLUME = ABS(VOLUME)/6
C
C  Initialise BASVAL and NULL
C
      DO 10 J = 1,NUMFUN
          BASVAL(J) = 0.
          DO 20 I = 1,6
              NULL(J,I) = 0.
   20     CONTINUE
   10 CONTINUE
      P = 1
C
C  Compute contributions from orbits with 1, 4, 6 and 12  points
C
      DO 4 TYPE = 0,3
          DO 3 NR = 1,K(TYPE)
              IF (TYPE.EQ.1) THEN
C                 Generator ( z(1) , z(2), z(2) ; z(2) )
                  Z(1) = TYPE1(NR)
                  Z(2) = (1-Z(1))/3
              ELSE IF (TYPE.EQ.2) THEN
C                 Generator ( z(1) , z(1), z(2) ; z(2) )
                  Z(1) = TYPE2(NR)
                  Z(2) = (1-2*Z(1))/2
              ELSE
C                 Generator ( z(1) , z(2), z(3) ; z(3) )
                  Z(1) = TYPE3(1,NR)
                  Z(2) = TYPE3(2,NR)
                  Z(3) = (1-Z(1)-Z(2))/2
              END IF
C                                RGNERR is work array for DORTHE
              CALL SORTET(TYPE,Z,VER,NUMFUN,FUNSUB,SUMVAL,RGNERR)
              DO 2 J = 1,NUMFUN
                  BASVAL(J) = BASVAL(J) + WEIGHT(1,P)*SUMVAL(J)
                  DO 1 I = 1,6
                      NULL(J,I) = NULL(J,I) + WEIGHT(I+1,P)*SUMVAL(J)
    1             CONTINUE
    2         CONTINUE
              P = P + 1
    3     CONTINUE
    4 CONTINUE
C
C  Compute error estimates
C
      GREATE = 0
      DO 30 J = 1,NUMFUN
          NOISE = ABS(BASVAL(J))*TRES
          DEG4 = SQRT(NULL(J,1)**2+NULL(J,2)**2)
          DEG3 = SQRT(NULL(J,3)**2+NULL(J,4)**2)
          IF (DEG4.LE.NOISE) THEN
              RGNERR(J) = NOISE
          ELSE
              DEG1 = SQRT(NULL(J,5)**2+NULL(J,6)**2)
              IF (DEG3.NE.0) THEN
                  R1 = (DEG4/DEG3)**2
              ELSE
                  R1 = 1
              END IF
              IF (DEG1.NE.0) THEN
                  R2 = DEG3/DEG1
              ELSE
                  R2 = 1
              END IF
              R = MAX(R1,R2)
              IF (R.GE.CRIVAL) THEN
                  RGNERR(J) = FACMED*R*DEG4
              ELSE
                  RGNERR(J) = FACOPT* (R**2)*DEG4
              END IF
              RGNERR(J) = MAX(NOISE,RGNERR(J))
          END IF
          RGNERR(J) = VOLUME*RGNERR(J)
          BASVAL(J) = VOLUME*BASVAL(J)
          GREATE = MAX(GREATE,RGNERR(J))
   30 CONTINUE
C
C***END SRLTET
C
      RETURN
      END
      PROGRAM STEST1
      REAL VER(3,4,50),EPSABS,EPSREL,RESULT(1),ABSERR(1),WORK(1000)
      INTEGER NUMFUN,NUMTET,MINPTS,MAXPTS,LENVER,NW,RESTAR,NEVAL,IFAIL,
     +        IWORK(51)
      EXTERNAL F

      VER(1,1,1) = 0.
      VER(1,2,1) = 1.
      VER(1,3,1) = 0.
      VER(1,4,1) = 0.
      VER(2,1,1) = 0.
      VER(2,2,1) = 0.
      VER(2,3,1) = 1.
      VER(2,4,1) = 0.
      VER(3,1,1) = 0.
      VER(3,2,1) = 0.
      VER(3,3,1) = 0.
      VER(3,4,1) = 1.
      EPSABS = 0.
      EPSREL = 1.E-5
      NUMFUN = 1
      NUMTET = 1
      MINPTS = 0
      MAXPTS = 1000
      LENVER = 50
      NW = 1000
      RESTAR = 0
      CALL SCUTET(F,NUMFUN,VER,NUMTET,MINPTS,MAXPTS,EPSABS,EPSREL,
     +            LENVER,NW,RESTAR,RESULT,ABSERR,NEVAL,IFAIL,WORK,IWORK)
      WRITE (*,FMT=*) 'RESULT=',RESULT(1)
      WRITE (*,FMT=*) 'ERROR =',ABSERR(1)
      WRITE (*,FMT=*) 'NEVAL =',NEVAL
      WRITE (*,FMT=*) 'IFAIL =',IFAIL
      END
      SUBROUTINE F(X,NUMFUN,FUNVLS)
      INTEGER NUMFUN
      REAL X(3),FUNVLS(NUMFUN)

      FUNVLS(1) = EXP(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
      RETURN
      END
      PROGRAM STEST2
C
C   STEST2 checks that SCUTET integrates to machine precision
C   all monomials of degree less than or equal to 8.
C   STEST2 checks that SCUTET integrates to machine precision
C   all invariant polynomials that were used to construct the
C   cubature formula of degree 8 used by SRLTET.
C   STEST2 checks that the restart feature of SCUTET works and
C   that the routine works with vectors of functions.
C
C   Output produced on a SUN 3/50
C
C        
C       INTEGRATE MONOMIALS OF DEGREE 0
C            SUBROUTINE CALLS =     43, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.00E+00       0.60E-05   1.0000000000000000
C        
C       INTEGRATE MONOMIALS OF DEGREE 1
C            SUBROUTINE CALLS =     43, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.00E+00       0.60E-05   1.0000000000000000
C               2       0.60E-07       0.60E-05   0.9999999403953552
C               3       0.60E-07       0.60E-05   0.9999999403953552
C        
C       INTEGRATE MONOMIALS OF DEGREE 2
C            SUBROUTINE CALLS =     43, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.00E+00       0.60E-05   1.0000000000000000
C               2       0.12E-06       0.60E-05   0.9999998807907104
C               3       0.00E+00       0.60E-05   1.0000000000000000
C               4       0.12E-06       0.60E-05   0.9999998807907104
C               5       0.60E-07       0.60E-05   0.9999999403953552
C               6       0.60E-07       0.60E-05   0.9999999403953552
C        
C       INTEGRATE MONOMIALS OF DEGREE 3
C            SUBROUTINE CALLS =     43, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.00E+00       0.60E-05   1.0000000000000000
C               2       0.60E-07       0.60E-05   0.9999999403953552
C               3       0.00E+00       0.60E-05   1.0000000000000000
C               4       0.00E+00       0.60E-05   1.0000000000000000
C               5       0.00E+00       0.60E-05   1.0000000000000000
C               6       0.00E+00       0.60E-05   1.0000000000000000
C               7       0.00E+00       0.60E-05   1.0000000000000000
C               8       0.24E-06       0.60E-05   1.0000002384185791
C               9       0.00E+00       0.60E-05   1.0000000000000000
C              10       0.00E+00       0.60E-05   1.0000000000000000
C        
C       INTEGRATE MONOMIALS OF DEGREE 4
C            SUBROUTINE CALLS =     43, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.00E+00       0.60E-05   1.0000000000000000
C               2       0.00E+00       0.60E-05   1.0000000000000000
C               3       0.00E+00       0.60E-05   1.0000000000000000
C               4       0.00E+00       0.60E-05   1.0000000000000000
C               5       0.00E+00       0.60E-05   1.0000000000000000
C               6       0.00E+00       0.60E-05   1.0000000000000000
C               7       0.12E-06       0.60E-05   1.0000001192092896
C               8       0.12E-06       0.60E-05   1.0000001192092896
C               9       0.12E-06       0.60E-05   1.0000001192092896
C              10       0.00E+00       0.60E-05   1.0000000000000000
C              11       0.12E-06       0.60E-05   1.0000001192092896
C              12       0.00E+00       0.60E-05   1.0000000000000000
C              13       0.60E-07       0.60E-05   0.9999999403953552
C              14       0.12E-06       0.60E-05   1.0000001192092896
C              15       0.00E+00       0.60E-05   1.0000000000000000
C        
C       INTEGRATE MONOMIALS OF DEGREE 5
C            SUBROUTINE CALLS =    387, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.00E+00       0.38E-04   1.0000000000000000
C               2       0.00E+00       0.18E-02   1.0000000000000000
C               3       0.00E+00       0.96E-03   1.0000000000000000
C               4       0.00E+00       0.96E-03   1.0000000000000000
C               5       0.00E+00       0.18E-02   1.0000000000000000
C               6       0.00E+00       0.38E-04   1.0000000000000000
C               7       0.00E+00       0.91E-03   1.0000000000000000
C               8       0.00E+00       0.27E-02   1.0000000000000000
C               9       0.00E+00       0.11E-02   1.0000000000000000
C              10       0.12E-06       0.11E-02   1.0000001192092896
C              11       0.00E+00       0.18E-02   1.0000000000000000
C              12       0.00E+00       0.27E-03   1.0000000000000000
C              13       0.00E+00       0.41E-03   1.0000000000000000
C              14       0.00E+00       0.11E-02   1.0000000000000000
C              15       0.00E+00       0.96E-03   1.0000000000000000
C              16       0.00E+00       0.27E-03   1.0000000000000000
C              17       0.00E+00       0.27E-02   1.0000000000000000
C              18       0.00E+00       0.96E-03   1.0000000000000000
C              19       0.00E+00       0.91E-03   1.0000000000000000
C              20       0.00E+00       0.18E-02   1.0000000000000000
C              21       0.00E+00       0.38E-04   1.0000000000000000
C        
C       INTEGRATE MONOMIALS OF DEGREE 6
C            SUBROUTINE CALLS =    731, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.24E-06       0.16E-03   0.9999997615814209
C               2       0.00E+00       0.28E-02   1.0000000000000000
C               3       0.12E-06       0.46E-03   1.0000001192092896
C               4       0.12E-06       0.34E-02   1.0000001192092896
C               5       0.12E-06       0.94E-03   1.0000001192092896
C               6       0.00E+00       0.18E-02   1.0000000000000000
C               7       0.24E-06       0.13E-03   0.9999997615814209
C               8       0.00E+00       0.16E-02   1.0000000000000000
C               9       0.12E-06       0.94E-02   1.0000001192092896
C              10       0.12E-06       0.24E-02   1.0000001192092896
C              11       0.12E-06       0.23E-02   1.0000001192092896
C              12       0.12E-06       0.42E-02   1.0000001192092896
C              13       0.12E-06       0.18E-02   1.0000001192092896
C              14       0.12E-06       0.49E-03   1.0000001192092896
C              15       0.12E-06       0.50E-02   1.0000001192092896
C              16       0.24E-06       0.78E-02   1.0000002384185791
C              17       0.00E+00       0.23E-02   1.0000000000000000
C              18       0.12E-06       0.94E-03   1.0000001192092896
C              19       0.00E+00       0.20E-02   1.0000000000000000
C              20       0.12E-06       0.50E-02   1.0000001192092896
C              21       0.24E-06       0.24E-02   1.0000002384185791
C              22       0.00E+00       0.34E-02   1.0000000000000000
C              23       0.12E-06       0.49E-03   1.0000001192092896
C              24       0.12E-06       0.94E-02   1.0000001192092896
C              25       0.00E+00       0.46E-03   1.0000000000000000
C              26       0.00E+00       0.16E-02   1.0000000000000000
C              27       0.00E+00       0.28E-02   1.0000000000000000
C              28       0.24E-06       0.16E-03   0.9999997615814209
C        
C       INTEGRATE MONOMIALS OF DEGREE 7
C            SUBROUTINE CALLS =   2107, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.00E+00       0.19E-03   1.0000000000000000
C               2       0.00E+00       0.12E-02   1.0000000000000000
C               3       0.12E-06       0.11E-02   1.0000001192092896
C               4       0.12E-06       0.15E-02   1.0000001192092896
C               5       0.12E-06       0.12E-02   1.0000001192092896
C               6       0.00E+00       0.15E-02   1.0000000000000000
C               7       0.60E-07       0.23E-02   0.9999999403953552
C               8       0.12E-06       0.16E-03   0.9999998807907104
C               9       0.00E+00       0.12E-02   1.0000000000000000
C              10       0.12E-06       0.65E-02   1.0000001192092896
C              11       0.00E+00       0.18E-02   1.0000000000000000
C              12       0.24E-06       0.79E-02   1.0000002384185791
C              13       0.12E-06       0.76E-03   1.0000001192092896
C              14       0.00E+00       0.38E-02   1.0000000000000000
C              15       0.00E+00       0.19E-02   1.0000000000000000
C              16       0.00E+00       0.99E-03   1.0000000000000000
C              17       0.00E+00       0.19E-02   1.0000000000000000
C              18       0.00E+00       0.24E-02   1.0000000000000000
C              19       0.00E+00       0.60E-02   1.0000000000000000
C              20       0.00E+00       0.16E-02   1.0000000000000000
C              21       0.00E+00       0.10E-02   1.0000000000000000
C              22       0.00E+00       0.15E-02   1.0000000000000000
C              23       0.12E-06       0.80E-02   1.0000001192092896
C              24       0.00E+00       0.58E-02   1.0000000000000000
C              25       0.24E-06       0.56E-02   1.0000002384185791
C              26       0.12E-06       0.97E-03   1.0000001192092896
C              27       0.12E-06       0.12E-02   1.0000001192092896
C              28       0.12E-06       0.91E-03   1.0000001192092896
C              29       0.00E+00       0.16E-02   1.0000000000000000
C              30       0.12E-06       0.97E-03   1.0000001192092896
C              31       0.00E+00       0.15E-02   1.0000000000000000
C              32       0.00E+00       0.39E-02   1.0000000000000000
C              33       0.00E+00       0.10E-02   1.0000000000000000
C              34       0.60E-07       0.23E-02   0.9999999403953552
C              35       0.00E+00       0.19E-02   1.0000000000000000
C              36       0.12E-06       0.16E-03   0.9999998807907104
C        
C       INTEGRATE MONOMIALS OF DEGREE 8
C            SUBROUTINE CALLS =   2795, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.00E+00       0.85E-04   1.0000000000000000
C               2       0.00E+00       0.40E-03   1.0000000000000000
C               3       0.12E-06       0.72E-03   1.0000001192092896
C               4       0.12E-06       0.62E-04   1.0000001192092896
C               5       0.00E+00       0.88E-03   1.0000000000000000
C               6       0.12E-06       0.62E-04   1.0000001192092896
C               7       0.24E-06       0.72E-03   1.0000002384185791
C               8       0.00E+00       0.40E-03   1.0000000000000000
C               9       0.00E+00       0.85E-04   1.0000000000000000
C              10       0.00E+00       0.32E-03   1.0000000000000000
C              11       0.24E-06       0.87E-03   1.0000002384185791
C              12       0.00E+00       0.22E-03   1.0000000000000000
C              13       0.00E+00       0.56E-03   1.0000000000000000
C              14       0.00E+00       0.46E-03   1.0000000000000000
C              15       0.00E+00       0.23E-03   1.0000000000000000
C              16       0.12E-06       0.82E-03   1.0000001192092896
C              17       0.00E+00       0.41E-03   1.0000000000000000
C              18       0.12E-06       0.53E-03   1.0000001192092896
C              19       0.00E+00       0.34E-03   1.0000000000000000
C              20       0.12E-06       0.89E-03   1.0000001192092896
C              21       0.00E+00       0.28E-02   1.0000000000000000
C              22       0.12E-06       0.10E-02   1.0000001192092896
C              23       0.00E+00       0.30E-03   1.0000000000000000
C              24       0.24E-06       0.73E-03   1.0000002384185791
C              25       0.12E-06       0.61E-04   1.0000001192092896
C              26       0.00E+00       0.11E-02   1.0000000000000000
C              27       0.00E+00       0.88E-03   1.0000000000000000
C              28       0.00E+00       0.27E-02   1.0000000000000000
C              29       0.00E+00       0.89E-03   1.0000000000000000
C              30       0.12E-06       0.74E-04   1.0000001192092896
C              31       0.00E+00       0.88E-03   1.0000000000000000
C              32       0.00E+00       0.11E-02   1.0000000000000000
C              33       0.12E-06       0.84E-03   1.0000001192092896
C              34       0.00E+00       0.93E-03   1.0000000000000000
C              35       0.00E+00       0.88E-03   1.0000000000000000
C              36       0.12E-06       0.71E-04   1.0000001192092896
C              37       0.00E+00       0.42E-03   1.0000000000000000
C              38       0.00E+00       0.21E-03   1.0000000000000000
C              39       0.12E-06       0.60E-04   1.0000001192092896
C              40       0.24E-06       0.72E-03   1.0000002384185791
C              41       0.24E-06       0.90E-03   1.0000002384185791
C              42       0.12E-06       0.55E-03   1.0000001192092896
C              43       0.00E+00       0.38E-03   1.0000000000000000
C              44       0.00E+00       0.34E-03   1.0000000000000000
C              45       0.00E+00       0.86E-04   1.0000000000000000
C        
C       INTEGRATE INVARIANT POLYNOMIALS
C            SUBROUTINE CALLS =   1763, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.00E+00       0.60E-05   1.0000000000000000
C               2       0.00E+00       0.60E-05   1.0000000000000000
C               3       0.00E+00       0.60E-05   1.0000000000000000
C               4       0.12E-06       0.60E-05   1.0000001192092896
C               5       0.00E+00       0.60E-05   1.0000000000000000
C               6       0.12E-06       0.60E-05   1.0000001192092896
C               7       0.00E+00       0.38E-03   1.0000000000000000
C               8       0.24E-06       0.64E-02   1.0000002384185791
C               9       0.24E-06       0.32E-01   1.0000002384185791
C              10       0.12E-06       0.83E-05   1.0000001192092896
C              11       0.24E-06       0.36E-04   1.0000002384185791
C              12       0.12E-06       0.92E-03   1.0000001192092896
C              13       0.48E-06       0.20E-01   1.0000004768371582
C              14       0.12E-06       0.35E-01   1.0000001192092896
C              15       0.60E-06       0.11E-01   1.0000005960464478
C            SUBROUTINE CALLS =   2795, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.60E-07       0.60E-05   0.9999999403953552
C               2       0.60E-07       0.60E-05   0.9999999403953552
C               3       0.00E+00       0.62E-05   1.0000000000000000
C               4       0.00E+00       0.60E-05   1.0000000000000000
C               5       0.12E-06       0.60E-05   1.0000001192092896
C               6       0.12E-06       0.11E-04   1.0000001192092896
C               7       0.00E+00       0.10E-03   1.0000000000000000
C               8       0.24E-06       0.16E-02   1.0000002384185791
C               9       0.24E-06       0.80E-02   1.0000002384185791
C              10       0.12E-06       0.27E-04   1.0000001192092896
C              11       0.24E-06       0.48E-03   1.0000002384185791
C              12       0.24E-06       0.25E-03   1.0000002384185791
C              13       0.48E-06       0.52E-02   1.0000004768371582
C              14       0.12E-06       0.91E-02   1.0000001192092896
C              15       0.60E-06       0.33E-02   1.0000005960464478
C            SUBROUTINE CALLS =   3139, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.00E+00       0.60E-05   1.0000000000000000
C               2       0.00E+00       0.60E-05   1.0000000000000000
C               3       0.00E+00       0.63E-05   1.0000000000000000
C               4       0.00E+00       0.60E-05   1.0000000000000000
C               5       0.12E-06       0.60E-05   1.0000001192092896
C               6       0.00E+00       0.12E-04   1.0000000000000000
C               7       0.00E+00       0.13E-04   1.0000000000000000
C               8       0.24E-06       0.41E-04   1.0000002384185791
C               9       0.24E-06       0.62E-04   1.0000002384185791
C              10       0.12E-06       0.32E-04   1.0000001192092896
C              11       0.24E-06       0.52E-03   1.0000002384185791
C              12       0.24E-06       0.26E-04   1.0000002384185791
C              13       0.60E-06       0.20E-03   1.0000005960464478
C              14       0.12E-06       0.42E-03   1.0000001192092896
C              15       0.60E-06       0.64E-03   1.0000005960464478
C            SUBROUTINE CALLS =  17587, IFAIL =  0
C              N   ABSOLUTE ERROR  ESTIMATED ERROR    INTEGRAL
C               1       0.00E+00       0.60E-05   1.0000000000000000
C               2       0.60E-07       0.60E-05   0.9999999403953552
C               3       0.00E+00       0.63E-05   1.0000000000000000
C               4       0.00E+00       0.60E-05   1.0000000000000000
C               5       0.12E-06       0.60E-05   1.0000001192092896
C               6       0.00E+00       0.65E-05   1.0000000000000000
C               7       0.00E+00       0.60E-05   1.0000000000000000
C               8       0.12E-06       0.87E-05   1.0000001192092896
C               9       0.24E-06       0.12E-04   1.0000002384185791
C              10       0.60E-07       0.94E-05   0.9999999403953552
C              11       0.12E-06       0.20E-04   1.0000001192092896
C              12       0.60E-07       0.10E-04   0.9999999403953552
C              13       0.12E-06       0.62E-04   1.0000001192092896
C              14       0.12E-06       0.61E-04   1.0000001192092896
C              15       0.48E-06       0.10E-03   1.0000004768371582
C
C   End of output
C
      INTEGER LENVER,NW,DEGREE
      PARAMETER (LENVER=1300,NW=200000)
      REAL VER(3,4,LENVER),EPSREL,WORK(NW),FACTOR
      INTEGER NUMFUN,MAXPTS,RESTAR,IWORK(LENVER+1),I
      COMMON /INFOFC/DEGREE
      EXTERNAL FMONOM,FINVAR
C
C  Test for integrating monomials of degree <= 8
C
      DO 10 DEGREE = 0,8
          WRITE (*,FMT=90) DEGREE
   90     FORMAT (' ',/,'INTEGRATE MONOMIALS OF DEGREE ',I1)
          VER(1,1,1) = 0.
          VER(1,2,1) = 1.
          VER(1,3,1) = 0.
          VER(1,4,1) = 0.
          VER(2,1,1) = 0.
          VER(2,2,1) = 0.
          VER(2,3,1) = 1.
          VER(2,4,1) = 0.
          VER(3,1,1) = 0.
          VER(3,2,1) = 0.
          VER(3,3,1) = 0.
          VER(3,4,1) = 1.
          EPSREL = 1.E-2
          NUMFUN = FACTOR(DEGREE+2)/ (FACTOR(DEGREE)*2)
          MAXPTS = 50000
          RESTAR = 0
          CALL ATEST(FMONOM,NUMFUN,VER,1,MAXPTS,EPSREL,LENVER,NW,RESTAR,
     +               WORK,IWORK)
   10 CONTINUE
c
c Test restart feature on the special basis of invariant polynomials
c
      WRITE (*,FMT=91)
   91 FORMAT (' ',/,'INTEGRATE INVARIANT POLYNOMIALS')
      EPSREL = 1
      NUMFUN = 15
      MAXPTS = 50000
      RESTAR = 0
      VER(1,1,1) = 0.
      VER(1,2,1) = SQRT(2.)/3.
      VER(1,3,1) = SQRT(2.)/3.
      VER(1,4,1) = -2.*SQRT(2.)/3.
      VER(2,1,1) = 0.
      VER(2,2,1) = SQRT(6.)/3.
      VER(2,3,1) = -SQRT(6.)/3.
      VER(2,4,1) = 0.
      VER(3,1,1) = 1.
      VER(3,2,1) = -1./3.
      VER(3,3,1) = -1./3.
      VER(3,4,1) = -1./3.
      DO 20 I = 1,4
          EPSREL = EPSREL/10
          CALL ATEST(FINVAR,NUMFUN,VER,1,MAXPTS,EPSREL,LENVER,NW,RESTAR,
     +               WORK,IWORK)
          RESTAR = 1
   20 CONTINUE
      STOP

      END
      SUBROUTINE ATEST(F,NUMFUN,VER,NUMTET,MAXPTS,EPSREL,LENVER,NW,
     +                 RESTAR,WORK,IWORK)
      INTEGER NUMFUN,NUMTET,MAXPTS,LENVER,NW,RESTAR,IWORK(LENVER+1),
     +        NEVAL,IFAIL,N
      REAL VER(3,4,LENVER),EPSREL,WORK(NW),RESULT(165),ABSERR(165)
      EXTERNAL F

      CALL SCUTET(F,NUMFUN,VER,NUMTET,0,MAXPTS,0.,EPSREL,LENVER,NW,
     +            RESTAR,RESULT,ABSERR,NEVAL,IFAIL,WORK,IWORK)
      WRITE (*,FMT=90) NEVAL,IFAIL
   90 FORMAT (5X,'SUBROUTINE CALLS = ',I6,', IFAIL = ',I2)
      WRITE (*,FMT=91)
   91 FORMAT (7X,'N',3X,'ABSOLUTE ERROR',2X,'ESTIMATED ERROR',4X,
     +       'INTEGRAL')
      DO 10 N = 1,NUMFUN
          WRITE (*,FMT=92) N,ABS(RESULT(N)-1),ABSERR(N),RESULT(N)
   92     FORMAT (6X,I3,2E15.2,F21.16)
   10 CONTINUE
      END
      SUBROUTINE FMONOM(X,NUMFUN,VALUES)
C
C   Evaluate the monomials of degree DEGREE
C
      INTEGER NUMFUN,COUNT,DEGREE,I,J
      REAL X(3),VALUES(NUMFUN),MOMENT
      COMMON /INFOFC/DEGREE

      COUNT = 1
      DO 20 J = 0,DEGREE
          DO 10 I = 0,DEGREE - J
              VALUES(COUNT) = X(1)**J*X(2)**I*X(3)** (DEGREE-J-I)
              VALUES(COUNT) = VALUES(COUNT)/MOMENT(J,I,DEGREE-J-I)
              COUNT = COUNT + 1
              IF (COUNT.GT.NUMFUN) RETURN
   10     CONTINUE
   20 CONTINUE
      RETURN
      END
      REAL FUNCTION FACTOR(L)
C
C  Computes the factorial of L
C
      INTEGER L,K

      FACTOR = 1
      DO 10 K = 2,L
          FACTOR = FACTOR*K
   10 CONTINUE
      RETURN
      END
      REAL FUNCTION MOMENT(I,J,K)
C
C  Computes the integral of x(1)**I * x(2)**J * x(3)**K
C
      INTEGER I,J,K
      REAL FACTOR

      MOMENT = FACTOR(I)*FACTOR(J)*FACTOR(K)/FACTOR(3+I+J+K)
      RETURN
      END
      SUBROUTINE FINVAR(X,NUMFUN,VALUES)
C
C  Evaluates the invariant polynomials used to construct
C  a cubature formula.
C
      INTEGER NUMFUN,I
      REAL X(3),VALUES(NUMFUN),P2,P3,P4,INVAR(15)
c
c  Basic invariant polynomials
c
      P2 = X(1)**2 + X(2)**2 + X(3)**2
      P3 = - (SQRT(2.)* (X(1)**2-3.*X(2)**2)*X(1)+3.*X(1)**2*X(3)+
     +     3.*X(2)**2*X(3)-2.*X(3)**3)/2.
      P4 = (4.*SQRT(2.)* (X(1)**2-3.*X(2)**2)*X(1)*X(3)+
     +     6.* (2.*X(2)**2+X(3)**2)*X(1)**2+6.*X(1)**4+6.*X(2)**4+
     +     6.*X(2)**2*X(3)**2+7.*X(3)**4)/7.
c
c  Basis used to construct cubature formula
c
      INVAR(1) = 1
      INVAR(2) = 5*P2
      INVAR(3) = 15*P3
      INVAR(4) = 315*P2**2/19
      INVAR(5) = (945./2.)* (7./24.)* (P2**2-P4)
      INVAR(6) = 35*P2*P3
      INVAR(7) = (1701./41.)*P2**3
      INVAR(8) = (3645./2.)* (7./24.)* (P2**3-P2*P4)
      INVAR(9) = (-8505./8.)* (3*P2**3+4*P3**2-7*P2*P4)/12
      INVAR(10) = (18900./67.)*P2**2*P3/4
      INVAR(11) = (14175./8.)* (7./6.)* (P2**2-P4)*P3
      INVAR(12) = (623700./1829.)*P2**4/4
      INVAR(13) = (155925./26.)* (7./24.)* (P2**4-P2**2*P4)
      INVAR(14) = (-467775./152.)* (1./12.)*
     +            (3*P2**4+4*P2*P3**2-7*P2**2*P4)
      INVAR(15) = (155925./232.)* (7./12.)*
     +            (7*P4**2-4*P2*P3**2-3*P2**2*P4)
      DO 10 I = 1,NUMFUN
          VALUES(I) = (9*SQRT(3.)/8)*INVAR(I)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE STRTET(DVFLAG,SBRGNS,GREATE,LIST,NEW)
C***BEGIN PROLOGUE STRTET
C***REFER TO SCUTET
C***PURPOSE  STRTET maintains a heap of subregions.
C***DESCRIPTION STRTET maintains a heap of subregions.
C            The subregions are stored in a partially sorted
C            binary tree, ordered according to the size of the
C            greatest error estimates of each subregion(GREATE).
C            The subregion with greatest error estimate is in the
C            first position of the heap.
C
C   PARAMETERS
C
C     DVFLAG Integer.
C            If DVFLAG = 1, we remove the subregion with
C            greatest error from the heap.
C            If DVFLAG = 2, we insert a new subregion in the heap.
C     SBRGNS Integer.
C            Number of subregions in the heap.
C     GREATE Real array of dimension SBRGNS.
C            Used to store the greatest estimated errors in
C            all subregions.
C     LIST   Integer array of dimension SBRGNS.
C            Used as a partially ordered list of pointers to the
C            different subregions. This list is a heap where the
C            element on top of the list is the subregion with the
C            greatest error estimate.
C     NEW    Integer.
C            Index to the new region to be inserted in the heap.
C***ROUTINES CALLED-NONE
C***END PROLOGUE STRTET
C
C   Global variables.
C
      INTEGER DVFLAG,NEW,SBRGNS,LIST(*)
      REAL GREATE(*)
C
C   Local variables.
C
C   GREAT  is used as intermediate storage for the greatest error of a
C          subregion.
C   SUBRGN Position of child/parent subregion in the heap.
C   SUBTMP Position of parent/child subregion in the heap.
      INTEGER SUBRGN,SUBTMP
      REAL GREAT
C
C***FIRST EXECUTABLE STATEMENT STRTET
C
C    If DVFLAG = 1, we will reduce the partial ordered list by the
C    element with greatest estimated error. Thus the element in
C    in the heap with index LIST(1) is vacant and can be used later.
C    Reducing the heap by one element implies that the last element
C    should be re-positioned.
C
      IF (DVFLAG.EQ.1) THEN
          GREAT = GREATE(LIST(SBRGNS))
          SBRGNS = SBRGNS - 1
          SUBRGN = 1
   20     SUBTMP = 2*SUBRGN
          IF (SUBTMP.LE.SBRGNS) THEN
              IF (SUBTMP.NE.SBRGNS) THEN
C
C   Find max. of left and right child.
C
                  IF (GREATE(LIST(SUBTMP)).LT.
     +                GREATE(LIST(SUBTMP+1))) THEN
                      SUBTMP = SUBTMP + 1
                  END IF
              END IF
C
C   Compare max.child with parent.
C   If parent is max., then done.
C
              IF (GREAT.LT.GREATE(LIST(SUBTMP))) THEN
C
C   Move the pointer at position subtmp up the heap.
C
                  LIST(SUBRGN) = LIST(SUBTMP)
                  SUBRGN = SUBTMP
                  GO TO 20
              END IF
          END IF
C
C   Update the pointer.
C
          IF (SBRGNS.GT.0) THEN
              LIST(SUBRGN) = LIST(SBRGNS+1)
          END IF
      ELSE IF (DVFLAG.EQ.2) THEN
C
C   If DVFLAG = 2, find the position for the NEW region in the heap.
C
          GREAT = GREATE(NEW)
          SUBRGN = SBRGNS
   40     SUBTMP = SUBRGN/2
          IF (SUBTMP.GE.1) THEN
C
C   Compare max.child with parent.
C   If parent is max, then done.
C
              IF (GREAT.GT.GREATE(LIST(SUBTMP))) THEN
C
C   Move the pointer at position subtmp down the heap.
C
                  LIST(SUBRGN) = LIST(SUBTMP)
                  SUBRGN = SUBTMP
                  GO TO 40
              END IF
          END IF
C
C    Set the pointer to the new region in the heap.
C
          LIST(SUBRGN) = NEW
      END IF
C
C***END STRTET
C
      RETURN
      END
