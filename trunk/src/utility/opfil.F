c----------------------------------------------------------------------
c
c   routine name       - opfil
c
c----------------------------------------------------------------------
c
c   latest revision    - Dec 07
c
c   purpose            - routine opens up files
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine opfil
c
      use control
      use parallelism
c
#include "syscom.blk"
#include "cinout.blk"
c      
      parameter(MAXSIZE=100)
      character(MAXSIZE) :: File_name
      character(MAXSIZE) :: filenameout
      character*9 :: pathin, pathout

#if PARALLEL_MODE
c  ...open files on scratch or local directory
      ionscratch = 0
c
      if (ionscratch .eq. 0) then
        pathin = ''
      else
        pathin = '/scratch/'
      endif
      pathout = pathin
c
      NIN  = 10
      filenameout = pathin // 'files/input'
      open(unit=NIN,file=filenameout,
     .     form='formatted',access='sequential',status='unknown')
c
      NOUT = 11
      filenameout = pathin // 'files/output' // PRINTRANK
      open(unit=NOUT,file=filenameout,form='formatted',
     .     access='sequential',status='unknown')
c
      NDUMP=12
c
      NHIST=13
      filenameout = pathin // 'files/history' // PRINTRANK
      open(unit=NHIST,file=filenameout,form='formatted',
     .     access='sequential',status='unknown')
c
c  ...file 'files/control' is open from main/read_control...
c
      NWE = 14
      filenameout = pathin // 'files/input_new' 
      open(unit=NWE,file=filenameout,form='formatted',
     .     access='sequential',status='unknown')
c
      NPROFB = 15
      filenameout = pathin // 'files/profilb'
      open(unit=NPROFB,file=filenameout,form='formatted',
     .     access='sequential',status='unknown')
c
      KIN = 16
      select case(INPUT_FILE)
      case(3)
        filenameout = pathin // 'files/input_compact' 
      case(4)
        filenameout = pathin // 'files/input_rec' 
      case(5)
        filenameout = pathin // 'files/input_tets' 
      end select
      open(unit=KIN,file=filenameout,form='formatted',
     .     access='sequential',status='unknown')
c
#else
      NIN  = 10
      open(unit=NIN,file='files/input',
     .     form='formatted',access='sequential',status='unknown')
c
      NOUT = 11
      open(unit=NOUT,file='files/output',form='formatted',
     .     access='sequential',status='unknown')
c
      NDUMP=12
c
      NHIST=13
      open(unit=NHIST,file='files/history',form='formatted',
     .     access='sequential',status='unknown')
c
c  ...file 'files/control' is open from main/read_control...
c
      NWE = 14
      open(unit=NWE,file='files/input_new',form='formatted',
     .     access='sequential',status='unknown')
c
      NPROFB = 15
      open(unit=NPROFB,file='files/profilb',form='formatted',
     .     access='sequential',status='unknown')
c
      KIN = 16
      select case(INPUT_FILE)
      case(3)
        open(unit=KIN,file='files/input_compact',form='formatted',
     .       access='sequential',status='unknown')
      case(4)
        open(unit=KIN,file='files/input_rec',form='formatted',
     .       access='sequential',status='unknown')
      case(5)
        open(unit=KIN,file='files/input_tets',form='formatted',
     .       access='sequential',status='unknown')
      end select
c-----------------------------------------------------------------------
c  ...geometry customization file
      N_GEO = 17
      open(unit = N_GEO, file = 'files/geometry', form = 'formatted',
     .     access = 'sequential', status = 'unknown')
c-----------------------------------------------------------------------
c
      NHIST=20
      open(unit=NHIST,file='files/history',
     .     form='formatted',access='sequential',status='unknown')
#endif
c
      return
      end
