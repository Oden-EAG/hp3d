#if HP3D_USE_X11

/*
  FORTRAN INTERFACE TO GRAPHICS (X11 Windows style):
  =======================================================
  XGOPW(narg,SD0,SD1,SD2,SD3,SD4,SD5)    init + open window
  XGDRAWL(ix0,iy0,ix1,iy1)   draw line
  XGSETCM(index,ir,ig,ib)    set color palette  !!! (not in X-window style !!)
  XGSETCOL(index)            select color nr
  XGSENDCOL(index,ir,ig,ib)  send the rgb values based on the index..
  XGLENGTH()                 length of event queue
  XGNXTEV(ev_type,ev_val,ev_x,ev_y)   event types:
  0  -  key press (ev_val = key code < 256+ 256*shifts
  1  -              ( shifts:Spec,Shift,Ctrl,Alt,Lock
  2  -  button press (         16   8    4    2   1
  3  -  button release ( ev_val = 1,2,3 + 256 *Shifts
  4  -  configure event (ev_x,ev_y - new size)
  5  -  expose event ( repetitions filtered out )
  key release - event 1 - is trapped inside
  see below for map of keyboard values
  XGDRRECT(xl,yl,xu,yu)      draw rectangle
  XGDRFRECT(xl,yl,xu,yu)     draw filled rectangle
  XGDRPOLY(npts,ix,iy)       draw polygon
  XGDRFPOLY(npts,ix,iy)      draw filled polygon
  XGSETTL(x,y)               set text position
  XGDRTXT(SD)                draw text
  XGCLOSW()                  close window
  XGCLEAR()                  clear window
  XGMOUSE(x,y)               get mouse position
  XMOUSEPOS(x,y)             get the position of the mouse key
  XGSETCLIP(xl,yl,xu,yu)     set clip region on screen
  XGINFO(width,height,maxcol,mxw,mxh)    - sideffect to XSync
  get info on window and screen (wind size, max colors, screen size)
  x,y coordinates start in lower left corner of window
  SD - text items - are assumed to end with \0 (C strings)
  xgopw strings for geom should not have blanks
  ========================================================
  XGPROGRESS () - not X windows related - displays turning symbol
  in the terminal. Need to be called repeatedly, e.g
  after processing of each element.
  ##############################################################
  HANDLING OF STRINGS IS WRONG !!! (correct only for Alliant)
  they are modified without check for room (added char(0) !)
  #############################################################
  =========================================================
  Keyboard map: All standard ASCII codes are obtained as MOD(ev_val,256)
  Capitals and special symbols are set properly (VT100 key setup assumed!)
  To check for Ctrl, Shift, Alt with special setup - check value of
  ev_val/256, eg. to check for Ctrl-Shift-Left ev_val = 7217 | 7633 (lock!)
  If result is gt 127 - function keys are pressed:
  f1 - f20 -> 190 - 209 (not mapped - f11=esc, f13=nl, f15,f16=Help,Do)
  left,up,right,down -> 209 - 212, PF1-PF4 -> 145-148, left==f20!
  Help,Menu,Find,Select,Insert,PrevScreen,NextScreen,Remove,Delete ->
  234, 231, 232, 224,   227,   213,       214,       128,   255
  to differenciate for num keypad: check Spec bit (ev_val.AND.0x1000)
  i.e. standard '0' is 48, '0' from num_pad is 48+0x1000 = 4096+48)
  ==========================================================
  Possible define's for machine dependent code:
  ARDENT
  SGI
  Alliant - defines SGI - only difference in handling strings
  IBM
  ==========================================================

*/

/* Fortan to C interfacing portability issues (not all, only
 *    those relevant here):
 *
 * 1. Function names are translated as:
 *    a)  name -> name_   (IRIX, Linux, SunOS, DEC, ....
 *    b)  name -> name    (AIX,HPUX,hpf90,WIN32, ...
 *    c)  name -> NAME    (CRAY, ...
 *   Some compilers allow to choose any of the above methods or even
 *   choose yet different variant, but it is always best to use default
 *   method, which avoids possible problems when linking with system and/or
 *   third party libraries.
 *
 * 2. Most standard arguments are passed as pointers
 *    a)  integer          -> int *
 *    b)  double precision -> double *
 *    c)  real*4           -> float *
 *   some compilers, notably f90/f95 allow more combinations, but the
 *   above ones are sufficient for most interfacing and we will limit
 *   the interfacing to such.
 *
 * 3. String arguments are most difficult to handle (the number variations
 *   defies common sense: In all known cases each string is passed as double
 *   argument: char *, int where the first one defines address of the string
 *   and the other one: its length. Note also that the string in FORTRAN
 *   needs to be augmented by NUL character (char(0)) which is automatically
 *   done for all arguments in symbol.f and openwind.f.
 *    The variation comes in the decision where to place the hidden (length)
 *   argument, esp. when there is more than a single character argument:
 *    a) some compilers add this argument immediately after each string pointer
 *    b) some add all of these arguments after the last 'real' argument
 *    c) some compilers additionally move all character arguments after all
 *       numeric variable arguments, and then use method a) or b) above
 *    d) there are more exotic combinations, e.g. FORTRAN compiler for ARDENT
 *       (this was quite fast computer some time ago ....) makes a struct
 *       combining {char * Addr; int Length; } and passes pointer to such
 *       struct in the argument list.
 *   We simplify the problem by allowing string arguments only at the end
 *   the argument list of a FORTRAN subroutine (this eliminates option c) and
 *   limiting string arguments to only two functions: xgopw and xgdrtxt
 *
 * 4. Most C compilers provide automatic #defines allowing to recognize
 *   the operating system, although sometimes it is necessary to use
 *   'manual' define on the compile command, esp. to differentiate between
 *   different platforms/versions of operating system (e.g. SGI on MIPS or
 *   Itanium architectures, HP on RISK or Itanium, Sun on Sparc or Intel x86,
 *   etc. For some machines there exist multiple compilers with potentially
 *   different methods (e.g on Linux: GNU g77, Intel IPF, Portland Group pgf,
 *   Lahey, LCC, etc.)
 *
 * 5. We do not handle Windows compilers, because Microsoft graphics
 *   differs significantly from X Windows style. If you need to run
 *   this code on Windows you need to use Unix emulators, e.g. excellent
 *   (and free) CygWin (which will fit Linux methods for our purpose)
 */

/* To solve issue 1, we define macro NEED_WRAPPER and WRAPPER_FUNC(name),
 *
 * when NEED_WRAPPER is defined as 1 then:
 * method 1a is handled with:
 *  #define WRAPPER_FUNC( lowerCaseName )  lowerCaseName##_
 *
 * and method 1b with:
 *  #define WRAPPER_FUNC( lowerCaseName )  lowerCaseName
 *
 * when NEED_WRAPPER is defined as zero then WRAPPER_FUNC is not used and
 * method 1c is ineffect.
 */

#define NEED_WRAPPER 1

#if defined(ibm) || (defined(hpux) && ! defined(hpf90)) || defined (_WIN32)

#  define WRAPPER_FUNC( lowerCaseName )  lowerCaseName

#elif (defined(CRAY))

/* For Cray, our wrapper function shouldn't do anything
 * so as to preserve the uppercase function declarations.
 * Just to make sure that this is indeed the case, we define
 * the wrapper to produce illegal names which would be caught
 * by the compiler.
 */

#   define WRAPPER_FUNC( lowerCaseName )  +=*lowerCaseName
#   undef  NEED_WRAPPER
#   define NEED_WRAPPER 0

#else

/* This case handles SGI IRIX, SunOS and Solaris, DEC,
   Linux x86, and ConvexOS.  */
#   define WRAPPER_FUNC( lowerCaseName )  lowerCaseName##_

#endif

/* To handle string length ordering we define two macro's:
 * #define STRING_LENGTHS_AFTER 0  --- to handle case 3a
 * #define STRING_LENGTHS_AFTER 1  ---     3b
 * #define STRINGS_AS_STRUCTS 1    ---     3c
 */

/* all our machines use this so far */
#define STRING_LENGTHS_AFTER 1
#define STRINGS_AS_STRUCTS 0
/* ################################################################## */

#if NEED_WRAPPER
#define XGOPW      WRAPPER_FUNC(xgopw)
#define XGDRAWL    WRAPPER_FUNC(xgdrawl)
#define XGSETCM    WRAPPER_FUNC(xgsetcm)
#define XGSETCOL   WRAPPER_FUNC(xgsetcol)
#define XGSENDCOL  WRAPPER_FUNC(xgsendcol)
#define XGLENGTH   WRAPPER_FUNC(xglength)
#define XGNXTEV    WRAPPER_FUNC(xgnxtev)
#define XGDRRECT   WRAPPER_FUNC(xgdrrect)
#define XGDRFRECT  WRAPPER_FUNC(xgdrfrect)
#define XGDRPOLY   WRAPPER_FUNC(xgdrpoly)
#define XGDRFPOLY  WRAPPER_FUNC(xgdrfpoly)
#define XGSETTL    WRAPPER_FUNC(xgsettl)
#define XGDRTXT    WRAPPER_FUNC(xgdrtxt)
#define XGCLOSW    WRAPPER_FUNC(xgclosw)
#define XGCLEAR    WRAPPER_FUNC(xgclear)
#define XGMOUSE    WRAPPER_FUNC(xgmouse)
#define XMOUSEPOS  WRAPPER_FUNC(xmousepos)
#define XGSETCLIP  WRAPPER_FUNC(xgsetclip)
#define XGINFO     WRAPPER_FUNC(xginfo)
#define XGPROGRESS WRAPPER_FUNC(xgprogress)
#define XGDOTS     WRAPPER_FUNC(xgdots)
#endif

/* Actual prototypes: */
/* extern "C" { */
void XGDRAWL   (int *ix0, int *iy0, int *ix1, int *iy1);
void XGSETCM   (int *index,int *ir,int *ig,int *ib);
void XGSETCOL  (int *index);
void XGSENDCOL (int *index,int *ir,int *ib,int *ig);
void XMOUSEPOS (int *imx,int *imy);
int  XGLENGTH  ( void );
void XGNXTEV   (int *ev_type,int *ev_val,int *ev_x,int *ev_y);
void XGSETCLIP (int *xl,int *yl,int *xu,int *yu);
void XGDRRECT  (int *xl,int *yl,int *xu,int *yu);
void XGDRFRECT (int *xl,int *yl,int *xu,int *yu);
void XGDRPOLY  (int *npts,int ix[],int iy[]);
void XGDRFPOLY (int *npts,int ix[],int iy[]);
void XGSETTL   (int *x,int *y);
void XGCLOSW   ( void );
void XGCLEAR   ( void );
void XGMOUSE   (int *x,int *y);
void XGINFO    (int *x,int *y,int *mxcol,int *mxx,int *mxy);
void XGPROGRESS( void );
void XGDOTS    ( void );

#if STRINGS_AS_STRUCTS
typedef struct {char * Addr; int Length; } Str_Desc;
void XGDRTXT(Str_Desc *SD);
void XGOPW(int *narg,Str_Desc *SD0,Str_Desc *SD1,Str_Desc *SD2,
           Str_Desc *SD3,Str_Desc *SD4,Str_Desc *SD5);
#else
void XGDRTXT(char *SD,int length);
#if STRING_LENGTHS_AFTER
void XGOPW(int *narg,char *SD0,char *SD1,char *SD2,char *SD3,char *SD4,
           char *SD5,int l0,int l1,int l2,int l3,int l4,int l5);
#else
void XGOPW(int *narg,char *SD0,int l0,char *SD1,int l1,char *SD2,
           int l2,char *SD3,int l3,char *SD4,int l4,char *SD5,int l5);
#endif
#endif


/* end of F to C interface magic */
/* ################################################################## */

#include <stdlib.h>           /* bogus include */
#include <stdio.h>            /* bogus include */
#include <ctype.h>            /* bogus include */
#include <X11/Xos.h>          /* bogus include */
#include <X11/Xlib.h>         /* bogus include */
#include <X11/Xutil.h>        /* bogus include */
#define XK_LATIN1
#include <X11/keysymdef.h>    /* bogus include */

typedef unsigned long Pixel;

static char *ProgramName;
static Display *dpy;
static Window xwindow;
static int screen;
static int window_height,window_width;
static GC gc;
static Colormap colormap;
static int textloc_x,textloc_y;
static int shift;

#define SHIFT (!(!(shift&ShiftMask)))<<3 | (!(!(shift&ControlMask)))<<2 \
  | (!(!(shift&Mod1Mask)))<<1 | !(!(shift&LockMask))

#define CHGY(A) (window_height - (A))
#define GABS(A) ((A)>0?(A):(-(A)))
#define MAXA(A,B) ((A) > (B) ? (A) : (B))
#define MINA(A,B) ((A) < (B) ? (A) : (B))

#define ProgramTitle "MeshGen"

#define MAXCOLS 4100
XColor colors[MAXCOLS];
/* .................................................................. */
static void xg_open_window (int argc,char ** argv);
static void xg_draw_fill_rectangle(int xl,int yl,int xu,int yu);
static void xg_draw_text(char *string);
static void set_sizehints (XSizeHints *hintp, int min_width, int min_height,
			   int defwidth, int defheight, int defx, int defy,
			   char *geom);
static void usage ( void );
/* .................................................................. */

/* FORTRAN interface ================ open display and window */

#if STRINGS_AS_STRUCTS
void XGOPW(int *narg,Str_Desc *SD0,Str_Desc *SD1,Str_Desc *SD2,
	   Str_Desc *SD3,Str_Desc *SD4,Str_Desc *SD5)
#else
#if STRING_LENGTHS_AFTER
  void XGOPW(int *narg,char *SD0,char *SD1,char *SD2,char *SD3,char *SD4,
             char *SD5,int l0,int l1,int l2,int l3,int l4,int l5)
#else
  void XGOPW(int *narg,char *SD0,int l0,char *SD1,int l1,char *SD2,
             int l2,char *SD3,int l3,char *SD4,int l4,char *SD5,int l5)
#endif
#endif

{

  int i,j;
  int leng[6];

  char *argv[6];

  printf(" in xgopw_, narg: %d\n",*narg);

#if STRINGS_AS_STRUCTS
  argv[0] = SD0->Addr;
  argv[1] = SD1->Addr;
  argv[2] = SD2->Addr;
  argv[3] = SD3->Addr;
  argv[4] = SD4->Addr;
  argv[5] = SD5->Addr;

#define MAX(A,B) ((A)>(B)?(A):(B))

#else

  argv[0] = SD0;
  argv[1] = SD1;
  argv[2] = SD2;
  argv[3] = SD3;
  argv[4] = SD4;
  argv[5] = SD5;

#endif


  /*  decode text parameters and lengths, remove trailing blanks  */

  for ( i=0 ; i < *narg ; i++) {
    leng[i]=strlen(argv[i]);
    for ( j=leng[i]-1; j>1; j-- ) {
      if(argv[i][j] != ' ') {
        //	argv[i][j+1] = '\0';
	break;
      }
    }
    printf("argv[%d] = %s\n",i,argv[i]);
  }
  fflush(stdout);

  xg_open_window(*narg,argv);

}
/* end of fortran - C routine */

static void xg_open_window (int argc,char ** argv)
{
  char *displayname = NULL;
  char *geom = NULL;
  int i;
  XSizeHints hints;
  Window w;
  Pixel fore, back, border;
  int done;
  XGCValues xgcv;
  int geom_result;
  int x,y,h,wd;


  ProgramName = argv[0];
  for (i = 1; i < argc; i++) {
    char *arg = argv[i];

    if (arg[0] == '-') {
      switch (arg[1]) {
      case 'd':			/* -display host:dpy */
	if (++i >= argc) usage ();
	displayname = argv[i];
	continue;
      case 'g':			/* -geometry geom */
	if (++i >= argc) usage ();
	geom = argv[i];
	continue;
      default:
	usage ();
	/* doesn't return */
      }
    } else
      usage ();
  }
  fprintf (stderr, "%s: Opening display '%s'\n",
	   ProgramName, XDisplayName (displayname));

  dpy = XOpenDisplay (displayname);
  if (!dpy) {
    fprintf (stderr, "%s:  unable to open display '%s'\n",
	     ProgramName, XDisplayName (displayname));
    exit (1);
  }

  set_sizehints (&hints, 10, 10, 1000, 750, 10, 10, geom);


  /*    printf (" Size after Parse: %d  %d  %d  %d : 0x%x  \n",
	hints.x, hints.y,
	hints.width, hints.height,
	geom_result);

	set_sizehints (hintp, min_width, min_height,
	defwidth, defheight, defx, defy, geom)
  */

  screen = DefaultScreen (dpy);
  fore = BlackPixel (dpy, screen);
  back = WhitePixel (dpy, screen);
  border = fore ;


  hints.min_width = 100;
  hints.max_width = 1500;
  hints.min_height= 100;
  hints.max_height= 1500;
  hints.width_inc = 1;
  hints.height_inc = 1;


  window_width  = hints.width;
  window_height = hints.height;


  w = XCreateSimpleWindow (dpy, RootWindow (dpy, screen), hints.x, hints.y,
			   hints.width, hints.height, 2, border, back);

  xwindow = w;

  XSetStandardProperties (dpy, w, ProgramName , NULL, (Pixmap) 0,
			  argv, argc, &hints);

  XSelectInput (dpy, w, (KeyPressMask |
			 ExposureMask | StructureNotifyMask |
			 ButtonPressMask | ButtonReleaseMask));


  XMapWindow (dpy, w);

  xgcv.function = GXcopy;

  xgcv.foreground = BlackPixel(dpy,screen);
  xgcv.background = WhitePixel(dpy,screen);

  xgcv.fill_style = FillSolid;

  gc = (GC)XCreateGC(dpy, w,
		     GCFunction | GCForeground | GCBackground
		     | GCFillStyle, &xgcv);

  colormap = DefaultColormap(dpy, screen);

}

/* FORTRAN interface =========================== draw line segment */

void XGDRAWL   (int *ix0, int *iy0, int *ix1, int *iy1)
{
  XDrawLine(dpy, xwindow, gc, *ix0, CHGY(*iy0), *ix1, CHGY(*iy1));
}


/* FORTRAN interface ============================== set color map */

void XGSETCM   (int *index,int *ir,int *ig,int *ib)
{

  colors[*index].red = ((*ir)<<8);
  colors[*index].green = ((*ig)<<8);
  colors[*index].blue = ((*ib)<<8);

  if ( XAllocColor(dpy,colormap,&colors[*index]) == 0 ) {
    printf("*** XAllocColor died\n");
    exit(1);
  }

}

/* FORTRAN interface ================================= set color index */
void XGSETCOL  (int *index)
{
  XGCValues xgcv;
  xgcv.foreground = colors[*index].pixel;
  XChangeGC(dpy, gc, GCForeground, &xgcv);

  /* SetForeground(dpy, gc, colors[index].pixel); */
}

/* To send back the rgb values of the present index value ================...*/

void XGSENDCOL (int *index,int *ir,int *ib,int *ig)
{
  XGCValues xgcv;
  xgcv.foreground = colors[*index].pixel;
  XChangeGC(dpy, gc, GCForeground, &xgcv);

  *ir = colors[*index].red;
  *ib = colors[*index].blue;
  *ig = colors[*index].green;

  return;
  /* SetForeground(dpy, gc, colors[index].pixel); */
}


/*  FORTRAN INTERFACE ======== Get the position of the mouse key */

void XMOUSEPOS (int *imx,int *imy)
{
  XEvent event;

  while ( XCheckWindowEvent( dpy, xwindow, ButtonPress, &event ) ){
    /* empty loop to discard old button events */
  }

  /* get next button event in the given window */
  XWindowEvent ( dpy, xwindow, ButtonPress, &event );

  *imx = event.xbutton.x;
  *imy = CHGY(event.xbutton.y);
}

/* FORTRAN interface ============ Get number of events in Queue,no wait */

int XGLENGTH  ( void )
{
  return ( QLength(dpy) );
}

/* FORTRAN interface ============= Get next event - wait if none */

void XGNXTEV   (int *ev_type,int *ev_val,int *ev_x,int *ev_y)
{

  XEvent event;
  XKeyEvent *kep;
  XButtonEvent *bep;
  XConfigureEvent *cep;
  XExposeEvent *exep;
  KeySym ks;
  char *ksname;
  int done;


  for (done = 0; !done; ) {
    XNextEvent (dpy, &event);
    switch (event.type) {
    case KeyPress:
      kep = (XKeyEvent *) &event;
      ks = XLookupKeysym (kep, 0);
      shift = kep->state &
	(ShiftMask | LockMask | ControlMask | Mod1Mask);
      switch (ks) {
      case 65507 :;
      case 65505 :;
      case 65509 :;
      case 65513 :;
	done = 0 ; continue;
      default : ;
      }
      /*      adjust for shift and shift lock - VT100 keyboard */
      if ((ks<='z' && ks >='a') &&
	  (!(shift&ShiftMask) != !(shift&LockMask)) )
	ks = toupper(ks);
      if (shift&ShiftMask) {
	switch ((char)ks) {
	case '1': ks='!'; break;
	case '2': ks='@'; break;
	case '3': ks='#'; break;
	case '4': ks='$'; break;
	case '5': ks='%'; break;
	case '6': ks='^'; break;
	case '7': ks='&'; break;
	case '8': ks='*'; break;
	case '9': ks='('; break;
	case '0': ks=')'; break;
	case '-': ks='_'; break;
	case '=': ks='+'; break;
	case '`': ks='~'; break;
	case '[': ks='{'; break;
	case ']': ks='}'; break;
	case ';': ks=':'; break;
	case '\'':ks='"'; break;
	case '\\':ks='|'; break;
	case '/': ks='?'; break;
	case '<': ks='>'; break;
	default:;
	}
      }

      done = 1;

      *ev_type = 0;
      *ev_val  = (int)((ks & 0x10ff)  |
		       ( SHIFT  ) << 8) | (ks>0xff)<<7 ;
      /* adjust for num_pad */
      ks &= 0xff;
      if (ks <= 0x1f && ks>0) *ev_val &= 0xfff ;
      ks &= 0x7f ;
      if ((ks<=0x0d && ks>0) || (ks<=0x2e && ks>=0x2c) ||
	  (ks<='9' && ks>='0'))
	*ev_val  &= 0xff7f;
      *ev_x    = (int)kep->x;
      *ev_y    = (int)CHGY(kep->y);

      continue;


      /*	  case KeyRelease:
                  kep = (XKeyEvent *) &event;
                  ks = XLookupKeysym (kep, 0);

                  printf ("%s\n\t(%d,%d) root:(%d,%d), state 0x%x, ",
                  "KeyRelease",
                  kep->x, kep->y, kep->x_root, kep->y_root,
                  kep->state);
                  printf ("keycode %d = 0x%x,\n", kep->keycode, kep->keycode);
                  ksname = XKeysymToString (ks);
                  printf ("\tkeysym %d = 0x%x (%s)\n",
                  ks, ks, ksname ? ksname : "?");
                  printf("integer of key: %d\n",(int)ks);
                  done = 0;
                  .........  only shift keys are tested for release ........
                  continue;       .... event locked out ....*/

    case ButtonPress:
      bep = (XButtonEvent *) &event;
      /*	    printf ("%s\n\t(%d,%d) root:(%d,%d), state 0x%x, ",
                    "ButtonPress",
                    bep->x, bep->y, bep->x_root, bep->y_root,
                    bep->state);
                    printf ("button %d = 0x%x\n", bep->button, bep->button);  */
      shift = bep->state &
	(ShiftMask | LockMask | ControlMask | Mod1Mask);

      done = 1;

      *ev_type = 2;
      *ev_val  = (int)bep->button |
	( SHIFT ) << 8;
      *ev_x    = (int)bep->x;
      *ev_y    = (int)CHGY(bep->y);

      continue;
    case ButtonRelease:
      bep = (XButtonEvent *) &event;
      shift = bep->state &
	(ShiftMask | LockMask | ControlMask | Mod1Mask);
      /*	    printf ("%s\n\t(%d,%d) root:(%d,%d), state 0x%x, ",
                    "ButtonRelease",
                    bep->x, bep->y, bep->x_root, bep->y_root,
                    bep->state);
                    printf ("button %d = 0x%x\n", bep->button, bep->button);  */

      done = 1;

      *ev_type = 3;
      *ev_val  = (int)bep->button |
	( SHIFT ) << 8;
      *ev_x    = (int)bep->x;
      *ev_y    = (int)CHGY(bep->y);

      continue;

    case ConfigureNotify:
      cep = (XConfigureEvent *) &event;
      /*	    printf ("ConfigureNotify\n");  */

      done = 1;

      *ev_type = 4;
      *ev_val  = 0;
      window_width =
	*ev_x    = (int)event.xconfigure.width;
      window_height =
	*ev_y    = (int)event.xconfigure.height;
      continue;

    case Expose:
    case GraphicsExpose:
    case NoExpose:
      /* get rid of all other Expose events on queue */
      while (XCheckTypedEvent (dpy, Expose, &event));

      exep = (XExposeEvent *) &event;
      /*          printf ("Expose  %d \n",event.type);  */

      done = 1;

      *ev_type = 5;
      *ev_val  = (int)exep->count;
      *ev_x    = 0;
      *ev_y    = 0;
      continue;

      /*	  case VisibilityNotify:
                  vep = (XVisibilityEvent *) &event;
                  printf ("Visibility  %d \n",event.type);

                  done = 1;

                  *ev_type = 6;
                  *ev_val  = (int)vep->state;
                  continue;
      */
    default:
      /*          printf ("Unhandled event type %d\n", event.type);  */
      continue;
    }
  }

}

/* FORTRAN interface ===================set clip region on screen */

void XGSETCLIP (int *xl,int *yl,int *xu,int *yu)
{

  XRectangle rects[1];

  rects[0].x = MINA(*xl,*xu);
  rects[0].y = CHGY(MAXA(*yl,*yu));
  rects[0].width = GABS(*xu-*xl);
  rects[0].height = GABS(*yu-*yl);

  XSetClipRectangles(dpy, gc, 0,0, rects, 1, YXBanded);
}

/* FORTRAN interface ============================= draw rectangle */

void XGDRRECT  (int *xl,int *yl,int *xu,int *yu)
{

  XRectangle rects[1];

  rects[0].x = MINA(*xl,*xu);
  rects[0].y = CHGY(MAXA(*yl,*yu));
  rects[0].width = GABS(*xu-*xl);
  rects[0].height = GABS(*yu-*yl);

  XDrawRectangles(dpy, xwindow, gc, rects, 1);
}

/* FORTRAN interface =============================== draw filled rectangle */

void XGDRFRECT (int *xl,int *yl,int *xu,int *yu)
{
  xg_draw_fill_rectangle(*xl,*yl,*xu,*yu);
}

static void xg_draw_fill_rectangle(int xl,int yl,int xu,int yu)
{
  XRectangle rects[1];

  rects[0].x = MINA(xl,xu);
  rects[0].y = CHGY(MAXA(yl,yu));
  rects[0].width = GABS(xu-xl);
  rects[0].height = GABS(yu-yl);

  XFillRectangles(dpy, xwindow, gc, rects, 1);
}

/* FORTRAN interface ============================ draw polygon */

void XGDRPOLY  (int *npts,int ix[],int iy[])
{
  XPoint points[256];
  int i;

  for ( i=0; i<MINA(*npts,256); i++ ) {
    points[i].x = (short)ix[i];
    points[i].y = (short)CHGY(iy[i]);
  }

  XDrawLines(dpy, xwindow, gc, points, MINA(*npts,256), CoordModeOrigin);
}


/* FORTRAN interface ================================== draw filled poly */

void XGDRFPOLY (int *npts,int ix[],int iy[])
{
  XPoint points[256];
  int i;

  for ( i=0; i<MINA(*npts,256); i++ ) {
    points[i].x = (short)ix[i];
    points[i].y = (short)CHGY(iy[i]);
  }

  XFillPolygon(dpy, xwindow, gc, points, MINA(*npts,256),
	       Nonconvex, CoordModeOrigin);

}

/* FORTRAN interface ============================= set text location */

void XGSETTL   (int *x,int *y)
{

  textloc_x = *x;
  textloc_y = CHGY(*y);
}

/* FORTRAN interface ================================== draw text */

#if STRINGS_AS_STRUCTS
void XGDRTXT(Str_Desc *SD)
{
  xg_draw_text(SD->Addr);
}
#else
void XGDRTXT(char *SD,int length)
{
  xg_draw_text(SD);
}
#endif

static void xg_draw_text(char *string)
{
  XDrawString(dpy, xwindow, gc, textloc_x, textloc_y, string, strlen(string) );
}


/* FORTRAN interface ================================== close window */

void XGCLOSW   ( void )
{
  XCloseDisplay (dpy);
}

/* ################################################################## */

static void set_sizehints (XSizeHints *hintp, int min_width, int min_height,
                           int defwidth, int defheight,int defx,int defy,char *geom)
{
  int geom_result;

  /* set the size hints, algorithm from xlib xbiff */

  hintp->width = hintp->min_width = min_width;
  hintp->height = hintp->min_height = min_height;
  hintp->flags = PMinSize;
  hintp->x = hintp->y = 0;
  geom_result = NoValue;
  if (geom != NULL) {

    geom_result = XParseGeometry (geom, &hintp->x, &hintp->y,
				  (unsigned int *)&hintp->width,
				  (unsigned int *)&hintp->height);

    if ((geom_result & WidthValue) && (geom_result & HeightValue)) {
      hintp->width = MAXA (hintp->width, hintp->min_width);
      hintp->height = MAXA (hintp->height, hintp->min_height);
      hintp->flags |= USSize;
    }
    if ((geom_result & XValue) && (geom_result & YValue)) {
      hintp->flags += USPosition;
    }
  }
  if (!(hintp->flags & USSize)) {
    hintp->width = defwidth;
    hintp->height = defheight;
    hintp->flags |= PSize;
  }
  if (!(hintp->flags & USPosition)) {
    hintp->x = defx;
    hintp->y = defy;
    hintp->flags |= PPosition;
  }

  if (geom_result & XNegative) {
    hintp->x = DisplayWidth (dpy, DefaultScreen (dpy)) + hintp->x -
      hintp->width;
  }
  if (geom_result & YNegative) {
    hintp->y = DisplayHeight (dpy, DefaultScreen (dpy)) + hintp->y -
      hintp->height;
  }
  return;
}

static void usage ( void )
{
  fprintf (stderr,
	   "usage:  %s [-display host:dpy] [-geometry geom]\n", ProgramName);
  exit (1);
}

/* FORTRAN interface ============================= clear window */

void XGCLEAR   ( void )
{
  /*	XClearWindow(dpy,xwindow);
   */
  xg_draw_fill_rectangle(1,1,window_width,window_height);
}

/* FORTRAN interface =============================== get mouse position */

void XGMOUSE   (int *x,int *y)
{
  Window root,child;
  int rx,ry,wx,wy;
  unsigned int keys_buttons;
  XQueryPointer(dpy,xwindow,&root,&child,&rx,&ry,&wx,&wy,
		&keys_buttons);
  *x = wx;
  *y = CHGY(wy);
}

/* FORTRAN interface ================================= get screen info */

void XGINFO    (int *x,int *y,int *mxcol,int *mxx,int *mxy)
{
  int scr_num,i1,i2,i3,i4;

  *x = window_width;
  *y = window_height;

  scr_num = DefaultScreen(dpy);

  *mxx = DisplayWidth(dpy,scr_num);
  *mxy = DisplayHeight(dpy,scr_num);

  *mxcol = DisplayCells(dpy,scr_num);

  XSync(dpy,0);     /* extra service : sync and delete queued events */
  while (QLength(dpy)) XGNXTEV ( &i1, &i2, &i3, &i4 );
  /* they are processed to get all Configure events */

}

void XGPROGRESS ( void )
{
  static int i = 0;
  static char rot [] = { "|/-\\" };

  fprintf(stderr, "%c%c%c", 13, rot[i], 13 );
  i = (i+1)%4;
}
void XGDOTS ( void )
{

  fprintf(stderr, "." );
}

#endif
