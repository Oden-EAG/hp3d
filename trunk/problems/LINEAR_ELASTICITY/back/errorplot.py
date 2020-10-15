"""Generate errorplots.

Useful for postprocessing hp3d output.

Created by B. Keith, 04/15
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def main():
    """Define user-controlled parameters."""

    # ------------------------------------------------------------------
    #                             PARAMETERS
    h0=float(1)  # initial h
    init=2       # first data point to measure rate of convergence
    nlFlag=0     # least squares curve fitting flag

    # define error file and line to read
    errorfile="files/errorlogs/errorlog.txt"
    line=4      # REMEMBER: start counting at 0

    # ------------------------------------------------------------------

    # Extract error as array
    # errorarray = ReadError(errorfile, line)
    errorarray = np.array([546,256,127,38.8])

    # Plot error points
    fig, ax = PlotError(h0, errorarray)

    # Calculate rate of convergence and plot
    if nlFlag == 1:
        print('-- Nonlinear fitting --')
        NonlinearFit(fig, ax, init)
    else:
        print('-- Linear fitting --')
        ax.set_xscale("log")
        ax.set_yscale("log")
        LinearFit(fig,ax,init)

    # Format plot
    PlotOptions()

    # Title
    ax.set_title('Error for L-shaped domain, $p=3$', fontsize='x-large')
    ax.set_xlabel('$h$', fontsize='large')
    ax.set_ylabel('$b(u-u_h,u-u_h) = \sqrt{2(J(u_h)-J(u))}$', fontsize='large')

    # Legend
    # plt.legend()

    # Display figure
    # plt.savefig('files/figures/Energy_plot_L_shape_p_5')
    plt.show()


#-----------------------------------------------------------------------

def PlotOptions():
    # use LaTeX
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', serif='')


#-----------------------------------------------------------------------
#                         Error calculations
#-----------------------------------------------------------------------

def ReadError(errorfile,line):
    # open file
    ef = open(errorfile, "r+")
    print('Reading the error file "%s".\n' % ef.name)
    print('Extracting the following data...\n')

    # read line and save as numpy array
    linelist = ef.readlines()
    print(linelist[line-1]+linelist[line])
    err = str.split(linelist[line])
    err = np.array([float(data) for data in err])

    # close file
    ef.close()
    return err

def PlotError(h0,error):
    # create plot environment
    fig, ax = plt.subplots()

    # generate h array assuming uniform h refinements
    L = error.size
    h = h0*np.power(0.5,np.arange(L),dtype=np.float)
    print(' h = '+str(h)[1:-1]+'\n')

    ax.plot(h,error,'o-')


    # plt.show()

    return fig, ax


#-----------------------------------------------------------------------
#                         Nonlinear fitting
#-----------------------------------------------------------------------

def NonlinearFit(fig,ax,init):

    # extract h, error
    line = ax.get_lines()[0]
    xdata = line.get_xdata()[init:]
    ydata = line.get_ydata()[init:]

    # nonlinear least squares fit
    popt, pcov = curve_fit(_func, xdata, ydata, p0=[1,1,ydata[-1]])

    print('Looking for a fit to C*h^s + J0 ...')
    print('C = %2.5f, s = %2.5f, J0 = %2.6f' % (popt[0],popt[1],popt[2]))

    # # set the plotting to 'hold' so the plotting commands won't
    # #overwrite one another
    # ax.hold(True)

    L = xdata.size
    x = xdata[0]*np.power(0.5,np.arange(0,L,.1),dtype=np.float)

    ax.plot(x,_func(x,popt[0],popt[1],popt[2]),'-')


    ax.hold(False)

    ### SACRED NUMBERS FOR L-SHAPED DOMAIN PROBLEM ###
    #  Dirichlet + Neumann BC's
    J0 =-0.05713414
    # #  Dirichlet BC's
    # J0 = 0.04103559

    ax.plot(xdata,np.sqrt(2*(ydata-popt[2])),'o-')
    # ax.plot(xdata,np.sqrt(2*(ydata-J0)),'o-')
    ax.set_xscale("log")
    ax.set_yscale("log")
    LinearFit(fig,ax,init)


def _func(h,C,s,J0):
    return C*np.power(h,s) + J0


#-----------------------------------------------------------------------
#                          Linear fitting
#-----------------------------------------------------------------------

def LinearFit(fig,ax,init=0):

    # extract h, error
    line = ax.get_lines()[0]
    xdata = line.get_xdata()[init:]
    ydata = line.get_ydata()[init:]

    # linear least squares fit (of log data)
    coef = np.polyfit(np.log(xdata),np.log(ydata),1)

    print('Rate of convergence : %2.5f' % coef[0])
    # print('y-intercept of fitting curve : %2.5f' % coef[1])

    # set the plotting to 'hold' so the plotting commands won't
    # overwrite one another
    ax.hold(True)

    epsX = (xdata[-2]-xdata[-1])/2
    epsY = (ydata[-2]-ydata[-1])/10
    SlopeMarker((xdata[-2]-epsX,ydata[-1]+epsY), (round(coef[0]*100)/100, 1), size_frac=0.15, ax=ax)

    # L = xdata.size
    # x = xdata[0]*np.power(0.5,np.arange(0,L,.1),dtype=np.float)
    # plt.plot(x,np.exp(coef[1])*np.power(x,coef[0]),'-')

def SlopeMarker(origin, slope, size_frac=0.1, pad_frac=0.1, ax=None,
                 invert=False):
    """Plot triangular slope marker labeled with slope.

    Parameters
    ----------
    origin : (x, y)
        tuple of x, y coordinates for the slope
    slope : float or (rise, run)
        the length of the slope triangle
    size_frac : float
        the fraction of the xaxis length used to determine the size of the slope
        marker. Should be less than 1.
    pad_frac : float
        the fraction of the slope marker used to pad text labels. Should be less
        than 1.
    invert : bool
        Normally, the slope marker is below a line for positive slopes and above
        a line for negative slopes; `invert` flips the marker.
    """
    if ax is None:
        ax = plt.gca()

    if np.iterable(slope):
        rise, run = slope
        slope = float(rise) / run
    else:
        rise = run = None

    x0, y0 = origin
    xlim = ax.get_xlim()
    dx_linear = size_frac * (xlim[1] - xlim[0])
    dx_decades = size_frac * (np.log10(xlim[1]) - np.log10(xlim[0]))

    if invert:
        dx_linear = -dx_linear
        dx_decades = -dx_decades

    if ax.get_xscale() == 'log':
        log_size = dx_decades
        dx = _LogDistance(x0, log_size)
        x_run = _TextPosition(x0, log_size/2., scale='log')
        x_rise = _TextPosition(x0+dx, dx_decades*pad_frac, scale='log')
    else:
        dx = dx_linear
        x_run = _TextPosition(x0, dx/2.)
        x_rise = _TextPosition(x0+dx, pad_frac * dx)

    if ax.get_yscale() == 'log':
        log_size = dx_decades * slope
        dy = _LogDistance(y0, log_size)
        y_run = _TextPosition(y0, -dx_decades*slope*pad_frac, scale='log')
        y_rise = _TextPosition(y0, log_size/2., scale='log')
    else:
        dy = dx_linear * slope
        y_run = _TextPosition(y0, -(pad_frac * dy))
        y_rise = _TextPosition(y0, dy/2.)

    x_pad = pad_frac * dx
    y_pad = pad_frac * dy

    va = 'top' if y_pad > 0 else 'bottom'
    ha = 'left' if x_pad > 0 else 'right'
    if rise is not None:
        ax.text(x_run, y_run, str(run), va=va, ha='center')
        ax.text(x_rise, y_rise, str(rise), ha=ha, va='center')
    else:
        ax.text(x_rise, y_rise, str(slope), ha=ha, va='center')

    ax.add_patch(_SlopeTriangle(origin, dx, dy))


def LogDisplace(x0, dx_log=None, x1=None, frac=None):
    """Return point displaced by a logarithmic value.

    For example, if you want to move 1 decade away from `x0`, set `dx_log` = 1,
    such that for `x0` = 10, we have `displace(10, 1)` = 100

    Parameters
    ----------
    x0 : float
        reference point
    dx_log : float
        displacement in decades.
    x1 : float
        end point
    frac : float
        fraction of line (on logarithmic scale) between x0 and x1
    """
    if dx_log is not None:
        return 10**(np.log10(x0) + dx_log)
    elif x1 is not None and frac is not None:
        return 10**(np.log10(x0) + frac * np.log10(float(x1)/x0))
    else:
        raise ValueError('Specify `dx_log` or both `x1` and `frac`.')


def _LogDistance(x0, dx_decades):
    return LogDisplace(x0, dx_decades) - x0

def _TextPosition(x0, dx, scale='linear'):
    if scale == 'linear':
        return x0 + dx
    elif scale == 'log':
        return LogDisplace(x0, dx)
    else:
        raise ValueError('Unknown value for `scale`: %s' % scale)


def _SlopeTriangle(origin, dx, dy, ec='none', fc='0.8', **poly_kwargs):
    """Return Polygon representing slope.
          /|
         / | dy
        /__|
         dx
    """
    verts = [np.asarray(origin)]
    verts.append(verts[0] + (dx, 0))
    verts.append(verts[0] + (dx, dy))
    return plt.Polygon(verts, ec=ec, fc=fc, **poly_kwargs)

main()
